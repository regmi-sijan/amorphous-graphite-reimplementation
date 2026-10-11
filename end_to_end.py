"""Run the complete Python-led workflow, or resume completed stages safely.

Usage: .venv/bin/python end_to_end.py --config configs/end_to_end_smoke.json --out results/end_to_end_smoke
Stages: random carbon -> real GAP MD -> relaxation -> layers/rings/RDF -> periodic
PBE/HSE06 -> density/orbitals/DOS -> KG/SPC -> graphite comparison -> report.
Optional zero-pressure PBE+D3 relaxation is implemented, but costly on a laptop.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil
import time
import uuid
import numpy as np
from ase.io import read,write,iread
from ase.optimize import FIRE
from reproduce import ROOT,dump,prepare,simulate,make_calculator,geometry,layering,rdf_3d
from agworkflow.structure import graphite,identify_layers,in_plane_rdf,periodic_rings
from agworkflow.electronic import run_scf,analyze_electrons
from agworkflow.relaxation import dispersion,zero_pressure_relax
from agworkflow.plotting import electronic_plots,structural_plots,evolution_plot


def structural_analysis(atoms,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    layers=identify_layers(atoms)
    r,g=rdf_3d(atoms)
    np.savetxt(out/'rdf_3d.csv',np.column_stack([r,g]),delimiter=',',header='r_A,g_3d',comments='')
    rings=periodic_rings(atoms)
    summary={'layers':layers,'rings_chordless_periodic':rings,'in_plane_rdf_available':layers['detected']}
    if layers['detected']:
        r,g=in_plane_rdf(atoms,layers)
        np.savetxt(out/'rdf_in_plane.csv',np.column_stack([r,g]),delimiter=',',header='r_A,g_2d',comments='')
    else:
        summary['in_plane_rdf_reason']='No validated layers; a 2D sheet normalization is not physically justified'
    dump(out/'summary.json',summary)
    write(out/'structure.extxyz',atoms)
    structural_plots(out)
    return summary


def gap_relax(atoms,cfg,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    a=atoms.copy();a.set_momenta(np.zeros((len(a),3)));a.calc,_=make_calculator(cfg)
    optimizer=FIRE(a,logfile=str(out/'relax.log'),trajectory=str(out/'relax.traj'))
    converged=optimizer.run(fmax=cfg['relax_fmax_eV_A'],steps=cfg['relax_steps'])
    write(out/'final.extxyz',a)
    result={'converged':bool(converged),'optimizer':'FIRE; differs from paper CG',
            'max_force_eV_A':float(np.linalg.norm(a.get_forces(),axis=1).max()),
            'steps':optimizer.nsteps,'cell':'fixed','force_method':'2017 GAP'}
    dump(out/'summary.json',result)
    if not converged and not cfg.get('allow_unconverged_geometry',False):
        raise RuntimeError('Geometry relaxation did not converge; increase the iteration budget')
    return result


class Workflow:
    def __init__(self,config,out,resume=False):
        self.cfg=json.loads(Path(config).read_text());self.out=Path(out);self.resume=resume
        sources=[ROOT/'end_to_end.py',ROOT/'reproduce.py',*sorted((ROOT/'agworkflow').glob('*.py'))]
        current_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
        digest=hashlib.sha256(json.dumps(self.cfg,sort_keys=True).encode()).hexdigest()
        if self.out.exists():
            if not resume: raise FileExistsError('Output exists; use --resume with the same config or a new directory')
            old=json.loads((self.out/'manifest.json').read_text())
            if digest!=old['config_sha256']: raise ValueError('Config changed; use a fresh output directory')
            self.manifest=old
        else:
            self.out.mkdir(parents=True)
            self.manifest={'config_sha256':digest,'status':'running','completed':{},
                'versions':{m:importlib.metadata.version(m) for m in ['ase','quippy-ase','pyscf','dftd3','numpy','scipy']},
                'source_sha256':current_hashes}
            dump(self.out/'config.json',self.cfg)
        invocations=self.manifest.setdefault('invocations',[])
        snapshot=self.out/'source_snapshots'/f'invocation_{len(invocations):03d}'
        for source in sources:
            dest=snapshot/source.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(source,dest)
        invocations.append({'unix_time':time.time(),'source_sha256':current_hashes,
                            'snapshot':str(snapshot.relative_to(self.out))})
        self.invocation=len(invocations)-1
        self.save()
    def save(self): dump(self.out/'manifest.json',self.manifest)
    def stage(self,name,fn):
        path=self.out/name
        if name in self.manifest['completed']:
            print(f'Reusing completed stage: {name}',flush=True)
            return self.manifest['completed'][name]['result']
        if path.exists():
            if not self.resume: raise FileExistsError(path)
            backup=self.out/'failed_attempts';backup.mkdir(exist_ok=True)
            path.rename(backup/(name.replace('/','_')+'_'+uuid.uuid4().hex[:8]))
        path.parent.mkdir(parents=True,exist_ok=True)
        self.manifest.update(status='running',active_stage=name);self.save()
        print(f'Starting {name}',flush=True);start=time.monotonic()
        try: result=fn(path)
        except BaseException as e:
            self.manifest.update(status='interrupted_or_failed',error=f'{type(e).__name__}: {e}')
            self.save();raise
        self.manifest['completed'][name]={'seconds':time.monotonic()-start,'result':result,
                                          'invocation':self.invocation}
        self.manifest.pop('error',None);self.save()
        return result

    def run(self):
        c=self.cfg;s=c['structural'];e=c['electronic'];t=c['transport']
        if max(t['photon_eV'])+6*max(t['broadening_eV'])>t['energy_window_eV']:
            raise ValueError('Transport energy window is too small for the requested spectrum/broadening')
        def md(path):
            prepare(s,path);simulate(path)
            rows=[]
            for a in iread(path/'trajectory.traj'):
                rows.append({'time_ps':a.info['time_ps'],**geometry(a),**layering(a)})
            dump(path/'structural_time_series.json',rows)
            evolution_plot(path)
            # Include momenta and physical coordinates for independent analysis.
            write(path/'trajectory.extxyz',list(iread(path/'trajectory.traj')))
            return json.loads((path/'provenance.json').read_text())
        self.stage('sample_md',md)
        self.stage('sample_relax',lambda p:gap_relax(read(self.out/'sample_md/final.extxyz'),s,p))
        def reference(path):
            path.mkdir();a=graphite(c['reference_repeat'],s['density_g_cm3'])
            write(path/'final.extxyz',a)
            return {'n_atoms':len(a),'type':'AB graphite reference at requested density'}
        self.stage('reference',reference)
        cases={'sample':self.out/'sample_relax/final.extxyz','reference':self.out/'reference/final.extxyz'}
        if c.get('defect_reference',False):
            def defect(path):
                a=read(cases['reference']);del a[[0,1]]
                return gap_relax(a,s,path)
            self.stage('defect_relax',defect)
            cases['defect']=self.out/'defect_relax/final.extxyz'
        summaries={}
        for case,path in cases.items():
            atoms=read(path)
            zp=None
            if c['zero_pressure']['enabled']:
                zp=self.stage(f'{case}_zero_pressure',lambda p,a=atoms:zero_pressure_relax(a,e,c['d3'],c['zero_pressure'],p))
                if not zp['converged'] and not c.get('smoke_test',False):
                    raise RuntimeError('PBE+D3 zero-pressure relaxation is not converged')
                # Preserve the fixed-density comparison. The zero-pressure result
                # is a separate spacing/geometry calculation, as in the paper.
                relaxed=read(self.out/f'{case}_zero_pressure/final.extxyz')
                self.stage(f'{case}_zero_pressure_structure',lambda p,a=relaxed:structural_analysis(a,p))
            structure=self.stage(f'{case}_structure',lambda p,a=atoms:structural_analysis(a,p))
            pbe=self.stage(f'{case}_pbe',lambda p,a=atoms:run_scf(a,e,p,'PBE'))
            def post(path,case=case,layers=structure['layers']):
                path.mkdir()
                q=analyze_electrons(self.out/f'{case}_pbe',t,do_transport=True)
                electronic_plots(self.out/f'{case}_pbe',layers)
                tolerance=e['field_relative_tolerance']
                checks={k:q[k]<=tolerance for k in ['density_relative_integration_error','momentum_grid_relative_error',
                                                   'spc_sum_relative_error','active_orbital_overlap_grid_max_error']}
                q['field_accuracy_checks']=checks
                dump(path/'summary.json',q)
                if e['require_field_accuracy'] and not all(checks.values()):
                    raise RuntimeError('Field quadrature checks failed; increase field_mesh and rerun in a fresh directory')
                return q
            quality=self.stage(f'{case}_postprocess',post)
            if e.get('hse06',True):
                self.stage(f'{case}_hse06',lambda p,a=atoms:run_scf(a,e,p,'HSE06'))
                def hsepost(p,case=case,layers=structure['layers']):
                    p.mkdir();q=analyze_electrons(self.out/f'{case}_hse06',t,do_transport=False)
                    electronic_plots(self.out/f'{case}_hse06',layers);dump(p/'summary.json',q)
                    if e['require_field_accuracy'] and q['density_relative_integration_error']>e['field_relative_tolerance']:
                        raise RuntimeError('HSE06 density grid accuracy failed; increase field_mesh')
                    return q
                self.stage(f'{case}_hse06_postprocess',hsepost)
            d3=dispersion(atoms,c['d3'])
            summaries[case]={'structure':structure,'pbe':pbe,'grid_quality':quality,
                            'zero_pressure_relaxation':zp,
                            'd3_correction_eV':d3,'pbe_d3_energy_eV_atom':pbe['energy_eV_atom']+d3/len(atoms)}
        self.report(summaries)
        self.manifest.update(status='completed',active_stage=None);self.save()
        print(f'Completed pipeline. Report: {self.out / "REPORT.md"}',flush=True)

    def report(self,summaries):
        base=summaries['reference'];base_sigma=base['grid_quality']['sigma_iso_at_spc_frequency_S_m']
        comparisons={}
        lines=[]
        for name,s in summaries.items():
            sigma=s['grid_quality']['sigma_iso_at_spc_frequency_S_m']
            ratio=sigma/base_sigma if abs(base_sigma)>1e-8 else None
            comparisons[name]={'delta_pbe_energy_eV_atom':s['pbe']['energy_eV_atom']-base['pbe']['energy_eV_atom'],
                'isotropic_sigma_ratio':ratio,'ratio_reason':None if ratio is not None else 'Reference conductivity is numerically near zero',
                'in_plane_ratio':None,'in_plane_ratio_reason':'See conductivity tensor; evaluated below only if both layer normals are established'}
            data=np.load(self.out/f'{name}_pbe/conductivity.npz')
            reference=np.load(self.out/'reference_pbe/conductivity.npz')
            if s['structure']['layers']['detected'] and base['structure']['layers']['detected']:
                n=np.array(s['structure']['layers']['normal']);nr=np.array(base['structure']['layers']['normal'])
                sig=data['tensor_S_m'];sr=reference['tensor_S_m']
                parallel=(np.trace(sig,axis1=2,axis2=3)-np.einsum('i,abij,j->ab',n,sig,n))/2
                pref=(np.trace(sr,axis1=2,axis2=3)-np.einsum('i,abij,j->ab',nr,sr,nr))/2
                ratios=np.divide(parallel,pref,out=np.full_like(parallel,np.nan),where=abs(pref)>1e-8)
                np.savez_compressed(self.out/f'{name}_in_plane_comparison.npz',photon_eV=data['photon_eV'],
                    broadening_eV=data['broadening_eV'],sigma_parallel_S_m=parallel,reference_S_m=pref,ratio=ratios)
                comparisons[name]['in_plane_ratio_reason']='Frequency/broadening-resolved results in in_plane_comparison.npz; unresolved denominators are NaN'
            lines.append(f'| {name} | {s["structure"]["layers"]["detected"]} | {s["pbe"]["energy_eV_atom"]:.6f} | {sigma:.6g} | {ratio if ratio is not None else "undefined"} |')
        dump(self.out/'comparison.json',comparisons);dump(self.out/'summary.json',summaries)
        warnings=[]
        if self.cfg.get('smoke_test'):
            warnings.append('This is a small software test, not a reproduction-quality simulation. Numerical conductivity ratios here have no validated physical interpretation.')
        if not summaries['sample']['structure']['layers']['detected']:
            warnings.append('The simulated sample has no validated graphitic layers. Its in-plane RDF and in-plane conductivity ratio are intentionally unavailable.')
        if not all(all(s['grid_quality']['field_accuracy_checks'].values()) for s in summaries.values()):
            warnings.append('At least one grid-quadrature check failed the configured tolerance. Increase the field grid before scientific interpretation.')
        if not self.cfg['zero_pressure']['enabled']:
            warnings.append('Zero-pressure PBE+D3 relaxation was disabled in this config; the complete implementation is available as an optional expensive stage.')
        self.out.joinpath('REPORT.md').write_text('''# End-to-end workflow results

This report contains outputs of actual GAP and periodic all-electron DFT calculations.
Software execution is distinct from reproducing the paper's numerical findings.

'''+ '\n'.join('- '+w for w in warnings)+'''

| Case | Layers detected | PBE energy (eV/atom) | Isotropic conductivity at SPC frequency (S/m) | Ratio to reference |
|---|---|---|---|---|
'''+ '\n'.join(lines)+'''

## Outputs

- `sample_md/`: initial configuration, actual trajectory, thermodynamics and structural time series.
- `sample_relax/`: fixed-cell relaxation, convergence status and final geometry.
- `*_structure/`: periodic rings, layer identification, 3D RDF, and a 2D RDF when justified.
- `*_pbe/`: converged SCF logs, reusable states, density/orbital cubes, DOS, KG tensors,
  signed and positive SPC cubes, and plots.
- `*_hse06/`: independent HSE06 electronic calculation and density/orbital exports.
- `*_zero_pressure/`: PBE+D3 full-cell relaxation if enabled, with convergence status.
- `comparison.json`: consistently computed energy differences and conductivity comparisons.
- `summary.json`: all numerical results and quadrature quality checks.
- `manifest.json`: stage completion, configuration hash, code hashes, and versions.

## Interpretation and limitations

The DFT backend is PySCF with all-electron Gaussian orbitals, not VASP/PAW. PBE
supplies conductivity states; HSE06 is a separate charge/orbital calculation.
The minimally sized smoke basis and k mesh are intentionally inexpensive and
cannot resolve graphite's quantitative band structure or conductivity.

The KG implementation uses a Gaussian-broadened, divided-difference Greenwood
form with a finite degenerate limit, including diagonal matrix elements. At finite
broadening this is a declared regularization of the delta-function expression.
DC is a finite-broadening estimate; no scattering lifetime has been measured.
The tensor, broadening sweep and grid sum checks must converge before a physical
DC value or the paper's approximately 10^-2 in-plane ratio can be inferred.
The absolute SPC field is the modulus of a complex row sum. Its integral does
not generally equal conductivity; the signed/complex row sum is used for that check.

The RDF and layer tests use documented independent definitions; the original
paper's exact planarity statistic is unavailable. No density peak alone is
treated as proof of layers. Orbitals near the chemical potential are exported
without automatically calling them pi or pi-star. Vacancy structures, when enabled,
are not assumed to reconstruct into a specific 5-8-5 topology.

GAP model identity and some original numerical settings remain unconfirmed.
Missing supplementary information cannot be recovered by writing more code.
The code contains every computational stage, but production-duration trajectories,
independent seeds, basis/k-point/grid/broadening tests and comparison with the
original data are still needed to establish reproduction.
''')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',default=str(ROOT/'configs/end_to_end_smoke.json'))
    p.add_argument('--out',default=str(ROOT/'results/end_to_end_smoke'))
    p.add_argument('--resume',action='store_true')
    a=p.parse_args();Workflow(a.config,a.out,a.resume).run()

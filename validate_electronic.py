"""Actual small DFT integration checks for optional and non-Gamma branches."""
import argparse
import json
from pathlib import Path
from ase import Atoms
from ase.build import bulk
from reproduce import ROOT,dump
from agworkflow.relaxation import zero_pressure_relax
from agworkflow.electronic import run_scf,analyze_electrons

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['cell','kpoints','exchange'])
    p.add_argument('--out',required=True)
    a=p.parse_args();out=Path(a.out)
    if out.exists(): raise FileExistsError('Use a fresh output directory')
    cfg=json.loads((ROOT/'configs/end_to_end_smoke.json').read_text())
    e=cfg['electronic'].copy()
    e.update(precision=1e-6,grid_level=0,scf_mesh=[15]*3,conv_tol_Ha=1e-7,field_mesh=[48]*3)
    atoms=bulk('C','diamond',a=3.567)
    if a.stage=='cell':
        settings=cfg['zero_pressure'].copy();settings.update(steps=1)
        result=zero_pressure_relax(atoms,e,cfg['d3'],settings,out)
        # An optimization step and finite derivatives exercise the whole branch.
        # A one-step calculation is explicitly not a converged relaxation.
        result['test_purpose']='One real PBE+D3 cell/atom optimization step; no convergence claim'
        dump(out/'summary.json',result)
    elif a.stage=='kpoints':
        e['kmesh']=[3,1,1]
        result=run_scf(atoms,e,out,'PBE')
        result['field_checks']=analyze_electrons(out,cfg['transport'])
        dump(out/'integration_test.json',result)
    else:
        out.mkdir(parents=True)
        atoms=Atoms('C2',positions=[[4.29,5,5],[5.71,5,5]],cell=[10]*3,pbc=True)
        result={}
        for method in ['full-minus-long-range','direct']:
            e['screened_exchange']=method
            result[method]=run_scf(atoms,e,out/method,'HSE06')
        difference=abs(result['direct']['energy_eV']-result['full-minus-long-range']['energy_eV'])
        result['absolute_energy_difference_eV']=difference
        result['test_purpose']='Compare screened-exchange evaluations on the same small periodic cell'
        dump(out/'comparison.json',result)
        if difference>.01: raise RuntimeError('Screened-exchange methods differ by more than 0.01 eV; refine fitting precision')
    print(json.dumps(result,indent=2))

"""Periodic all-electron PySCF DFT, grid fields and canonical-momentum transport.

The Gaussian all-electron representation is an explicit alternative to VASP/PAW.
PBE supplies transport states; HSE06 supplies a separate optional charge analysis.
No pseudopotentials are allowed in this backend, so an omitted nonlocal
pseudopotential velocity correction cannot silently contaminate transport.
"""
from pathlib import Path
import time
import json
import numpy as np
from scipy.special import expit
from ase import units
from ase.io import write
from reproduce import dump
from .transport import HA_EV,fermi_level,spectra,transition_weights,spc_row_density,conductivity_tensor


def build_cell(atoms,cfg):
    from pyscf.pbc import gto
    if not np.all(atoms.pbc): raise ValueError('This backend requires 3D periodicity')
    if cfg.get('pseudo') is not None: raise ValueError('Transport backend requires all-electron DFT')
    cell=gto.Cell()
    cell.a=atoms.cell.array
    cell.atom=list(zip(atoms.get_chemical_symbols(),atoms.positions))
    cell.unit='Angstrom'; cell.basis=cfg['basis']; cell.spin=0
    cell.precision=cfg.get('precision',1e-7)
    cell.mesh=cfg.get('scf_mesh',[21,21,21])
    cell.max_memory=cfg.get('max_memory_MB',2000)
    cell.verbose=0; cell.build()
    nk=int(np.prod(cfg['kmesh']))
    matrix_memory=cell.nao_nr()**2*nk*16*12/1e6
    if matrix_memory>cell.max_memory:
        raise MemoryError(f'Estimated dense-matrix working space {matrix_memory:.0f} MB exceeds configured {cell.max_memory} MB')
    return cell


def solve(atoms,cfg,xc='PBE',logfile=None):
    from pyscf.pbc import dft,scf
    cell=build_cell(atoms,cfg)
    if logfile:
        cell.stdout=open(logfile,'w'); cell.verbose=4
    kpts=cell.make_kpts(cfg['kmesh'])
    mf=dft.KRKS(cell,kpts=kpts).density_fit()
    mf.xc=xc
    mf.grids=dft.gen_grid.BeckeGrids(cell)
    mf.grids.level=cfg.get('grid_level',2)
    mf.conv_tol=cfg.get('conv_tol_Ha',1e-8)
    mf.max_cycle=cfg.get('max_cycle',120)
    mf=scf.addons.smearing_(mf,sigma=cfg['smearing_eV']/HA_EV,method='fermi')
    exchange_mode=cfg.get('screened_exchange','full-minus-long-range')
    if exchange_mode not in ('full-minus-long-range','direct'):
        raise ValueError('Unknown screened_exchange evaluation')
    if xc.upper()=='HSE06' and exchange_mode=='full-minus-long-range':
        # erfc(w*r)/r = 1/r - erf(w*r)/r. Reuse full-Coulomb GDF and
        # PySCF's analytic Fourier long-range integrator instead of rebuilding
        # very expensive small-omega short-range three-center integrals.
        # Their finite fitting errors must still be converged independently.
        original_get_k=mf.get_k
        def screened_k(cell=None,dm_kpts=None,hermi=1,kpts=None,kpts_band=None,omega=None):
            args=(cell,dm_kpts,hermi,kpts,kpts_band)
            if omega is not None and omega<0:
                return original_get_k(*args,omega=None)-original_get_k(*args,omega=-omega)
            return original_get_k(*args,omega=omega)
        mf.get_k=screened_k
    begun=time.monotonic()
    mf.kernel()
    if not mf.converged:
        dm=mf.make_rdm1(); mf.level_shift=.2; mf.damp=.25
        mf.kernel(dm0=dm)
    if not mf.converged:
        if logfile: cell.stdout.close()
        raise RuntimeError(f'{xc} SCF did not converge; inspect {logfile}')
    elapsed=time.monotonic()-begun
    return cell,mf,kpts,elapsed


def run_scf(atoms,cfg,out,xc='PBE'):
    out=Path(out); out.mkdir(parents=True,exist_ok=False)
    cell,mf,kpts,elapsed=solve(atoms,cfg,xc,out/'scf.log')
    energy=np.array(mf.mo_energy); coeff=np.array(mf.mo_coeff); occ=np.array(mf.mo_occ)
    mu=fermi_level(energy,cell.nelectron,cfg['smearing_eV']/HA_EV)
    predicted=2*expit((mu-energy)/(cfg['smearing_eV']/HA_EV))
    if np.max(abs(predicted-occ))>1e-5: raise RuntimeError('SCF occupations do not match transport Fermi function')
    overlap=np.asarray(mf.get_ovlp())
    norm_error=max(np.max(abs(c.conj().T@s@c-np.eye(c.shape[1]))) for c,s in zip(coeff,overlap))
    np.savez_compressed(out/'states.npz',energies_Ha=energy,coefficients=coeff,occupations=occ,kpts_bohr_inv=kpts,mu_Ha=mu)
    (out/'cell.json').write_text(cell.dumps())
    dump(out/'settings.json',cfg)
    summary={'engine':'PySCF periodic all-electron Gaussian KRKS/GDF','xc':xc,'converged':bool(mf.converged),
             'energy_eV':float(mf.e_tot*HA_EV),'free_energy_eV':float(mf.e_free*HA_EV),
             'energy_eV_atom':float(mf.e_tot*HA_EV/len(atoms)),
             'mu_eV':float(mu*HA_EV),'nelectrons':cell.nelectron,'nao':cell.nao_nr(),
             'kmesh':cfg['kmesh'],'basis':cfg['basis'],'wall_seconds':elapsed,
             'screened_exchange':cfg.get('screened_exchange','full-minus-long-range') if xc.upper()=='HSE06' else None,
             'orbital_overlap_max_error':float(norm_error),'smearing_eV':cfg['smearing_eV'],
             'note':'SCF convergence is not basis, k-point, grid or thermodynamic convergence'}
    dump(out/'summary.json',summary)
    write(out/'structure.extxyz',atoms)
    # Make the state bundle independent of the current open text log.
    cell.stdout.close()
    return summary


def load_state(out):
    from pyscf.pbc import gto
    out=Path(out)
    cell=gto.loads((out/'cell.json').read_text())
    z=np.load(out/'states.npz')
    cfg=json.loads((out/'settings.json').read_text())
    return cell,z,cfg


def write_cube(path,cell,data,comment):
    """Periodic cube in atomic units; no duplicated endpoint grid planes."""
    mesh=np.array(data.shape); lattice=cell.lattice_vectors()
    with Path(path).open('w') as f:
        f.write(comment+'\nPeriodic grid; origin zero; coordinates and voxel vectors in bohr\n')
        f.write(f'{cell.natm:5d} 0.0 0.0 0.0\n')
        for n,v in zip(mesh,lattice/mesh[:,None]): f.write(f'{n:5d} {v[0]:.10f} {v[1]:.10f} {v[2]:.10f}\n')
        for charge,r in zip(cell.atom_charges(),cell.atom_coords()):
            f.write(f'{charge:5d} {float(charge):.5f} {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}\n')
        flat=np.asarray(data).ravel()
        for offset in range(0,len(flat),6): f.write(' '.join(f'{v:.8e}' for v in flat[offset:offset+6])+'\n')


def analyze_electrons(out,transport_cfg,do_transport=True):
    out=Path(out); cell,z,cfg=load_state(out)
    eigen=z['energies_Ha']; coefficients=z['coefficients']; kpts=z['kpts_bohr_inv']; mu=float(z['mu_Ha'])
    occupations=z['occupations']; nk=len(kpts)
    xc=json.loads((out/'summary.json').read_text())['xc']
    if do_transport and xc.upper()!='PBE':
        raise ValueError('Transport is implemented for local PBE; hybrid nonlocal-velocity terms are not included')
    mesh=np.array(cfg['field_mesh'],int)
    coords=cell.gen_uniform_grids(mesh)
    dv=cell.vol/len(coords)
    rho=np.zeros(len(coords)); selected_fields=np.zeros((4,len(coords)))
    occupied=np.flatnonzero(occupations[0]>1)
    empty=np.flatnonzero(occupations[0]<=1)
    chosen=[int(i) for i in list(occupied[-2:])+list(empty[:2])]
    if len(chosen)!=4: raise ValueError('Need at least two occupied and two empty bands for orbital export')
    band_lists=[]; selected_coeff=[]; selected_energies=[]; momenta=[]
    window=transport_cfg['energy_window_eV']/HA_EV
    for ik,ek in enumerate(eigen):
        ids=np.flatnonzero(abs(ek-mu)<=window)
        if not len(ids): raise ValueError('Transport window contains no states')
        band_lists.append(ids); selected_coeff.append(coefficients[ik][:,ids]); selected_energies.append(ek[ids])
    if do_transport:
        # Integral is derivative of the bra; integration by parts gives p=+i*ipovlp.
        ip=np.asarray(cell.pbc_intor('int1e_ipovlp',comp=3,kpts=kpts))
        for c,b in zip(selected_coeff,ip):
            p=1j*np.einsum('mi,amn,nj->aij',c.conj(),b,c,optimize=True)
            error=np.max(abs(p-p.conj().transpose(0,2,1)))
            if error>1e-6: raise RuntimeError(f'Momentum is not Hermitian: {error}')
            momenta.append(p)
        photon=np.asarray(transport_cfg['photon_eV'])
        eta=np.asarray(transport_cfg['broadening_eV'])
        tensor=spectra(selected_energies,momenta,mu,cfg['smearing_eV']/HA_EV,photon,eta,cell.vol)
        np.savez_compressed(out/'conductivity.npz',photon_eV=photon,broadening_eV=eta,tensor_S_m=tensor)
        rows=[]
        for ib,b in enumerate(eta):
            for iw,w in enumerate(photon): rows.append([b,w,*np.diag(tensor[ib,iw]),np.trace(tensor[ib,iw])/3])
        np.savetxt(out/'conductivity.csv',rows,delimiter=',',header='eta_eV,photon_eV,sigma_xx_S_m,sigma_yy_S_m,sigma_zz_S_m,sigma_iso_S_m',comments='')
    spc_complex=np.zeros(len(coords),complex)
    grid_p=[np.zeros((3,len(ids),len(ids)),complex) for ids in band_lists]
    grid_overlap=[np.zeros((len(ids),len(ids)),complex) for ids in band_lists]
    weights=[]
    if do_transport:
        for ek in selected_energies:
            weights.append(transition_weights(ek,mu,cfg['smearing_eV']/HA_EV,
                transport_cfg['spc_photon_eV']/HA_EV,transport_cfg['spc_broadening_eV']/HA_EV,cell.vol,nk))
    for begin in range(0,len(coords),cfg.get('grid_block',1024)):
        end=min(len(coords),begin+cfg.get('grid_block',1024))
        for ik,kpt in enumerate(kpts):
            ao=cell.pbc_eval_gto('GTOval_sph_deriv1',coords[begin:end],kpt=kpt)
            allpsi=ao[0]@coefficients[ik]
            rho[begin:end]+=np.sum(abs(allpsi)**2*occupations[ik],axis=1)/nk
            if ik==0:
                selected_fields[:,begin:end]=abs(allpsi[:,chosen].T)**2
            if do_transport:
                psi=allpsi[:,band_lists[ik]]
                grad=np.einsum('agm,mi->agi',ao[1:4],selected_coeff[ik],optimize=True)
                spc_complex[begin:end]+=spc_row_density(psi,grad,momenta[ik],weights[ik])
                grid_p[ik]+=(-1j)*np.einsum('gi,agj->aij',psi.conj(),grad)*dv
                grid_overlap[ik]+=psi.conj().T@psi*dv
    cube_rho=rho.reshape(mesh)
    write_cube(out/'density.cube',cell,cube_rho,f'{xc} all-electron density, electrons/bohr^3')
    for j,band in enumerate(chosen):
        write_cube(out/f'orbital_k0_band{band}.cube',cell,selected_fields[j].reshape(mesh),
                   f'{xc} |psi|^2, band {band}, k point 0; not automatically classified as pi or pi-star')
    np.savez_compressed(out/'fields.npz',density_e_bohr3=cube_rho,
        orbitals_e_bohr3=selected_fields.reshape((4,*mesh)),orbital_band_indices=chosen,
        lattice_bohr=cell.lattice_vectors())
    # DOS sums both spins and equally weighted k points. Integral over all energies
    # tends to twice the number of bands as the plotting interval expands.
    emin,emax=transport_cfg.get('dos_range_eV',[-15,10])
    axis=np.linspace(emin,emax,600); width=transport_cfg.get('dos_broadening_eV',.15)
    dos=np.zeros(len(axis))
    for ek in eigen:
        dos+=2*np.exp(-.5*((axis[:,None]-(ek-mu)*HA_EV)/width)**2).sum(axis=1)/(nk*np.sqrt(2*np.pi)*width)
    np.savetxt(out/'dos.csv',np.column_stack([axis,dos]),delimiter=',',header='energy_minus_mu_eV,states_per_eV_cell',comments='')
    quality={'integrated_density_electrons':float(rho.sum()*dv),'expected_electrons':cell.nelectron,
             'density_relative_integration_error':float(abs(rho.sum()*dv/cell.nelectron-1)),
             'field_mesh':mesh.tolist(),'orbital_kpoint_index':0,
             'orbital_band_indices':chosen,'transport_computed':do_transport}
    if do_transport:
        direct=sum(np.trace(conductivity_tensor(p,w))/3 for p,w in zip(momenta,weights))
        integrated=spc_complex.sum()*dv
        quality.update(sigma_iso_at_spc_frequency_S_m=float(direct),
            spc_integral_real_S_m=float(integrated.real),spc_integral_imag_S_m=float(integrated.imag),
            spc_sum_relative_error=float(abs(integrated-direct)/max(abs(direct),1e-12)),
            momentum_grid_relative_error=float(max(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-12) for a,b in zip(grid_p,momenta))),
            active_orbital_overlap_grid_max_error=float(max(np.max(abs(s-np.eye(len(s)))) for s in grid_overlap)),
            transport_energy_window_eV=transport_cfg['energy_window_eV'],
            momentum_operator='canonical all-electron momentum, local PBE',
            dc_interpretation='finite-broadening independent-particle estimate; requires eta/k/basis/size/temperature convergence')
        write_cube(out/'spc_absolute.cube',cell,abs(spc_complex).reshape(mesh),'Positive SPC = absolute complex row-sum density; (S/m)/bohr^3')
        write_cube(out/'spc_signed.cube',cell,spc_complex.real.reshape(mesh),'Signed SPC row-sum density; integral should recover isotropic conductivity')
        np.savez_compressed(out/'spc.npz',row_density=spc_complex.reshape(mesh),dv_bohr3=dv)
    dump(out/'quality.json',quality)
    return quality

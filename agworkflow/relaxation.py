"""Dispersion-inclusive DFT geometry/cell relaxation.

Finite differences of converged periodic PBE+D3 free energies provide forces and
stress. This is deliberately transparent and portable, but expensive. It avoids
pretending GAP forces are derivatives of a PBE+D3 electronic energy.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
from ase.calculators.calculator import Calculator,all_changes
from ase.calculators.fd import FiniteDifferenceCalculator
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from ase.io import write
from reproduce import dump,AMU_A3_TO_G_CM3
from .electronic import solve,HA_EV
from .structure import identify_layers


def dispersion(atoms,cfg):
    from dftd3.ase import DFTD3
    a=atoms.copy(); a.calc=DFTD3(method=cfg['method'],damping=cfg['damping'])
    return float(a.get_potential_energy())


class PBEPlusD3(Calculator):
    implemented_properties=['energy','free_energy']
    def __init__(self,electronic,d3,**kwargs):
        super().__init__(**kwargs)
        self.electronic=electronic; self.d3=d3; self.cache={}; self.evaluations=0

    def calculate(self,atoms=None,properties=('energy',),system_changes=all_changes):
        super().calculate(atoms,properties,system_changes)
        key=hashlib.sha256(atoms.numbers.tobytes()+atoms.positions.tobytes()+atoms.cell.array.tobytes()).hexdigest()
        if key not in self.cache:
            cell,mf,kpts,elapsed=solve(atoms,self.electronic,'PBE')
            correction=dispersion(atoms,self.d3)
            self.cache[key]={'energy':float(mf.e_tot*HA_EV+correction),
                             'free_energy':float(mf.e_free*HA_EV+correction)}
            self.evaluations+=1
        self.results=dict(self.cache[key])


def zero_pressure_relax(atoms,electronic,d3,settings,out):
    out=Path(out); out.mkdir(parents=True,exist_ok=False)
    a=atoms.copy()
    base=PBEPlusD3(electronic,d3)
    a.calc=FiniteDifferenceCalculator(base,eps_disp=settings['displacement_A'],
                     eps_strain=settings['strain_step'],force_consistent=True)
    filt=FrechetCellFilter(a,scalar_pressure=0)
    opt=BFGS(filt,logfile=str(out/'optimizer.log'),trajectory=str(out/'optimizer.traj'),maxstep=.05)
    converged=opt.run(fmax=settings['fmax_eV_A'],steps=settings['steps'])
    write(out/'final.extxyz',a)
    summary={'converged':bool(converged),'method':'PBE+D3; numerical free-energy derivatives; full cell and atoms',
        'damping':d3['damping'],'density_g_cm3':float(sum(a.get_masses())*AMU_A3_TO_G_CM3/a.get_volume()),
        'volume_A3':float(a.get_volume()),'n_scf_evaluations':base.evaluations,
        'max_atomic_force_eV_A':float(np.linalg.norm(a.get_forces(),axis=1).max()),
        'stress_eV_A3':a.get_stress().tolist(),'layers':identify_layers(a)}
    dump(out/'summary.json',summary)
    return summary

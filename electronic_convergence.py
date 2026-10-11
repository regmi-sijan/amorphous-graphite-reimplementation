"""Converge electronic settings on one unchanged structure, using real DFT.

Plan format: a JSON list of {"name":"...", "electronic":{...}, "transport":{...}}.
This avoids rerunning MD when changing basis, k mesh or integration grids.
"""
import argparse
import copy
import json
from pathlib import Path
import numpy as np
from ase.io import read
from reproduce import dump
from agworkflow.electronic import run_scf,analyze_electrons
from agworkflow.structure import identify_layers
from agworkflow.plotting import electronic_plots


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--structure',required=True)
    p.add_argument('--config',required=True)
    p.add_argument('--plan',required=True)
    p.add_argument('--out',required=True)
    a=p.parse_args();out=Path(a.out)
    out.mkdir(parents=True,exist_ok=False)
    atoms=read(a.structure);cfg=json.loads(Path(a.config).read_text())
    plan=json.loads(Path(a.plan).read_text());results=[];names=set()
    for step in plan:
        name=step['name']
        if name in names or Path(name).name!=name or name in ('.','..'):
            raise ValueError('Each plan name must be a unique plain directory name')
        names.add(name)
        e=copy.deepcopy(cfg['electronic']);e.update(step.get('electronic',{}))
        t=copy.deepcopy(cfg['transport']);t.update(step.get('transport',{}))
        folder=out/name
        print(f'Calculating {name}',flush=True)
        summary=run_scf(atoms,e,folder,'PBE');quality=analyze_electrons(folder,t)
        electronic_plots(folder,identify_layers(atoms))
        result={'name':name,'electronic':e,'transport':t,'scf':summary,'quality':quality}
        if results:
            result['energy_change_eV_atom']=summary['energy_eV_atom']-results[-1]['scf']['energy_eV_atom']
        results.append(result);dump(out/'convergence.json',results)
    print(f'Finished {len(results)} actual calculations. Inspect convergence.json; no automatic physical convergence claim.')

if __name__=='__main__': main()

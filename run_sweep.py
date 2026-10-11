"""Generate and optionally execute sequential parameter studies.

Nested configuration paths use dots. Default is preparation only; --run actually
launches every job, which can take days for production settings. Each output is
isolated; --resume reuses finished stages via end_to_end.py.
"""
import argparse
import copy
import itertools
import json
from pathlib import Path
import subprocess
import sys
from reproduce import ROOT,dump


def assign(cfg,path,value):
    parts=path.split('.');target=cfg
    for part in parts[:-1]: target=target[part]
    if parts[-1] not in target: raise KeyError(path)
    target[parts[-1]]=value


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',required=True)
    p.add_argument('--grid',required=True,help='JSON mapping dotted config paths to value lists')
    p.add_argument('--out',required=True)
    p.add_argument('--run',action='store_true')
    p.add_argument('--resume',action='store_true')
    args=p.parse_args()
    base=json.loads(Path(args.config).read_text());grid=json.loads(Path(args.grid).read_text())
    if not grid or any(not isinstance(v,list) or not v for v in grid.values()):
        raise ValueError('Each sweep axis must be a nonempty JSON list')
    out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=True)
    jobs=[]
    for index,values in enumerate(itertools.product(*grid.values())):
        config=copy.deepcopy(base)
        for key,value in zip(grid,values): assign(config,key,value)
        cfgpath=out/f'job_{index:04d}.json';result=out/f'job_{index:04d}'
        if cfgpath.exists() and json.loads(cfgpath.read_text())!=config:
            raise ValueError('Existing sweep differs; choose a fresh output directory')
        dump(cfgpath,config)
        command=[sys.executable,str(ROOT/'end_to_end.py'),'--config',str(cfgpath),'--out',str(result)]
        if args.resume: command.append('--resume')
        jobs.append({'parameters':dict(zip(grid,values)),'command':command,'output':str(result)})
    dump(out/'jobs.json',jobs)
    print(f'Prepared {len(jobs)} jobs in {out}',flush=True)
    if args.run:
        for job in jobs: subprocess.run(job['command'],check=True,cwd=ROOT)
        summaries=[]
        for job in jobs:
            summaries.append({'parameters':job['parameters'],
                'comparison':json.loads((Path(job['output'])/'comparison.json').read_text())})
        dump(out/'sweep_results.json',summaries)
    else: print('No simulations launched. Add --run to execute.')

if __name__=='__main__': main()

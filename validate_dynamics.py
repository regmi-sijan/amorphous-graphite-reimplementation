"""Short, real-GAP numerical stability and thermal-settling checks.

This is a screening study on one 80-atom random start, not a proof of equilibrium
or of the paper's structural transition. All durations are simulated times.
"""
from __future__ import annotations

from reproduce import ROOT, load_config, make_calculator, dump, geometry, layering
import argparse
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import time

import numpy as np
from ase import units
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.md.verlet import VelocityVerlet
from ase.md.nose_hoover_chain import NoseHooverChainNVT
from ase.md.velocitydistribution import thermalize_momenta


def run_segment(initial, config, out, ensemble, dt_fs, duration_fs):
    out.mkdir(parents=True, exist_ok=False)
    atoms = initial.copy()
    atoms.calc, _ = make_calculator(config)
    if ensemble == "NVE":
        dyn = VelocityVerlet(atoms, timestep=dt_fs * units.fs)
    else:
        dyn = NoseHooverChainNVT(atoms, timestep=dt_fs * units.fs,
            temperature_K=config["temperature_K"],
            tdamp=config["thermostat_damping_fs"] * units.fs,
            tchain=config["thermostat_chain_length"])
    steps = round(duration_fs / dt_fs)
    if not np.isclose(steps * dt_fs, duration_fs):
        raise ValueError("Duration must be divisible by timestep")
    interval = max(1, round(1 / dt_fs))  # output every approximately 1 fs
    begun = time.monotonic()
    status = {"ensemble": ensemble, "timestep_fs": dt_fs, "duration_fs": duration_fs,
              "thermostat_damping_fs": config["thermostat_damping_fs"] if ensemble == "NVT" else None,
              "status": "running", "completed_steps": 0}
    dump(out / "status.json", status)
    records = []
    with (out / "thermo.csv").open("w") as stream, Trajectory(str(out / "trajectory.traj"), "w", atoms) as trajectory:
        fields = ["time_fs", "temperature_K", "potential_eV_atom", "physical_total_eV_atom",
                  "conserved_eV_atom", "max_force_eV_A", "wall_seconds"]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        def record():
            t = dyn.nsteps * dt_fs
            if records and t == records[-1]["time_fs"]:
                return
            ep = atoms.get_potential_energy() / len(atoms)
            total = ep + atoms.get_kinetic_energy() / len(atoms)
            conserved = total if ensemble == "NVE" else dyn.get_conserved_energy() / len(atoms)
            row = dict(zip(fields, [t, atoms.get_temperature(), ep, total, conserved,
                float(np.linalg.norm(atoms.get_forces(), axis=1).max()), time.monotonic()-begun]))
            if not np.isfinite(list(row.values())).all() or row["temperature_K"] > 50000:
                raise RuntimeError("Unstable trajectory; inspect saved records")
            atoms.info.update(time_ps=t/1000, ensemble=ensemble)
            trajectory.write(atoms)
            writer.writerow(row); stream.flush()
            records.append(row)
            status.update(completed_steps=dyn.nsteps, last_temperature_K=row["temperature_K"],
                          wall_seconds=row["wall_seconds"])
            dump(out / "status.json", status)
        dyn.attach(record, interval=interval)
        try:
            dyn.run(steps)
            record()
        except BaseException:
            status["status"] = "failed_or_interrupted"
            dump(out / "status.json", status)
            raise
    status.update(status="completed", wall_seconds=time.monotonic()-begun)
    dump(out / "status.json", status)
    write(out / "final.extxyz", atoms)
    energies = np.array([r["conserved_eV_atom"] for r in records])
    times = np.array([r["time_fs"] for r in records])
    summary = {**status,
        "max_abs_conserved_error_meV_atom": float(np.max(abs(energies-energies[0])) * 1000),
        "endpoint_conserved_error_meV_atom": float((energies[-1]-energies[0]) * 1000),
        "linear_conserved_drift_meV_atom_ps": float(np.polyfit(times/1000, energies-energies[0], 1)[0]*1000),
        "final_temperature_K": records[-1]["temperature_K"],
        "final_geometry": geometry(atoms), "final_layering_proxy": layering(atoms)}
    if ensemble == "NVT":
        blocks = []
        # Non-overlapping 25 fs windows; means are correlated, not independent CIs.
        for start in np.arange(max(0, duration_fs-100), duration_fs, 25):
            block = [r for r in records if start < r["time_fs"] <= start+25]
            if block:
                blocks.append({"start_fs": float(start), "end_fs": float(start+25),
                    "temperature_mean_K": float(np.mean([r["temperature_K"] for r in block])),
                    "potential_mean_eV_atom": float(np.mean([r["potential_eV_atom"] for r in block]))})
        summary["last_100fs_blocks"] = blocks
        summary["tail_mean_temperature_K"] = float(np.mean([b["temperature_mean_K"] for b in blocks]))
        summary["all_tail_blocks_within_10pct_target"] = all(
            abs(b["temperature_mean_K"]/config["temperature_K"]-1) < 0.10 for b in blocks)
        summary["equilibrium_claim"] = False
    dump(out / "summary.json", summary)
    print(json.dumps(summary), flush=True)
    return summary


def study(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    c = load_config(ROOT / "configs/mac_pilot.json")
    atoms = read(ROOT / "results/mac_pilot_verified/initial.extxyz")
    thermalize_momenta(atoms, temperature_K=c["temperature_K"], rng=np.random.default_rng(c["seed"]+1))
    write(out / "common_initial.extxyz", atoms)
    _, xml = make_calculator(c)
    dump(out / "study.json", {
        "config": c, "nve_duration_fs": 20, "nvt_duration_fs": 300,
        "nve_timesteps_fs": [0.5, 0.25, 0.125],
        "acceptance_screen_max_energy_error_meV_atom": 1.0,
        "scope": "One 80-atom random start; screening criterion chosen here, not from paper",
        "same_positions_and_momenta_for_all_segments": True,
        "versions": {k: importlib.metadata.version(k) for k in ["ase", "quippy-ase", "numpy"]},
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "model_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in xml.parent.glob("carbon.xml*")},
    })
    nve = [run_segment(atoms, c, out / f"nve_dt_{dt}", "NVE", dt, 20) for dt in [0.5,0.25,0.125]]
    acceptable = [x["timestep_fs"] for x in nve if x["max_abs_conserved_error_meV_atom"] <= 1.0]
    if not acceptable:
        dump(out / "decision.json", {"proceed_to_nvt": False, "reason": "No timestep passed 1 meV/atom short NVE screen"})
        return
    dt = max(acceptable)
    dump(out / "decision.json", {"proceed_to_nvt": True, "chosen_timestep_fs": dt,
        "reason": "Largest tested dt passing <=1 meV/atom error over 20 fs; provisional choice only"})
    run_segment(atoms, c, out / "nvt", "NVT", dt, 300)


def report(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = Path(out)
    fig, axes = plt.subplots(2, 2, figsize=(11,7), constrained_layout=True)
    for folder in sorted(out.glob("nve_dt_*")):
        a = np.genfromtxt(folder / "thermo.csv", delimiter=",", names=True)
        axes[0,0].plot(a["time_fs"], (a["conserved_eV_atom"]-a["conserved_eV_atom"][0])*1000, label=folder.name.replace("nve_dt_", "dt = ")+" fs")
    axes[0,0].set(xlabel="Time (fs)", ylabel="NVE energy error (meV/atom)")
    axes[0,0].legend()
    nvt = out / "nvt/thermo.csv"
    if nvt.exists():
        a = np.genfromtxt(nvt, delimiter=",", names=True)
        axes[0,1].plot(a["time_fs"], a["temperature_K"])
        axes[0,1].axhline(2700, ls="--", color="gray", label="Target: 2700 K")
        axes[0,1].legend()
        axes[0,1].set(xlabel="Time (fs)", ylabel="NVT temperature (K)")
        axes[1,0].plot(a["time_fs"], a["potential_eV_atom"])
        axes[1,0].set(xlabel="Time (fs)", ylabel="NVT potential energy (eV/atom)")
        axes[1,1].plot(a["time_fs"], (a["conserved_eV_atom"]-a["conserved_eV_atom"][0])*1000)
        axes[1,1].set(xlabel="Time (fs)", ylabel="NVT extended-energy error (meV/atom)")
    fig.suptitle("80-atom numerical screening; not a layer-formation result")
    fig.savefig(out / "validation.png", dpi=160)
    plt.close(fig)
    nve_summaries = [json.loads(p.read_text()) for p in sorted(out.glob("nve_dt_*/summary.json"))]
    nvt_summary_path = out / "nvt/summary.json"
    if len(nve_summaries) != 3 or not nvt_summary_path.exists():
        raise RuntimeError("Plot saved, but the completed study is required for the final report")
    s = json.loads(nvt_summary_path.read_text())
    if any(x["status"] != "completed" for x in [*nve_summaries, s]):
        raise RuntimeError("Cannot issue a completed report for an unfinished run")
    nve_table = "\n".join(
        f'| {x["timestep_fs"]} | {x["max_abs_conserved_error_meV_atom"]:.3f} | '
        f'{"Pass" if x["max_abs_conserved_error_meV_atom"] <= 1 else "Fail"} |'
        for x in sorted(nve_summaries, key=lambda x: -x["timestep_fs"]))
    temperature_table = "\n".join(
        f'| {b["start_fs"]:.0f}-{b["end_fs"]:.0f} | {b["temperature_mean_K"]:.1f} |'
        for b in s["last_100fs_blocks"])
    elapsed = sum(x["wall_seconds"] for x in [*nve_summaries, s])
    (out / "REPORT.md").write_text(f'''# Completed dynamics validation

The 80-atom model passed preliminary numerical-stability and temperature-settling
checks on this Mac. It has **not** demonstrated structural equilibrium or amorphous
graphite formation. The threefold-coordination fraction is
{s["final_geometry"]["threefold_fraction"]:.0%}, which is only a geometric proxy for sp2.

## Timestep test

Three unthermostatted (NVE) trajectories ran for 20 fs each, starting from exactly
the same positions and momenta. The acceptance screen was a maximum sampled total
energy error of 1 meV/atom, chosen for this project rather than taken from the paper.

| Timestep (fs) | Maximum energy error (meV/atom) | Screen |
|---|---|---|
{nve_table}

Reducing the timestep consistently reduces the error. Retain **0.25 fs** for now:
it is the largest tested timestep that passed this short screen. This does not
replace longer energy-drift and observable-convergence studies. A fitted slope
over 20 fs should not be extrapolated to a long production trajectory.

## Thermostat test

The NVT trajectory ran for 300 fs (0.3 ps), with a 2700 K target, 25 fs damping,
a three-variable Nose-Hoover chain, and the selected {s["timestep_fs"]} fs timestep.
The common random starting structure was used without pre-relaxation or imposed
layering. Initial heating was followed by cooling toward the thermostat target.

| Final time window (fs) | Mean temperature (K) |
|---|---|
{temperature_table}

The last 100 fs averaged **{s["tail_mean_temperature_K"]:.1f} K**. All four
25 fs block averages were within 10% of 2700 K: **{s["all_tail_blocks_within_10pct_target"]}**.
These blocks are correlated diagnostics, not independent samples or confidence
intervals. The final instantaneous temperature was {s["final_temperature_K"]:.1f} K;
instantaneous fluctuations around the target are expected.

The maximum error in the thermostat's conserved extended energy was
**{s["max_abs_conserved_error_meV_atom"]:.3f} meV/atom**. Physical total energy alone
is not conserved in NVT because the thermostat exchanges energy with the atoms.

## What is complete and what remains

The three NVE runs and the NVT run all completed. Combined measured simulation
time was **{elapsed/60:.1f} minutes** on this Mac. No simulation from this study
remains running. Raw records, trajectories, final structures and JSON summaries
are retained in the corresponding subdirectories. `validation.png` contains the
four-panel diagnostic plot. `study.json` records settings, versions and model
hashes; `executed_driver.py` preserves the source used for the runs.

This supports the chosen numerical settings and preliminary thermal settling for
one 80-atom configuration. It does not establish long-time equilibrium, independent
seed reproducibility, finite-size convergence, or the paper's >95% sp2-plus-layering
criterion. Charge density and conductivity were not calculated. Exact identity
of the recovered 2017 GAP file with the one used in the 2022 PRL remains unconfirmed.

The Mac is sufficient for these checks and further short runs. The earlier
400-atom benchmark estimated about 10 days for one 100 ps trajectory at the
measured settings. Larger production runs and multiple seeds would benefit from
additional computing resources, after benchmarking the intended backend.

The next scientific stage is longer structural relaxation/annealing, followed by
layer identification and comparisons across seeds and sizes. Preserve the present
results as a preliminary check rather than treating them as a completed paper
reproduction. The current trajectory files omit full thermostat restart state;
starting from the endpoint would require documented thermostat re-equilibration.
''')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["run", "report"])
    parser.add_argument("--out", default=str(ROOT / "results/dynamics_validation"))
    args = parser.parse_args()
    if args.action == "run": study(args.out)
    else: report(args.out)

"""Structural reproduction workflow for PRL 128, 236402 (2022).

Real forces come from the published carbon GAP through QUIP. No force-field
fallback or imposed layering is used. This does not calculate electrons/SPC.
Run `python reproduce.py --help` for the individual workflow stages.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import itertools
import json
import os
from pathlib import Path
import platform
import time

# Set before importing numerical libraries. Users can override these explicitly.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).parent / ".mplconfig"))

import numpy as np
from ase import Atoms, units
from ase.io import read, write, iread
from ase.neighborlist import neighbor_list

ROOT = Path(__file__).resolve().parent
AMU_A3_TO_G_CM3 = 1.66053906660


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def load_config(path):
    c = json.loads(Path(path).read_text())
    for k in ("n_atoms", "density_g_cm3", "temperature_K", "duration_ps",
              "timestep_fs", "thermostat_damping_fs", "save_interval_steps",
              "minimum_initial_distance_A", "bond_cutoff_A", "thermostat_chain_length"):
        if not np.isfinite(c[k]) or c[k] <= 0:
            raise ValueError(f"{k} must be finite and positive")
    for k in ("n_atoms", "save_interval_steps", "thermostat_chain_length"):
        if int(c[k]) != c[k]:
            raise ValueError(f"{k} must be an integer")
    if c["n_atoms"] < 8:
        raise ValueError("Use at least 8 atoms")
    steps = c["duration_ps"] * 1000 / c["timestep_fs"]
    if steps < 1 or not np.isclose(steps, round(steps), rtol=0, atol=1e-7):
        raise ValueError("duration_ps must be an integer number of timesteps")
    return c


def random_carbon(n, density, seed, min_distance=1.2):
    """Independent random start with periodic hard-core rejection; no sheets."""
    length = (n * 12.011 * AMU_A3_TO_G_CM3 / density) ** (1 / 3)
    rng = np.random.default_rng(seed)
    positions = []
    for _ in range(n):
        for attempt in range(100000):
            candidate = rng.uniform(0, length, 3)
            if positions:
                delta = candidate - np.array(positions)
                delta -= length * np.rint(delta / length)
                if np.min(np.linalg.norm(delta, axis=1)) < min_distance:
                    continue
            positions.append(candidate)
            break
        else:
            raise RuntimeError("Random placement failed; reduce the minimum distance")
    return Atoms("C" * n, positions=positions, cell=[length] * 3, pbc=True)


def make_calculator(c):
    xml = Path(c["gap_xml"])
    if not xml.is_absolute():
        xml = ROOT / xml
    if not xml.is_file():
        raise FileNotFoundError(f"Missing GAP model: {xml}. See README.md; no substitute is used.")
    from quippy.potential import Potential
    calc = Potential(c["gap_args"], param_filename=str(xml.resolve()))
    # QUIP 0.10.3's name getter expects a setter that newer ASE no longer calls.
    # Set the public property so ASE trajectory serialization works.
    calc.name = "quippy"
    return calc, xml


def geometry(atoms, cutoff=1.85):
    """Coordination and a local geometric proxy for threefold planar carbon.

    Coordination alone is NOT an electronic hybridization measurement. Local
    planarity uses the smallest/total eigenvalue of neighbor-vector covariance.
    """
    i, j, vectors = neighbor_list("ijD", atoms, cutoff)
    coord = np.bincount(i, minlength=len(atoms))
    normals, planar_ids = [], []
    for index in np.flatnonzero(coord == 3):
        v = vectors[i == index]
        eigenvalues, eigenvectors = np.linalg.eigh(v.T @ v)
        if eigenvalues[0] / max(eigenvalues.sum(), 1e-12) < 0.02:
            normals.append(eigenvectors[:, 0])
            planar_ids.append(int(index))
    alignment = 0.0
    if normals:
        nn = np.array(normals)
        q = 1.5 * (nn.T @ nn / len(nn)) - 0.5 * np.eye(3)
        alignment = float(np.linalg.eigvalsh(q)[-1])
    return {
        "twofold_fraction": float(np.mean(coord == 2)),
        "threefold_fraction": float(np.mean(coord == 3)),
        "fourfold_fraction": float(np.mean(coord == 4)),
        "planar_threefold_fraction_proxy": len(normals) / len(atoms),
        "normal_alignment_proxy": alignment,
    }


def layering(atoms, min_spacing=2.5, max_spacing=4.0):
    """Search commensurate reciprocal-cell density modes for possible layering.

    |mean exp(2*pi*i*h.s)| ranges 0..1. It is a diagnostic, not the paper's
    planarity fraction. Maximizing over modes biases the random baseline upward.
    A peak alone cannot establish graphite; inspect bonding, topology, snapshots.
    """
    reciprocal = np.asarray(atoms.cell.reciprocal())
    bound = int(np.ceil(max(atoms.cell.lengths()) / min_spacing))
    hs = []
    for h in itertools.product(range(-bound, bound + 1), repeat=3):
        first = next((x for x in h if x), 0)
        if first <= 0:
            continue
        q = np.array(h) @ reciprocal
        d = 1 / np.linalg.norm(q)
        if min_spacing <= d <= max_spacing:
            hs.append(h)
    if not hs:
        return {"layering_amplitude_proxy": 0.0, "candidate_spacing_A": None,
                "candidate_hkl": None}
    hs = np.array(hs)
    amplitudes = abs(np.exp(2j * np.pi * (atoms.get_scaled_positions() @ hs.T)).mean(axis=0))
    best = int(np.argmax(amplitudes))
    return {"layering_amplitude_proxy": float(amplitudes[best]),
            "candidate_spacing_A": float(1 / np.linalg.norm(hs[best] @ reciprocal)),
            "candidate_hkl": hs[best].tolist()}


def rdf_3d(atoms, bins=120):
    """Isotropic 3D g(r), not the in-plane RDF shown in Figure 1."""
    heights = 1 / np.linalg.norm(atoms.cell.reciprocal(), axis=1)
    rmax = min(5.0, 0.49 * min(heights))
    distances = neighbor_list("d", atoms, rmax)
    edges = np.linspace(0, rmax, bins + 1)
    counts, _ = np.histogram(distances, edges)
    shell_volume = 4 * np.pi / 3 * np.diff(edges**3)
    density = len(atoms) / atoms.get_volume()
    g = counts / (len(atoms) * density * shell_volume)
    return (edges[:-1] + edges[1:]) / 2, g


def ring_counts(atoms, cutoff=1.85, max_size=8):
    """Count unique chordless, non-winding graph cycles of sizes 3..max_size.

    This declared ring convention is not asserted to equal the authors' method.
    Periodic winding loops are excluded. Require a cell larger than 2*cutoff so
    each bonded pair has an unambiguous periodic image.
    """
    heights = 1 / np.linalg.norm(atoms.cell.reciprocal(), axis=1)
    if min(heights) <= 2 * cutoff:
        raise ValueError("Cell too small for this ring implementation")
    i, j, shifts = neighbor_list("ijS", atoms, cutoff)
    adj = [set() for _ in atoms]
    shift = {}
    for a, b, s in zip(i, j, shifts):
        adj[int(a)].add(int(b))
        shift[int(a), int(b)] = s
    cycles = set()
    for start in range(len(atoms)):
        def visit(path):
            last = path[-1]
            for nxt in sorted(adj[last]):
                if nxt == start:
                    if len(path) >= 3:
                        winding = sum((shift[a, b] for a, b in zip(path, path[1:] + path[:1])), np.zeros(3, int))
                        if np.all(winding == 0):
                            cycles.add(min(tuple(path), (start, *reversed(path[1:]))))
                    continue
                if nxt <= start or nxt in path or len(path) == max_size:
                    continue
                # A connection to start closes the cycle on the next visit;
                # connections to any other earlier vertex introduce chords.
                if any(nxt in adj[p] for p in path[1:-1]):
                    continue
                if len(path) > 2 and last in adj[start]:
                    continue
                visit(path + [nxt])
        visit([start])
    return {str(n): sum(len(c) == n for c in cycles) for n in range(3, max_size + 1)}


def prepare(c, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    atoms = random_carbon(c["n_atoms"], c["density_g_cm3"], c["seed"], c["minimum_initial_distance_A"])
    atoms.info.update(time_ps=0.0, origin="independent random configuration")
    write(out / "initial.extxyz", atoms)
    dump(out / "config.json", c)
    print(f"Prepared {len(atoms)} atoms, box {atoms.cell.lengths()[0]:.4f} A in {out}")


def simulate(out, steps_override=None):
    from ase.io.trajectory import Trajectory
    from ase.md.nose_hoover_chain import NoseHooverChainNVT
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    out = Path(out)
    c = load_config(out / "config.json")
    if (out / "trajectory.traj").exists():
        raise FileExistsError("Trajectory already exists. Prepare a new run directory.")
    atoms = read(out / "initial.extxyz")
    atoms.calc, xml = make_calculator(c)
    steps = round(c["duration_ps"] * 1000 / c["timestep_fs"])
    if steps_override is not None:
        if steps_override < 1:
            raise ValueError("steps must be positive")
        steps = steps_override
    # No COM constraint: the thermostat uses 3N degrees of freedom consistently.
    MaxwellBoltzmannDistribution(atoms, temperature_K=c["temperature_K"], rng=np.random.default_rng(c["seed"] + 1))
    md = NoseHooverChainNVT(atoms, timestep=c["timestep_fs"] * units.fs,
                          temperature_K=c["temperature_K"],
                          tdamp=c["thermostat_damping_fs"] * units.fs,
                          tchain=c["thermostat_chain_length"])
    provenance = {
        "platform": platform.platform(), "python": platform.python_version(),
        "versions": {k: importlib.metadata.version(k) for k in ("ase", "quippy-ase", "numpy", "scipy", "matplotlib")},
        "gap_xml": str(xml), "potential_files_sha256": {
            str(p.relative_to(xml.parent)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(xml.parent.rglob("*")) if p.is_file()},
        "planned_steps": steps, "actual_duration_ps": steps*c["timestep_fs"]/1000,
        "notes": "Independent GAP simulation. No claim of paper reproduction from a short pilot.",
        "status": "running",
    }
    dump(out / "provenance.json", provenance)
    started = time.monotonic()
    with Trajectory(str(out / "trajectory.traj"), "w", atoms) as trajectory, (out / "thermo.csv").open("w") as f:
        writer = csv.writer(f)
        writer.writerow(["step", "time_ps", "temperature_K", "potential_eV_per_atom", "total_eV_per_atom", "wall_seconds"])
        last_saved = [-1]
        def record():
            step = md.nsteps
            if step == last_saved[0]:
                return
            energy = atoms.get_potential_energy()
            temperature = atoms.get_temperature()
            if not np.isfinite(energy) or not np.isfinite(temperature) or temperature > 50000:
                raise RuntimeError("MD became unstable. Inspect the start, timestep, and forces.")
            atoms.info["time_ps"] = step * c["timestep_fs"] / 1000
            trajectory.write(atoms)
            writer.writerow([step, atoms.info["time_ps"], temperature, energy/len(atoms),
                             (energy+atoms.get_kinetic_energy())/len(atoms), time.monotonic()-started])
            f.flush()
            last_saved[0] = step
        md.attach(record, interval=c["save_interval_steps"])
        try:
            md.run(steps)
            record()
        except BaseException:
            provenance.update(status="interrupted_or_failed", completed_steps=md.nsteps)
            dump(out / "provenance.json", provenance)
            raise
    write(out / "final.extxyz", atoms)
    elapsed = time.monotonic() - started
    provenance.update(status="completed", completed_steps=steps, wall_seconds=elapsed,
                      estimated_hours_per_100ps=elapsed/steps*(100000/c["timestep_fs"])/3600)
    dump(out / "provenance.json", provenance)
    print(json.dumps(provenance, indent=2))


def analyze(out, stride=1):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = Path(out)
    if stride < 1:
        raise ValueError("stride must be positive")
    c = load_config(out / "config.json")
    rows = []
    last = None
    for frame in iread(out / "trajectory.traj", index=f"::{stride}"):
        metrics = geometry(frame, c["bond_cutoff_A"])
        layer = layering(frame)
        rows.append({"time_ps": frame.info["time_ps"], **metrics,
                     **{k:v for k,v in layer.items() if k != "candidate_hkl"}})
        last = frame
    if last is None:
        raise ValueError("No frames to analyze")
    with (out / "structure.csv").open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    # Analyze the actual final saved frame even when stride omits it.
    last = read(out / "trajectory.traj", index=-1)
    r, g = rdf_3d(last)
    np.savetxt(out / "rdf_3d.csv", np.column_stack([r, g]), delimiter=",", header="r_A,g_3d", comments="")
    summary = {"final_time_ps": last.info["time_ps"], **geometry(last, c["bond_cutoff_A"]),
               **layering(last), "ring_counts_chordless_nonwinding": ring_counts(last, c["bond_cutoff_A"]),
               "interpretation": "Geometric diagnostics only; threefold != measured sp2. Candidate spacing is not a confirmed layer spacing. No electrons or conductivity computed."}
    dump(out / "analysis.json", summary)
    thermo = np.atleast_1d(np.genfromtxt(out / "thermo.csv", delimiter=",", names=True))
    fig, ax = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    ax[0,0].plot(thermo["time_ps"], thermo["potential_eV_per_atom"])
    ax[0,0].set(xlabel="Time (ps)", ylabel="GAP potential energy (eV/atom)")
    ax[0,1].plot(thermo["time_ps"], thermo["temperature_K"])
    ax[0,1].axhline(c["temperature_K"], ls="--", color="gray")
    ax[0,1].set(xlabel="Time (ps)", ylabel="Temperature (K)")
    for key, label in [("threefold_fraction", "Threefold coordination"),
                       ("planar_threefold_fraction_proxy", "Planar threefold (proxy)"),
                       ("layering_amplitude_proxy", "Density modulation (proxy)")]:
        ax[1,0].plot([x["time_ps"] for x in rows], [x[key] for x in rows], label=label)
    ax[1,0].set(xlabel="Time (ps)", ylabel="Fraction / amplitude", ylim=(0,1))
    ax[1,0].legend(fontsize=8)
    ax[1,1].plot(r, g)
    ax[1,1].set(xlabel="Distance (A)", ylabel="3D g(r); not in-plane RDF")
    fig.suptitle(f"Independent GAP run: {len(last)} C atoms, {last.info['time_ps']:.4g} ps")
    fig.savefig(out / "diagnostics.png", dpi=160)
    plt.close(fig)
    # Calculator-free trajectory export for OVITO or other viewers.
    def frames():
        for frame in iread(out / "trajectory.traj"):
            frame.calc = None
            yield frame
    write(out / "trajectory.extxyz", frames(), format="extxyz")
    print(json.dumps(summary, indent=2))


def relax(out, steps=500):
    from ase.optimize import FIRE
    out = Path(out)
    if (out / "relaxed.extxyz").exists():
        raise FileExistsError("A relaxed structure already exists")
    c = load_config(out / "config.json")
    atoms = read(out / "final.extxyz")
    atoms.set_momenta(np.zeros((len(atoms),3)))
    atoms.calc, _ = make_calculator(c)
    opt = FIRE(atoms, logfile=str(out / "relax.log"), trajectory=str(out / "relax.traj"))
    converged = opt.run(fmax=0.03, steps=steps)
    write(out / "relaxed.extxyz", atoms)
    dump(out / "relaxation.json", {"optimizer": "FIRE (differs from paper's CG)",
          "cell": "fixed; does not reproduce zero-pressure vdW relaxation", "converged": bool(converged),
          "max_force_eV_A": float(np.max(np.linalg.norm(atoms.get_forces(), axis=1))),
          **geometry(atoms, c["bond_cutoff_A"]), **layering(atoms)})
    print(f"Relaxation converged: {converged}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    a = sub.add_parser("prepare", help="Generate an unbiased random carbon structure")
    a.add_argument("--config", default=str(ROOT / "configs/mac_pilot.json"))
    a.add_argument("--out", required=True)
    for command in ("run", "analyze", "relax"):
        a = sub.add_parser(command)
        a.add_argument("--out", required=True)
        if command in ("run", "relax"):
            a.add_argument("--steps", type=int, default=None if command == "run" else 500)
        else:
            a.add_argument("--stride", type=int, default=1)
    args = p.parse_args()
    if args.command == "prepare": prepare(load_config(args.config), args.out)
    elif args.command == "run": simulate(args.out, args.steps)
    elif args.command == "analyze": analyze(args.out, args.stride)
    elif args.command == "relax": relax(args.out, args.steps)


if __name__ == "__main__":
    main()

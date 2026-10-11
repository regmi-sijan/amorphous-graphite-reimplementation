# Structural runner: detailed usage and historical context

For the main entry point and current scope, start with [README.md](../README.md).

Target: Thapa et al., *Ab Initio Simulation of Amorphous Graphite*,
Physical Review Letters **128**, 236402 (2022),
https://doi.org/10.1103/PhysRevLett.128.236402.

## What this project does

Python generates independent random carbon configurations at a specified density,
runs real constant-volume GAP molecular dynamics through ASE/QUIP, saves the
trajectory and thermodynamic data, and analyzes structural changes. No layer
formation is prescribed. A simulation can legitimately fail to layer.

The Python-led end-to-end workflow now also includes periodic PBE/HSE06, charge
and orbital exports, DOS, Kubo-Greenwood/SPC transport, graphite comparisons and
optional PBE+D3 cell relaxation. See [END_TO_END.md](../END_TO_END.md) for the complete
workflow, method definitions and scientific limitations. These are actual
calculations, not mock data. Execution does not establish paper reproduction.

```sh
.venv/bin/python end_to_end.py --out results/my_complete_test
```

The sections below document the original structural runner and model provenance.

Run commands from the repository root. Raw local results and the paper PDF are
not included in the source repository.

See `VALIDATION.md` for the completed tests, measured Mac performance, and the
limits of the short pilot. A complete local run is in `results/mac_pilot_verified`.

The subsequent numerical-stability and thermal-settling study is also complete:
see `results/dynamics_validation/REPORT.md`. The 0.25 fs timestep passed the short
energy-conservation screen, and the final 100 fs of the 300 fs thermostat run
averaged 2724 K at a 2700 K target. This is not a layer-formation result.

To regenerate its plots and report:

```sh
.venv/bin/python validate_dynamics.py report --out results/dynamics_validation
```

The validation driver uses the original pilot's saved starting geometry. To
repeat the numerical study, retain `results/mac_pilot_verified/initial.extxyz`
and run `validate_dynamics.py run --out results/dynamics_validation_repeat`.
Always use a fresh output directory; no existing trajectory is overwritten.

## Installation

Create an isolated Python 3.12 environment on a fresh machine:

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

`requirements-lock.txt` records the exact package versions used in the local tests.

QUIP/GAP supplies the force engine; the surrounding workflow is Python. It is not
a pure-Python reimplementation of DFT. QUIP publishes Apple Silicon wheels:
https://github.com/libAtoms/QUIP.

## Potential and provenance

Use the **2017 Deringer-Csanyi amorphous-carbon GAP**, the potential cited in
reference 28 of the PRL, with all companion sparse-data files. Do not substitute
GAP-20 or an arbitrary Tersoff/ReaxFF potential and call it an exact reproduction.

Sources:

- https://libatoms.github.io/GAP/data.html
- https://doi.org/10.1103/PhysRevB.95.094203
- https://doi.org/10.17863/CAM.7453

The historical standalone `aC_GAP.tar.gz` link returned 404 during setup.
The model was recovered from `training/full_das` in the Cambridge archive:
`carbon.xml` and its three `carbon.xml.sparseX.*` companions. Its label is
`GAP_2016_10_18_60_23_11_22_108`. This is the full two-body/three-body/SOAP fit in
the 2017 dataset; **exact file identity with the model used by the 2022 PRL
authors remains unconfirmed**. The archive includes several other fitted variants.
Its dataset
license is CC BY-NC-SA 4.0; see the repository's license for the original files.
QUIP and GAP have their own software licenses. Preserve the original provenance.

The model is not included in the repository. To obtain it, download the
1.19 GiB main archive from the repository and run the checksum-verifying extractor:

```sh
curl -L --fail 'https://api.repository.cam.ac.uk/server/api/core/bitstreams/8834f03c-6fc9-46af-aeff-1b2517bf8f2d/content' -o aC_GAP_data_main.tar.gz
python3 extract_model.py aC_GAP_data_main.tar.gz --out potentials
```

The extractor refuses to overwrite an existing nonempty potential directory.

Set `gap_xml` in the configuration to the recovered XML, and `gap_args` to the
model's initialization string. Relative XML paths resolve from this project
directory. Runs fail explicitly if the model is absent. A `provenance.json` file
records software versions and potential-file SHA256 hashes for each run.

## Run a small pilot

Run these commands from this directory. Use a **new output directory** each time;
the program refuses to overwrite a trajectory.

```sh
.venv/bin/python reproduce.py prepare --config configs/mac_pilot.json --out results/my_pilot
.venv/bin/python reproduce.py run --out results/my_pilot
.venv/bin/python reproduce.py analyze --out results/my_pilot
```

The pilot has 80 atoms at 2.44 g/cm3, target temperature 2700 K, and only 0.01 ps
of dynamics (40 steps). Its purpose is to verify installation and data flow.
The paper's layer formation takes tens to hundreds of ps. A short pilot cannot
demonstrate that result, establish temperature equilibration, or measure transport.

Optional fixed-cell relaxation of the final structure:

```sh
.venv/bin/python reproduce.py relax --out results/my_pilot --steps 500
```

This uses FIRE rather than the paper's conjugate-gradient optimizer. Check
`relaxation.json` for convergence; reaching a step limit is not convergence.

## Longer runs and a fair comparison

`configs/m5_like.json` supplies 400 atoms, 2700 K, 2.44 g/cm3 and 100 ps, resembling
the paper's GAP M5 case. It generates a new random configuration; it is not the
authors' M5 structure. At 0.25 fs this is **400,000 MD steps**. Benchmark a short
400-atom run before scheduling this on a laptop:

```sh
.venv/bin/python reproduce.py prepare --config configs/m5_like.json --out results/m5_benchmark
.venv/bin/python reproduce.py run --out results/m5_benchmark --steps 40
.venv/bin/python reproduce.py analyze --out results/m5_benchmark
```

For the full duration, prepare a fresh directory and run without `--steps`.
The provenance includes a rough wall-time estimate for 100 ps at the same atom
count. Startup and structural changes affect this estimate.

Before interpreting production results:

1. Check equilibration, forces, numerical stability and timestep convergence.
2. Repeat independent random seeds; a single trajectory is not a statistical study.
3. Compare sizes (400 and 1000 atoms), temperatures (2500, 2700, 3000, 3300 K), and
   densities. Use different directories/configurations for each case.
4. Inspect trajectories and relaxed snapshots alongside the numerical diagnostics.
5. Record failure to layer as a result; do not modify the starting structure into
   sheets to force agreement.

The initial hard-core random placement, 0.25 fs timestep, 25 fs thermostat damping,
three-variable Nose-Hoover chain, 1.85 A bond cutoff and planarity tolerance are
**our explicit choices**, not values recovered from the main PRL. They require
validation/sensitivity checks. There is no hidden pre-relaxation or melt-quench.
The high-energy random start may need a carefully documented preparation protocol.

The current runner saves frames periodically but does not checkpoint the complete
thermostat state. An interrupted run retains its data and is marked interrupted;
it cannot currently be resumed as an identical dynamical trajectory.

## Outputs and their meanings

| File | Contents |
|---|---|
| `initial.extxyz` | Random starting structure, cell, periodicity |
| `config.json` | Chosen simulation settings |
| `trajectory.traj` | ASE trajectory with energies, forces, momenta and times |
| `trajectory.extxyz` | Portable trajectory for visual inspection |
| `thermo.csv` | Time, temperature, potential/total energy, elapsed wall time |
| `final.extxyz` | Last MD structure |
| `structure.csv` | Coordination and structural proxy time series |
| `rdf_3d.csv` | Final-frame **3D** pair distribution function |
| `analysis.json` | Final geometry, candidate layering mode, ring counts |
| `diagnostics.png` | Energy, temperature, structural proxies, 3D RDF |
| `provenance.json` | Versions, model hashes, actual duration, completion status |

Three neighbors are only a coordination-based proxy for sp2 hybridization. The
planar-threefold fraction additionally checks neighbor geometry. It is not the
paper's odds-ratio planarity fraction. Local normal alignment measures orientation
of planar neighborhoods and is unreliable when few neighborhoods qualify.

The layering diagnostic searches reciprocal-cell density modes with wavelengths
2.5-4.0 A. It is rotation/translation invariant, but scanning several modes raises
the random baseline. Its best wavelength is only a **candidate** spacing; it is
not evidence of actual layers by itself. Its restricted wavelength range is an
analysis choice that can miss layering at other spacings.

Ring analysis counts chordless, non-winding cycles of size 3-8 in a distance-based
bond graph; it excludes periodic loops that do not close in real space. This is a
declared convention, not a verified match to the paper's ring algorithm. Hot-frame
ring counts fluctuate; analyze relaxed structures for a less thermally distorted
picture. Coordination and rings depend on the bond cutoff.

The original structural runner uses spherical 3D RDF normalization. The new
`end_to_end.py` also implements a 2D in-plane RDF when validated layers are found.

## Full workflow and remaining scientific validation

[END_TO_END.md](../END_TO_END.md) explains all implemented stages, commands and
outputs. [END_TO_END_VALIDATION.md](../END_TO_END_VALIDATION.md) records completed
small tests. PySCF and DFT-D3 are installed and electronic calculations have been
run on this Mac. The alternative all-electron backend requires convergence and
cross-validation against the paper's VASP/PAW results.

Production-duration dynamics, independent seeds, electronic convergence and
matching the authors' unavailable settings/data remain necessary. Exact model
identity is unconfirmed. No short software test establishes the paper's claims.

## Validate structural analysis

```sh
.venv/bin/python -m unittest -v test_structure.py
```

Checks cover density/minimum-image distances, graphene coordination and hexagons,
diamond coordination, periodic translation and rotation, known layer spacing,
ideal-gas RDF normalization, and a ring crossing a periodic boundary.

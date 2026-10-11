# Reproducibility guide

## Fresh checkout

Work in the repository root. Install Python 3.12 dependencies and obtain the
external model as described in [README.md](../README.md). The recorded platform
is Apple Silicon macOS with Python 3.12.12; a new machine still needs validation.
No existing virtual environment or local result directory is assumed.

Use `requirements-lock.txt` to match the recorded dependency versions, subject
to wheel availability. It is an environment snapshot, not a cross-platform
package or binary guarantee. Scientific engines are external dependencies.

## Level 1: fast analytic and geometry checks

```sh
.venv/bin/python -m unittest -v test_structure.py test_workflow.py
```

These tests require ASE, NumPy, SciPy and Matplotlib, but do not execute GAP or
DFT. They validate selected geometry, normalization, symmetry and transport
identities, not the material model's accuracy.

## Level 2: actual small calculation

```sh
.venv/bin/python end_to_end.py --config configs/end_to_end_extended_smoke.json --out results/verification
```

This runs sample, graphite and vacancy-reference branches. All completed stages
and their timings are recorded in `manifest.json`. Check `REPORT.md`, individual
SCF logs, relaxation convergence and field-accuracy checks. Defaults are small
software-test settings. A run can finish without forming graphitic layers.

Use a fresh output directory for new settings. For an interrupted run with the
same settings, append `--resume`. Source snapshots distinguish invocations when
code changes between resumed stages. The entire thermostat state is not saved;
an incomplete MD stage restarts rather than continuing the identical trajectory.

## Level 3: optional branches and numerical convergence

```sh
.venv/bin/python validate_electronic.py cell --out results/cell_test
.venv/bin/python validate_electronic.py kpoints --out results/kpoint_test
.venv/bin/python validate_electronic.py exchange --out results/exchange_test
```

The cell test deliberately limits optimization to one step and does not claim
converged equilibrium. The other tests exercise complex Bloch states and compare
screened-exchange evaluations. Check their saved numeric results.

For an unchanged structure, vary electronic settings without rerunning MD:

```sh
.venv/bin/python electronic_convergence.py --structure results/verification/reference/final.extxyz --config configs/end_to_end_smoke.json --plan configs/electronic_convergence.json --out results/electronic_convergence
```

This sample plan is a starting sequence, not a convergence certificate. Broadenings
are already swept in the transport settings. Test basis, k mesh, grids, active
window, smearing, finite-difference steps and observables relevant to the study.

To regenerate the historical dynamics-screening input first:

```sh
.venv/bin/python reproduce.py prepare --config configs/mac_pilot.json --out results/mac_pilot_verified
.venv/bin/python validate_dynamics.py run --out results/dynamics_validation_repeat
```

The dynamics driver expects that specific saved initial geometry path. It does
not need the original raw run to be distributed. If that directory already
exists, inspect and reuse the intended input rather than overwriting it.

## Level 4: scientific reproduction

Prepare a seed study without starting simulations:

```sh
.venv/bin/python run_sweep.py --config configs/research_start.json --grid configs/seed_sweep.json --out results/seed_study
```

Adding `--run` starts the jobs sequentially. Large runs can take days or longer;
benchmark the intended size before committing to a study. Study temperature,
density, size and independent initial conditions, and report negative results.
The research template enables costly cell relaxation and is not a validated
production prescription.

For each scientific result retain the code commit, configuration, potential
hashes, software versions, actual duration, convergence evidence and data used
for the figure. Publish matched comparisons and numerical uncertainty. See
[method differences](METHODS_AND_LIMITATIONS.md) before describing agreement
with the original paper.

## Recorded status

[END_TO_END_VALIDATION.md](../END_TO_END_VALIDATION.md) records the completed
small checks. [VALIDATION.md](../VALIDATION.md) records the earlier dynamics
screening. Raw runs remain local and are excluded from the source bundle.

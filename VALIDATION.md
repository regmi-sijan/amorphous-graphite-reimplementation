# Validation and current status

Validated locally on 2026-10-10, Apple Silicon macOS, Python 3.12.12.
The exact installed package set is recorded in `requirements-lock.txt`.

## Subsequent dynamics validation: completed

The next-stage study is now complete; see
`results/dynamics_validation/REPORT.md` and `validation.png` in that directory.
Three 20 fs NVE runs gave maximum energy errors of 1.380, 0.352 and 0.095 meV/atom
at 0.5, 0.25 and 0.125 fs timesteps respectively. The 0.25 fs timestep passed the
chosen 1 meV/atom short-run screen.

A 300 fs NVT run then completed at a target of 2700 K. The last 100 fs averaged
2724.4 K; all four final 25 fs block averages were within 10% of the target. The
maximum extended-energy error was 0.665 meV/atom. This establishes preliminary
numerical stability and thermal settling for one 80-atom configuration, not
structural equilibrium or layer formation. Final threefold coordination was 70%.

The sections below describe the earlier pilot and its original limitations.

## What ran successfully

- Seven structural tests: periodic placement/density, graphene coordination and
  hexagonal rings, diamond coordination, rotation/translation invariance,
  a known 3.1 A layer spacing, ideal-gas 3D RDF normalization, and a hexagon
  crossing a periodic boundary.
- Loaded the recovered 2017 full two-/three-body/SOAP carbon GAP model.
- Energy-gradient check on one component of a random configuration:
  analytic force 4.5170250342 eV/A, central difference 4.5170264548 eV/A;
  absolute error 1.42e-6 eV/A at displacement 1e-4 A.
- Full prepare -> MD -> analysis -> trajectory export for 80 atoms and 0.01 ps.
  Outputs are in `results/mac_pilot_verified`.
- Ten-step, 400-atom performance benchmark in `results/m5_benchmark`.
- A three-step relaxation I/O test. As expected for this short limit, it is
  **not converged**; `relaxed.extxyz` is not a validated relaxed endpoint.
- Checksum-based extraction of the four potential files from the research archive.

An initial trial in `results/mac_pilot` failed before completing a saved frame
because QUIP 0.10.3 does not initialize its calculator name with newer ASE.
The workflow now sets the calculator's public `name` property. The subsequent
`mac_pilot_verified` run completed. The failed trial is retained as provenance.

## Measured performance

| System | Steps / simulated duration | Wall time | Linear estimate for 100 ps |
|---|---|---|---|
| 80 atoms | 40 / 0.01 ps | 16.16 s | 44.9 h |
| 400 atoms | 10 / 0.0025 ps | 21.82 s | 242.5 h (10.1 days) |

These estimates use 0.25 fs steps and one OpenMP thread. Short benchmarks include
startup effects; changing structures, temperature and hardware load can change
performance. A separate four-thread force benchmark gave approximately 2.09 s per
force call for 400 atoms and did not demonstrate useful acceleration. No full
100 ps run has been launched or completed.

## What these checks do not establish

The 80-atom pilot is only 10 femtoseconds long. Temperature rises from roughly
3000 K to 7600 K as the random starting geometry releases energy. It has **not
equilibrated** to the 2700 K thermostat target. Its approximately 31% threefold
coordination does not meet the paper's >95% sp2-plus-layering criterion.
There is no demonstrated amorphous-graphite formation yet.

At the time of that pilot, no timestep convergence, NVE energy-drift study or
thermostat-settling study had been completed; the subsequent study above now
provides preliminary checks. An independent-seed ensemble and long-time
finite-size comparison remain outstanding.
The numerical force check validates a software interface, not model accuracy.

The potential comes from `training/full_das` in the 2017 data archive. The PRL
cites that work, but does not identify the exact XML file. Matching this specific
file to the 2022 authors' model remains unconfirmed.

The paper's charge densities, orbitals, DOS and conductivity are not reproduced.
The in-plane RDF and the original planarity statistic are not reproduced.
The generated plot is a diagnostic of this pilot, not a recreation of a paper figure.

## Next scientific milestones

1. Establish a stable, equilibrated protocol with timestep and thermostat checks.
2. Benchmark any speed improvements against energies/forces and trajectory
   statistics; do not silently change the potential.
3. Run several sufficiently long trajectories and inspect whether layers emerge.
4. Identify actual layers and calculate their in-plane RDF and spacing.
5. Obtain missing supplemental/settings information and verify the potential file.
6. Converge the implemented electronic workflow and validate physical transport;
   see END_TO_END_VALIDATION.md for the completed small electronic checks.

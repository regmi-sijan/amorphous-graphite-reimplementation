# End-to-end validation on this Mac

Completed October 10, 2026. These are software and numerical-consistency tests,
not a reproduction of the paper's scientific results.

## Completed calculations

- Main run: `results/end_to_end_verified/REPORT.md`.
- Configuration: `configs/end_to_end_extended_smoke.json`.
- All 19 stages completed, about 129.6 seconds of summed stage time.
- Eight random carbon atoms, 2.44 g/cm³, eight MD steps (0.002 ps), then GAP/FIRE
  relaxation. An eight-atom AB graphite reference and six-atom vacancy reference
  were also calculated. Both GAP relaxations converged to their smoke tolerances.
- PBE and HSE06 converged for all three structures. Actual charge/orbital cubes,
  DOS, KG tensors, SPC fields, comparisons and plots were produced.
- Field mesh 64³; all configured field checks passed a 1% tolerance. This tests
  quadrature consistency, not basis/k-point/transport convergence.
- A completed-run `--resume` reused all 19 completed stages successfully. Earlier,
  a stopped HSE06 stage was also backed up and restarted while preserving prior
  stages; its slower direct-integral attempt remains under the earlier run's
  `failed_attempts/` directory.

| Structure | PBE charge integral error | HSE06 charge integral error | Momentum quadrature error | SPC integral error |
|---|---|---|---|---|
| sample | 0.00848% | 0.00848% | 0.04134% | 0.0015143% |
| reference | 0.49891% | 0.49925% | 0.00427% | 2.4569e-09% |
| defect | 0.48979% | 0.49011% | 0.00704% | 0.35268% |

The random sample did not form validated graphitic layers. Its in-plane RDF and
in-plane conductivity ratio are correctly absent. The crystalline reference has
layers and four periodic six-membered rings. The smoke basis and Gamma sampling
make quantitative conductivity ratios physically unreliable. No claim of the
paper's approximately 10^-2 ratio is made.

## Independent tests

- `python -m unittest -v test_structure.py test_workflow.py`: **14 tests passed**.
  They cover density, periodic distances, coordination, translation/rotation,
  RDF normalization, layer rejection for diamond, periodic covering-graph rings,
  Fermi electron count, degenerate transitions, gauge and rotational invariance,
  tensor positivity, and SPC agreement with a dense kernel and analytic
  plane-wave integrals. Analytic test fixtures are not research simulation data.
- `results/cell_relaxation_smoke/`: an actual two-carbon PBE+D3 full-cell/atom
  optimization step, requiring 50 SCF evaluations. Forces/stress and output
  writing worked. The deliberately one-step relaxation **did not converge**;
  this does not establish an equilibrium layer spacing.
- `results/complex_kpoints_verified/`: actual PBE on a 3×1×1 mesh, including
  non-time-reversal-invariant k points. Charge integral error about 1.86e-6,
  momentum error 1.35e-9 and SPC sum relative error 5.53e-8.
- `results/exchange_comparison/`: direct screened exchange and the full-minus-
  long-range HSE06 evaluation both converged on the same periodic C2 test cell.
  Their total energies differed by 2.07e-06 eV.
  This is a small comparison, not a convergence guarantee for every geometry.
- `results/sweep_preparation_test/`: three independent-seed job configurations
  and exact commands generated; production jobs were not launched.
- `results/convergence_driver_verified/`: the convergence driver completed two
  actual PBE calculations on an unchanged two-carbon structure using 32³ and 48³
  field grids, saving the comparison and plots.
- Python compilation checks passed. Structural and electronic plots were
  inspected. The original 3D RDF's half-cell-height cutoff lies below the first
  neighbor for the tiny graphite reference; its 2D RDF includes periodic images
  and correctly shows the 1.42 Å first-neighbor peak.

## Scope still requiring scientific work

No 100–500 ps production study, large-cell HSE06 calculation, converged
zero-pressure graphite optimization, statistical formation-kinetics study or
converged physical DC transport result has been completed. The new code enables
these studies, subject to computational resources and convergence testing.
PySCF all-electron Gaussian DFT substitutes for VASP/PAW; exact original model
identity and missing supplementary settings remain unresolved.

The earlier stability study remains at `results/dynamics_validation/REPORT.md`.
Read `END_TO_END.md` for commands, formulas, defaults and interpretation.

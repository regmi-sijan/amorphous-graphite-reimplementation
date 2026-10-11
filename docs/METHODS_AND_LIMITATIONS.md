# Relationship to the original paper

This is an independent reimplementation inspired by
[Thapa et al. (2022)](https://doi.org/10.1103/PhysRevLett.128.236402).
The paper established the scientific questions and reported findings. This project
provides separate code to explore those questions with accessible software.
It does not supply the authors' original source code, configurations or raw data.

## Method mapping

| Stage | Original paper | This implementation |
|---|---|---|
| Atomic motion | NVT with ab initio PBE, PBE with dispersion, and GAP studies | NVT using the recovered 2017 GAP through ASE/QUIP |
| Model identity | Cites the 2017 carbon GAP work | `training/full_das` variant from the associated archive; exact identity unconfirmed |
| Cell sizes and durations | 80–1000 atoms; total trajectories 100–500 ps | Default software test: eight atoms, 0.002 ps; larger runs configurable |
| Structural relaxation | Conjugate gradient | GAP/FIRE at fixed volume |
| Planarity | Reported odds-ratio-based statistic | Independent local-planarity and reciprocal-density proxies |
| Ring topology | Discusses nonhexagonal ring disorder | Distance-based chordless, non-winding periodic rings; not a verified identical algorithm |
| Electronic engine | VASP with PAW | PySCF periodic all-electron Gaussian orbitals |
| Charge and orbitals | HSE06; interpretation of pi and pi-star states | HSE06 charge and near-Fermi orbital exports; no automatic orbital assignment |
| Transport | Kubo–Greenwood/SPC analysis | Local-PBE canonical momentum, finite Gaussian broadening and a documented divided-difference regularization |
| Zero-pressure structure | Relaxation including van der Waals interactions | Optional numerical PBE+D3 free-energy derivatives and cell/atom optimization |
| Defect comparison | Discussion of a 5-8-5 defect | Adjacent atoms removed and structure relaxed; resulting topology must be inspected |

## Interpreting the outputs

Threefold coordination is a geometric proxy for sp² bonding. Layering requires
multiple checks; a density peak alone is insufficient. An in-plane RDF is omitted
when layers are not established. The paper's numerical order parameter cannot be
reconstructed exactly from the available description.

The all-electron density includes core charge. Its values and visual color scales
cannot be compared uncritically to a PAW-derived density representation. Band
indices and near-Fermi probability densities are not, by themselves, evidence
that a selected state has a particular bonding character.

PBE is used for transport because this implementation does not include the
nonlocal velocity corrections needed for hybrid or pseudopotential Hamiltonians.
The finite-broadening DC estimate is not a scattering-limited experimental DC
conductivity. The complex/signed SPC row integrates to the computed conductivity;
the integral of its absolute value generally does not. See
[the equations and units](../END_TO_END.md).

SCF convergence, field quadrature checks, finite-size convergence, statistical
sampling and agreement with published observables are different validation tasks.
A pass in one does not establish the others. Failure to form layers is a valid
simulation outcome and must not be hidden by constructing a layered sample.

## Known gaps before a quantitative reproduction claim

- Verify the exact original GAP variant and obtain missing numerical settings.
- Establish equilibrated, production-duration trajectories at relevant conditions.
- Test multiple initial states, seeds, cell sizes and analysis thresholds.
- Converge basis, k mesh, SCF and field grids, smearing and transport broadening.
- Complete a converged dispersion-inclusive zero-pressure study.
- Compare matched observables and density definitions with the paper's results.

The code enables an independent study. It does not establish the original paper's
formation kinetics, reported spacing, energy difference or conductivity ratio.

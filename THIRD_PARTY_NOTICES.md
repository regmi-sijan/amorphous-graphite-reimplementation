# Sources, attribution and third-party materials

## Scientific origin

The scientific motivation and published findings belong to R. Thapa,
C. Ugwumadu, K. Nepal, J. Trembly and D. A. Drabold:
[Ab Initio Simulation of Amorphous Graphite, PRL 128, 236402 (2022)](https://doi.org/10.1103/PhysRevLett.128.236402).
This is a separate implementation, with no claimed endorsement or affiliation.
The paper PDF and its published figures are not distributed in this repository.
The workflow diagram is a newly created explanatory diagram, with source credit.

## Carbon potential and data

The force model is an external dependency from the
[Cambridge research-data record](https://www.repository.cam.ac.uk/items/575ab523-3e3a-469b-a1b7-ac22cd4a65b8),
DOI [10.17863/CAM.7453](https://doi.org/10.17863/CAM.7453), associated with
[Deringer and Csányi (2017)](https://doi.org/10.1103/PhysRevB.95.094203).
The source record states CC BY-NC-SA 4.0 except where otherwise noted. Consult
that record and the accompanying files for their applicable terms. This project
neither redistributes nor relicenses the model/training files in its source bundle.
Its provenance and selected-file hashes are documented separately.

## Conductivity method

The SPC analysis is informed by K. N. Subedi, K. Prasai and D. A. Drabold,
[Space-Projected Conductivity and Spectral Properties of the Conduction Matrix](https://doi.org/10.1002/pssb.202000438),
physica status solidi (b) 258, 2000438 (2021). The implementation's broadening,
velocity-operator scope and numerical checks are described in END_TO_END.md.

## External software

Dependencies are installed separately and retain their own licenses and citation
requirements. This document is an attribution inventory, not a relicensing of
their code. Consult the release actually installed:

| Dependency | Role | Upstream source/documentation |
|---|---|---|
| ASE | Structures, dynamics, optimizers and I/O | [ASE](https://ase-lib.org/) |
| QUIP / quippy-ase / GAP | Carbon energies and forces | [QUIP](https://github.com/libAtoms/QUIP), [GAP](https://libatoms.github.io/GAP/) |
| PySCF | Periodic electronic calculations | [PySCF](https://pyscf.org/) |
| simple-dftd3 / dftd3 | Dispersion correction | [simple-dftd3](https://github.com/dftd3/simple-dftd3) |
| NumPy / SciPy | Arrays and numerical algorithms | [NumPy](https://numpy.org/), [SciPy](https://scipy.org/) |
| Matplotlib | Scientific plots | [Matplotlib](https://matplotlib.org/) |
| Graphviz | Optional regeneration of the workflow image | [Graphviz](https://graphviz.org/) |

Other transitive packages are recorded in requirements-lock.txt. No proprietary
VASP executable, PAW dataset, downloaded paper or external model is bundled.

## Original project code

The maintainer has chosen to leave the original workflow code and documentation
unlicensed for now. Citation metadata identifies the collective project contributors;
it does not list the original paper's authors as authors of this software.
The original code and explanatory diagram were developed with AI coding assistance.

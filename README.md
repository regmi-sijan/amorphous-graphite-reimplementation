# Amorphous graphite — an independent Python reimplementation

A Python-led computational workflow inspired by **“Ab Initio Simulation of
Amorphous Graphite”** by R. Thapa, C. Ugwumadu, K. Nepal, J. Trembly and
D. A. Drabold, *Physical Review Letters* **128**, 236402 (2022).
[Read the original paper](https://doi.org/10.1103/PhysRevLett.128.236402).

This is an **independent reimplementation and reproduction study inspired by the
paper**. The original authors receive credit for the scientific ideas and
published findings. This repository supplies separate workflow and analysis
code; it is not the authors' official code or an endorsed implementation.

**Status:** a small end-to-end software test has completed. The paper's full
scientific results have **not** been reproduced. Production trajectories,
statistical sampling and electronic convergence studies remain necessary.

## What the project does

- Generates periodic random carbon and runs real GAP molecular dynamics.
- Relaxes structures and measures coordination, planar geometry, layering,
  periodic rings and pair distributions.
- Computes periodic PBE and HSE06 electronic states using PySCF.
- Exports charge density, orbital densities, DOS, plots and cube files.
- Computes broadened Kubo–Greenwood conductivity tensors and space-projected
  conductivity (SPC), with numerical consistency checks.
- Compares sample, graphite and optional vacancy-reference structures.
- Supports optional PBE+D3 cell relaxation, parameter studies and stage recovery.

The workflow is mostly Python; QUIP, PySCF and other compiled libraries perform
the force and electronic calculations. It does not manufacture layers or fit
outputs to the paper's reported values.

## Workflow

![Scientific workflow of the original paper](paper_workflow.png)

This original explanatory diagram summarizes the paper's workflow. See
[what is measured and how](docs/PAPER_WORKFLOW.md) and
[the method differences](docs/METHODS_AND_LIMITATIONS.md) for the relationship to
this implementation. The editable diagram source is [paper_workflow.dot](paper_workflow.dot).

## Quick start

Run commands from the repository root. Python 3.12 on Apple Silicon macOS was
tested; other platforms and wheel availability have not been verified here.

### 1. Install dependencies

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

For the exact package versions recorded in the local validation, use
`requirements-lock.txt` instead. A virtual environment is not included in Git.

### 2. Obtain the external GAP model

The potential is an external dependency and is **not included** in this
repository. Download the main archive from the
[Cambridge dataset](https://doi.org/10.17863/CAM.7453), then extract the selected
model with checksum verification:

```sh
curl -L --fail 'https://api.repository.cam.ac.uk/server/api/core/bitstreams/8834f03c-6fc9-46af-aeff-1b2517bf8f2d/content' -o aC_GAP_data_main.tar.gz
.venv/bin/python extract_model.py aC_GAP_data_main.tar.gz --out potentials
```

The main archive is approximately 1.19 GB. See [data provenance and retrieval](docs/DATA_AND_MODELS.md)
for the selected variant, checksums and source terms. Its exact identity with the
model used in the 2022 paper is **unconfirmed**.

### 3. Run the small end-to-end test

```sh
.venv/bin/python end_to_end.py --config configs/end_to_end_smoke.json --out results/my_complete_test
```

This uses eight atoms and only 0.002 ps of dynamics, followed by electronic
calculations. It tests software execution; it is far too short and small to
establish amorphous-graphite formation or quantitative transport.

The output includes `REPORT.md`, configuration/provenance files, structures,
electronic states, charge/orbital cubes, conductivity arrays and plots.

### 4. Resume an interrupted workflow

```sh
.venv/bin/python end_to_end.py --config configs/end_to_end_smoke.json --out results/my_complete_test --resume
```

Completed stages are reused. An incomplete stage is backed up and restarted.
An interrupted MD stage restarts from its seed; this is not an exact continuation
of the complete thermostat state. Use the same configuration for a resume.

## Validation and limitations

The recorded extended small test completed **19 stages**, including converged
PBE and HSE06 calculations for sample, graphite and vacancy-reference structures.
Fourteen analytic/geometry tests passed. Separate small calculations exercised
complex k-points, screened-exchange evaluation and one PBE+D3 cell-relaxation step.

Run the fast tests with:

```sh
.venv/bin/python -m unittest -v test_structure.py test_workflow.py
```

See [the validation record](END_TO_END_VALIDATION.md). Raw local run directories
are excluded from the repository; they can be regenerated. Reported past checks
are not a claim that a new installation has passed them.

Material differences from the paper include PySCF all-electron Gaussian DFT in
place of VASP/PAW, independently chosen structural proxies, FIRE in place of
conjugate-gradient relaxation, and a declared finite-broadening transport
regularization. HSE06 supplies charge/orbital analysis; transport uses local PBE.
See [methods and limitations](docs/METHODS_AND_LIMITATIONS.md).

## Documentation and files

| File or directory | Purpose |
|---|---|
| [END_TO_END.md](END_TO_END.md) | Detailed commands, equations, units and outputs |
| [docs/PAPER_WORKFLOW.md](docs/PAPER_WORKFLOW.md) | Scientific questions, measurements and figure mapping |
| [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) | Reproduction from a fresh checkout and validation levels |
| [docs/DATA_AND_MODELS.md](docs/DATA_AND_MODELS.md) | External model provenance, retrieval and hashes |
| [docs/GITHUB_UPLOAD.md](docs/GITHUB_UPLOAD.md) | Uploading the source and documentation |
| [REFERENCES.bib](REFERENCES.bib), [CITATION.cff](CITATION.cff) | Scientific citations and GitHub citation metadata |
| [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) | External software, data and publication attribution |
| `end_to_end.py`, `agworkflow/` | Complete Python-led workflow and analysis |
| `reproduce.py`, `validate_dynamics.py` | Structural MD and stability study |
| `electronic_convergence.py`, `run_sweep.py` | Electronic convergence and parameter-study drivers |
| `validate_electronic.py`, `test_*.py` | Small integration checks and analytic/geometry tests |
| `configs/` | Small-test and research starting configurations |

Large configurations are research starting points, not validated production
recipes. No long production study runs automatically during installation.

## Attribution and citation

If this project informs your work, cite the original scientific paper:

> R. Thapa, C. Ugwumadu, K. Nepal, J. Trembly, and D. A. Drabold,
> “Ab Initio Simulation of Amorphous Graphite,” Physical Review Letters 128,
> 236402 (2022). DOI: [10.1103/PhysRevLett.128.236402](https://doi.org/10.1103/PhysRevLett.128.236402).

Also cite the [carbon GAP model](https://doi.org/10.1103/PhysRevB.95.094203),
[its dataset](https://doi.org/10.17863/CAM.7453), and the
[SPC method](https://doi.org/10.1002/pssb.202000438) when you use those components.
The BibTeX file includes these references. A scientific use of this implementation
should additionally identify this repository and the exact commit or release used.
The source repository is [regmi-sijan/amorphous-graphite-reimplementation](https://github.com/regmi-sijan/amorphous-graphite-reimplementation).
No project DOI has been assigned; the paper's DOI must not be used as the software's DOI.

The workflow code and explanatory diagram were developed with AI coding
assistance. The saved tests document what was checked; scientific conclusions
still require independent scrutiny and convergence testing.

## Licensing and contributions

The maintainer has chosen to leave the original project code unlicensed for now.
External software, model data and publications retain their own terms; see
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). The paper PDF, downloaded model,
training data and generated run directories are not part of the source bundle.

See [CONTRIBUTING.md](CONTRIBUTING.md) for reporting problems and proposing changes.

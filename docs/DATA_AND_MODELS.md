# External data and model provenance

## Original research

The original PRL paper is linked by DOI rather than redistributed:
[10.1103/PhysRevLett.128.236402](https://doi.org/10.1103/PhysRevLett.128.236402).
The authors' original configurations and complete numerical settings have not
been recovered. Project-generated random configurations are independent inputs.

## Selected GAP model

- Model paper: [Deringer and Csányi, PRB 95, 094203 (2017)](https://doi.org/10.1103/PhysRevB.95.094203).
- Data record: [Cambridge, DOI 10.17863/CAM.7453](https://doi.org/10.17863/CAM.7453).
- Archive: `aC_GAP_data_main.tar.gz`, approximately 1.19 GB.
- Subdirectory: `training/full_das`.
- Model label: `GAP_2016_10_18_60_23_11_22_108`.
- Source record terms: CC BY-NC-SA 4.0 except where otherwise noted.
- Exact identity with the model files used for the 2022 PRL: **unconfirmed**.

The historic standalone download was unavailable during initial setup, so the
selected model was recovered from the associated training archive. The archive
contains multiple fits; belonging to this dataset does not prove exact identity
with the later authors' files.

## Retrieval

From the repository root, with dependencies installed:

```sh
curl -L --fail 'https://api.repository.cam.ac.uk/server/api/core/bitstreams/8834f03c-6fc9-46af-aeff-1b2517bf8f2d/content' -o aC_GAP_data_main.tar.gz
.venv/bin/python extract_model.py aC_GAP_data_main.tar.gz --out potentials
```

If the download endpoint changes, use the DOI landing page to locate the same
main archive. The extractor checks the four selected members against known
SHA256 values. It does not verify a checksum of the entire archive. It refuses
a nonempty destination, reads only explicitly named regular files, and executes
no archive scripts. Do not bypass a checksum mismatch without investigating it.

[model-provenance.json](model-provenance.json) preserves the source and expected
member hashes without including the model data. The extractor also creates a
local `potentials/SOURCE.json`. Simulation provenance records the files used.

## Files to retain locally

`.gitignore` excludes `potentials/`, `results/`, downloaded archives, virtual
environments and common scientific output formats. Source documentation and the
new explanatory workflow image are included. The archive builder uses an explicit
source allowlist, so local results, PDFs and external data cannot enter its ZIP
merely because they exist next to the code.

Local validation reports summarize completed checks. The raw trajectories,
wavefunctions and volumes are not distributed by default. If publishing research
data later, use a separately documented data deposit and preserve configuration,
seed, source version, dependency versions, units and provenance. Respect the
applicable external source terms when including third-party material.

# Uploading this project to GitHub

Public source repository: [regmi-sijan/amorphous-graphite-reimplementation](https://github.com/regmi-sijan/amorphous-graphite-reimplementation).
The instructions below describe source bundles and Git uploads for this project.

Suggested About description:

> Independent Python reimplementation inspired by Thapa et al., PRL 128, 236402
> (2022): GAP dynamics, PySCF electronic analysis and Kubo–Greenwood/SPC transport.
> Small tests validated; full paper reproduction remains a research task.

## Source-only bundle

From the project root:

```sh
python3 scripts/make_release_bundle.py
```

This creates `dist/amorphous-graphite-reimplementation.zip`. Extract it and upload
the **contents of the extracted project folder**, not the ZIP itself, if using
GitHub's file-upload interface. Include `.gitignore` and `.gitattributes` (hidden
files on some systems). Preserve the folder structure. The PNG is an explanatory
project diagram; downloaded publications and model data are excluded.

`BUNDLE_MANIFEST.json` inside the ZIP records SHA256 hashes of included files.
The bundle has an explicit allowlist and contains no Git history, local results,
virtual environment, external potential or training archive. Regenerate it after
changing documentation or source files.

## Upload with Git

Alternatively, create an empty GitHub repository under your own account and use
its actual URL. Run the following from the project root. `git init` is needed
only if the folder is not already a Git repository:

```sh
git init
git add README.md CITATION.cff REFERENCES.bib THIRD_PARTY_NOTICES.md CONTRIBUTING.md CHANGELOG.md
git add END_TO_END.md END_TO_END_VALIDATION.md VALIDATION.md .gitignore .gitattributes
git add requirements.txt requirements-lock.txt paper_workflow.dot paper_workflow.png
git add end_to_end.py reproduce.py extract_model.py validate_dynamics.py validate_electronic.py
git add electronic_convergence.py run_sweep.py test_structure.py test_workflow.py
git add agworkflow configs docs scripts
git diff --cached --stat
git status --short
git commit -m "Add independent amorphous graphite reimplementation and documentation"
```

Then connect the actual repository URL and push the current branch:

```sh
git remote add origin https://github.com/regmi-sijan/amorphous-graphite-reimplementation.git
git push -u origin HEAD
```

The URL above is this project's repository. If `origin` already exists, inspect
`git remote -v` before changing it. Set the pushed branch as the
repository's default if needed. The `.gitignore` rules do not remove previously
tracked files; inspect the staged file list before publishing.

## Attribution and repository metadata

- Keep the paper citation and independent-reimplementation statement prominent.
- `CITATION.cff` supplies the original paper as the preferred scientific citation.
  Its software contributor name is a collective label, not a verified personal
  author list. Replace it with the desired actual contributor names before a
  formal software release.
- `repository-code` in CITATION.cff identifies the public source repository.
  Add a software version or DOI only when one has actually been issued.
- Cite the original paper for its science and this repository's commit/release
  for the particular code used. The paper DOI is not this project's software DOI.
- The maintainer chose to leave the original code **unlicensed for now**. No
  LICENSE file or open-source license badge is included. External components
  retain their own terms, as documented in THIRD_PARTY_NOTICES.md.
- The recorded local validation is not a claim of CI coverage or complete
  scientific reproduction. Do not add passing CI or reproduction badges without
  corresponding evidence.

You can upload the prepared source as it stands. Longer simulations and external
model files do not need to be uploaded to make the code repository usable.

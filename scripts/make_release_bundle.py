"""Build a source/documentation ZIP without local models or generated runs.

Uses an explicit allowlist; includes no .git history, environments, results,
external models or paper PDFs. Run from any directory with Python 3.12.
"""
from pathlib import Path
import hashlib
import json
from zipfile import ZipFile, ZIP_DEFLATED

ROOT = Path(__file__).resolve().parents[1]
TOP_LEVEL = [
    '.gitignore', '.gitattributes', 'README.md', 'CITATION.cff', 'REFERENCES.bib',
    'THIRD_PARTY_NOTICES.md', 'CONTRIBUTING.md', 'CHANGELOG.md',
    'END_TO_END.md', 'END_TO_END_VALIDATION.md', 'VALIDATION.md',
    'requirements.txt', 'requirements-lock.txt', 'paper_workflow.dot', 'paper_workflow.png',
    'end_to_end.py', 'reproduce.py', 'extract_model.py', 'validate_dynamics.py',
    'validate_electronic.py', 'electronic_convergence.py', 'run_sweep.py',
    'test_structure.py', 'test_workflow.py',
]
TREES = {'agworkflow': {'.py'}, 'configs': {'.json'}, 'docs': {'.md', '.json'}, 'scripts': {'.py'}}
PREFIX = 'amorphous-graphite-reimplementation'


def main():
    files = [ROOT / name for name in TOP_LEVEL]
    # Include an original-code license only if the maintainer adds one later.
    if (ROOT / 'LICENSE').is_file():
        files.append(ROOT / 'LICENSE')
    for directory, suffixes in TREES.items():
        files.extend(p for p in (ROOT / directory).rglob('*')
                     if p.is_file() and p.suffix in suffixes and '__pycache__' not in p.parts)
    files = sorted(set(files))
    for path in files:
        if path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(ROOT):
            raise ValueError(f'Missing or unexpected source file: {path}')
    contents = {str(p.relative_to(ROOT)): p.read_bytes() for p in files}
    manifest = {name: hashlib.sha256(data).hexdigest() for name, data in contents.items()}
    destination = ROOT / 'dist' / (PREFIX + '.zip')
    destination.parent.mkdir(exist_ok=True)
    with ZipFile(destination, 'w', compression=ZIP_DEFLATED) as archive:
        for name, data in contents.items():
            archive.writestr(f'{PREFIX}/{name}', data)
        archive.writestr(f'{PREFIX}/BUNDLE_MANIFEST.json', json.dumps(manifest, indent=2) + '\n')
    print(f'Created {destination} ({len(files)} files plus SHA256 manifest)')
    print('Extract before uploading to GitHub; no models, results or paper PDF are included.')


if __name__ == '__main__':
    main()

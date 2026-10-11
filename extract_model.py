"""Recover checksum-verified 2017 carbon GAP files from the Cambridge archive.

Usage: python extract_model.py /path/to/aC_GAP_data_main.tar.gz --out potentials
No archive scripts are executed. Only four explicitly named files are read.
"""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile

EXPECTED = {
    "carbon.xml": "03158254caf24f358cfcfc36f26416ac2e87032cf7f6d73e2e2b059554b50a9d",
    "carbon.xml.sparseX.GAP_2016_10_18_60_23_11_22_1081": "75f5ceee8a18a6ad0f5f7360c2af5d0697ca0d44bf73abcaa64604933d3096b1",
    "carbon.xml.sparseX.GAP_2016_10_18_60_23_11_22_1082": "b9780079852a89d231035feabb1fd8bdac35f26087f90d279877f0bd58274a9b",
    "carbon.xml.sparseX.GAP_2016_10_18_60_23_11_22_1083": "93aaf5b32576e8629c90aead381dcdc711b9a0161da9cf61479bf463d99b786a",
}


def extract(archive, out):
    out = Path(out)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty {out}")
    found = {}
    with tarfile.open(archive, "r|gz") as tar:
        for member in tar:
            base = member.name.removeprefix("training/full_das/")
            if member.name != "training/full_das/" + base or base not in EXPECTED:
                continue
            if not member.isfile() or member.size > 40_000_000:
                raise ValueError(f"Unexpected archive member {member.name}")
            data = tar.extractfile(member).read()
            if hashlib.sha256(data).hexdigest() != EXPECTED[base]:
                raise ValueError(f"Checksum mismatch: {member.name}")
            found[base] = data
    if set(found) != set(EXPECTED):
        raise ValueError(f"Archive is missing {set(EXPECTED) - set(found)}")
    out.mkdir(parents=True, exist_ok=True)
    for name, data in found.items():
        (out / name).write_bytes(data)
    (out / "SOURCE.json").write_text(json.dumps({
        "source": "https://doi.org/10.17863/CAM.7453",
        "archive_subdirectory": "training/full_das",
        "license": "CC BY-NC-SA 4.0; see repository record",
        "scope": "Full fit in the 2017 dataset; exact file identity with PRL 2022 is unconfirmed",
        "sha256": EXPECTED,
    }, indent=2) + "\n")
    print(f"Recovered and verified {len(found)} model files in {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive")
    parser.add_argument("--out", default="potentials")
    args = parser.parse_args()
    extract(args.archive, args.out)

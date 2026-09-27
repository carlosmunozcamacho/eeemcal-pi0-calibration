# SPDX-License-Identifier: LGPL-3.0-or-later
"""Check built release contents; no runtime datasets or scientific calculation."""
import argparse
from email.parser import BytesParser
from pathlib import Path, PurePosixPath
import tarfile
import zipfile

PDF = "docs/reports/EEEMCal_pi0_calibration_DIS_SIDIS_v0.1.0.pdf"
FIGURES = ["docs/figures/validation/" + name for name in (
    "calibration_coefficient_distribution_DIS_SIDIS.png",
    "calibration_coefficients_by_crystal_DIS_SIDIS.png",
    "calibration_coefficient_maps_DIS_SIDIS.png",
    "photon_energy_DIS_SIDIS.png",
)]
NOTICES = ("LICENSE", "COPYING", "NOTICE.md", "THIRD_PARTY_NOTICES.md")
FORBIDDEN = {".git", ".ssh", ".venv", "venv", "env", "__pycache__",
             "build", "dist", "cache", "work", "logs", "results",
             ".pytest_cache", ".mypy_cache", ".ruff_cache"}


def check_paths(names):
    for name in names:
        path = PurePosixPath(name)
        assert not path.is_absolute() and ".." not in path.parts, name
        assert not (FORBIDDEN & set(path.parts)), name
        assert path.suffix not in {".pyc", ".root", ".npz", ".log", ".pem", ".key"}, name
        assert not path.name.startswith(".env"), name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dist", type=Path)
    args = parser.parse_args()
    wheels = list(args.dist.glob("*.whl"))
    sources = list(args.dist.glob("*.tar.gz"))
    assert len(wheels) == len(sources) == 1, "Expected one wheel and one sdist"
    with zipfile.ZipFile(wheels[0]) as archive:
        names = archive.namelist()
        check_paths(names)
        metadata_name, = [name for name in names if name.endswith(".dist-info/METADATA")]
        metadata = BytesParser().parsebytes(archive.read(metadata_name))
        assert metadata["Name"] == "eeemcal-pi0"
        assert metadata["Version"] == "0.1.0"
        assert metadata["License-Expression"] == "LGPL-3.0-or-later"
        license_root = str(PurePosixPath(metadata_name).parent) + "/licenses/"
        assert all(license_root + name in names for name in NOTICES)
        assert not any(name.endswith((".pdf", ".png")) for name in names)
        assert all(name.startswith(("eeemcal_pi0/", "eeemcal_pi0-0.1.0.dist-info/")) for name in names)
        print(f"Wheel: {len(names)} entries; software + all four license/notices; no report/figures")
    with tarfile.open(sources[0]) as archive:
        members = archive.getmembers()
        assert all(member.isfile() or member.isdir() for member in members), "Unexpected archive link/type"
        names = [member.name for member in members]
        check_paths(names)
        assert len({PurePosixPath(name).parts[0] for name in names}) == 1
        relative = {str(PurePosixPath(*PurePosixPath(name).parts[1:])) for name in names}
        required = set(NOTICES) | set(FIGURES) | {PDF, "CITATION.cff", "CHANGELOG.md",
            "RELEASE_NOTES_v0.1.0.md", ".github/workflows/ci.yml", "tests/test_model.py"}
        assert required <= relative, sorted(required - relative)
        print(f"Sdist: {len(names)} entries; notices, citation, CI, tests, report and four figures present")
    print("Release contents: PASS")


if __name__ == "__main__":
    main()

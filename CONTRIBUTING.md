# Contributing

Use Python 3.10+ in a virtual environment. From the checkout:

```bash
python -m pip install -e .
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest discover -s tests -p 'test_*.py' -v
```

The focused tests generate small synthetic fixtures in temporary directories;
they require neither central production files nor Slurm. Scheduler retry tests
use mocks. Keep scientific changes separate from packaging/documentation work,
preserve uint64 cellIDs and assigned-energy semantics, and add focused tests
for changed behavior. Do not retune defaults merely to match a saved result.

Describe input/configuration changes, validation performed and remaining limits
in contributions. Preserve LGPL/SPDX notices, Axel's pinned source attribution
and the Jones citation. Use generic paths in examples and never commit data,
credentials, logs, virtual environments or build products.

The saved full compact regression is evidence for the documented implementation;
do not run `tests/full_regression.py` as part of routine unit tests or replace
its evidence with a new study. A future scientific change needs an explicitly
planned regression with the appropriate compact inputs. No production data are
included here.

For EICRecon integration, discuss the native cellID payload, baseline/geometry
validity and reconstruction requirements with the relevant ePIC developers;
see [integration notes](docs/eicrecon_integration.md). This package derives
candidate constants and does not implement a conditions service.

## Portable CI and distribution checks

[CI](.github/workflows/ci.yml) tests Python 3.10 (the declared minimum), 3.12
and 3.13 on Ubuntu. These are representative versions satisfying the current
`requires-python >=3.10` metadata, not an upper-bound declaration. It installs
with `python -m pip install . build twine cffconvert`, runs the portable suite
above, checks import/CLI and citation metadata, then runs:

```bash
cffconvert --validate
python -m build
python -m twine check --strict dist/*
python tests/check_release.py dist
```

Start with an empty `dist/` directory when rebuilding a release. The workflow
installs the fresh wheel in a separate clean environment and checks it from
outside the checkout with no `PYTHONPATH`. The wheel excludes documentation
assets; the sdist deliberately includes the authoritative PDF and four PNGs.
The workflow follows the official [checkout](https://github.com/actions/checkout)
and [setup-python](https://github.com/actions/setup-python) interfaces;
[cffconvert](https://github.com/citation-file-format/cffconvert) validates CFF.

## Optional full compact regression (site allocation only)

This is deliberately **excluded** from portable CI: it needs saved large
compact inputs and reference products. The existing PASS remains the release
evidence. Do not rerun it for documentation/packaging edits. For a separately
approved future numerical change, from a disposable source copy inside a site
allocation, configure `examples/dis_regression_jlab.yaml` with your saved
manifest, then use:

```bash
python tests/full_regression.py \
  --compact /shared/pi0/saved-broad-compact \
  --reference /shared/pi0/saved-relaxed-B-min10-reference \
  --config examples/dis_regression_jlab.yaml \
  --output /shared/pi0/new-regression \
  --work /shared/pi0/new-regression-work \
  --source "$PWD" --archive /shared/pi0/new-regression-source.tar.gz
```

The helper overwrites `docs/full_regression.md`/`.json` in the supplied source
copy and compares with the reference; preserve the released evidence separately.
Paths are explicit site placeholders. This command does not submit jobs.

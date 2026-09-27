# Producer validation

**Full compact DIS regression: PASS.** Slurm **11839626** completed with
ExitCode **0:0**. The [report](full_regression.md) and
[numerical comparison](full_regression.json) retain the validated results.
The JSON is byte-for-byte unchanged; only non-numerical report wording was
cleaned up for the public source package. **The full regression was not rerun.**

The comparison covered 11,379 files / 4,967,376 reconstructed events against
the saved relaxed-B/min10 result: 45,385 training-window pairs, 2,740
active/retained parameters/modes, median gain 1.05040330855, and 2,678 channels
half-stable below 10%. Held-out peak: 128.442359 → 134.646174 MeV; relative
width: 3.87069966% → 3.92835328%. Event/pair IDs, training flags, cellIDs,
sparse structure, per-cell support and window-pair counts matched exactly.

Saved tolerances: exact discrete counts/IDs/rank; sample floating arrays
rtol=3e-14, atol=1e-18; gains and fitted mass/width/window rtol=2e-5, atol=2e-7
(GeV for masses); half-change diagnostics rtol=2e-4, atol=2e-5. The JSON keeps
the actual differences. No selection or solver retuning was performed.

## Local release-candidate checks — 2026-09-27

**11 portable tests passed; 0 failures; 0 skips**, in 18.423 seconds on Python
3.12.3. A clean virtual environment installed this source with pip, without
an editable/source-path import. Command from the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1 \
  python -m unittest discover -s tests -p 'test_*.py' -v
```

The environment used NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2, Uproot 5.7.6,
Awkward 2.14.0 and PyYAML 6.0.3. No scientific source, configuration defaults,
existing unit assertions, full-regression evidence or final PDF/PNG was changed.
The archive helper's file allowlist was extended for release metadata only.

CITATION.cff validates against CFF 1.2.0. Internal documentation links resolve.
The [CI workflow](../.github/workflows/ci.yml) will test Python 3.10, 3.12 and
3.13 after publication; no hosted CI run is claimed by this local check.
The full compact regression helper is deliberately excluded from portable
unittest discovery and was not executed. No central production files,
reconstruction or calibration jobs were run for release preparation.

## Focused public-package checks — 2026-09-26

**11 tests passed in 14.243 seconds**, Python 3.12, against a non-editable
wheel installation. The suite uses generated toy fixtures and temporary
paths; it does not need central production files or a scheduler.

```bash
python -m pip install --force-reinstall --no-deps /path/to/eeemcal_pi0-0.1.0-py3-none-any.whl
cd "$(mktemp -d)"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest discover \
  -s /path/to/producer/tests -v
```

Coverage includes derivatives with asymmetric/shared showers, ideal gain
recovery, unsupported-channel handling, unity response, uint64 IDs, additive
shard equations, configuration/split/ownership, truth-free synthetic ROOT
extraction, missing track relations, cache compatibility, normalization,
mocked scheduler retries, installed CLI outputs, cache reuse, plot regeneration
and no-peak reporting.

Additional path checks exercise local and Slurm dry runs with relative paths
and spaces: both worker scripts receive the resolved work plan, the same
physics configuration and the configured output/cache/work locations. Logs
and task manifests derive from output. Input/output and symlink collisions
are rejected before writing; compact reuse leaves every input-cache file
unchanged. No batch-processing logic was changed.

Both `pi0-calibrate --help` and `python -m eeemcal_pi0 --help` ran successfully.
A separate 2,200-event generated fixture exercised `--stage extract` followed
by `--stage analyze`, producing eight active gains and numerical PASS. Status
with no recorded jobs returned successfully. These are focused interface
checks, not production extraction or a substitute for the saved regression.

Live Slurm submission/status/retry, site containers and XRootD access were not
exercised during this cleanup. Scheduler behavior was checked with mocks and
dry runs. Generic dataset/site examples require real paths and resources;
they were not submitted. Source-checkout installation instructions use the
standard pip interface; the focused suite here used the rebuilt wheel.

All 17 runtime module ASTs, scientific defaults and calibration parameters
remain unchanged. Only source attribution comments, report-generation wording,
packaging, examples, documentation and focused tests changed. The full
regression helper was not executed. LICENSE/COPYING texts are unchanged.

## Interpretation and integration boundaries

The package is **LGPL-3.0-or-later**; see [release notes](release.md) and the
license/attribution files. Candidate coefficients still need detector review.

**Full reconstruction feedback-loop validation: not yet performed. Native
EICRecon cellID-to-gain consumer: still to be integrated.** The temporary
samplingFraction-expression replay was a bounded intermediate check, not the
production gain-loading mechanism. This numerical PASS does not establish
real-data accuracy, detector-systematic uncertainties, significance of the
slight relative-width increase, or approval to apply every candidate gain.
General non-HepMC weighted normalization and arbitrary detector schemas need
explicit supported adapters.

# eeemcal-pi0-calibration v0.1.0 — local release candidate

| Item | Status |
|---|---|
| Producer numerical regression | **PASS**, Slurm **11839626**, COMPLETED, ExitCode **0:0** |
| License | **LGPL-3.0-or-later** |
| Redistribution clearance | **Resolved** |
| Full reconstruction feedback-loop validation | **Not yet performed** |
| Native EICRecon cellID-to-gain consumer | **Still to be integrated** |
| GitHub publication | **Not performed**; local release candidate awaiting publication approval |

The package derives candidate relative gains with fixed shower membership and
directions. The previous samplingFraction-expression replay was a bounded
intermediate check, not the production gain-loading mechanism. Numerical
regression does not establish detector-systematic uncertainties or approve
constants for production.

## Numerical evidence

The existing [full regression report](full_regression.md) and unchanged
[numerical JSON](full_regression.json) retain these results:

| Quantity | Result |
|---|---:|
| Files / reconstructed events | 11,379 / 4,967,376 |
| Training-window pairs | 45,385 |
| Active parameters / retained modes | 2,740 / 2,740 |
| Median gain | 1.05040330855 |
| Half-stable below 10% | 2,678 |
| Held-out peak, before → after | 128.442359 → 134.646174 MeV |
| Relative width, before → after | 3.87069966% → 3.92835328% |

Event/pair IDs, training flags, cellIDs, sparse structure, per-cell support and
training-window pair counts matched exactly. Public documentation edits do not
change these measurements. **The full regression has not been rerun.**

## Public package preparation

The standalone repository boundary is this producer directory. The public
project name is `eeemcal-pi0-calibration`; distribution `eeemcal-pi0`, import
`eeemcal_pi0` and version 0.1.0 remain unchanged. No remote, tag or public
release is created during local preparation.

The [README](../README.md) introduces installation and the synthetic workflow.
Detailed [workflow](workflows.md) and [configuration](configuration.md) guides
retain site/path/cache contracts. [CITATION.cff](../CITATION.cff) records the
confirmed author; repository URL, DOI and publication date are omitted until
established. [Changelog](../CHANGELOG.md) and
[release notes](../RELEASE_NOTES_v0.1.0.md) describe the delivered scope.

The authoritative [validation PDF](reports/EEEMCal_pi0_calibration_DIS_SIDIS_v0.1.0.pdf)
and [four final figures](figures/validation/README.md) are preserved unchanged
and included as repository documentation. The software source distribution
includes them; the wheel remains software-only apart from license notices.
Attach the PDF separately when a GitHub Release is authorized.

[LICENSE](../LICENSE), [COPYING](../COPYING), [NOTICE.md](../NOTICE.md) and
[THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md) preserve LGPL-3.0-or-later,
Axel's pinned upstream provenance and the Jones citation. All four accompany
both distributions; wheel notices reside under `.dist-info/licenses/`.
Neither upstream article/notebook nor dependency implementation is bundled.

Fresh standard-build artifacts reside in ignored `dist/`:

```text
eeemcal_pi0-0.1.0.tar.gz
eeemcal_pi0-0.1.0-py3-none-any.whl
```

The source archive contains package code, tests, examples, schemas, CI,
notices and deliberate documentation. It excludes Git internals, nested build
products, environments, caches, credentials and scientific/runtime datasets.
The [portable checks](validation.md) preserve the distinction between
numerical regression and full detector validation.

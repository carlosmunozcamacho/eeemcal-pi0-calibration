# eeemcal-pi0-calibration

Derive **candidate per-crystal relative gains** for the ePIC backward
electromagnetic calorimeter (EEEMCal) from reconstructed π0 → γγ showers.
The Python distribution is `eeemcal-pi0`; the import package is `eeemcal_pi0`.
Version **0.1.0**, licensed **LGPL-3.0-or-later**.

The two-photon mass constrains the energy response. This producer retains each
crystal's assigned shower energy and solves the constrained Jones-style gain
system, with separate normalization of the two showers and shared-crystal
handling. Gains multiply the existing reconstruction baseline. Memberships,
directions and shower-response factors remain fixed during this calculation.

## Install

Python **3.10+** and a Unix-like environment are required. From this checkout:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install .
pi0-calibrate --help
```

Pip installs the numerical, plotting and ROOT-reading Python dependencies.
Portable analysis does not require EICRecon, CVMFS or Slurm. Remote XRootD
inputs additionally need a compatible XRootD environment supplied by the site.

## Try a synthetic example

This small, truth-free fixture is generated locally; no production data are
needed. Run from the repository root using a new output directory:

```bash
python tests/fixture.py results/demo_input
pi0-calibrate --input results/demo_input/manifest.txt \
  --config examples/calibration.yaml --output results/demo --executor local
```

Inspect `results/demo/report.md`, `gains.csv`, `summary.json` and `plots/`.
Missing luminosity or geometry identifiers in this toy example produce QA
warnings; they do not make its coefficients detector constants.

## Use reconstructed showers

The input must provide cluster energies/xyz, exact uint64 hit cellIDs,
cluster-to-hit relations and assigned hit-energy contributions. Neutral
selection needs charged-track matches. A reconstructed vertex is preferred;
a configured nominal IP is the documented fallback. MC truth is optional.

The default relaxed Selection B requires neutral photons with E ≥0.2 GeV,
asymmetry <0.8, overlap <0.2 and different seeds. One eligible pair per event
is chosen by largest E1+E2, without using mass or truth. A deterministic 80/20
split separates training and held-out events. The training peak defines the
frozen calibration window.

The package provides compact extraction, sparse/chunked equations, iterative
positive gains, support/rank checks, training-half stability, optional
fixed-shower injected-gain closure, held-out diagnostics and exposure helpers.
Unsupported channels remain unchanged and flags identify candidate limitations.
Defaults include minimum support 10 and `rcond=0.001`; support is not a
statistical precision guarantee.

- [Configuration examples](examples/calibration.yaml), including a
  [minimum-support-20 override](examples/min20.yaml).
- [Workflow guide](docs/workflows.md): local execution, compact reuse,
  outputs, Slurm arrays, resume and external orchestration.
- [Configuration and path contracts](docs/configuration.md),
  [gain fields](schemas/gains.schema.json) and [normalization](docs/normalization.md).
- [Site-specific Slurm example](examples/jlab.yaml): replace site placeholders;
  portable local execution needs none of those facilities.

## Validation

The authoritative [v0.1.0 DIS and SIDIS/photoproduction report](docs/reports/EEEMCal_pi0_calibration_DIS_SIDIS_v0.1.0.pdf)
and [four final comparison figures](docs/figures/validation/README.md) document
18×275 GeV studies using the common relaxed selection, per-crystal gains,
training-half stability and independent held-out response.

**Producer full compact DIS regression: PASS.** The preserved
[regression evidence](docs/full_regression.md) covers 11,379 files,
4,967,376 reconstructed events and 45,385 training-window pairs, with
2,740 active parameters and retained modes. It reproduces event identities,
per-cell support and the saved numerical solution without retuning. The
held-out peak improves in scale; a width improvement is not universal.

[Portable tests](CONTRIBUTING.md) cover contribution extraction, deterministic
selection/splitting, separate-shower derivatives, shared cells, rank/support
handling, ideal gain recovery and CLI workflows. The physics report includes
fixed-shower injected-gain closure; this does not establish detector accuracy.

The full **gains → EICRecon → reclustering/reselection → validation** feedback
loop remains a later detector-level validation step. A native cellID-to-gain
consumer still needs integration. The previous samplingFraction-expression
replay was a bounded intermediate check, not a production gain-loading
mechanism. See [EICRecon handoff](docs/eicrecon_integration.md).
These candidate gains are **not official production detector constants**;
real-data calibration is not yet validated.

## License, provenance and citation

See [LICENSE](LICENSE), incorporated GPLv3 terms in [COPYING](COPYING),
[NOTICE.md](NOTICE.md) and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
The adapted solver retains Axel Perez Ruiz / `master_thesis` provenance at
commit `c36e96ce463b3bfc29abac40243cf5825400c5d1`. The negative sign of D was
already present upstream.

The method is R. T. Jones et al., NIM A 566 (2006) 366–374,
[DOI 10.1016/j.nima.2006.07.061](https://doi.org/10.1016/j.nima.2006.07.061),
Section 3, equations (1)–(9). The article and upstream notebook are not bundled.
Software citation metadata: [CITATION.cff](CITATION.cff).
See the [changelog](CHANGELOG.md) and [v0.1.0 release notes](RELEASE_NOTES_v0.1.0.md).

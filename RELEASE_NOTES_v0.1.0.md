# eeemcal-pi0-calibration v0.1.0

First public release candidate of the EEEMCal π0 calibration producer. The
`eeemcal-pi0` distribution derives candidate per-crystal relative gains from
reconstructed showers using a constrained Jones-style solve. It preserves
assigned cell energies and supplies support/rank, stability and held-out QA.

- Compact extraction and mass-independent relaxed Selection B.
- Deterministic training/held-out split, frozen training mass window, sparse
  equations and positive multiplicative gains with numerical safeguards.
- Split-sample stability, optional fixed-shower closure, plots, gain CSV and
  normalization helpers.
- **Producer full compact DIS regression: PASS**, with exact identity/support
  checks and preserved [numerical evidence](docs/full_regression.md).
- DIS and low-Q² SIDIS/photoproduction validation in the
  [v0.1.0 report](docs/reports/EEEMCal_pi0_calibration_DIS_SIDIS_v0.1.0.pdf)
  and [final figures](docs/figures/validation/README.md).

Python 3.10+ and pip-installed dependencies suffice for portable local/toy
workflows. Slurm, shared storage, containers and remote XRootD access require
site configuration; CI needs no JLab or physics datasets.

License: **LGPL-3.0-or-later**. Redistribution clearance is resolved; provenance
and notices accompany the software. Coefficients remain candidates: the native
EICRecon cellID gain consumer and full reconstruction/reclustering/reselection
feedback validation are not yet complete. The temporary expression replay is
not a production interface. Real-data calibration is not yet validated.

At publication, attach the unchanged PDF separately alongside the wheel and
source distribution. The source distribution also includes the report and
figures; the wheel contains software and license notices only.

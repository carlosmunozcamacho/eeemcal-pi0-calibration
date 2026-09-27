# Changelog

## v0.1.0 — first public release candidate

Publication date and repository URL will be recorded at publication.

- Portable `pi0-calibrate` CLI, configuration examples and local/Slurm workflows
  with scalable compact extraction, restartable shards and one global solve.
- Persisted per-hit/per-crystal assigned cluster-energy contributions, uint64
  cellIDs, xyz/vertex geometry and original shower response.
- Relaxed Selection B: neutral photons ≥0.2 GeV, asymmetry <0.8, overlap <0.2,
  different seeds and one pair/event ranked by energy sum without mass or truth.
- Deterministic training/held-out partition, training-only peak fit and frozen
  calibration window; independent held-out mass response.
- Constrained Jones-style gain determination with separate two-shower Jacobian
  normalization and contributions from crystals shared by both showers.
- Sparse/chunked equations, eigenmode/rank handling, positive multiplicative
  updates, step/convergence safeguards and configurable minimum event support.
- Split-sample stability, fixed-shower injected-gain closure and explicit
  unsupported/unstable/rank-sensitive candidate flags.
- Optional normalization/exposure helpers, separate from gain fitting; bounded
  compatible HepMC metadata reading and explicit supplied data normalization.
- Relative gain CSV/schema and EICRecon handoff documentation. Native gain
  consumption and the complete reconstruction feedback loop remain integration
  tasks; the earlier expression replay is not the production interface.
- Full compact DIS producer regression PASS, preserving the saved comparison
  without retuning. DIS and low-Q² SIDIS/photoproduction fixed-shower validation
  are documented in the [authoritative report](docs/reports/EEEMCal_pi0_calibration_DIS_SIDIS_v0.1.0.pdf).
- Inclusion of that report and [four curated figures](docs/figures/validation/README.md),
  portable CI, citation metadata, LGPL-3.0-or-later texts and source provenance.

Candidate coefficients are not official production detector constants.

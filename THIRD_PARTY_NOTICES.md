# Third-party notices and implementation provenance

Package: **eeemcal-pi0 0.1.0**. Package license:
**LGPL-3.0-or-later**; see [LICENSE](LICENSE), [COPYING](COPYING), and
[NOTICE.md](NOTICE.md). These notices accompany both source and wheel releases.
The package license also covers its project-authored tests, documentation,
examples and schemas; JSON files remain valid JSON without comment headers.
Existing copyright in adapted material is retained, not reassigned by this
notice. No unsupported copyright year or institutional ownership is asserted.

## Axel Perez Ruiz / master_thesis

Upstream repository: https://github.com/AxelPerezR/master_thesis

Pinned source revision: **c36e96ce463b3bfc29abac40243cf5825400c5d1**.

- [pi0_calibration.ipynb](https://github.com/AxelPerezR/master_thesis/blob/c36e96ce463b3bfc29abac40243cf5825400c5d1/eeemcal_simulation/pi0_energy_calibration/pi0_calibration.ipynb): zero-based notebook cell 3, `minimize_eps`, especially lines 19–38 (C, D, L, B and the constrained solve), and cell 4, especially lines 36–115 (iteration/cumulative gains).
- [reconstruct_pi0_mass_optimized.C](https://github.com/AxelPerezR/master_thesis/blob/c36e96ce463b3bfc29abac40243cf5825400c5d1/eeemcal_simulation/pi0_energy_calibration/reconstruct_pi0_mass_optimized.C): historical topology/direction/mass and extractor reference; no substantial copied C++ block was identified in the current Python extractor by the existing audit.

The directly documented notebook-derived implementation is
`src/eeemcal_pi0/model.py`, through its immediate EEEMCal parent
`calibration/pi0/src/model.py`. `equations`/`constrained_step` and the iterative
`fit` have that lineage. Refactoring and mathematical corrections do not erase
it; this is not claimed to be a clean-room implementation.

Substantive project changes include actual assigned per-cell shower
contributions instead of nearest-crystal cluster assignment; separately
normalized shower derivatives; preserved xyz and vertex geometry; matched
held-out comparisons; sparse/chunked equations; frozen unsupported channels
with their energy kept in shower totals; eigenmode truncation; positive
multiplicative updates and convergence/constraint-merit checks. **D already
had its negative sign in Axel's `minimize_eps` (cell 3, line 25)**. Correcting
that attribution changes documentation only. The producer retains the tested
EEEMCal numerical implementation.

## EEEMCal project contributions

Immediate parent paths below are relative to the EEEMCal project root, not
additional bundled third-party packages. They identify the project
implementations adapted into this producer. New packaging, orchestration
and interface code was developed within the EEEMCal project.

| Producer module / region | Immediate project origin |
|---|---|
| `model.py`, all numerical definitions | `calibration/pi0/src/model.py`; notebook ancestry described above |
| `analysis.py`, matrix/run_fit/training/support/halves | `src/selection_B_calibration.py`, `src/dis_analysis.py`, `src/support_threshold_study.py` under `calibration/pi0/` |
| `candidates.py`, CandidateSink/load/hit/truth_proxy | `calibration/pi0/src/broad_candidates.py` |
| `diagnostics.py`, pair_features/event_ids/meaningful_usage | `calibration/pi0/src/diagnostics.py` |
| `extract.py`, relations/vertex/contribution extraction | `calibration/pi0/src/extract.py`; earlier prototype inspiration distinguished above |
| `mc_diagnostics.py`, fields/EventTruth | `calibration/pi0/src/mc_diagnostics.py` |
| `peak_diagnostics.py`, peak/background/window routines | `calibration/pi0/src/peak_diagnostics.py`, earlier `pilot.find_peak` |
| `selection.py`, choose/pair_shard/build_sample | `calibration/pi0/src/broad_candidates.py`, `sidis_compact.py`, `selection_B_calibration.py` |
| `truth.py`, track checks / optional setup | `calibration/pi0/src/mc_diagnostics.py` |
| `config.py`, deterministic partition only | `calibration/pi0/src/extract.py` |
| `io.py`, small clean/JSON helper pattern | `calibration/pi0/src/dis_analysis.py` |
| `normalization.py`, metadata interpretation | `calibration/pi0/src/sample_luminosity.py`; bounded reader is new |
| `tests/test_model.py`, toy/PhysicsTests | `calibration/pi0/tests/test_model.py` |

`batch.py`, `cli.py`, the module entry/version declarations, most configuration
and ownership handling, reporting implementation and other test/helpers are
new project work. Existing project plot/scheduling concepts and input/API
contracts informed them. No additional third-party implementation was found
in these files. This does not attribute all project modules to Axel merely
because they call the adapted solver.

## Published method and data-format references

R. T. Jones et al., *A bootstrap method for gain calibration and resolution
determination of a lead-glass calorimeter*, Nuclear Instruments and Methods A
566 (2006) 366–374, [DOI 10.1016/j.nima.2006.07.061](https://doi.org/10.1016/j.nima.2006.07.061),
Section 3, equations (1)–(9). This is the mathematical reference, not a software
license. The Elsevier article/PDF, text and figures are not redistributed.
Section 4's resolution-extraction study is outside the implementation.

[HepMC3 GenCrossSection](https://eic-code-browser.sdcc.bnl.gov/lxr/source/include/HepMC3/GenCrossSection.h)
is the format reference for cross-section metadata in pb; no HepMC/Pythia
implementation is copied into this package. Interoperation with EICRecon/ROOT
likewise does not bundle their implementation or replace their own licenses.

## Separately installed dependencies

The source distribution and wheel do not vendor NumPy, SciPy, Matplotlib,
Uproot, Awkward or PyYAML. They retain their own copyright and license notices:

- [NumPy](https://github.com/numpy/numpy/blob/main/LICENSE.txt): BSD terms, NumPy Developers.
- [SciPy](https://github.com/scipy/scipy/blob/main/LICENSE.txt): BSD terms, Enthought/SciPy Developers and other notices for included components.
- [Matplotlib](https://github.com/matplotlib/matplotlib/blob/main/LICENSE/LICENSE): Matplotlib license and historical PSF-based terms.
- [Uproot](https://github.com/scikit-hep/uproot5/blob/main/LICENSE): BSD-3-Clause, Scikit-HEP Project.
- [Awkward](https://github.com/scikit-hep/awkward/blob/main/LICENSE): BSD-3-Clause, Jim Pivarski and listed contributors.
- [PyYAML](https://github.com/yaml/pyyaml/blob/main/LICENSE): MIT, the copyright holders identified in the installed version.

Those links identify upstream notices; the actual installed versions' complete
license files govern those dependencies. A future bundled environment must
retain their required notices, including transitive and bundled components.
The LGPL license texts included here are verbatim standard texts, copyright
2007 Free Software Foundation, Inc.; their own copying notice is retained.

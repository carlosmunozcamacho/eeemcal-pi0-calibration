# Configuration reference

All keys below are accepted by `--config YAML` and `--set dotted.key=YAML`.
Unknown keys fail. One configuration file overlays built-in defaults; repeated
`--set` overrides are applied afterward, in order. `examples/min20.yaml`
overlays the defaults, not a second user configuration file. This table lists
the complete defaults.

## Path resolution and output ownership

- CLI input/compact/output/cache/work/config paths resolve against the invoking
  working directory. Manifest entries resolve against the manifest's directory;
  `root://` entries remain URLs. Blank lines and whole-line comments are allowed.
  Duplicate normalized strings are rejected; distinct aliases for identical
  events must be resolved by the caller. Local manifest identities preserve
  symlink spelling; use `input.identity_map` for deliberate relocation.
- Physics `input.inventory`, `input.expected_manifest`, `input.identity_map`
  and `normalization.metadata` paths resolve against the physics YAML's
  directory (or the working directory without `--config`), including values
  supplied through `--set`. These fields and local manifest entries expand `~`.
  No YAML field performs `$VAR` or `${VAR}` expansion. The shell may expand
  variables in CLI arguments before invocation.
- Site YAML `launcher` and `python` are literal command arguments, not paths
  resolved against the YAML location. Use explicit absolute shared paths or
  executables on the compute-node PATH. `~` and shell variables are not expanded
  there. If `python` is null, the submission interpreter is used and must be
  accessible from the batch node/container.
- `--output` owns results, `logs/` and `batch/`; these two subdirectories have no
  independent relocation flag. The work plan is `batch/work_plan.json`, scripts
  are `batch/extract.sh` and `batch/analyze.sh`, manifests are under
  `batch/task_manifests/`, and submission records are `batch/submissions.json`.
  `--plan` is an internal worker argument, not an output-location setting.
- `--cache-dir` sets newly extracted broad-cache storage (default: `OUTPUT/cache`).
  Existing `--compact` inputs are read-only. `--work-dir` sets selected-sample
  buffers (default: `OUTPUT/work`); work products are not automatically deleted.
  Farm output/cache/work/interpreter paths must be shared across nodes; local
  execution can use local scratch. The global plan stores resolved absolute
  paths and the same physics configuration used by local workers.
- Output/cache/work cannot overlap protected inputs. Output cannot equal,
  contain or lie inside the compact input; symlink-resolved overlaps are
  rejected. Cache and work cannot overlap each other or be ancestors of output.
  Directories are created with ownership records; unowned existing output,
  cache or work directories are refused.
- Existing output requires `--resume` or `--overwrite` with the same run identity
  and configuration. Resume reuses completed work. Overwrite replaces known
  analysis products and reuses signature-compatible completed extraction shards;
  it retains caches and unknown files, and refuses queued/running workflow jobs.
  Use a new output directory for a different selection or configuration.

## Physics defaults

| Key | Default |
|---|---|
| `input.clusters` | `"EcalEndcapNClusters"` |
| `input.hits` | `"EcalEndcapNRecHits"` |
| `input.event_header` | `"EventHeader"` |
| `input.primary_vertices` | `"PrimaryVertices"` |
| `input.track_matches` | `"EcalEndcapNTrackClusterMatches"` |
| `input.tracks` | `"CentralCKFTracks"` |
| `input.mc_particles` | `"MCParticles"` |
| `input.associations` | `"EcalEndcapNClusterAssociations"` |
| `input.chunk` | `128` |
| `input.candidate_floor_GeV` | `0.15` |
| `input.mode` | `"auto"` |
| `input.geometry_id` | `null` |
| `input.cellid_encoding` | `null` |
| `input.baseline_id` | `null` |
| `input.inventory` | `null` |
| `input.expected_manifest` | `null` |
| `input.identity_map` | `null` |
| `input.legacy_compatibility_acknowledged` | `false` |
| `vertex.allow_nominal` | `true` |
| `vertex.nominal_mm` | `[0.0, 0.0, 0.0]` |
| `vertex.max_xy_mm` | `20.0` |
| `vertex.max_z_mm` | `300.0` |
| `split.training_fraction` | `0.8` |
| `split.seed` | `"pi0-v1"` |
| `split.convention` | `"blake2b-mod5"` |
| `selection.energy_GeV` | `0.2` |
| `selection.neutral` | `true` |
| `selection.asymmetry` | `0.8` |
| `selection.overlap` | `0.2` |
| `selection.different_seeds` | `true` |
| `selection.ranking` | `"largest_energy_sum"` |
| `peak.fit_range` | `[0.04, 0.24]` |
| `peak.bins` | `50` |
| `peak.background_form` | `"exponential"` |
| `peak.criteria.entries` | `60` |
| `peak.criteria.signal` | `15` |
| `peak.criteria.significance` | `4.0` |
| `peak.criteria.mean_bounds` | `[0.076, 0.194]` |
| `peak.criteria.sigma_bounds` | `[0.0041, 0.0279]` |
| `peak.window_sigma` | `2.0` |
| `solver.min_usage` | `10` |
| `solver.meaningful_fraction` | `0.05` |
| `solver.rcond` | `0.001` |
| `solver.iterations` | `24` |
| `solver.max_step` | `0.1` |
| `solver.step_tolerance` | `0.0001` |
| `solver.bias_tolerance` | `1e-06` |
| `solver.chunk` | `2048` |
| `qa.halves` | `true` |
| `qa.half_seed` | `20260924` |
| `qa.stability_warn` | `0.05` |
| `qa.stability_fail` | `0.1` |
| `qa.low_support` | `20` |
| `qa.projection_min` | `0.99` |
| `qa.width_warning` | `0.05` |
| `qa.closure` | `false` |
| `qa.closure_seed` | `41256` |
| `qa.detailed_events` | `false` |
| `qa.formats` | `["png", "pdf"]` |
| `qa.histogram_bins` | `60` |
| `normalization.metadata` | `null` |
| `normalization.instantaneous_cm2_s` | `[1e+33, 1e+34]` |

`input.mode=data` skips optional MC branches; auto/simulation read supported truth only if present. No mode requires MC for gains. Track-match existence defines a charged match even for zero-weight matches. Missing required contributions/relations or neutral-veto track matches fails extraction.

`geometry_id`, `cellid_encoding`, `baseline_id` are compatibility assertions by the production owner, not autodetected conditions. New caches record them; reuse compares them. Unspecified identifiers produce QA warnings. Legacy v1 caches require all three, explicit compatibility acknowledgement and an expected complete manifest. Channel inventory CSV columns: cellID,x_mm,y_mm,z_mm. Decimal cellIDs are parsed directly as integers; missing xyz may be blank/NaN. An observed channel outside the inventory fails.

`input.identity_map` is an optional JSON object mapping physical input paths/PFNs to original canonical strings. Those identities are stored in caches and hashed. Existing compact train assignments are never regenerated. The legacy split is blake2b64 of canonical:entry:pi0-v1, train when digest mod5 !=0. blake2b-fraction supports other training fractions for new extraction.

`vertex`: use valid PrimaryVertices links, ndf>0, finite chi2>=0, configured xyz bounds; highest ndf then lowest chi2/ndf. Otherwise use nominal_mm only if allow_nominal. Extraction stores full xyz and vertex source. Changing vertex/split settings is incompatible with an existing cache.

`selection`: energy floor may not be below the stored cache floor. Neutral veto, asymmetry (null disables), overlap and seed test are applied before ranking the eligible pairs. Ranking is largest E1+E2, tie by original cluster indices. No mass or truth is used. All original contributions remain in the solver.

`peak`: configurable histogram range/binning, background form, identification criteria, and window sigma. The preserved optimizer searches mean0.075–0.195 GeV, sigma0.004–0.028 GeV using original initializations; identification interior bounds can be stricter but do not expand that search. `window_sigma` works for any positive width; global training fit occurs once.

`solver`: meaningful_fraction defines independent event support from either original shower. min_usage sets eligibility only; all channels remain in energy denominators. rcond truncates positive matrix eigenvalues below rcond*maximum. Convergence requires maximum multiplicative step <step_tolerance, absolute squared-mass bias <bias_tolerance, and squared-mass residual-spread change <bias_tolerance. Equations accumulate in chunks; no dense event-by-detector matrix is formed.

`qa`: halves use the full-training active set and deterministic permutation seed. Half convergence or rank failure disables stability acceptance. stability_warn/fail mark channels, projection_min marks rank-sensitive channels, low_support labels exploratory statistics. width_warning is a descriptive fractional increase in sigma/peak, not a significance test. Any reason produces physics WARNING; no-peak/numerical exceptions produce numerical FAIL. PASS does not approve constants. closure is optional algebraic 5% injection compared to nominal fitted gains. Detailed pair CSVs are off by default.

`qa.formats` controls png/pdf/svg; histogram_bins controls coefficient/support histograms. Mass counts use the configured peak binning/range for identical comparisons. Original broad counts, all gain values and resolved settings suffice to regenerate plots.

Normalization never enters selection or fitting. See normalization.md. All runtime outputs, cache paths, work paths and batch logs derive from the CLI output/cache/work arguments. Batch resources and launcher are specified separately in site YAML.

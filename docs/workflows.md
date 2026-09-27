# Workflow guide

Practical local, compact-cache and Slurm workflows for `eeemcal-pi0-calibration`.
Start with [installation and the synthetic example](../README.md). Commands
below run from the repository root after installation.

## Input manifest

Edit [examples/manifest.txt](../examples/manifest.txt) to list your RECO inputs:

```text
/path/to/reco/file_0001.eicrecon.edm4eic.root
/path/to/reco/file_0002.eicrecon.edm4eic.root
# root://your-xrootd-host//path/to/reco/file_0003.eicrecon.edm4eic.root
```

One file per line; blank lines and whole-line `#` comments are ignored.
Local relative paths resolve against the manifest directory; explicit
`root://` URLs are supported. Duplicate normalized entries are rejected.
Distinct aliases/copies of the same generated events need caller bookkeeping;
the producer does not scan central files to identify duplicates.

Required reconstructed information includes cluster energy/xyz, RecHit
cellID/energy/xyz, cluster-to-hit ObjectIDs and **assigned hit-energy
contributions**, plus event identity. Defaults are `EcalEndcapNClusters`,
`EcalEndcapNRecHits`, `EventHeader`, `EcalEndcapNTrackClusterMatches` and
`CentralCKFTracks` with their persisted relations. The default neutral veto
requires track-match information. A usable `PrimaryVertices` entry is preferred;
the configured nominal IP is the documented fallback. MC particles and
cluster/MC associations are optional. See [configuration](configuration.md)
for collection names, vertex rules and cache contracts. Schema compatibility
must be checked for a new production; familiar collection names alone do not
establish usable relations.

## Basic local workflow

After replacing the manifest placeholders, run:

```bash
pi0-calibrate \
  --input examples/manifest.txt \
  --config examples/calibration.yaml \
  --output results/run001 \
  --executor local
```

This extracts broad candidate shards, selects one eligible pair per event,
fits the training π0 peak, solves channel gains and produces held-out QA and
reports. Defaults place the broad cache in `results/run001/cache/` and selected
sample buffers in `results/run001/work/`. To place them separately, add:

```text
--cache-dir /shared/pi0/cache/run001 --work-dir /scratch/pi0/run001
```

Use new directories: the producer creates them with ownership records.
Declare your geometry, cellID encoding and existing calibration baseline in
the configuration; unspecified identifiers generate QA warnings. Optionally
supply a `cellID,x_mm,y_mm,z_mm` inventory CSV so unsupported crystals also
appear in the output. Otherwise the channel inventory is the observed cells.

## Reusing a compact cache

```bash
pi0-calibrate \
  --compact /shared/pi0/broad_compact \
  --config examples/calibration.yaml \
  --output results/selection_test \
  --executor local
```

This reads the existing broad cache without rereading RECO ROOT files or
modifying the input cache. Selection and support studies use new output roots.
`--input` and `--compact` are mutually exclusive; `--cache-dir` controls new
extraction storage, not the location of an existing compact input.

The cache must contain the needed candidates and relations: the requested
photon threshold cannot be below its extraction floor. Geometry, encoding,
baseline, split and vertex settings must be compatible. Existing training
assignments are preserved. Legacy caches additionally require explicit
compatibility declarations and a complete expected manifest. See
[configuration](configuration.md) and [compact format](../schemas/compact.schema.json).

## Changing calibration parameters

[examples/calibration.yaml](../examples/calibration.yaml) lists the complete
current defaults. Copy it to a new configuration file and edit the fields you
need. For example, these keys control photon energy, asymmetry, support and
rank truncation:

```yaml
selection:
  energy_GeV: 0.2
  asymmetry: 0.8
solver:
  min_usage: 20
  rcond: 0.001
```

A single YAML file overlays built-in defaults. Repeated `--set` overrides are
applied after that file, in order:

```bash
pi0-calibrate --compact /shared/pi0/broad_compact \
  --config examples/calibration.yaml --set solver.min_usage=20 \
  --output results/min20
```

[examples/min20.yaml](../examples/min20.yaml) is a minimal override of the built-in
defaults. It does not automatically load another configuration file. The
shipped defaults remain minimum support 10, meaningful contribution ≥5% of
either shower, `rcond=0.001` and 24 maximum iterations. Selection defaults are
neutral photons with E≥0.2 GeV, asymmetry<0.8, overlap<0.2 and different seeds;
there is no exactly-two-total-clusters requirement.

The default split is the original deterministic 80/20 `blake2b-mod5` partition.
For **new extraction**, another fraction requires both
`split.convention: blake2b-fraction` and `split.training_fraction`, for example
0.75. An existing cache cannot be repartitioned. Changing thresholds is an
explicit new study, not an improvement guaranteed by a larger active set.

## Output products

```text
results/run001/
├── gains.csv
├── summary.json
├── provenance.json
├── config_resolved.yaml
├── inputs.txt
├── report.md
├── plots/
├── diagnostics/
├── logs/
├── batch/
├── cache/                 # default for new extraction
└── work/                  # default selected-sample buffers
```

`gains.csv` contains exact decimal `cellID`, positive `relative_gain`,
`meaningful_support` (distinct training-window events), `active`, `half_change`,
`retained_projection`, detector xyz, a display `crystal_number`, and status
flags. Flags distinguish unsupported/unchanged, low support, instability,
rank sensitivity and numerical/QA failures. Unity for an unsupported channel
means **unchanged**, not calibrated. See [gain fields](../schemas/gains.schema.json).

Coefficients multiply the declared existing hit-energy baseline. They are
**candidates for review**, not automatically accepted conditions constants.
Training-half reproducibility and retained eigenmodes are not uncertainty
estimates. `summary.json` separates execution, numerical, physics-QA and
normalization status; `report.md` explains findings. Provenance records input
identity, settings and software versions; `diagnostics/` retains plot inputs.

Plots include training mass, before/after selection, held-out mass before/after
calibration, gain histogram, gain versus crystal number, 2D gain map, support
map/histogram, stability, convergence and eigenvalues. PNG/PDF are the defaults;
SVG is configurable. Before-selection plots count all eligible cached pairs,
so pair counts and distinct-event counts are reported separately. Rebuild plots
without reading ROOT/cache inputs or fitting again:

```bash
pi0-calibrate --stage plots --output results/run001
```

Exit 0 means completed or submitted, possibly with QA warnings; exit 2 denotes
configuration/input/execution failure; exit 3 denotes a diagnostic completion
without a usable peak/numerical solution. Exit 0 is not approval to apply gains.

## Simulation luminosity and equivalent beam time

Normalization is metadata, separate from gain fitting. For compatible,
confirmed unit-weight simulation:

```text
L_int = N_generated_represented / sigma_generator
t_equivalent = L_int / L_instantaneous
```

Reporting defaults are **10³³** and **10³⁴ cm⁻² s⁻¹**. These are ideal
luminosity-equivalent times, not accelerator/DAQ live time. The generated
numerator describes the input sample, not selected pairs or the training split.
Missing or ambiguous metadata produces unavailable/warning normalization
without preventing the calibration.

Set `normalization.metadata` to your own YAML/JSON record. The
[simulation](../examples/simulation_normalization.yaml) and
[data](../examples/data_normalization.yaml) examples contain illustrative values
that must be replaced. The optional explicit HepMC3 metadata reader is bounded;
it does not infer EVGEN mappings or normalize arbitrary weighted samples.
Details are in [normalization](normalization.md).

## Running on a Slurm farm

**One user-facing command can submit many jobs:**

```text
submission command
    ├── extraction array task 0 ──┐
    ├── extraction array task 1 ──┤
    └── ... extraction task N ────┤
                                 ↓
                       broad compact shards
                                 ↓
                  one dependent global calibration / QA job
```

Extraction is parallel over input files. **Shard gain maps are not separately
fitted and averaged.** One global analysis combines the compact showers with
a common channel mapping and accumulates C,D,L,B for each solve. That fit runs
on one node reading the compact shards. `--executor slurm` submits and returns;
no persistent foreground session is required.

Edit [examples/jlab.yaml](../examples/jlab.yaml) for your site before submitting:

```bash
pi0-calibrate \
  --input manifest.txt \
  --config examples/calibration.yaml \
  --output /shared/pi0/run001 \
  --executor slurm \
  --batch-config examples/jlab.yaml
```

The site YAML controls scheduling and the container launcher, separately from
physics settings. Its actual keys are:

```yaml
account: <ACCOUNT>
partition: <PARTITION>
qos: null
launcher: [/path/to/eic-shell, --]
python: /shared/pi0/venv/bin/python
files_per_shard: 25
files_per_task: 50
concurrency: 16
extraction: {cpus: 2, memory: 8G, time: '04:00:00'}
analysis: {cpus: 2, memory: 16G, time: '01:00:00'}
```

Replace literal placeholders with site values; use null for scheduler options
that should be omitted. `launcher` is an argv prefix that must accept the Python
command appended to it. Use `launcher: []` without a container. No YAML shell
variables are expanded. Use absolute shared paths for the interpreter and
launcher, or commands available on each compute node's PATH. If `python` is
null, it uses the submitting interpreter, which must also be accessible there.

The demonstrated settings of roughly **25 files/shard**, **50 files/task** and
**16 concurrent extraction tasks** are examples, not intrinsic limits.
`files_per_task` must be a positive multiple of `files_per_shard`; the last
shard/task can be smaller. BLAS/OpenMP thread counts follow allocated CPUs.
Output, cache, work buffers, interpreter and configuration payload must be
accessible across the farm. Local scratch is suitable for local/external-stage
jobs only when their inputs and final products are managed appropriately.

Locations derive from `--output`; no independent log/work-plan relocation
option exists:

| Product | Location under output |
|---|---|
| Combined Slurm stdout/stderr | `logs/extract-%A_%a.log`, `logs/analyze-%A_%a.log` |
| Resolved paths and common physics configuration | `batch/work_plan.json` |
| Worker scripts | `batch/extract.sh`, `batch/analyze.sh` |
| Task manifests | `batch/task_manifests/` |
| Submission IDs/commands | `batch/submissions.json` |
| Proposed submission commands | `batch/proposed_submissions.json` |

Compact-only farm execution uses the same command with `--compact` instead of
`--input`; it submits one analysis job and no extraction array.

## Batch dry run

```bash
pi0-calibrate --input manifest.txt --config examples/calibration.yaml \
  --output /shared/pi0/preview --executor slurm \
  --batch-config examples/jlab.yaml --dry-run
```

This creates the owned output/work/cache directories and writes the resolved
plan, worker scripts, task manifests and proposed commands. It reads the
manifest/configuration but does not open ROOT/cache payloads or submit jobs.
To submit that same planned run later, repeat the command without `--dry-run`
and add `--resume`. Do not treat a dry-run directory as an unrelated empty run.

## Batch status, resume and failed tasks

Query the recorded jobs once, then exit:

```bash
pi0-calibrate --stage status --output /shared/pi0/run001
```

This uses `squeue` for current jobs; it does not poll continuously or provide a
full `sacct` history. For completed work, inspect `summary.json`, `report.md`
and completion records. Retry/resume the same command with `--resume`:

```bash
pi0-calibrate --input manifest.txt --config examples/calibration.yaml \
  --output /shared/pi0/run001 --executor slurm \
  --batch-config examples/jlab.yaml --resume
```

Signature-compatible completed shards are reused. Recorded active tasks stay
running; missing/failed tasks are resubmitted, and workflow-owned pending
analysis dependencies are repaired. There is no separate retry subcommand.
The global fit refuses incomplete extraction; it is not allowed to silently
fit only successful jobs. A completed analysis is reused without refitting.

An existing output is refused without `--resume` or `--overwrite`.
`--overwrite` recomputes standard analysis results for the **same identity and
configuration**, preserving completed caches and unknown files; it refuses
queued/running workflow jobs. Changed inputs, physics or site settings require
a new output root. Output/cache/work cannot overlap protected input data.
Details and exact relative-path rules are in [configuration](configuration.md).

## Stages under an external production system

A site production manager can schedule stages itself; the Slurm adapter is
optional. Inside its allocation, use `--executor local` for each stage:

```bash
pi0-calibrate --stage extract --input task_manifest.txt \
  --config examples/calibration.yaml --output /shared/pi0/task001 \
  --cache-dir /shared/pi0/task001_cache --executor local

pi0-calibrate --stage analyze --compact /shared/pi0/combined_cache \
  --config examples/calibration.yaml --output /shared/pi0/global \
  --work-dir /shared/pi0/global_work --executor local
```

Extraction writes `candidates.npz`, its companion JSON, shard completion
records and a complete cache index. For multiple external tasks, assemble a
complete `index.json` in the combined-cache directory, listing paths to their
candidate shards:

```json
{"schema_version": 1, "complete": true, "shards": ["shards/shard_000/candidates.npz"]}
```

Paths resolve against the index directory; companion JSON files must also be
present. The external scheduler must establish that **all expected tasks**
finished before marking this index complete. Set `input.expected_manifest` to
the complete disjoint file set for a further completeness check. It is required
for legacy caches. Keep canonical event identities and compatible settings;
do not concatenate independently fitted coefficient maps.

## Real-data usage

The calibration is truth-independent. `input.mode: data` skips MC branches;
real data do not provide `GenCrossSection`. Recorded/delivered luminosity can
be supplied separately through the metadata interface, or remain unavailable.
A lack of luminosity does not prevent fitting. This interface has been tested
with truth-free synthetic fixtures; **real-data calibration has not yet been
validated**.


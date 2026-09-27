# Full compact regression: PASS

Compared against the saved detailed relaxed-B/min10 result. All event/pair IDs,
train flags, cellIDs, sparse structure, per-cell support and training-window
pair counts match exactly. Numerical tolerances and differences are recorded
in full_regression.json; no selection or solver retuning was performed.

Files: 11379; reconstructed events: 4967376.
Training-window pairs: 45385; active/retained: 2740/2740.
Median gain: 1.05040330855; half-stable below10%: 2678.
Held-out peak: 128.442359 → 134.646174 MeV;
relative width: 3.87069966% → 3.92835328%.

Only existing broad compact inputs and saved reference products were read.
This reproduces the fixed-shower calculation; it does not validate real data,
a full reconstruction feedback loop or detector systematic errors.
Normalization was supplied as existing approximate
regression metadata, independently of the gain calculation.

Recorded Slurm job: **11839626**, COMPLETED, ExitCode **0:0**.
The original numerical JSON is preserved; this public documentation cleanup
did not rerun the regression or alter its numerical results.

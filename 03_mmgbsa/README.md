# Stage 3 — MM-GBSA physics-based scoring

This is the longest and most elaborate stage of the pipeline. It takes the full-atom Boltz-2 predictions from Stage 2, runs each one through a molecular dynamics simulation, and computes a binding free energy using MM-GBSA via AmberTools' MMPBSA.py. The output is a physics-informed ΔG_bind for each candidate that can be compared against the ML-based rankings from ODesign and Boltz-2.

## Pipeline overview

The stage consists of four sub-stages chained together by SLURM dependencies:

1. **Amber parameter generation** (CPU, ~30 s) — `prepare_amber.py` converts the Boltz CIF into Amber topology and coordinate files using tLeap, with explicit handling of the head-to-tail peptide cyclization bond.

2. **Energy minimization** (GPU, ~1 min) — `minimize.py` relaxes the tLeap-built structure using OpenMM with GB-OBC2 implicit solvent. This resolves any geometric strain introduced by automated hydrogen placement around the cyclization.

3. **Molecular dynamics** (GPU, ~15 min) — `run_md.py` runs 100 ps of heating (50→300 K) followed by 10 ns of production at 300 K, saving a frame every 20 ps. Implicit solvent is used throughout for consistency with the scoring step.

4. **MM-GBSA scoring** (CPU, ~10 min) — `run_mmgbsa.py` invokes AmberTools' MMPBSA.py in single-trajectory mode to compute ΔG_bind from the 500-frame trajectory.

Each sub-stage has both a Python script that does the actual work and a SLURM wrapper (`run_<stage>.sh`) that manages the job submission. The wrappers are idempotent: if the expected outputs already exist, the job exits successfully without re-running. This means you can freely re-submit after a failure without worrying about wasted compute.

## Running the pipeline

You have three options depending on how much control you want.

### Option A — run everything for all designs at once (recommended)

This is the normal path for students:

```bash
bash submit_all_designs.sh \
    $SCRATCH/boltz/inputs/odesign_eval_manifest.tsv \
    $SCRATCH/mmgbsa
```

The script reads the Stage 2 manifest and submits the four-sub-stage pipeline for each of the 5 designs, with SLURM dependencies between sub-stages so that nothing runs until its predecessor succeeds. You will have 20 jobs queued total (4 sub-stages × 5 designs), though at any given time only a few will actually be running.

Each design requires about 20 minutes compute time with an expected total wall time: 30-90 minutes depending on queue load. The MMPBSA.py sub-stages run in parallel on CPU once all the MD jobs finish.

### Option B — run the pipeline for one design

If you want to test the pipeline on a single design before running all five, or if one design failed and you want to retry just that one:

```bash
bash submit_design.sh \
    $SCRATCH/boltz/outputs/odesign_eval/boltz_results_odesign_eval/predictions/design_01_seed123/design_01_seed123_model_0.cif \
    $SCRATCH/mmgbsa/design_01_seed123
```

This submits the four sub-stages as a dependency chain for that single design. The idempotence checks run at submission time, so if some sub-stages are already done, only the remaining ones will be submitted — no wasted jobs.

### Option C — run each sub-stage manually

If you want to inspect the output of each sub-stage before proceeding to the next, submit them one at a time:

```bash
sbatch run_prepare_amber.sh <boltz_cif> <design_dir>
# wait for it to finish, inspect outputs, then:
sbatch run_minimize.sh <design_dir>
# wait, inspect, then:
sbatch run_md.sh <design_dir>
# wait, inspect, then:
sbatch run_mmgbsa.sh <design_dir>
```

Each wrapper checks that its upstream inputs exist before running, so you cannot accidentally skip a sub-stage.

## Monitoring progress

Use the included status reporter to see which sub-stages of which designs are done, running, queued, or failed:

```bash
bash check_status.sh $SCRATCH/mmgbsa
```

Output looks like:

```
Design                  Prep    Min     MD      MMGB    DELTA TOTAL (kcal/mol)
----------------------  ------  ------  ------  ------  ----------------------
design_01_seed123        done    done    done    done       -37.79 +/- 0.34
design_02_seed123        done    done    done    RUN
design_03_seed123        done    done    done    RUN
design_04_seed123        done    done    RUN     PEND
design_05_seed123        done    done    PEND    PEND
```

You can also pass a single design directory instead of the root:

```bash
bash check_status.sh $SCRATCH/mmgbsa/design_01_seed123
```

## Outputs

Each per-design directory (`$SCRATCH/mmgbsa/<design_id>/`) accumulates outputs from all four sub-stages:

**Stage 1 (Amber prep):**
- `receptor.prmtop`, `receptor.inpcrd` — MDM2 alone
- `ligand.prmtop`, `ligand.inpcrd` — cyclic peptide alone (with cyclization bond)
- `complex.prmtop`, `complex.inpcrd` — full complex
- `complex.pdb` — complex structure for visualization
- `tleap.log` — tLeap run log

**Stage 2 (minimization):**
- `complex_min.pdb` — minimized structure for visualization
- `complex_min.rst7` — minimized coordinates (Amber restart format)
- `minimize.log` — pre/post energies, forces, and sanity check output

**Stage 3 (MD):**
- `production.nc` — 10 ns trajectory with 500 frames (NetCDF format)
- `production.rst7` — final coordinates
- `md.log` — per-segment temperature and energy report from heating, plus production state reporter output

**Stage 4 (MM-GBSA):**
- `mmgbsa.dat` — human-readable summary with average energy components
- `mmgbsa.csv` — per-frame energies for downstream analysis
- `mmgbsa.log` — run log from MMPBSA.py
- `mmgbsa.in` — the input file we generated for MMPBSA.py

The `mmgbsa.csv` file is consumed by `04_analysis/analyze_results.py` to generate the per-frame ΔG time trace plot.

## Key design decisions

**Implicit solvent throughout.** Both the MD and the MMPBSA.py scoring use the GB-OBC2 model (`igb=5` in Amber's terminology) with `mbondi2` radii and 0.15 M salt. This is faster than explicit solvent, consistent with the MM-GBSA scoring model (no impedance mismatch between simulation and analysis), and appropriate for a teaching exercise where the goal is relative comparison among candidates rather than publication-quality absolute numbers.

**Single-trajectory MM-GBSA.** We only simulate the complex, and MMPBSA.py slices each complex frame into receptor-only and ligand-only snapshots for the decomposition. Internal energies cancel exactly, which dramatically reduces noise compared to a dual-trajectory approach. The downside is that we ignore any conformational reorganization that the peptide or MDM2 undergo upon binding, but for a small series of related binders this is probably fine.

**No entropy term.** Normal-mode entropy calculations would multiply the MM-GBSA runtime by roughly 50× and rarely improve rankings for closely related compounds. The literature consensus is that MM-GBSA without entropy is the right default for relative ranking among similar binders.

**10 ns production trajectories.** This is short by MD standards but adequate for the comparison. Real publication-quality MM-GBSA on similar systems uses 100 ns to 1 μs. We're staying short to fit in a single SLURM job within the course wall-time budget and to produce a deliverable in a week. Students should acknowledge this limitation and reason about its implications (e.g., insufficient sampling of alternative binding poses).

## Troubleshooting specific to this stage

**tLeap reports "N3-C bond parameter missing"** — this means the sequence-based construction didn't apply. Check that the sequence in `prepare_amber.py` is using the internal residue template and that the explicit bond was created.

**OpenMM fails with a CUDA error** — usually means another job grabbed your GPU first. Re-submit the minimize or MD sub-stage.

**MMPBSA.py fails with `ModuleNotFoundError: numpy.compat`** — this is the module-load-order conflict documented in `docs/troubleshooting.md`. The fix is already baked into `run_mmgbsa.sh` (we call MMPBSA.py via its absolute path without loading the ambertools module in the conda shell), so this should not happen on a fresh checkout.

**Cyclic peptide design has ΔG near zero** — this is real, not a bug. The peptide may have dissociated from MDM2 during MD. Verify by loading the first and last frames of `production.nc` in PyMOL. A dissociated peptide is actually an informative negative result and should be discussed in your report.

For other issues, see `docs/troubleshooting.md` at the top level.

## What's next

Once all 5 designs have completed the full pipeline, move to `04_analysis/` to compare the three methods and produce the final report.

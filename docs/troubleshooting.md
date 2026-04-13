# Troubleshooting guide

This document collects errors you are likely to encounter at least once while running this pipeline, along with the cause and fix for each. If you hit something not listed here, the next places to look are the SLURM `.err` and `.out` files from the failing job, and the stage-specific log file (`tleap.log`, `minimize.log`, `md.log`, `mmgbsa.log`) in the per-design directory.

## Setup and environment issues

### `setup_openmm.sh` fails with "conda not found"

You need to load the `conda` module first:

```bash
ml conda
bash 00_setup/setup_openmm.sh
```

The setup script does `ml conda` internally, but some older Scholar configurations require it to be loaded in the parent shell as well.

### `setup_openmm.sh` runs but CUDA platform is not detected

The last step of setup verifies that OpenMM can see the CUDA platform by listing available platforms. If it reports only `['Reference', 'CPU']` without CUDA, the GPU will not be accelerated but the environment still works for CPU runs. Most Stage 3 sub-stages are fast enough on CPU to still complete within the SLURM time limits, though MD will take substantially longer. To get CUDA support, re-run the setup on an A40 node — the verification runs on whatever node the script is launched from, and login nodes sometimes don't have the CUDA libraries visible.

### OpenMM environment already exists and I want to rebuild it

```bash
rm -rf $SCRATCH/conda/envs/openmm
bash 00_setup/setup_openmm.sh
```

The setup script refuses to overwrite an existing environment to prevent accidents.

## Stage 1 (ODesign) issues

### ODesign container fails with "fastfold_layer_norm_cuda JIT compile error"

The SIF file should have this extension pre-compiled during container build. If you get this error, the container was rebuilt without pre-compilation. Workaround: provide a writable torch extensions cache directory via bind mount. Earlier versions of `run_odesign.sh` did this; see git history if you need the older code.

### ODesign output is empty or CIF files are not in `predictions/` subdirectories

ODesign's output path convention includes timestamps, so the specific directory nesting can vary. The downstream `prepare_boltz_inputs.py` script handles this by walking the whole tree and finding any `predictions/*.cif` files. If no CIFs are found, check the ODesign stderr log — the run probably crashed.

### "gcc version is too new" or CUDA compiler errors in ODesign

The container runs with `--cleanenv` precisely to isolate it from Scholar's host compilers. If you see this error, `--cleanenv` was removed from the SLURM script — restore it.

## Stage 2 (Boltz-2) issues

### Boltz fails with "Unable to parse filetype .tsv"

Boltz processes every file in its input directory and tried to parse the manifest TSV as a prediction input. The manifest should live in the parent directory, not alongside the YAMLs. This is fixed in the current version of `prepare_boltz_inputs.py` — it writes the manifest to `<parent>/<dirname>_manifest.tsv` rather than `<dirname>/manifest.tsv`. If you hit this error with the current code, the fix was applied but you have a stale manifest from an older run; delete it manually:

```bash
rm $SCRATCH/boltz/inputs/odesign_eval/manifest.tsv
```

### Boltz fails with OOM (out of memory) during structure prediction

The DataLoader workers can be memory-hungry. The SLURM script sets `--cpus-per-task=8 --mem=96G` and passes `--num_workers 0` to boltz to avoid this, but if you modified those settings you may see OOM kills. Restore the original values.

### Boltz reports very low ipTM values for all designs

This is unlikely but possible if the MSA lookup for MDM2 fails. Boltz uses the mmseqs2 server to find homologs for the target protein; if the server is unreachable, Boltz falls back to no-MSA mode which produces lower-quality predictions. Check the Boltz log for MSA warnings.

### `conda activate boltz2` fails

The boltz2 environment is separate from the openmm environment used for Stages 3-4. Setting it up is outside the scope of this document — see the Boltz-2 documentation or ask your instructor.

## Stage 3 (MM-GBSA) issues

### tLeap reports "Could not find bond parameter for atom types: N3 - C"

This was the symptom of the early cyclization bug. It means the first residue of the cyclic peptide was assigned an N-terminal patch (`N3` atom type) instead of being treated as an internal residue. The fix is in the current `prepare_amber.py`: we use `sequence { ... }` followed by `loadpdbusingseq` to force internal residue templates, then add an explicit bond. If you hit this error, the fix was not applied — verify that `prepare_amber.py` has the `write_tleap_script` function using `sequence { }` construction.

### tLeap reports "The value must be of the type: Atom"

This was the symptom of an earlier attempted fix (`set pep head none`). The keyword `none` is not a valid Atom reference and tLeap silently ignored the command, then processed the peptide as linear. The current fix avoids the `set head/tail` commands entirely.

### OpenMM minimization: "max force after minimization is very large"

The minimization did not fully converge. Usually means the starting geometry has clashing atoms that could not be resolved within `max_iterations`. Try increasing the iteration limit:

```bash
python minimize.py --input-dir $SCRATCH/mmgbsa/design_XX --max-iterations 20000
```

If it still doesn't converge, there is probably a topology problem — check `tleap.log` for warnings about close contacts or missing parameters.

### OpenMM MD fails with "particle coordinate is NaN"

Some atom ended up in an unphysical position during MD, typically because of a force explosion from bad starting geometry. The usual cause is that minimization didn't run first (or was skipped because a stale `complex_min.rst7` was sitting around). Solution: delete the stale minimization output and re-run both stages:

```bash
rm $SCRATCH/mmgbsa/design_XX/complex_min.*
bash submit_design.sh <cif> $SCRATCH/mmgbsa/design_XX
```

### MMPBSA.py fails with `ModuleNotFoundError: No module named 'numpy.compat'`

This is the classic module-load-order conflict. When you `ml ambertools/25` in a shell that also has the openmm conda env active, AmberTools' bundled Python (with an outdated parmed) ends up on the path and shadows the conda env's modern parmed. The fix is to not load ambertools/25 in the same shell as the conda env, and instead call MMPBSA.py via its absolute path. The current `run_mmgbsa.sh` already does this — if you're hitting this error, you may have modified the script to add `ml ambertools/25`. Remove it.

### MMPBSA.py fails with "Could not import Amber Python modules"

MMPBSA.py's internal imports need `AMBERHOME` set so it can find its `MMPBSA_mods` package. The current `run_mmgbsa.py` sets `AMBERHOME` and `LD_LIBRARY_PATH` in the subprocess environment before calling MMPBSA.py. If you hit this error, either those environment variable settings were lost, or the AmberTools installation path has changed from `/apps/external/ambertools/25/ambertools25/`.

### `run_mmgbsa.py` fails with `No module named 'netCDF4'`

The earlier version used netCDF4 for frame counting. The current version uses `parmed.amber.NetCDFTraj` instead, which is already in the openmm environment. If you hit this error, update `run_mmgbsa.py` from git.

### Status script shows "FAIL" for a stage that is actually running

The status script consults SLURM's queue before checking file-based state, but if the SLURM snapshot is stale (for example, because squeue momentarily returned empty), you may briefly see a stage marked FAIL when it is really still running. Just re-run `check_status.sh` a few seconds later.

### A design reports ΔG_bind near zero

This is almost certainly not a bug — the peptide actually dissociated from MDM2 during the MD simulation. Verify by loading the last frame of `production.nc` in PyMOL (see the "Visualizing trajectories" note in `03_mmgbsa/README.md`) and checking whether the peptide is still in contact with the protein. A dissociated peptide is an informative negative result and should be discussed in your report rather than suppressed.

### Design_XX's prep_amber stage fails with "Expected exactly 2 chains, found 1"

The Boltz-2 CIF has both the target and the peptide on the same chain ID. This is unusual but can happen if Boltz processed the input with an unexpected chain assignment. Inspect the CIF file and look at the `label_asym_id` column — there should be two distinct values, one for MDM2 and one for the peptide. If there isn't, regenerate that design's Boltz input with a fresh YAML.

## Stage 4 (Analysis) issues

### "Per-frame plot is empty"

The script could not parse `mmgbsa.csv`. MMPBSA.py's CSV format has evolved across AmberTools versions. Check that the CSV has a section labeled "DELTA Energy Terms" followed by a header row containing "TOTAL". If the format is different, the parser in `analyze_results.py` needs updating.

### Scatter plot shows no labels

The script strips `_seed\d+` suffixes from design IDs for cleaner labels. If your designs are named differently, the labels may be hidden or misplaced. Edit the regex in `plot_scatter_iptm_vs_ddg()` if needed.

### Correlations all show "None" in the markdown table

Happens when either Boltz-2 confidence or MM-GBSA values are missing for most designs. Run `bash ../03_mmgbsa/check_status.sh $SCRATCH/mmgbsa` and confirm everything completed.

### matplotlib import error

The openmm conda environment doesn't include matplotlib by default. Install it once:

```bash
ml conda
source activate $SCRATCH/conda/envs/openmm
conda install -y -c conda-forge matplotlib
```

The current `setup_openmm.sh` includes matplotlib in the default package list, so this should only happen on environments created with an older version of the setup script.

## General SLURM issues

### My jobs are stuck in "PD" (pending) for a long time

Check `squeue -u $USER --start` to see the estimated start time. The `gpu` queue on Scholar can have long waits during peak hours. If the estimated start is unreasonably far out, you may be hitting a QoS or account limit. Check with `sacctmgr show assoc user=$USER` for your limits.

### A dependency-chained job was cancelled with reason "DependencyNeverSatisfied"

An upstream job in the chain failed, which caused SLURM to cancel everything downstream. This is working as intended — check the upstream job's `.err` file to find out why it failed, fix the issue, and re-run `submit_design.sh` or `submit_all_designs.sh` (which will skip already-done stages).

### "Submitted batch job" but the job never shows up in squeue

Rarely, SLURM accepts a submission but then silently drops it. Run `sacct -j <jobid>` to check the job's state. If it's COMPLETED in sacct but never appeared in squeue, the job ran and finished before you checked. If it's CANCELLED with no other jobs of yours running, something is wrong with your account — contact Scholar support.

### Email notifications are not arriving

The SLURM scripts use `--mail-user=bsdrown@purdue.edu` as a placeholder from the original development. If you haven't changed it, emails are going to the wrong person. Change the lines in each `run_*.sh` script to `--mail-user=${USER}@purdue.edu`.

# Stage 2 — Boltz-2 structure refinement

This stage takes the backbone-only ODesign outputs from Stage 1 and re-predicts each candidate's full-atom structure using [Boltz-2](https://github.com/jwohlwend/boltz), an open-source AlphaFold3-style co-folding model. Boltz-2 adds the side chain atoms that ODesign didn't place, gives us an interface confidence score (ipTM), and provides a properly-structured input for the physics-based Stage 3 pipeline.

## Why re-predict instead of using ODesign's output directly

ODesign is a backbone-only generative model. Its output CIF files contain only Cα, N, C, and O atoms — no side chains, no hydrogens. The MM-GBSA pipeline in Stage 3 cannot work without full atom coordinates. We could in principle ask a lightweight tool like Rosetta to pack side chains onto the ODesign backbone, but Boltz-2 is a better choice for two reasons: it produces a calibrated structural confidence score (ipTM) that gives us a second opinion on binding quality, and its joint prediction of backbone and side chains is likely to be self-consistent in a way that side chain packing onto a fixed backbone is not.

Boltz-2 also has an affinity prediction head, but we deliberately do not use it here. The Boltz-2 authors have stated that the affinity head is trained on drug-like small molecules (<50 heavy atoms) and is not reliable for peptidic binders with ~100 heavy atoms. We rely on MM-GBSA in Stage 3 for binding energy estimates instead.

## Inputs

This stage expects ODesign output CIF files produced by Stage 1 at `$SCRATCH/odesign/outputs/mdm2_cyclic/`. No manual input files are required.

## Running the stage

Stage 2 has two steps: preparing Boltz-2 YAML inputs from the ODesign CIFs, and then running Boltz-2 itself.

### Step 2a — Generate Boltz YAML inputs (login node, ~5 seconds)

```bash
python prepare_boltz_inputs.py \
    --odesign-output $SCRATCH/odesign/outputs/mdm2_cyclic \
    --output-dir $SCRATCH/boltz/inputs/odesign_eval \
    --n-designs 5 \
    --min-length 7
```

This script walks the ODesign output tree, picks the 5 most recent prediction CIF files, extracts each peptide's sequence, and writes one Boltz-2 input YAML per design. Each YAML specifies the MDM2 target sequence plus the cyclic peptide with `cyclic: true` to enable head-to-tail cyclization in the prediction.

The `--min-length` argument skips designs shorter than the threshold (default 7 residues) for synthetic accessibility reasons. The script walks through candidates in order until it has collected `--n-designs` that pass the filter, or until it runs out of ODesign outputs.

A manifest TSV file is written to `$SCRATCH/boltz/inputs/odesign_eval_manifest.tsv`. This manifest lists each selected design along with its source CIF, peptide sequence, and the path to the generated Boltz YAML. It is used by Stage 3 (`submit_all_designs.sh`) and Stage 4 (`analyze_results.py`) to track which designs belong to the evaluation set.

**Note**: The manifest is deliberately placed outside the YAML input directory. Boltz-2 processes every file in its input directory and would fail trying to parse the TSV as a prediction input if we placed it alongside the YAMLs.

### Step 2b — Run Boltz-2 (GPU job, ~20-30 minutes)

```bash
sbatch run_boltz_evaluation.sh $SCRATCH/boltz/inputs/odesign_eval
```

This submits a SLURM job that runs Boltz-2 on all 5 YAMLs in a single batch. Boltz-2 generates 5 diffusion samples per YAML (so 25 total structure predictions across the 5 designs), runs an MSA search for MDM2 via the mmseqs2 server, and produces CIF files ranked by ML confidence.

The job requests an A40 GPU with 96 GB of RAM. The memory is needed because Boltz-2's DataLoader can be memory-hungry during inference — we've set `--num_workers 0` to avoid DataLoader workers eating memory but still need generous RAM for the main process.

You will receive an email when the job completes.

## Outputs

Boltz-2 writes results to `$SCRATCH/boltz/outputs/odesign_eval/boltz_results_odesign_eval/predictions/<design_id>/`. For each design, you will find:

- `<design_id>_model_0.cif` — top-ranked predicted structure (used in Stage 3)
- `<design_id>_model_1.cif` through `_model_4.cif` — alternative samples
- `confidence_<design_id>_model_0.json` — confidence metrics (ipTM, pTM, pLDDT)
- `pae_<design_id>_model_0.json` — predicted aligned error matrix
- `plddt_<design_id>_model_0.json` — per-residue pLDDT

The `model_0` files are the ones Stage 3 consumes. Stage 4 reads `confidence_<design_id>_model_0.json` to get the ipTM values for the comparison table.

## Verifying the output

Three quick checks before moving on:

1. **All 5 designs have output directories**. If any are missing, Boltz-2 may have failed on that design specifically — check the SLURM error file.

2. **ipTM values look reasonable**. Open one of the confidence JSON files and check that `iptm` is a number between 0 and 1 (higher is better). Typical values for successfully predicted cyclic peptide-protein complexes are 0.7-0.95. Very low values (<0.5) suggest Boltz-2 is uncertain about the prediction.

3. **Visual spot-check**. Load one of the `model_0.cif` files in PyMOL and confirm the peptide is docked in the MDM2 cleft. Side chains should be present (this is the main reason we ran Boltz-2), and the peptide should look cyclic.

## What's next

Once Stage 2 completes, move to `03_mmgbsa/` to run the physics-based evaluation pipeline on the refined structures.

# Cyclic Peptide Binder Evaluation Pipeline

A multi-method computational pipeline for generating and evaluating cyclic peptide binders against a protein target, developed for CHM696 Graduate Chemical Biology at Purdue University. This pipeline uses MDM2 as the example target and compares three different scoring approaches: generative-model ranking (ODesign), ML-based structure prediction with confidence scores (Boltz-2), and physics-based binding free energy from molecular dynamics (MM-GBSA).

The core pedagogical goal is to let students see firsthand how these three methods — each representing a different philosophical approach to the binding problem — agree and disagree on a common set of candidate peptides, and to develop critical thinking about when each method can and cannot be trusted.

## What this pipeline does

Given a target protein (MDM2 in the example) and a set of binding-site hotspot residues, the pipeline:

1. Generates candidate cyclic peptide binders using ODesign, a flexible-receptor generative model
2. Predicts full-atom structures of the top candidates in complex with the target using Boltz-2, and extracts interface confidence scores (ipTM)
3. Runs short (10 ns) molecular dynamics simulations of each complex in implicit solvent using OpenMM
4. Computes binding free energies from the trajectories using MM-GBSA via AmberTools' MMPBSA.py
5. Compares the three methods side by side in a summary report with plots

The complete pipeline runs in roughly two hours of wall time on Purdue's Scholar HPC cluster, most of which is spent waiting in queues rather than computing.

## Prerequisites

You need an account on Scholar with access to the `gpu` partition. The pipeline has been tested with:

- Scholar's `ambertools/25` module (for tLeap, MMPBSA.py, pdb4amber)
- Scholar's `conda` module (for creating the OpenMM environment)
- An A40 GPU (`--constraint=J`) for Boltz-2 and OpenMM jobs
- The ODesign container at `/class/bsdrown/apps/odesign/odesign.sif`

You do not need to install anything manually. The one-time setup script in `setup/` creates a conda environment with OpenMM and its dependencies; everything else is loaded via Scholar modules or already present on the cluster.

## Quick start

Clone the repository into your home directory on Scholar:

```bash
cd ~
git clone https://github.com/bdrown/chm696-cyclic-peptide-design.git
cd chm696-cyclic-peptide-design
```

All pipeline outputs land under your `$SCRATCH` space in conventional subdirectories (`$SCRATCH/odesign`, `$SCRATCH/boltz`, `$SCRATCH/mmgbsa`, `$SCRATCH/results`). You don't need to create these manually; each stage creates its own directories as needed.

Create the two conda environments (one-time, total about 15-20 minutes). See `setup/README.md` for details on what each environment contains and why they are kept separate:

```bash
bash 00_setup/setup_boltz.sh
bash 00_setup/setup_openmm.sh
```

Then walk through the four stages in order. Each stage has its own README with detailed instructions, but the high-level sequence is:

```bash
# Stage 1 — ODesign generation (GPU, ~15 min)
cd 01_odesign
sbatch run_odesign.sh inputs/mdm2_cyclic.json
cd ..

# Stage 2 — Boltz-2 structure refinement (GPU, ~30 min)
cd 02_boltz_refinement
python prepare_boltz_inputs.py \
    --odesign-output $SCRATCH/odesign/outputs/mdm2_cyclic \
    --output-dir $SCRATCH/boltz/inputs/odesign_eval
sbatch run_boltz_evaluation.sh $SCRATCH/boltz/inputs/odesign_eval
cd ..

# Stage 3 — MM-GBSA pipeline (GPU + CPU, ~1 hr)
cd 03_mmgbsa
bash submit_all_designs.sh \
    $SCRATCH/boltz/inputs/odesign_eval_manifest.tsv \
    $SCRATCH/mmgbsa
bash check_status.sh $SCRATCH/mmgbsa  # monitor progress
cd ..

# Stage 4 — Comparison analysis (CPU, ~1 min)
cd 04_analysis
python analyze_results.py \
    --manifest $SCRATCH/boltz/inputs/odesign_eval_manifest.tsv \
    --mmgbsa-root $SCRATCH/mmgbsa \
    --boltz-root $SCRATCH/boltz/outputs/odesign_eval/boltz_results_odesign_eval/predictions \
    --output-dir $SCRATCH/results
cd ..
```

Your final report and figures will be in `$SCRATCH/results/`, including a markdown file (`results.md`) and figures that you can convert to answer questions.

## Repository layout

```
cyclic-peptide-mmgbsa/
├── README.md                  # this file
├── 00_setup/                     # one-time conda environment creation
├── 01_odesign/                # cyclic peptide generation
├── 02_boltz_refinement/       # ML-based structure + confidence
├── 03_mmgbsa/                 # physics-based binding energy
├── 04_analysis/               # cross-method comparison and plots
├── docs/
│   ├── pipeline_overview.md   # detailed explanation of what each stage does
│   ├── problem_set.md         # student-facing questions for the assignment
│   └── troubleshooting.md     # common errors and how to fix them
└── example_results/           # reference outputs for comparison
```

Each stage directory contains its own README with stage-specific details, inputs, expected outputs, and runtime estimates. Start with the stage READMEs when you're running a particular step; come back to this top-level document when you need to understand how stages fit together.

## What to do when things go wrong

Three resources, in order of how often you'll use them:

1. **`bash 03_mmgbsa/check_status.sh $SCRATCH/mmgbsa`** — the status reporter tells you which stages of the MM-GBSA pipeline have completed, which are running, which are queued, and which have failed. This is almost always the first thing to consult.

2. **`docs/troubleshooting.md`** — a list of errors you will likely encounter at least once and the fix for each. Topics include SLURM module conflicts, tLeap cyclization errors, OpenMM CUDA failures, and MMPBSA.py import issues.

3. **SLURM job output files** (`*-JOBID-$USER.err` and `*-JOBID-$USER.out`) — these are created in whichever directory you ran `sbatch` from and contain the raw output of each job. Look here if the troubleshooting guide doesn't cover your issue.

## Customizing for a different target

The default pipeline uses MDM2 with the 1YCR-based construct as the target. To run against a different protein:

1. Update the target sequence in `01_odesign/inputs/mdm2_cyclic.json`
2. Update the `MDM2_SEQUENCE` constant in `02_boltz_refinement/prepare_boltz_inputs.py`

Nothing else in the pipeline is target-specific; the MM-GBSA stage and all plotting is generic.

## Citation and acknowledgements

If you use this pipeline in a publication or course derivative, please cite the underlying tools:

- **ODesign**: Zhang et al., *ODesign: A World Model for Biomolecular Interaction Design*, 2025
- **Boltz-2**: Wohlwend et al., *Boltz-2: A biomolecular foundation model*, 2025
- **AmberTools**: Case et al., *The Amber biomolecular simulation programs*, J. Comput. Chem., 2005
- **OpenMM**: Eastman et al., *OpenMM 8*, J. Phys. Chem. B, 2023

The course materials and pipeline organization were developed by Prof. Bryon Drown for CHM696 - Chemical Biology at Purdue University.

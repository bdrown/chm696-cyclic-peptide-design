# Stage 1 — ODesign cyclic peptide generation

This stage uses [ODesign](https://odesign1.github.io/) to generate candidate cyclic peptide binders targeted at a specified binding site on a protein. ODesign is a flexible-receptor diffusion-based generative model that can design cyclic peptides by jointly sampling sequence, backbone geometry, and binder-receptor pose.

## Inputs

**`1ycr.pdb`** - Structure of MDM2 in complex with transactivation domain of p53

**`mdm2_cyclic.json`** — ODesign input specification. The file defines:

- The target protein (MDM2 residues 17-125, matching the 1YCR crystal structure construct)
- The cyclic peptide to be designed (length range 8-12 residues, `if_cyc: true`)
- Hotspot residues on MDM2 defining the binding cleft (the p53 binding site with key residues Leu54, Gly58, Ile61, Met62, Tyr67, Gln72, His73, Val75, Phe91, Val93, His96, Ile99)
- The `hotspot_center` specification that tells ODesign to focus sampling near these residues

Change these files if you want to design against a different target or change the peptide length range. The valid length range is tool-dependent; for reference, 8-12 residues is a typical starting point for synthetically accessible head-to-tail cyclic peptides.

## Running the stage

### Copy input files to Scratch

The ODesign script expects that all input files are located on scratch in the `$SCRATCH/odesign/inputs/` directory. Create this directory and copy files over:

```bash
mkdir -p $SCRATCH/odesign/inputs/
cp mdm2_cyclic.json $SCRATCH/odesign/inputs/
cp 1ycr.pdb $SCRATCH/odesign/inputs/
```

### Run ODesign (GPU queue, 10 minutes)

The ODesign SLURM wrapper takes the input JSON filename (relative to the `inputs/` directory) as its sole argument:

```bash
sbatch run_odesign.sh mdm2_cyclic.json
```

This submits a GPU job that runs the ODesign container on Scholar. Expected wall time is 10 minutes on an A40 GPU. The job runs with 3 seeds × 5 samples per seed = 15 total candidates generated. You can adjust the `SEEDS` and `N_SAMPLE` variables in the SLURM script if you want more or fewer candidates. While the required runtime is reasonably fast, there are limited compute nodes equipped with A40 GPU's, so you may need to wait in line depending on time of day.


## Outputs

ODesign writes its results to `$SCRATCH/odesign/outputs/<input_name>/` where `<input_name>` is the basename of your input JSON without the extension (so `mdm2_cyclic.json` produces `mdm2_cyclic/`). The output directory structure is:

```
$SCRATCH/odesign/outputs/mdm2_cyclic/
└── <timestamp>/
    └── <sample_name>/
        └── seed_<seed>/
            └── predictions/
                ├── design_01.cif
                ├── design_02.cif
                └── ...
```

Each `.cif` file contains a predicted cyclic peptide in complex with MDM2. The peptide's backbone coordinates and sequence are specified, but side chain atoms are not placed — ODesign only predicts Cα/N/C/O backbone positions. Adding side chains is the job of Stage 2 (Boltz-2 refinement).

The ODesign output uses its own CIF format, which has standard atom_site records but does not include bond records for the head-to-tail cyclization. The downstream stages (both `prepare_boltz_inputs.py` in Stage 2 and `prepare_amber.py` in Stage 3) handle cyclization explicitly rather than reading it from the CIF.

## Verifying the output

Load one of the generated CIF files in PyMOL to spot-check that the peptide is actually docked into the MDM2 cleft and that its length is in the expected range:

```
load $SCRATCH/odesign/outputs/mdm2_cyclic/<timestamp>/.../design_01.cif
```

If the peptide is floating away from MDM2 or has an unusual backbone geometry, this is a failed ODesign run — try different seeds or adjust the input JSON.

## Infrastructure notes

The ODesign container sits at `/class/bsdrown/apps/odesign/odesign.sif` on Scholar. The SLURM script binds the container's checkpoint and data directories read-only and its output directory writable. It uses `--cleanenv` and `--writable-tmpfs` to isolate the container from the host environment, which is necessary because host compiler versions conflict with the CUDA 12.1 toolchain baked into the container.

## What's next

Once Stage 1 completes, move to `02_boltz_refinement/` to re-evaluate the top candidates with Boltz-2 and get full-atom structures with side chains.

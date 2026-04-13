# One-time environment setup

This directory contains two setup scripts that create the conda environments the pipeline depends on. You run these exactly once, before starting any of the pipeline stages. Both scripts are idempotent — if the environment already exists they refuse to overwrite it.

## Why two environments

The pipeline depends on two independent software stacks that don't coexist peacefully in a single conda environment:

**`boltz2`** — used in Stage 2. Contains Boltz-2 and its PyTorch + cuEquivariance stack, built against CUDA 12.6 to satisfy both torch and cuEquivariance version constraints.

**`openmm`** — used in Stages 3 and 4. Contains OpenMM (built against CUDA 12) for molecular dynamics, plus parmed, mdtraj, numpy, scipy, and matplotlib for topology handling and analysis.

Keeping these separate means the PyTorch CUDA stack used by Boltz doesn't collide with OpenMM's CUDA stack, and updating one doesn't risk breaking the other. Scripts in each stage load the environment they need. You never need both active at the same time.

Stage 3 (MM-GBSA) also uses AmberTools for tLeap, pdb4amber, and MMPBSA.py, but we deliberately do *not* put AmberTools in a conda environment. Scholar already has a well-tested `ambertools/25` module that we load on demand. Conda-installing AmberTools would duplicate it and has been observed to cause numpy compatibility problems because AmberTools' bundled Python shadows the conda env's Python when both are active simultaneously. If you're curious about the gory details see `docs/troubleshooting.md`.

## Prerequisites

Before running either setup script:

- You need a Scholar account with access to your `$SCRATCH` space
- You need the `conda` module available (`ml conda` should succeed)
- You should run these on a login node, not inside a SLURM job
- Your home directory should have some free space — the scripts redirect caches to scratch, but conda itself puts a small amount of state in `$HOME/.conda`

## Setting up Boltz-2

From the top-level repository directory:

```bash
bash 00_setup/setup_boltz.sh
```

The script takes 10-15 minutes to download and install everything. It does the following:

1. Cleans up any stale pip cache from your home directory (earlier failed installs often leave megabytes of broken cache files behind)
2. Redirects the conda package cache, pip cache, and XDG cache to `$SCRATCH/conda/` to preserve your home directory quota
3. Creates the environment at `$SCRATCH/conda/envs/boltz2` using Python 3.11
4. Installs PyTorch built for CUDA 12.6 from PyTorch's custom index
5. Installs Boltz with CUDA extras and cuEquivariance kernels from NVIDIA's package index
6. Writes activate/deactivate hooks so the cache redirects persist whenever the environment is activated
7. Verifies the installation by checking that the installed torch has CUDA 12.6 and that `boltz predict --help` works
8. Creates `$SCRATCH/boltz/inputs` and `$SCRATCH/boltz/outputs` convenience directories

**Why CUDA 12.6 specifically.** Boltz-2 needs cuEquivariance for fast CUDA kernels, and cuEquivariance's CUDA ops require `nvidia-cublas-cu12 >= 12.5`. PyTorch's default cu124 build ships cublas 12.4, which is too old. PyTorch's cu126 build ships cublas 12.6+, which satisfies both. Scholar's A40 nodes run driver 575.57.08 which exposes CUDA 12.9, so cu126 is comfortably within the driver's compatibility range. Don't attempt to use PyTorch's default cu124 build or you'll hit a cryptic runtime error about cublas versions.

**First-run model download.** Boltz downloads its model weights (~5 GB) the first time it runs a prediction, not during setup. This happens inside your SLURM job on Stage 2, which is why `run_boltz_evaluation.sh` has a generous 1-hour wall time — most of that is headroom for the first-run download, and subsequent runs are much faster because the weights are cached in `$SCRATCH/conda/boltz-cache`.

## Setting up OpenMM

```bash
bash 00_setup/setup_openmm.sh
```

This is faster, about 5 minutes. It creates `$SCRATCH/conda/envs/openmm` with:

- OpenMM 8+ with CUDA 12 support
- parmed (reading/writing Amber topology files)
- mdtraj (trajectory analysis)
- numpy, scipy (general numerics)
- matplotlib (Stage 4 analysis plots)

At the end of the script the verification step prints the OpenMM version and lists available platforms. You should see `['Reference', 'CPU', 'CUDA', 'OpenCL']`, confirming CUDA acceleration is available. If CUDA is missing from the list, the environment still works but Stage 3 will run slowly on CPU. To fix this, re-run the setup on a compute node that has CUDA libraries visible.

## Activating the environments

Activation is the same pattern for both environments:

```bash
# For Boltz (Stage 2)
ml conda
source activate $SCRATCH/conda/envs/boltz2

# For OpenMM (Stages 3 and 4)
ml conda
source activate $SCRATCH/conda/envs/openmm
```

All SLURM wrappers in the pipeline handle activation automatically. You only need to activate manually if you want to run one of the Python scripts interactively for debugging.

Never activate both environments simultaneously in the same shell. Their Python interpreters and torch builds will conflict in subtle ways.

## Adding packages later

If you need additional packages for custom analysis, install them into the existing environment rather than creating a new one:

```bash
ml conda
source activate $SCRATCH/conda/envs/openmm
conda install -y -c conda-forge <package>
# or
pip install <package>
```

The same applies to the `boltz2` environment if you need additional dependencies for custom inference work.

## Rebuilding a broken environment

Both setup scripts refuse to overwrite an existing environment. If one gets into a broken state and you want to rebuild:

```bash
rm -rf $SCRATCH/conda/envs/boltz2
bash 00_setup/setup_boltz.sh
```

Or for OpenMM:

```bash
rm -rf $SCRATCH/conda/envs/openmm
bash 00_setup/setup_openmm.sh
```

## Troubleshooting

Most setup problems fall into three categories:

**Home directory quota exceeded** — conda's default package cache can overflow a small home quota. Both setup scripts redirect the cache to `$SCRATCH/conda/pkgs` to avoid this, but if you've already started a failed attempt with the default location, clean up with `rm -rf ~/.conda ~/.cache/pip` and try again.

**Network timeouts fetching packages** — the Boltz setup in particular has to download a lot of PyTorch + CUDA wheels from multiple package indices. If the network flakes, just re-run the script. It will resume from cached packages where possible.

**CUDA version mismatch after Boltz install** — if the verification step at the end of `setup_boltz.sh` reports a CUDA version other than 12.6, something in the dependency resolution pulled a different torch build. This usually means pip resolved a package that brought in an incompatible torch as a dependency. The script prints the exact fix command; run it to force-install the cu126 build.

For other problems, see `docs/troubleshooting.md` at the top level of the repository.

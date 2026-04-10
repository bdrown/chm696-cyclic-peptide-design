#!/usr/bin/env python3
from __future__ import annotations

"""
run_mmgbsa.py

Run MM-GBSA analysis on a cyclic peptide-MDM2 complex trajectory using
AmberTools' MMPBSA.py.

Strategy: single-trajectory approach. The complex trajectory is the
only input MD trajectory; receptor-only and ligand-only snapshots are
extracted by MMPBSA.py from the complex frames. This means internal
strain energies (bond, angle, torsion, internal vdW/electrostatic)
cancel exactly between bound and unbound states, leaving only the
interaction free energy components.

Inputs (in --input-dir):
    complex.prmtop, receptor.prmtop, ligand.prmtop  (from prepare_amber.py)
    production.nc                                    (from run_md.py)
Outputs:
    mmgbsa.in            -- MMPBSA.py input file
    mmgbsa.dat           -- per-frame and average energies
    mmgbsa.csv           -- machine-readable summary
    mmgbsa.log           -- run log

Usage:
    ml conda
    ml ambertools/25
    source activate $SCRATCH/conda/envs/openmm
    python run_mmgbsa.py --input-dir $SCRATCH/mmgbsa/design_01

The script expects both ambertools/25 (for MMPBSA.py and ante-MMPBSA.py)
and the openmm conda env (for parmed) to be loaded simultaneously.
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import parmed


# MMPBSA.py input template. Settings:
#   igb=5: GB-OBC2 model, matches our tLeap mbondi2 radii and OpenMM OBC2.
#   saltcon=0.15: 0.15 M salt for Debye-Huckel correction. Matches MD.
#   surften=0.0072: standard nonpolar surface tension (kcal/mol/A^2).
#   interval=1: process every frame from the trajectory. Adjust if too slow.
#   keep_files=0: clean up intermediate snapshot PDBs/inpcrds after run.
#   netcdf=1: input trajectory is NetCDF (Amber format).
#
# We deliberately do NOT request entropy (no nmode= block). Normal-mode
# entropy calculations are extremely expensive and rarely improve
# rankings for closely related compounds. The literature consensus is
# that MM-GBSA without entropy is the right choice for relative
# ranking among similar binders.
MMPBSA_INPUT_TEMPLATE = """MMPBSA input file for cyclic peptide -- MDM2 binding
&general
   startframe = 1,
   endframe   = {endframe},
   interval   = {interval},
   verbose    = 2,
   keep_files = 0,
   netcdf     = 1,
/
&gb
   igb        = 5,
   saltcon    = 0.15,
   surften    = 0.0072,
/
"""


def count_frames(traj_path: Path, prmtop_path: Path) -> int:
    """Count frames in the NetCDF trajectory using parmed."""
    from parmed.amber import NetCDFTraj
    traj = NetCDFTraj.open_old(str(traj_path))
    return traj.frame


def get_residue_ranges(complex_prmtop: Path,
                       receptor_prmtop: Path) -> tuple:
    """
    Determine the residue index ranges for receptor and ligand within
    the complex topology.

    The complex was built by tLeap with `combine { rec pep }`, so
    receptor residues come first followed by peptide residues.
    Returns (n_receptor_residues, n_total_residues).
    """
    rec = parmed.load_file(str(receptor_prmtop))
    comp = parmed.load_file(str(complex_prmtop))
    n_rec = len(rec.residues)
    n_total = len(comp.residues)
    return n_rec, n_total


def write_mmpbsa_input(input_path: Path,
                       n_frames: int,
                       interval: int = 1) -> None:
    """Write the MMPBSA.py input file."""
    content = MMPBSA_INPUT_TEMPLATE.format(
        endframe=n_frames,
        interval=interval,
    )
    input_path.write_text(content)


def run_mmpbsa(input_dir: Path,
               input_file: Path,
               complex_prmtop: Path,
               receptor_prmtop: Path,
               ligand_prmtop: Path,
               trajectory: Path,
               log_path: Path,
               mmpbsa_exe: str = "MMPBSA.py") -> None:
    """
    Invoke MMPBSA.py with the prepared inputs. Runs in input_dir so
    that intermediate files land there rather than wherever the user
    happened to be when they submitted.
    """
    cmd = [
        mmpbsa_exe,
        "-O",                              # overwrite outputs
        "-i", str(input_file),
        "-o", "mmgbsa.dat",
        "-eo", "mmgbsa.csv",               # per-frame data as CSV
        "-cp", str(complex_prmtop),
        "-rp", str(receptor_prmtop),
        "-lp", str(ligand_prmtop),
        "-y", str(trajectory),
    ]

    print(f"Running MMPBSA.py in {input_dir}", file=sys.stderr)
    print(f"Command: {' '.join(cmd)}", file=sys.stderr)

    # Build a subprocess environment so MMPBSA.py finds its bundled
    # Python modules. AMBERHOME tells it where its package tree lives;
    # we also have to put MMPBSA_mods on PYTHONPATH because MMPBSA.py
    # imports it as a top-level package rather than relative to amber.
    amber_bin = Path(mmpbsa_exe).parent
    amber_home = amber_bin.parent
    sub_env = os.environ.copy()
    sub_env["PATH"] = str(amber_bin) + os.pathsep + sub_env.get("PATH", "")
    sub_env["AMBERHOME"] = str(amber_home)
    sub_env["LD_LIBRARY_PATH"] = (
        str(amber_home / "lib") + os.pathsep
        + sub_env.get("LD_LIBRARY_PATH", "")
    )
    # MMPBSA.py expects MMPBSA_mods, parmed, etc. in AmberTools' own
    # site-packages, not our conda env's. Point PYTHONPATH there and
    # clear PYTHONHOME so the bundled Python uses its own stdlib.
    amber_site = amber_home / "lib" / "python3.12" / "site-packages"
    sub_env["PYTHONPATH"] = str(amber_site)
    sub_env.pop("PYTHONHOME", None)

    with log_path.open("w") as logf:
        logf.write(f"Command: {' '.join(cmd)}\n")
        logf.write(f"Working directory: {input_dir}\n\n")
        logf.flush()

        result = subprocess.run(
            cmd, cwd=str(input_dir), env=sub_env,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        )
        logf.write(result.stdout)

    if result.returncode != 0:
        print(f"MMPBSA.py failed with exit code {result.returncode}",
              file=sys.stderr)
        print(f"See log at {log_path}", file=sys.stderr)
        # Print last few lines of log to stderr for quick diagnosis
        print("\n--- Last 30 lines of log ---", file=sys.stderr)
        for line in result.stdout.splitlines()[-30:]:
            print(line, file=sys.stderr)
        raise RuntimeError(f"MMPBSA.py failed (exit {result.returncode})")


def parse_mmgbsa_output(dat_path: Path) -> dict:
    """
    Parse the MMPBSA.py text output to extract the binding energy
    components and final DELTA TOTAL.

    Returns a dict with keys for each energy term (vdW, EEL, EGB, ESURF,
    GGAS, GSOLV, TOTAL) plus their standard errors. All in kcal/mol.
    """
    if not dat_path.is_file():
        return {}

    text = dat_path.read_text()

    # MMPBSA.py output has a "DELTA Energy Decomposition" or
    # "Differences (Complex - Receptor - Ligand)" section near the end.
    # We look for lines like:
    #   VDWAALS         -52.3456         5.1234         0.7245
    # where columns are Average, StdDev, StdErrofMean.
    result = {}
    in_delta_section = False
    for line in text.splitlines():
        line = line.rstrip()
        if "Differences (Complex - Receptor - Ligand)" in line:
            in_delta_section = True
            continue
        if not in_delta_section:
            continue

        parts = line.split()
        if len(parts) < 4:
            continue

        # Try to parse "NAME  avg  std  sem" pattern
        name = parts[0]
        if name in ("BOND", "ANGLE", "DIHED", "VDWAALS", "EEL", "1-4",
                    "EGB", "ESURF", "GGAS", "GSOLV", "TOTAL", "DELTA"):
            try:
                avg = float(parts[-3])
                std = float(parts[-2])
                sem = float(parts[-1])
                result[name] = {"avg": avg, "std": std, "sem": sem}
            except (ValueError, IndexError):
                continue

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Run MM-GBSA on a cyclic peptide-MDM2 complex"
    )
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--interval", type=int, default=1,
                        help="Frame interval (1 = every frame)")
    parser.add_argument("--start-frame", type=int, default=1,
                        help="First frame to process (1-indexed, default: 1)")
    args = parser.parse_args()

    input_dir = args.input_dir
    complex_prmtop = input_dir / "complex.prmtop"
    receptor_prmtop = input_dir / "receptor.prmtop"
    ligand_prmtop = input_dir / "ligand.prmtop"
    trajectory = input_dir / "production.nc"

    for f in (complex_prmtop, receptor_prmtop, ligand_prmtop, trajectory):
        if not f.is_file():
            print(f"Error: required input not found: {f}", file=sys.stderr)
            sys.exit(1)

    # Resolve MMPBSA.py executable. Prefer the MMPBSA_PY env var (set by
    # the SLURM wrapper to an absolute path) so we avoid module-loading
    # ambertools/25 in the same shell as the conda env, which would
    # shadow our modern parmed with the AmberTools-bundled outdated one.
    mmpbsa_exe = os.environ.get("MMPBSA_PY") or shutil.which("MMPBSA.py")
    if not mmpbsa_exe or not os.path.isfile(mmpbsa_exe):
        print("Error: MMPBSA.py not found. Set MMPBSA_PY to its absolute "
              "path or load ambertools/25 in a non-conflicting shell.",
              file=sys.stderr)
        sys.exit(1)

    print(f"=== MM-GBSA Analysis ===", file=sys.stderr)
    print(f"Input directory: {input_dir}", file=sys.stderr)

    # Topology consistency check
    n_rec, n_total = get_residue_ranges(complex_prmtop, receptor_prmtop)
    n_lig = n_total - n_rec
    print(f"Receptor: {n_rec} residues", file=sys.stderr)
    print(f"Ligand:   {n_lig} residues", file=sys.stderr)
    print(f"Complex:  {n_total} residues", file=sys.stderr)

    # Frame count
    n_frames = count_frames(trajectory, complex_prmtop)
    print(f"Trajectory: {n_frames} frames", file=sys.stderr)
    n_processed = (n_frames - args.start_frame + 1) // args.interval
    print(f"Will process {n_processed} frames "
          f"(start={args.start_frame}, interval={args.interval})",
          file=sys.stderr)
    print("", file=sys.stderr)

    # Write MMPBSA input file
    input_file = input_dir / "mmgbsa.in"
    write_mmpbsa_input(input_file, n_frames=n_frames,
                       interval=args.interval)
    print(f"Wrote {input_file}", file=sys.stderr)

    # Run MMPBSA.py
    log_path = input_dir / "mmgbsa.log"
    run_mmpbsa(input_dir, input_file,
               complex_prmtop, receptor_prmtop, ligand_prmtop,
               trajectory, log_path, mmpbsa_exe=mmpbsa_exe)

    # Parse and report results
    dat_path = input_dir / "mmgbsa.dat"
    results = parse_mmgbsa_output(dat_path)
    print("", file=sys.stderr)
    print(f"=== MM-GBSA Results ({input_dir.name}) ===", file=sys.stderr)
    if not results:
        print("(could not parse results -- check mmgbsa.dat)",
              file=sys.stderr)
    else:
        print(f"{'Component':<12} {'Average':>12} {'StdDev':>10} "
              f"{'SEM':>10}  (kcal/mol)", file=sys.stderr)
        for term in ("VDWAALS", "EEL", "EGB", "ESURF",
                     "GGAS", "GSOLV", "DELTA TOTAL", "TOTAL"):
            if term in results:
                r = results[term]
                print(f"{term:<12} {r['avg']:>12.3f} {r['std']:>10.3f} "
                      f"{r['sem']:>10.3f}", file=sys.stderr)

    print("", file=sys.stderr)
    print(f"Full output:    {dat_path}", file=sys.stderr)
    print(f"Per-frame CSV:  {input_dir}/mmgbsa.csv", file=sys.stderr)
    print(f"Run log:        {log_path}", file=sys.stderr)


if __name__ == "__main__":
    main()

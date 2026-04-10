#!/usr/bin/env python3
from __future__ import annotations

"""
minimize.py

Energy-minimize an Amber-prepared MDM2-cyclic peptide complex using
OpenMM with the GB-OBC implicit solvent model (igb=5 equivalent).

The minimization serves two purposes:
  1. Relax any geometric strain from tLeap's automated hydrogen
     placement around the cyclization bond (the phantom clashes
     that PyMOL inferred are resolved here).
  2. Produce a sensible starting point for subsequent MD equilibration
     and production runs.

Usage:
    ml conda
    source activate $SCRATCH/conda/envs/openmm
    python minimize.py --input-dir $SCRATCH/mmgbsa/design_01

Inputs expected in --input-dir:
    complex.prmtop, complex.inpcrd
Outputs written to --input-dir:
    complex_min.pdb    -- minimized coordinates as PDB
    complex_min.rst7   -- minimized coordinates as Amber restart
    minimize.log       -- energy report before and after
"""

import argparse
import sys
import time
from pathlib import Path

import openmm as mm
import openmm.app as app
import openmm.unit as unit
import parmed


def minimize_complex(prmtop_path: Path,
                     inpcrd_path: Path,
                     output_dir: Path,
                     max_iterations: int = 5000,
                     tolerance_kj_mol_nm: float = 10.0,
                     use_gpu: bool = True) -> None:
    """
    Load an Amber topology/coordinate pair, set up GB-OBC implicit
    solvent, minimize energy, and write out the minimized structure.

    Parameters
    ----------
    prmtop_path : Path
        Amber parameter/topology file.
    inpcrd_path : Path
        Amber coordinate file.
    output_dir : Path
        Where to write minimized PDB, restart file, and log.
    max_iterations : int
        Maximum L-BFGS iterations. 0 means run until convergence.
    tolerance_kj_mol_nm : float
        Convergence criterion: RMS force in kJ/mol/nm. OpenMM default
        is 10 which corresponds to ~0.24 kcal/mol/Angstrom, a standard
        minimization tolerance.
    use_gpu : bool
        If True, use CUDA platform. Falls back to CPU if unavailable.
    """
    log_path = output_dir / "minimize.log"
    log_lines = []

    def log(msg: str) -> None:
        print(msg, file=sys.stderr)
        log_lines.append(msg)

    log(f"=== OpenMM minimization ===")
    log(f"Topology:    {prmtop_path}")
    log(f"Coordinates: {inpcrd_path}")
    log(f"Output dir:  {output_dir}")
    log("")

    # Load Amber files via parmed for a clean OpenMM handoff. parmed
    # reads the prmtop exactly as Amber understands it, which is safer
    # than OpenMM's native AmberPrmtopFile loader for unusual topologies
    # like our cyclic peptide.
    log("Loading Amber topology with parmed...")
    structure = parmed.load_file(str(prmtop_path), xyz=str(inpcrd_path))
    log(f"  {len(structure.atoms)} atoms, {len(structure.residues)} residues")
    log(f"  {len(structure.bonds)} bonds")

    # Check topology for the cyclization bond. Look for bonds between
    # the first and last residues of any chain.
    cyclic_bonds = []
    for bond in structure.bonds:
        r1 = bond.atom1.residue
        r2 = bond.atom2.residue
        if abs(r1.idx - r2.idx) > 1 and bond.atom1.name in ("N", "C") \
                and bond.atom2.name in ("N", "C"):
            cyclic_bonds.append((r1.idx, r2.idx, bond.atom1.name,
                                 bond.atom2.name))
    if cyclic_bonds:
        log(f"  Found {len(cyclic_bonds)} non-sequential N-C bond(s) "
            f"(likely cyclization):")
        for r1, r2, a1, a2 in cyclic_bonds:
            log(f"    residue {r1+1} {a1} -- residue {r2+1} {a2}")
    else:
        log("  WARNING: no cyclization bond detected. Expected one for "
            "a cyclic peptide system.")

    # Set up GB-OBC implicit solvent (matches Amber igb=5 / mbondi2).
    # HBonds constraints are standard for protein MD and allow the
    # integrator to use a longer timestep during MD (later), but they
    # also slightly accelerate minimization by removing stiff H bond
    # degrees of freedom.
    log("")
    log("Creating OpenMM system with GB-OBC implicit solvent...")
    system = structure.createSystem(
        implicitSolvent=app.OBC2,
        implicitSolventSaltConc=0.15 * unit.molar,
        constraints=app.HBonds,
        nonbondedMethod=app.NoCutoff,
        removeCMMotion=False,
    )
    log(f"  System has {system.getNumParticles()} particles, "
        f"{system.getNumForces()} force terms")

    # Dummy integrator — minimization doesn't actually integrate,
    # but OpenMM requires one to create a Context.
    integrator = mm.LangevinMiddleIntegrator(
        300 * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )

    # Platform selection
    if use_gpu:
        try:
            platform = mm.Platform.getPlatformByName("CUDA")
            properties = {"Precision": "mixed"}
            log("Using CUDA platform with mixed precision")
        except Exception as e:
            log(f"CUDA platform unavailable ({e}), falling back to CPU")
            platform = mm.Platform.getPlatformByName("CPU")
            properties = {}
    else:
        platform = mm.Platform.getPlatformByName("CPU")
        properties = {}
        log("Using CPU platform")

    simulation = app.Simulation(
        structure.topology, system, integrator, platform, properties,
    )
    simulation.context.setPositions(structure.positions)

    # Pre-minimization energy and forces
    log("")
    state_before = simulation.context.getState(
        getEnergy=True, getForces=True,
    )
    e_before = state_before.getPotentialEnergy().value_in_unit(
        unit.kilocalorie_per_mole
    )
    forces = state_before.getForces(asNumpy=True).value_in_unit(
        unit.kilocalorie_per_mole / unit.angstrom
    )
    max_force_before = float((forces ** 2).sum(axis=1).max() ** 0.5)
    log(f"Pre-minimization:")
    log(f"  Potential energy: {e_before:15.3f} kcal/mol")
    log(f"  Max force:        {max_force_before:15.3f} kcal/mol/A")

    # Run minimization
    log("")
    log(f"Minimizing (max {max_iterations} iterations, "
        f"tolerance {tolerance_kj_mol_nm} kJ/mol/nm)...")
    t_start = time.time()
    simulation.minimizeEnergy(
        tolerance=tolerance_kj_mol_nm * unit.kilojoule_per_mole / unit.nanometer,
        maxIterations=max_iterations,
    )
    t_elapsed = time.time() - t_start
    log(f"Minimization finished in {t_elapsed:.1f} seconds")

    # Post-minimization energy and forces
    state_after = simulation.context.getState(
        getPositions=True, getEnergy=True, getForces=True,
    )
    e_after = state_after.getPotentialEnergy().value_in_unit(
        unit.kilocalorie_per_mole
    )
    forces = state_after.getForces(asNumpy=True).value_in_unit(
        unit.kilocalorie_per_mole / unit.angstrom
    )
    max_force_after = float((forces ** 2).sum(axis=1).max() ** 0.5)
    log("")
    log(f"Post-minimization:")
    log(f"  Potential energy: {e_after:15.3f} kcal/mol")
    log(f"  Max force:        {max_force_after:15.3f} kcal/mol/A")
    log(f"  Energy change:    {e_after - e_before:+15.3f} kcal/mol")

    # Sanity check: energy should have decreased
    if e_after > e_before + 1.0:
        log("")
        log("WARNING: potential energy did not decrease. This usually "
            "indicates a topology error or extreme starting geometry.")
    # Sanity check: max force should be reasonable
    if max_force_after > 100.0:
        log("")
        log(f"WARNING: max force after minimization is {max_force_after:.1f} "
            f"kcal/mol/A, which is unusually large. Consider increasing "
            f"max_iterations or checking for topology errors.")

    # Write minimized coordinates to both PDB (for visualization) and
    # Amber restart format (for downstream MD / MMPBSA.py).
    log("")
    log("Writing minimized structure...")
    positions = state_after.getPositions()

    pdb_out = output_dir / "complex_min.pdb"
    with pdb_out.open("w") as f:
        app.PDBFile.writeFile(structure.topology, positions, f,
                              keepIds=True)
    log(f"  {pdb_out}")

    # Update parmed structure with minimized coordinates and save as
    # Amber restart. parmed handles unit conversion automatically.
    structure.positions = positions
    rst_out = output_dir / "complex_min.rst7"
    structure.save(str(rst_out), format="rst7", overwrite=True)
    log(f"  {rst_out}")

    log("")
    log("=== Minimization complete ===")

    log_path.write_text("\n".join(log_lines) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Minimize an Amber-prepared complex with OpenMM",
    )
    parser.add_argument(
        "--input-dir", type=Path, required=True,
        help="Directory containing complex.prmtop and complex.inpcrd",
    )
    parser.add_argument(
        "--max-iterations", type=int, default=5000,
        help="Maximum minimization iterations (default: 5000)",
    )
    parser.add_argument(
        "--tolerance", type=float, default=10.0,
        help="Convergence tolerance in kJ/mol/nm (default: 10.0)",
    )
    parser.add_argument(
        "--cpu", action="store_true",
        help="Force CPU platform (default: try CUDA first)",
    )
    args = parser.parse_args()

    prmtop = args.input_dir / "complex.prmtop"
    inpcrd = args.input_dir / "complex.inpcrd"

    if not prmtop.is_file():
        print(f"Error: {prmtop} not found", file=sys.stderr)
        sys.exit(1)
    if not inpcrd.is_file():
        print(f"Error: {inpcrd} not found", file=sys.stderr)
        sys.exit(1)

    try:
        minimize_complex(
            prmtop, inpcrd, args.input_dir,
            max_iterations=args.max_iterations,
            tolerance_kj_mol_nm=args.tolerance,
            use_gpu=not args.cpu,
        )
    except Exception as e:
        print(f"ERROR: minimization failed: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

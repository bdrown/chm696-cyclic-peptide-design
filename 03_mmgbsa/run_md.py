#!/usr/bin/env python3
from __future__ import annotations

"""
run_md.py

Heat and produce a short MD trajectory for an Amber-prepared
MDM2-cyclic peptide complex using OpenMM with GB-OBC implicit solvent.

Protocol
--------
1. Load the minimized complex (complex_min.rst7).
2. Heating: 100 ps from 50 K to 300 K using Langevin dynamics,
   2 fs timestep, HBonds constraints.
3. Production: 500 ps at 300 K, saving a frame every 10 ps.
   Output: 50 frames in NetCDF format for downstream MM-GBSA analysis.

The trajectory is written in Amber NetCDF format (.nc) rather than
DCD because MMPBSA.py natively reads NetCDF and the conversion
avoids a separate step.

Usage:
    ml conda
    source activate $SCRATCH/conda/envs/openmm
    python run_md.py --input-dir $SCRATCH/mmgbsa/design_01

Inputs expected in --input-dir:
    complex.prmtop, complex_min.rst7
Outputs written to --input-dir:
    production.nc        -- NetCDF trajectory, 50 frames
    production.rst7      -- final coordinates (Amber restart)
    md.log               -- energy/temperature report
"""

import argparse
import sys
import time
from pathlib import Path

import openmm as mm
import openmm.app as app
import openmm.unit as unit
import parmed


def setup_system(prmtop_path: Path, rst7_path: Path):
    """
    Load the Amber topology and restart, build an OpenMM system with
    GB-OBC2 implicit solvent. Returns (structure, system).
    """
    structure = parmed.load_file(str(prmtop_path), xyz=str(rst7_path))
    system = structure.createSystem(
        implicitSolvent=app.OBC2,
        implicitSolventSaltConc=0.15 * unit.molar,
        constraints=app.HBonds,
        nonbondedMethod=app.NoCutoff,
        removeCMMotion=True,
    )
    return structure, system


def get_platform(use_gpu: bool):
    """Select CUDA if available, else CPU."""
    if use_gpu:
        try:
            platform = mm.Platform.getPlatformByName("CUDA")
            return platform, {"Precision": "mixed"}
        except Exception:
            pass
    return mm.Platform.getPlatformByName("CPU"), {}


def heat(simulation: app.Simulation,
         total_ps: float,
         start_temp_k: float,
         end_temp_k: float,
         log_interval_ps: float,
         log_fn) -> None:
    """
    Linearly ramp the Langevin target temperature from start to end
    over total_ps picoseconds.
    """
    timestep_ps = simulation.integrator.getStepSize().value_in_unit(
        unit.picoseconds
    )
    total_steps = int(total_ps / timestep_ps)
    log_steps = max(1, int(log_interval_ps / timestep_ps))
    n_segments = 20  # update temperature in discrete steps
    steps_per_segment = max(1, total_steps // n_segments)

    log_fn(f"Heating: {start_temp_k:.0f} -> {end_temp_k:.0f} K over "
           f"{total_ps:.0f} ps ({total_steps} steps in {n_segments} segments)")

    t_start_wall = time.time()
    for seg in range(n_segments):
        frac = (seg + 1) / n_segments
        target_t = start_temp_k + frac * (end_temp_k - start_temp_k)
        simulation.integrator.setTemperature(target_t * unit.kelvin)
        simulation.step(steps_per_segment)

        state = simulation.context.getState(getEnergy=True)
        ke = state.getKineticEnergy().value_in_unit(unit.kilocalorie_per_mole)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        # Instantaneous temperature from KE
        # KE = (3/2) N kT  => T = 2 KE / (3 N k)
        n_dof = 3 * simulation.system.getNumParticles() \
            - simulation.system.getNumConstraints()
        kT = (2.0 * ke * 4.184 * 1000.0) / (n_dof * 8.314462618)  # K
        log_fn(f"  segment {seg+1:2d}/{n_segments} "
               f"target={target_t:6.1f} K  "
               f"T={kT:6.1f} K  PE={pe:10.1f}  KE={ke:8.1f} kcal/mol")

    t_wall = time.time() - t_start_wall
    log_fn(f"Heating completed in {t_wall:.1f} s "
           f"({total_steps * timestep_ps / t_wall:.1f} ps/s)")


def produce(simulation: app.Simulation,
            structure: parmed.Structure,
            total_ps: float,
            save_interval_ps: float,
            output_nc: Path,
            log_interval_ps: float,
            log_fn) -> None:
    """
    Run constant-temperature production MD and write a NetCDF
    trajectory at the specified frame interval.
    """
    timestep_ps = simulation.integrator.getStepSize().value_in_unit(
        unit.picoseconds
    )
    total_steps = int(total_ps / timestep_ps)
    save_steps = max(1, int(save_interval_ps / timestep_ps))
    log_steps = max(1, int(log_interval_ps / timestep_ps))
    n_frames = total_steps // save_steps

    log_fn(f"Production: {total_ps:.0f} ps at 300 K "
           f"({total_steps} steps, {n_frames} frames saved every "
           f"{save_interval_ps:.0f} ps)")
    log_fn(f"Trajectory: {output_nc}")

    # Use parmed's NetCDFReporter for Amber-compatible NetCDF output.
    # This is what MMPBSA.py expects natively.
    nc_reporter = parmed.openmm.NetCDFReporter(
        str(output_nc), save_steps, crds=True, vels=False, frcs=False,
    )
    simulation.reporters.append(nc_reporter)

    # Also attach a state data reporter for per-log-interval monitoring.
    state_reporter = app.StateDataReporter(
        sys.stdout, log_steps,
        step=True, time=True,
        potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
        temperature=True, speed=True, remainingTime=True,
        totalSteps=total_steps, separator="  ",
    )
    simulation.reporters.append(state_reporter)

    t_start_wall = time.time()
    simulation.step(total_steps)
    t_wall = time.time() - t_start_wall

    log_fn(f"Production completed in {t_wall:.1f} s "
           f"({total_steps * timestep_ps / t_wall:.2f} ps/s)")

    # Remove reporters so subsequent calls don't double-report
    simulation.reporters.remove(nc_reporter)
    simulation.reporters.remove(state_reporter)


def main():
    parser = argparse.ArgumentParser(
        description="Heat and produce MD for an Amber-prepared complex",
    )
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--heat-ps", type=float, default=100.0)
    parser.add_argument("--prod-ps", type=float, default=10000.0)
    parser.add_argument("--save-interval-ps", type=float, default=20.0)
    parser.add_argument("--timestep-fs", type=float, default=2.0)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--cpu", action="store_true")
    args = parser.parse_args()

    prmtop = args.input_dir / "complex.prmtop"
    rst7_min = args.input_dir / "complex_min.rst7"

    if not prmtop.is_file():
        print(f"Error: {prmtop} not found", file=sys.stderr)
        sys.exit(1)
    if not rst7_min.is_file():
        print(f"Error: {rst7_min} not found. Run minimize.py first.",
              file=sys.stderr)
        sys.exit(1)

    log_path = args.input_dir / "md.log"
    log_lines = []

    def log(msg: str) -> None:
        print(msg, file=sys.stderr)
        log_lines.append(msg)

    log("=== OpenMM MD (heating + production) ===")
    log(f"Input dir:     {args.input_dir}")
    log(f"Topology:      {prmtop}")
    log(f"Start coords:  {rst7_min}")
    log(f"Heat:          {args.heat_ps} ps (50 -> {args.temperature_k} K)")
    log(f"Production:    {args.prod_ps} ps @ {args.temperature_k} K")
    log(f"Save interval: {args.save_interval_ps} ps "
        f"({int(args.prod_ps / args.save_interval_ps)} frames)")
    log(f"Timestep:      {args.timestep_fs} fs")
    log("")

    try:
        # Set up system
        log("Loading system...")
        structure, system = setup_system(prmtop, rst7_min)
        log(f"  {len(structure.atoms)} atoms, "
            f"{len(structure.residues)} residues")

        # Langevin integrator
        integrator = mm.LangevinMiddleIntegrator(
            50.0 * unit.kelvin,  # Will be updated during heating
            1.0 / unit.picosecond,
            args.timestep_fs * unit.femtoseconds,
        )

        platform, properties = get_platform(not args.cpu)
        log(f"Using platform: {platform.getName()}")

        simulation = app.Simulation(
            structure.topology, system, integrator, platform, properties,
        )
        simulation.context.setPositions(structure.positions)
        simulation.context.setVelocitiesToTemperature(50 * unit.kelvin)

        # Heating
        log("")
        heat(simulation,
             total_ps=args.heat_ps,
             start_temp_k=50.0,
             end_temp_k=args.temperature_k,
             log_interval_ps=10.0,
             log_fn=log)

        # Production
        log("")
        produce(simulation, structure,
                total_ps=args.prod_ps,
                save_interval_ps=args.save_interval_ps,
                output_nc=args.input_dir / "production.nc",
                log_interval_ps=50.0,
                log_fn=log)

        # Save final restart
        log("")
        log("Saving final coordinates...")
        final_state = simulation.context.getState(
            getPositions=True, getVelocities=True,
        )
        structure.positions = final_state.getPositions()
        structure.velocities = final_state.getVelocities(asNumpy=True)
        final_rst = args.input_dir / "production.rst7"
        structure.save(str(final_rst), format="rst7", overwrite=True)
        log(f"  {final_rst}")

        log("")
        log("=== MD complete ===")

    except Exception as e:
        log("")
        log(f"ERROR: {e}")
        import traceback
        log(traceback.format_exc())
        log_path.write_text("\n".join(log_lines) + "\n")
        sys.exit(1)

    log_path.write_text("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    main()

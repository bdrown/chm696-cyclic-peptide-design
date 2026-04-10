#!/usr/bin/env python3
from __future__ import annotations

"""
prepare_amber.py

Convert a Boltz-2 CIF output into Amber topology/coordinate files for
MM-GBSA analysis. Generates three parameter sets per design:
  - complex.prmtop / complex.inpcrd (MDM2 + cyclic peptide)
  - receptor.prmtop / receptor.inpcrd (MDM2 alone)
  - ligand.prmtop / ligand.inpcrd (cyclic peptide alone)

Crucially, tLeap does not automatically detect head-to-tail peptide
cyclization. This script handles it explicitly by:
  1. Identifying the peptide chain by length
  2. Stripping terminal hydrogens from the peptide in the input PDB
  3. Generating a tLeap script that uses 'set pep head none / tail none'
     followed by an explicit 'bond' command to close the cycle

The resulting topology treats every peptide residue as an internal
(non-terminal) amino acid.

Usage:
    module load ambertools/25
    python prepare_amber.py \\
        --boltz-cif /path/to/design_01_model_0.cif \\
        --output-dir /path/to/amber/design_01

Compatible with Python 3.9+.
"""

import argparse
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple

# Three-letter to one-letter and back
AA_3LETTER = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # Alternative histidine forms Amber uses
    "HID", "HIE", "HIP",
}


def parse_cif_atoms(cif_path: Path) -> List[dict]:
    """
    Minimal CIF parser that returns a list of atom records as dicts.
    Each record has: group_PDB, atom_id, atom_name, resname, chain_id,
    seq_id, x, y, z, element.
    """
    atoms = []
    in_atom_site = False
    column_names: List[str] = []

    with cif_path.open() as f:
        for line in f:
            line = line.rstrip()

            if line.startswith("loop_"):
                column_names = []
                in_atom_site = False
                continue

            if line.startswith("_atom_site."):
                column_names.append(line.split(".", 1)[1].strip())
                in_atom_site = True
                continue

            if not in_atom_site:
                continue

            if line.startswith("#") or not line.strip():
                in_atom_site = False
                continue

            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                in_atom_site = False
                continue

            fields = line.split()
            if len(fields) < len(column_names):
                continue

            row = dict(zip(column_names, fields))
            try:
                atoms.append({
                    "group_PDB": row.get("group_PDB", "ATOM"),
                    "atom_id": int(row.get("id", "0")),
                    "atom_name": row.get("label_atom_id", "").strip('"'),
                    "resname": row.get("label_comp_id", "").upper(),
                    "chain_id": row.get("label_asym_id", ""),
                    "seq_id": int(row.get("label_seq_id", "0")),
                    "x": float(row.get("Cartn_x", "0")),
                    "y": float(row.get("Cartn_y", "0")),
                    "z": float(row.get("Cartn_z", "0")),
                    "element": row.get("type_symbol", "").strip(),
                })
            except (ValueError, KeyError):
                continue

    return atoms


def identify_chains(atoms: List[dict]) -> Tuple[str, str]:
    """
    Identify receptor and peptide chains by residue count.
    Returns (receptor_chain_id, peptide_chain_id).
    The receptor is the longer chain (MDM2, ~100 residues); the peptide
    is the shorter chain (cyclic binder, 8-12 residues).
    """
    chain_residues: dict = {}
    for atom in atoms:
        key = (atom["chain_id"], atom["seq_id"])
        chain_residues.setdefault(atom["chain_id"], set()).add(atom["seq_id"])

    chain_lengths = {c: len(r) for c, r in chain_residues.items()}
    if len(chain_lengths) != 2:
        raise ValueError(
            f"Expected exactly 2 chains, found {len(chain_lengths)}: "
            f"{chain_lengths}"
        )

    sorted_chains = sorted(chain_lengths.items(), key=lambda x: x[1])
    peptide_chain = sorted_chains[0][0]
    receptor_chain = sorted_chains[1][0]
    return receptor_chain, peptide_chain


def write_pdb(atoms: List[dict],
              output_path: Path,
              chain_override: str = "A",
              renumber: bool = True,
              strip_terminal_h: bool = False) -> List[str]:
    """
    Write atoms to a PDB file in standard Amber-compatible format.

    If strip_terminal_h is True, terminal hydrogens (H1, H2, H3, HXT,
    OXT) are omitted. This is required for cyclic peptides so tLeap
    does not add conflicting terminal patches.

    Returns the ordered list of residue names (3-letter codes).
    """
    # Build a renumbering map from old seq_id -> new 1-indexed residue number
    unique_residues = []
    seen = set()
    resname_by_seq = {}
    for atom in atoms:
        key = atom["seq_id"]
        if key not in seen:
            seen.add(key)
            unique_residues.append(atom["seq_id"])
            resname_by_seq[key] = atom["resname"]
    unique_residues.sort()
    resnum_map = {old: new for new, old in enumerate(unique_residues, start=1)}
    ordered_resnames = [resname_by_seq[s] for s in unique_residues]

    terminal_h_names = {"H1", "H2", "H3", "HXT", "OXT"}

    with output_path.open("w") as f:
        atom_serial = 1
        for atom in atoms:
            if strip_terminal_h and atom["atom_name"] in terminal_h_names:
                continue

            resnum = resnum_map[atom["seq_id"]] if renumber else atom["seq_id"]

            # PDB format is fixed-column
            name = atom["atom_name"]
            if len(name) < 4 and not name[0].isdigit():
                name = " " + name
            name = name.ljust(4)[:4]

            line = (
                f"ATOM  "
                f"{atom_serial:>5d} "
                f"{name}"
                f" "
                f"{atom['resname']:>3s} "
                f"{chain_override:1s}"
                f"{resnum:>4d}    "
                f"{atom['x']:>8.3f}"
                f"{atom['y']:>8.3f}"
                f"{atom['z']:>8.3f}"
                f"{1.00:>6.2f}"
                f"{0.00:>6.2f}"
                f"          "
                f"{atom['element']:>2s}"
            )
            f.write(line + "\n")
            atom_serial += 1
        f.write("TER\n")
        f.write("END\n")

    return ordered_resnames


def run_pdb4amber(input_pdb: Path, output_pdb: Path) -> None:
    """Run pdb4amber to clean up the PDB for tLeap."""
    result = subprocess.run(
        ["pdb4amber", "-i", str(input_pdb), "-o", str(output_pdb),
         "--nohyd", "--dry"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"pdb4amber stderr:\n{result.stderr}", file=sys.stderr)
        raise RuntimeError(f"pdb4amber failed on {input_pdb}")


def write_tleap_script(script_path: Path,
                       receptor_pdb: Path,
                       ligand_pdb: Path,
                       peptide_resnames: List[str],
                       output_dir: Path) -> None:
    """
    Write a tLeap script that loads both chains, creates three parameter
    sets, and closes the peptide cycle with an explicit bond command.

    Cyclic peptide construction strategy:
      - The receptor is loaded from PDB normally.
      - The peptide is built from sequence using 'sequence { ... }',
        which creates a chain where every residue uses its internal
        (non-terminal) residue template — no NGLY/CILE with wrong
        atom types. This is the critical fix: trying to load a cyclic
        peptide directly from PDB fails because tLeap auto-assigns
        N-terminal and C-terminal residue templates to the endpoints.
      - Coordinates are then loaded onto the sequence-built chain with
        'loadpdbusingseq', which snaps atom positions from a PDB file
        onto a pre-built sequence without re-assigning residue names.
      - An explicit 'bond pep.N.C pep.1.N' command closes the cycle.

    The GB radii (mbondi2) match the igb=5 GB model used later in MM-GBSA.
    """
    n = len(peptide_resnames)
    seq_str = " ".join(peptide_resnames)

    script = f"""# tLeap script for MDM2 + cyclic peptide parameter generation
# Generated by prepare_amber.py

source leaprc.protein.ff14SB
source leaprc.gaff2

# GB radii for MM-GBSA (igb=5 / mbondi2)
set default PBradii mbondi2

# Load receptor (MDM2) from PDB normally
rec = loadpdb {receptor_pdb}

# Build the cyclic peptide from sequence first so every residue uses
# its internal (non-terminal) template. Then snap coordinates from
# the cleaned PDB onto the sequence-built chain.
pep = sequence {{ {seq_str} }}
pep = loadpdbusingseq {ligand_pdb} {{ {seq_str} }}

# Create head-to-tail peptide bond
bond pep.{n}.C pep.1.N

# Build the complex by combining
comp = combine {{rec pep}}

# Check each unit for problems before saving
check rec
check pep
check comp

# Save Amber parameter and coordinate files
saveamberparm rec {output_dir}/receptor.prmtop {output_dir}/receptor.inpcrd
saveamberparm pep {output_dir}/ligand.prmtop {output_dir}/ligand.inpcrd
saveamberparm comp {output_dir}/complex.prmtop {output_dir}/complex.inpcrd

# Also save a PDB of the complex for visualization
savepdb comp {output_dir}/complex.pdb

quit
"""
    script_path.write_text(script)


def run_tleap(script_path: Path, log_path: Path) -> None:
    """Run tLeap and capture the log."""
    result = subprocess.run(
        ["tleap", "-f", str(script_path)],
        capture_output=True, text=True,
    )
    log_path.write_text(result.stdout + "\n--- STDERR ---\n" + result.stderr)
    if result.returncode != 0:
        print(f"tLeap log written to {log_path}", file=sys.stderr)
        raise RuntimeError(
            f"tLeap failed (exit {result.returncode}). See {log_path}"
        )
    # Even on success, check for FATAL errors in output
    if "FATAL" in result.stdout or "FATAL" in result.stderr:
        print(f"tLeap log written to {log_path}", file=sys.stderr)
        raise RuntimeError(f"tLeap reported FATAL errors. See {log_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Prepare Amber parameters from a Boltz-2 CIF output"
    )
    parser.add_argument("--boltz-cif", type=Path, required=True,
                        help="Path to Boltz model_0.cif file")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Directory to write Amber files into")
    args = parser.parse_args()

    if not args.boltz_cif.is_file():
        print(f"Error: {args.boltz_cif} not found", file=sys.stderr)
        sys.exit(1)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Parsing {args.boltz_cif}", file=sys.stderr)
    atoms = parse_cif_atoms(args.boltz_cif)
    if not atoms:
        print("Error: no atoms parsed from CIF", file=sys.stderr)
        sys.exit(1)

    receptor_chain, peptide_chain = identify_chains(atoms)
    print(f"Receptor chain: {receptor_chain}, peptide chain: {peptide_chain}",
          file=sys.stderr)

    receptor_atoms = [a for a in atoms if a["chain_id"] == receptor_chain]
    peptide_atoms = [a for a in atoms if a["chain_id"] == peptide_chain]

    # Write raw PDB files with chain IDs normalized to A (receptor) and
    # B (peptide), and strip terminal hydrogens from the peptide so
    # tLeap will not complain about the cyclization setup.
    raw_receptor_pdb = args.output_dir / "receptor_raw.pdb"
    raw_peptide_pdb = args.output_dir / "peptide_raw.pdb"

    receptor_resnames = write_pdb(
        receptor_atoms, raw_receptor_pdb, chain_override="A",
        renumber=True, strip_terminal_h=False,
    )
    peptide_resnames = write_pdb(
        peptide_atoms, raw_peptide_pdb, chain_override="B",
        renumber=True, strip_terminal_h=True,
    )
    print(f"Receptor: {len(receptor_resnames)} residues", file=sys.stderr)
    print(f"Peptide:  {len(peptide_resnames)} residues "
          f"({' '.join(peptide_resnames)})", file=sys.stderr)

    # Clean with pdb4amber
    print("Running pdb4amber...", file=sys.stderr)
    clean_receptor_pdb = args.output_dir / "receptor_clean.pdb"
    clean_peptide_pdb = args.output_dir / "peptide_clean.pdb"
    run_pdb4amber(raw_receptor_pdb, clean_receptor_pdb)
    run_pdb4amber(raw_peptide_pdb, clean_peptide_pdb)

    # Generate tLeap script and run it
    script_path = args.output_dir / "build.leap"
    log_path = args.output_dir / "tleap.log"
    write_tleap_script(
        script_path, clean_receptor_pdb, clean_peptide_pdb,
        peptide_resnames, args.output_dir,
    )
    print("Running tLeap...", file=sys.stderr)
    run_tleap(script_path, log_path)

    # Verify outputs
    expected = [
        "complex.prmtop", "complex.inpcrd", "complex.pdb",
        "receptor.prmtop", "receptor.inpcrd",
        "ligand.prmtop", "ligand.inpcrd",
    ]
    missing = [f for f in expected if not (args.output_dir / f).is_file()]
    if missing:
        print(f"ERROR: missing expected outputs: {missing}", file=sys.stderr)
        print(f"See tLeap log at {log_path}", file=sys.stderr)
        sys.exit(1)

    print(f"\nSuccess. Output files in {args.output_dir}:", file=sys.stderr)
    for f in expected:
        size = (args.output_dir / f).stat().st_size
        print(f"  {f} ({size} bytes)", file=sys.stderr)


if __name__ == "__main__":
    main()

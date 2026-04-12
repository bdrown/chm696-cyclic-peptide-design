#!/usr/bin/env python3
from __future__ import annotations

"""
analyze_results.py

Stage 5 of the MDM2 cyclic peptide evaluation pipeline. Collects
results from all three evaluation methods (ODesign ranking, Boltz-2
confidence, MM-GBSA binding energy) and produces a comparison report
with plots suitable for a problem set submission.

Outputs (written to --output-dir):
    results.tsv             -- machine-readable table, one row per design
    results.md              -- formatted markdown summary for reports
    bar_delta_total.png     -- ΔG_bind ± SEM bar chart across designs
    scatter_iptm_vs_ddg.png -- Boltz-2 ipTM vs MM-GBSA ΔG_bind
    components.png          -- MM-GBSA energy component breakdown
    per_frame_delta.png     -- per-frame ΔG time traces (diagnostic)

Usage:
    ml conda
    source activate $SCRATCH/conda/envs/openmm
    python analyze_results.py \\
        --manifest $SCRATCH/boltz/inputs/odesign_eval_manifest.tsv \\
        --mmgbsa-root $SCRATCH/mmgbsa \\
        --boltz-root $SCRATCH/boltz/outputs/odesign_eval/boltz_results_odesign_eval/predictions \\
        --output-dir $SCRATCH/results

Compatible with Python 3.9+. Requires numpy, matplotlib, scipy
(all in the openmm conda env).
"""

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless backend for SLURM/login nodes
import matplotlib.pyplot as plt
from scipy import stats


def read_manifest(manifest_path: Path) -> List[dict]:
    """
    Read the Boltz manifest TSV written by prepare_boltz_inputs.py.
    Returns a list of dicts with keys:
        design_id, source_cif, peptide_sequence, boltz_yaml
    The ordering is preserved and indicates ODesign's sample order.
    """
    designs = []
    with manifest_path.open() as f:
        for line in f:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            designs.append({
                "design_id": parts[0],
                "source_cif": parts[1],
                "peptide_sequence": parts[2],
                "boltz_yaml": parts[3],
            })
    return designs


def read_boltz_confidence(boltz_root: Path,
                          design_id: str) -> Optional[dict]:
    """
    Read Boltz-2 confidence JSON for a given design.

    Boltz writes one confidence JSON per diffusion sample:
        <boltz_root>/<design_id>/confidence_<design_id>_model_N.json
    We read model_0 (the top-ranked sample) by default.
    """
    pred_dir = boltz_root / design_id
    if not pred_dir.is_dir():
        return None

    # Look for confidence_*_model_0.json
    candidates = list(pred_dir.glob(f"confidence_{design_id}_model_0.json"))
    if not candidates:
        return None

    with candidates[0].open() as f:
        data = json.load(f)

    # Normalize the keys we care about. Boltz has evolved its output
    # schema over versions, so we try several possibilities.
    result = {
        "confidence_score": data.get("confidence_score",
                                     data.get("confidence", None)),
        "ptm": data.get("ptm", None),
        "iptm": data.get("iptm", None),
        "ligand_iptm": data.get("ligand_iptm", None),
        "complex_plddt": data.get("complex_plddt",
                                  data.get("plddt", None)),
    }
    return result


def read_mmgbsa_dat(dat_path: Path) -> Optional[dict]:
    """
    Parse mmgbsa.dat to extract the DELTA section energy components.
    Returns a dict mapping term name to (avg, std, sem) tuples, or
    None if the file is missing / incomplete.
    """
    if not dat_path.is_file():
        return None

    text = dat_path.read_text()
    if "DELTA TOTAL" not in text:
        return None

    result = {}
    in_delta = False
    known = {
        "VDWAALS", "EEL", "EGB", "ESURF",
        "DELTA G gas", "DELTA G solv", "DELTA TOTAL",
    }

    for line in text.splitlines():
        stripped = line.strip()
        if "Differences (Complex - Receptor - Ligand)" in stripped:
            in_delta = True
            continue
        if not in_delta or not stripped or stripped.startswith("-"):
            continue

        parts = stripped.split()
        if len(parts) < 4:
            continue
        try:
            avg = float(parts[-3])
            std = float(parts[-2])
            sem = float(parts[-1])
        except ValueError:
            continue

        name = " ".join(parts[:-3])
        if name in known:
            result[name] = (avg, std, sem)

    return result if result else None


def read_mmgbsa_per_frame(csv_path: Path) -> Optional[np.ndarray]:
    """
    Parse mmgbsa.csv to compute per-frame DELTA TOTAL values.

    MMPBSA.py writes the CSV in three sections — Complex Energy Terms,
    Receptor Energy Terms, Ligand Energy Terms — each preceded by a
    title line and a column header. There is no standalone DELTA
    section, so we compute ΔG per frame as:

        ΔG[i] = TOTAL_complex[i] − TOTAL_receptor[i] − TOTAL_ligand[i]

    The layout looks like:

        GENERALIZED BORN:
        Complex Energy Terms
        Frame #,BOND,ANGLE,...,TOTAL
        0,340.16,938.12,...,-2783.55
        1,400.72,939.13,...,-2759.49
        ...
        Receptor Energy Terms
        Frame #,BOND,ANGLE,...,TOTAL
        0,...
        ...
        Ligand Energy Terms
        Frame #,BOND,ANGLE,...,TOTAL
        0,...
        ...
    """
    if not csv_path.is_file():
        return None

    # Parse into three dicts: {frame_index: total_energy}
    section_totals = {"Complex": {}, "Receptor": {}, "Ligand": {}}
    current_section = None
    total_col_idx = None

    with csv_path.open() as f:
        for raw in f:
            line = raw.rstrip()
            if not line:
                continue

            # Section title lines
            if line.startswith("Complex Energy Terms"):
                current_section = "Complex"
                total_col_idx = None
                continue
            if line.startswith("Receptor Energy Terms"):
                current_section = "Receptor"
                total_col_idx = None
                continue
            if line.startswith("Ligand Energy Terms"):
                current_section = "Ligand"
                total_col_idx = None
                continue

            if current_section is None:
                # Still in preamble (e.g., "GENERALIZED BORN:")
                continue

            # Column header — first non-numeric line after a section title
            if line.startswith("Frame"):
                cols = [c.strip() for c in line.split(",")]
                try:
                    total_col_idx = cols.index("TOTAL")
                except ValueError:
                    total_col_idx = None
                continue

            if total_col_idx is None:
                continue

            # Data row — leading column is frame index, TOTAL is at total_col_idx
            parts = line.split(",")
            if len(parts) <= total_col_idx:
                continue
            try:
                frame_idx = int(parts[0])
                total = float(parts[total_col_idx])
            except ValueError:
                continue
            section_totals[current_section][frame_idx] = total

    if not (section_totals["Complex"]
            and section_totals["Receptor"]
            and section_totals["Ligand"]):
        return None

    # Compute ΔG per frame on the intersection of frames present in all three
    common_frames = sorted(
        set(section_totals["Complex"])
        & set(section_totals["Receptor"])
        & set(section_totals["Ligand"])
    )
    if not common_frames:
        return None

    deltas = np.array([
        section_totals["Complex"][i]
        - section_totals["Receptor"][i]
        - section_totals["Ligand"][i]
        for i in common_frames
    ])
    return deltas


def collect_results(manifest: List[dict],
                    mmgbsa_root: Path,
                    boltz_root: Path) -> List[dict]:
    """
    For each design in the manifest, gather ODesign rank, Boltz-2
    confidence metrics, and MM-GBSA energies into a single record.
    """
    records = []
    for i, design in enumerate(manifest, start=1):
        design_id = design["design_id"]
        rec = {
            "rank_odesign": i,
            "design_id": design_id,
            "peptide_sequence": design["peptide_sequence"],
            "peptide_length": len(design["peptide_sequence"]),
        }

        # Boltz-2 confidence
        boltz = read_boltz_confidence(boltz_root, design_id)
        if boltz is not None:
            rec.update({
                "boltz_confidence": boltz.get("confidence_score"),
                "boltz_ptm": boltz.get("ptm"),
                "boltz_iptm": boltz.get("iptm"),
                "boltz_complex_plddt": boltz.get("complex_plddt"),
            })
        else:
            rec.update({
                "boltz_confidence": None, "boltz_ptm": None,
                "boltz_iptm": None, "boltz_complex_plddt": None,
            })

        # MM-GBSA
        design_dir = mmgbsa_root / design_id
        dat_path = design_dir / "mmgbsa.dat"
        csv_path = design_dir / "mmgbsa.csv"

        mm = read_mmgbsa_dat(dat_path)
        if mm is not None:
            for term in ("VDWAALS", "EEL", "EGB", "ESURF",
                         "DELTA G gas", "DELTA G solv", "DELTA TOTAL"):
                key = term.lower().replace(" ", "_")
                if term in mm:
                    avg, std, sem = mm[term]
                    rec[f"{key}_avg"] = avg
                    rec[f"{key}_std"] = std
                    rec[f"{key}_sem"] = sem
                else:
                    rec[f"{key}_avg"] = None
                    rec[f"{key}_std"] = None
                    rec[f"{key}_sem"] = None

        per_frame = read_mmgbsa_per_frame(csv_path)
        rec["per_frame_total"] = per_frame  # numpy array or None

        records.append(rec)

    return records


def write_tsv(records: List[dict], output_path: Path) -> None:
    """Write records as a TSV with all scalar fields."""
    scalar_keys = [
        "rank_odesign", "design_id", "peptide_sequence", "peptide_length",
        "boltz_confidence", "boltz_ptm", "boltz_iptm", "boltz_complex_plddt",
        "vdwaals_avg", "vdwaals_sem",
        "eel_avg", "eel_sem",
        "egb_avg", "egb_sem",
        "esurf_avg", "esurf_sem",
        "delta_g_gas_avg", "delta_g_gas_sem",
        "delta_g_solv_avg", "delta_g_solv_sem",
        "delta_total_avg", "delta_total_sem",
    ]

    def fmt(v):
        if v is None:
            return ""
        if isinstance(v, float):
            return f"{v:.4f}"
        return str(v)

    with output_path.open("w") as f:
        f.write("\t".join(scalar_keys) + "\n")
        for rec in records:
            f.write("\t".join(fmt(rec.get(k)) for k in scalar_keys) + "\n")


def write_markdown(records: List[dict],
                   output_path: Path,
                   correlations: dict) -> None:
    """Write a markdown summary report."""
    lines = []
    lines.append("# MDM2 Cyclic Peptide Binder Evaluation")
    lines.append("")
    lines.append("Comparison of three methods for scoring ODesign-generated "
                 "cyclic peptide binders to MDM2:")
    lines.append("")
    lines.append("1. **ODesign** – generative ranking (order of generation)")
    lines.append("2. **Boltz-2** – interface confidence (ipTM) from the "
                 "predicted complex structure")
    lines.append("3. **MM-GBSA** – physics-based binding free energy from "
                 "a 10 ns implicit-solvent MD trajectory")
    lines.append("")
    lines.append("## Summary Table")
    lines.append("")
    lines.append("| Design | Sequence | Length | ODesign rank | Boltz-2 ipTM | "
                 "MM-GBSA ΔG_bind (kcal/mol) |")
    lines.append("|---|---|---|---|---|---|")

    def cell(v, fmt_spec="{:.3f}"):
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return "–"
        try:
            return fmt_spec.format(v)
        except (TypeError, ValueError):
            return str(v)

    for rec in records:
        dg_avg = rec.get("delta_total_avg")
        dg_sem = rec.get("delta_total_sem")
        if dg_avg is not None and dg_sem is not None:
            dg_str = f"{dg_avg:+.2f} ± {dg_sem:.2f}"
        else:
            dg_str = "–"

        lines.append(
            f"| {rec['design_id']} "
            f"| `{rec['peptide_sequence']}` "
            f"| {rec['peptide_length']} "
            f"| {rec['rank_odesign']} "
            f"| {cell(rec.get('boltz_iptm'))} "
            f"| {dg_str} |"
        )

    lines.append("")
    lines.append("## Rank Correlations")
    lines.append("")
    lines.append("Spearman rank correlations between the three ranking methods "
                 "(higher-magnitude values indicate stronger agreement). "
                 "Note that for MM-GBSA, more negative ΔG_bind means stronger "
                 "binding, so the rank is inverted for comparison purposes.")
    lines.append("")
    lines.append("| Comparison | Spearman ρ | p-value |")
    lines.append("|---|---|---|")
    for label, (rho, pval) in correlations.items():
        if rho is None:
            lines.append(f"| {label} | – | – |")
        else:
            lines.append(f"| {label} | {rho:+.3f} | {pval:.3f} |")

    lines.append("")
    lines.append("## MM-GBSA Component Decomposition")
    lines.append("")
    lines.append("Energy terms in kcal/mol (mean ± SEM over trajectory frames).")
    lines.append("")
    lines.append("| Design | VDW | EEL | EGB | ESURF | ΔG_gas | ΔG_solv | ΔG_total |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for rec in records:
        def tc(prefix):
            a = rec.get(f"{prefix}_avg")
            s = rec.get(f"{prefix}_sem")
            if a is None:
                return "–"
            return f"{a:+.2f}±{s:.2f}" if s is not None else f"{a:+.2f}"

        lines.append(
            f"| {rec['design_id']} "
            f"| {tc('vdwaals')} "
            f"| {tc('eel')} "
            f"| {tc('egb')} "
            f"| {tc('esurf')} "
            f"| {tc('delta_g_gas')} "
            f"| {tc('delta_g_solv')} "
            f"| **{tc('delta_total')}** |"
        )

    lines.append("")
    lines.append("## Figures")
    lines.append("")
    lines.append("- `bar_delta_total.png` – MM-GBSA ΔG_bind with error bars "
                 "across all designs")
    lines.append("- `scatter_iptm_vs_ddg.png` – Boltz-2 ipTM versus MM-GBSA "
                 "ΔG_bind, highlighting method agreement / disagreement")
    lines.append("- `components.png` – MM-GBSA energy component breakdown "
                 "(vdW, electrostatic, polar solvation, nonpolar)")
    lines.append("- `per_frame_delta.png` – per-frame ΔG trace for each "
                 "design, useful for diagnosing dissociation events")
    lines.append("")
    lines.append("## Methodological Notes")
    lines.append("")
    lines.append("- MD protocol: 100 ps heating (50→300 K) + 10 ns production "
                 "at 300 K, 2 fs timestep, Langevin thermostat, HBonds "
                 "constrained, GB-OBC2 implicit solvent, 0.15 M salt.")
    lines.append("- MM-GBSA: single-trajectory approach, 500 frames, "
                 "`igb=5` with `mbondi2` radii (consistent with the MD solvent "
                 "model), entropy term omitted.")
    lines.append("- Cyclic peptides: explicit head-to-tail N-C peptide bond "
                 "added to the ff14SB topology in tLeap via sequence-based "
                 "residue construction.")
    lines.append("- Boltz-2's affinity prediction head was NOT used because "
                 "its training data does not support reliable affinity "
                 "prediction for peptidic binders with >50 atoms.")

    output_path.write_text("\n".join(lines) + "\n")


def compute_correlations(records: List[dict]) -> dict:
    """
    Compute Spearman rank correlations between the three methods.
    Skips correlations where either vector has insufficient data or
    zero variance.
    """
    od_rank = [r["rank_odesign"] for r in records]
    iptm = [r.get("boltz_iptm") for r in records]
    ddg = [r.get("delta_total_avg") for r in records]

    def safe_spearman(a, b):
        pairs = [(x, y) for x, y in zip(a, b)
                 if x is not None and y is not None]
        if len(pairs) < 3:
            return (None, None)
        xs, ys = zip(*pairs)
        if len(set(xs)) < 2 or len(set(ys)) < 2:
            return (None, None)
        result = stats.spearmanr(xs, ys)
        return (float(result.statistic), float(result.pvalue))

    # Note: more negative ΔG = stronger binding. ODesign rank 1 = top.
    # Boltz-2 higher ipTM = more confident. Inverting ΔG sign so all
    # "higher is better" before correlating makes the sign intuition
    # consistent for the reader.
    ddg_neg = [-x if x is not None else None for x in ddg]
    od_rank_inv = [-x for x in od_rank]  # so rank 1 becomes highest

    return {
        "ODesign rank vs Boltz-2 ipTM": safe_spearman(od_rank_inv, iptm),
        "ODesign rank vs MM-GBSA (−ΔG)": safe_spearman(od_rank_inv, ddg_neg),
        "Boltz-2 ipTM vs MM-GBSA (−ΔG)": safe_spearman(iptm, ddg_neg),
    }


def plot_bar_delta_total(records: List[dict], output_path: Path) -> None:
    """Bar chart of MM-GBSA ΔG_bind with SEM error bars."""
    labels = [r["design_id"] for r in records]
    avgs = [r.get("delta_total_avg") for r in records]
    sems = [r.get("delta_total_sem") for r in records]

    # Replace None with 0 for plotting but mark visually
    plot_avgs = [a if a is not None else 0.0 for a in avgs]
    plot_sems = [s if s is not None else 0.0 for s in sems]

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ["#1f77b4" if a is not None and a < -5 else "#d62728"
              for a in avgs]
    bars = ax.bar(labels, plot_avgs, yerr=plot_sems,
                  color=colors, capsize=4, edgecolor="black")

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("MM-GBSA ΔG$_\\mathrm{bind}$ (kcal/mol)")
    ax.set_title("MM-GBSA binding free energy by design")
    ax.tick_params(axis="x", rotation=30)
    for lbl in ax.get_xticklabels():
        lbl.set_horizontalalignment("right")

    # Annotate each bar with its value
    for bar, a in zip(bars, avgs):
        if a is None:
            continue
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() - 1.5,
                f"{a:+.1f}",
                ha="center", va="top", fontsize=9, color="white",
                fontweight="bold")

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_scatter_iptm_vs_ddg(records: List[dict],
                             output_path: Path) -> None:
    """Scatter of Boltz-2 ipTM vs MM-GBSA ΔG_bind."""
    pts = [(r.get("boltz_iptm"), r.get("delta_total_avg"),
            r.get("delta_total_sem"), r["design_id"])
           for r in records]
    pts = [p for p in pts if p[0] is not None and p[1] is not None]
    if not pts:
        return

    iptm_vals, ddg_vals, ddg_sems, labels = zip(*pts)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.errorbar(iptm_vals, ddg_vals, yerr=ddg_sems,
                fmt="o", markersize=10, capsize=5,
                color="#1f77b4", ecolor="gray", markeredgecolor="black")

    for x, y, label in zip(iptm_vals, ddg_vals, labels):
        # Strip "_seed123" suffix for cleaner labels
        short = re.sub(r"_seed\d+", "", label)
        ax.annotate(short, (x, y), textcoords="offset points",
                    xytext=(8, 5), fontsize=9)

    ax.axhline(0, color="red", linewidth=0.8, linestyle="--",
               alpha=0.5, label="ΔG = 0 (no binding)")
    ax.set_xlabel("Boltz-2 ipTM (interface confidence)")
    ax.set_ylabel("MM-GBSA ΔG$_\\mathrm{bind}$ (kcal/mol)")
    ax.set_title("Boltz-2 confidence vs. MM-GBSA binding energy")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_components(records: List[dict], output_path: Path) -> None:
    """Stacked bar chart of MM-GBSA energy components."""
    labels = [r["design_id"] for r in records]
    n = len(records)

    vdw = np.array([r.get("vdwaals_avg") or 0 for r in records])
    eel = np.array([r.get("eel_avg") or 0 for r in records])
    egb = np.array([r.get("egb_avg") or 0 for r in records])
    esurf = np.array([r.get("esurf_avg") or 0 for r in records])
    total = np.array([r.get("delta_total_avg") or np.nan
                      for r in records])

    fig, ax = plt.subplots(figsize=(9, 5.5))
    x = np.arange(n)
    width = 0.17

    ax.bar(x - 1.5 * width, vdw, width, label="VDW", color="#2ca02c")
    ax.bar(x - 0.5 * width, eel, width, label="EEL (Coulomb)", color="#d62728")
    ax.bar(x + 0.5 * width, egb, width, label="EGB (polar solv)",
           color="#1f77b4")
    ax.bar(x + 1.5 * width, esurf, width, label="ESURF (nonpolar)",
           color="#ff7f0e")

    # Overlay total as markers
    ax.plot(x, total, "ko", markersize=10, label="ΔG total",
            markerfacecolor="black", zorder=10)

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_title("MM-GBSA component decomposition")
    ax.legend(loc="upper right", ncol=2, fontsize=9)
    ax.grid(True, alpha=0.3, axis="y")

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_per_frame_delta(records: List[dict], output_path: Path) -> None:
    """Per-frame ΔG time trace for each design on a shared axis."""
    fig, ax = plt.subplots(figsize=(9, 5.5))

    any_plotted = False
    for rec in records:
        trace = rec.get("per_frame_total")
        if trace is None:
            continue
        # X axis: frame index → time in ns, assuming 20 ps per frame
        time_ns = np.arange(len(trace)) * 0.020
        short_label = re.sub(r"_seed\d+", "", rec["design_id"])
        ax.plot(time_ns, trace, label=short_label, linewidth=1.5, alpha=0.85)
        any_plotted = True

    if not any_plotted:
        plt.close(fig)
        return

    ax.axhline(0, color="red", linewidth=0.8, linestyle="--", alpha=0.5)
    ax.set_xlabel("Simulation time (ns)")
    ax.set_ylabel("Per-frame ΔG (kcal/mol)")
    ax.set_title("Per-frame MM-GBSA binding energy over the trajectory")
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Compare ODesign / Boltz-2 / MM-GBSA evaluation results"
    )
    parser.add_argument("--manifest", type=Path, required=True,
                        help="Path to Boltz manifest TSV from prepare_boltz_inputs.py")
    parser.add_argument("--mmgbsa-root", type=Path, required=True,
                        help="Root containing per-design subdirectories")
    parser.add_argument("--boltz-root", type=Path, required=True,
                        help="Boltz predictions directory "
                             "(contains <design_id>/ subdirectories)")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Directory to write report and figures")
    args = parser.parse_args()

    if not args.manifest.is_file():
        print(f"Error: manifest {args.manifest} not found", file=sys.stderr)
        sys.exit(1)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading manifest: {args.manifest}", file=sys.stderr)
    manifest = read_manifest(args.manifest)
    print(f"  {len(manifest)} designs", file=sys.stderr)

    print(f"Collecting results...", file=sys.stderr)
    records = collect_results(manifest, args.mmgbsa_root, args.boltz_root)

    # Brief summary to stderr
    print("\nResults summary:", file=sys.stderr)
    print(f"{'Design':<25}  {'ipTM':>6}  {'ΔG_bind':>14}", file=sys.stderr)
    for r in records:
        iptm = r.get("boltz_iptm")
        dg = r.get("delta_total_avg")
        dgs = r.get("delta_total_sem")
        iptm_s = f"{iptm:.3f}" if iptm is not None else "–"
        if dg is not None and dgs is not None:
            dg_s = f"{dg:+7.2f}±{dgs:.2f}"
        else:
            dg_s = "–"
        print(f"{r['design_id']:<25}  {iptm_s:>6}  {dg_s:>14}",
              file=sys.stderr)

    # Correlations
    correlations = compute_correlations(records)

    # Write table and report
    tsv_path = args.output_dir / "results.tsv"
    write_tsv(records, tsv_path)
    print(f"\nWrote {tsv_path}", file=sys.stderr)

    md_path = args.output_dir / "results.md"
    write_markdown(records, md_path, correlations)
    print(f"Wrote {md_path}", file=sys.stderr)

    # Plots
    plot_bar_delta_total(records, args.output_dir / "bar_delta_total.png")
    print(f"Wrote bar_delta_total.png", file=sys.stderr)

    plot_scatter_iptm_vs_ddg(records,
                             args.output_dir / "scatter_iptm_vs_ddg.png")
    print(f"Wrote scatter_iptm_vs_ddg.png", file=sys.stderr)

    plot_components(records, args.output_dir / "components.png")
    print(f"Wrote components.png", file=sys.stderr)

    plot_per_frame_delta(records, args.output_dir / "per_frame_delta.png")
    print(f"Wrote per_frame_delta.png", file=sys.stderr)

    print(f"\nAll outputs written to {args.output_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()

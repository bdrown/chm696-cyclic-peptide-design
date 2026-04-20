# MDM2 Cyclic Peptide Binder Evaluation

Comparison of three methods for scoring ODesign-generated cyclic peptide binders to MDM2:

1. **ODesign** – generative ranking (order of generation)
2. **Boltz-2** – interface confidence (ipTM) from the predicted complex structure
3. **MM-GBSA** – physics-based binding free energy from a 10 ns implicit-solvent MD trajectory

## Summary Table

| Design | Sequence | Length | ODesign rank | Boltz-2 ipTM | MM-GBSA ΔG_bind (kcal/mol) |
|---|---|---|---|---|---|
| design_01_seed123 | `GENVIDGVKI` | 10 | 1 | 0.953 | -22.19 ± 0.31 |
| design_02_seed123 | `AGSEGSANGT` | 10 | 2 | 0.927 | -34.20 ± 0.28 |
| design_03_seed123 | `NPGNSAKAPG` | 10 | 3 | 0.941 | -55.77 ± 0.50 |
| design_04_seed123 | `WTPTCWFDNK` | 10 | 4 | 0.960 | -42.65 ± 0.31 |
| design_05_seed123 | `GKKIGKNIID` | 10 | 5 | 0.895 | -5.77 ± 0.33 |

## Rank Correlations

Spearman rank correlations between the three ranking methods (higher-magnitude values indicate stronger agreement). Note that for MM-GBSA, more negative ΔG_bind means stronger binding, so the rank is inverted for comparison purposes.

| Comparison | Spearman ρ | p-value |
|---|---|---|
| ODesign rank vs Boltz-2 ipTM | +0.300 | 0.624 |
| ODesign rank vs MM-GBSA (−ΔG) | +0.100 | 0.873 |
| Boltz-2 ipTM vs MM-GBSA (−ΔG) | +0.500 | 0.391 |

## MM-GBSA Component Decomposition

Energy terms in kcal/mol (mean ± SEM over trajectory frames).

| Design | VDW | EEL | EGB | ESURF | ΔG_gas | ΔG_solv | ΔG_total |
|---|---|---|---|---|---|---|---|
| design_01_seed123 | -11.63±0.42 | -71.05±1.18 | +62.96±1.30 | -2.46±0.05 | -82.68±1.46 | +60.49±1.25 | **-22.19±0.31** |
| design_02_seed123 | -46.43±0.24 | -103.80±1.04 | +122.14±1.04 | -6.11±0.03 | -150.23±1.14 | +116.03±1.03 | **-34.20±0.28** |
| design_03_seed123 | -54.47±0.33 | -181.80±1.53 | +188.21±1.35 | -7.71±0.04 | -236.27±1.67 | +180.50±1.32 | **-55.77±0.50** |
| design_04_seed123 | -60.96±0.32 | -37.93±0.87 | +63.55±0.82 | -7.31±0.04 | -98.89±0.99 | +56.24±0.80 | **-42.65±0.31** |
| design_05_seed123 | -7.42±0.41 | +0.77±0.93 | +2.01±0.94 | -1.14±0.06 | -6.64±1.14 | +0.87±0.90 | **-5.77±0.33** |

## Figures

- `bar_delta_total.png` – MM-GBSA ΔG_bind with error bars across all designs
- `scatter_iptm_vs_ddg.png` – Boltz-2 ipTM versus MM-GBSA ΔG_bind, highlighting method agreement / disagreement
- `components.png` – MM-GBSA energy component breakdown (vdW, electrostatic, polar solvation, nonpolar)
- `per_frame_delta.png` – per-frame ΔG trace for each design, useful for diagnosing dissociation events

## Methodological Notes

- MD protocol: 100 ps heating (50→300 K) + 10 ns production at 300 K, 2 fs timestep, Langevin thermostat, HBonds constrained, GB-OBC2 implicit solvent, 0.15 M salt.
- MM-GBSA: single-trajectory approach, 500 frames, `igb=5` with `mbondi2` radii (consistent with the MD solvent model), entropy term omitted.
- Cyclic peptides: explicit head-to-tail N-C peptide bond added to the ff14SB topology in tLeap via sequence-based residue construction.
- Boltz-2's affinity prediction head was NOT used because its training data does not support reliable affinity prediction for peptidic binders with >50 atoms.

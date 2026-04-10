# Stage 4 — Cross-method comparison and reporting

This stage pulls together the results from Stages 1-3 and produces a comparison report with summary table, rank correlations, and plots. It is the final step of the pipeline and the source of the figures and numbers you will include in your writeup.

## What this stage does

`analyze_results.py` reads:

- The Stage 2 manifest TSV for design identities, sequences, and ODesign sample order
- Boltz-2 confidence JSON files for ipTM, pTM, and related metrics
- `mmgbsa.dat` files from Stage 3 for MM-GBSA energy component averages
- `mmgbsa.csv` files from Stage 3 for per-frame ΔG time traces

and produces:

- `results.tsv` — machine-readable table with one row per design
- `results.md` — formatted markdown report with summary table, correlations, and component decomposition
- `bar_delta_total.png` — MM-GBSA ΔG_bind ± SEM bar chart
- `scatter_iptm_vs_ddg.png` — Boltz-2 ipTM versus MM-GBSA ΔG_bind scatter plot
- `components.png` — MM-GBSA energy component breakdown (vdW, EEL, EGB, ESURF)
- `per_frame_delta.png` — per-frame ΔG time traces for all designs (diagnostic plot)

## Running the analysis

Activate the conda environment and run the script with paths to the three input sources:

```bash
ml conda
source activate $SCRATCH/conda/envs/openmm

python analyze_results.py \
    --manifest $SCRATCH/boltz/inputs/odesign_eval_manifest.tsv \
    --mmgbsa-root $SCRATCH/mmgbsa \
    --boltz-root $SCRATCH/boltz/outputs/odesign_eval/boltz_results_odesign_eval/predictions \
    --output-dir $SCRATCH/results
```

The script runs in a few seconds and prints a summary table to the terminal as it goes. Check the terminal output for any warnings about missing data — these usually indicate that one of the upstream stages didn't complete for a particular design.

## The outputs

### `results.md`

The markdown report has three sections:

**Summary Table** — one row per design with ODesign rank, Boltz-2 ipTM, and MM-GBSA ΔG_bind ± SEM. This is the headline table for your writeup.

**Rank Correlations** — Spearman rank correlations between each pair of methods. We use Spearman rather than Pearson because rank-based statistics are more robust with small sample sizes and are not distorted by outliers. For n=5 designs the p-values will not be statistically significant even with perfect agreement, but the direction of the correlation is what matters.

The correlations normalize all three metrics so that "higher is better" — ODesign rank 1 becomes the highest score, Boltz-2 ipTM is unchanged, and MM-GBSA ΔG_bind is negated. This makes the sign of each correlation interpretable as "do the methods agree on which designs are best."

**MM-GBSA Component Decomposition** — per-design breakdown into vdW, Coulomb (EEL), polar solvation (EGB), nonpolar solvation (ESURF), and the total. This table is where the interesting chemistry shows up — you can often identify different binding modes (electrostatic-dominated vs hydrophobic-dominated) by looking at which terms dominate.

### The figures

**`bar_delta_total.png`** — the simplest visualization. Bar chart of ΔG_bind for each design with error bars. Strong binders (ΔG < −5 kcal/mol) are plotted in blue; anything weaker is plotted in red to draw visual attention to potential non-binders.

**`scatter_iptm_vs_ddg.png`** — this is the "method comparison" plot. X-axis is Boltz-2 ipTM, Y-axis is MM-GBSA ΔG_bind. A design that Boltz-2 is confident about but MM-GBSA reports zero binding energy will show up as an obvious outlier (high ipTM, Y near zero). This is the figure that most directly illustrates the pedagogical point that ML confidence and physics-based binding energy measure different things.

**`components.png`** — grouped bar chart showing vdW, EEL, EGB, and ESURF for each design side by side, with the total ΔG_bind overlaid as black dots. Useful for identifying binding modes (look at which designs have large EEL contributions vs which rely on vdW).

**`per_frame_delta.png`** — per-frame ΔG vs time, all designs on the same axes. This is the diagnostic plot that reveals trajectory stability problems. A design that dissociated mid-trajectory will show a visible walk from negative ΔG values toward zero; a design that was never bound will hover near zero throughout; stable binders will fluctuate around their mean value.

## Interpreting the results

The correlation structure typically tells an interesting story. With our example dataset:

- **ODesign vs anything**: near-zero correlation. ODesign's "rank" is just generation order, not a confidence ranking, so this is expected and should be discussed in your report — it illustrates that generative ranking and scoring ranking are different things.
- **Boltz-2 ipTM vs MM-GBSA**: moderate but noisy correlation. These methods often agree on the extremes (bad designs are recognized by both) but can disagree on the middle of the ranking.
- **The most informative result is disagreement.** If Boltz-2 is confident about a structure but MM-GBSA reports zero binding, that tells you something important about the difference between structural plausibility and thermodynamic stability.

Use the report as the starting point for answering the questions in `docs/problem_set.md`.

## Converting to PDF for submission

The markdown report is designed to be convertible to PDF via pandoc or any markdown-to-PDF tool. On Scholar:

```bash
# If pandoc is available via module
ml pandoc 2>/dev/null || true
pandoc $SCRATCH/results/results.md -o $SCRATCH/results/results.pdf

# Or copy to your local machine and convert with your preferred tool
```

Include the four PNG figures alongside the PDF when you submit to Gradescope. The figures are at 150 DPI which is appropriate for printing at reasonable sizes.

## Troubleshooting

**"Per-frame plot is empty"** — the script could not parse `mmgbsa.csv`. MMPBSA.py's CSV format has changed across AmberTools versions, so if you are running a different version than this pipeline was tested with, the column parsing may fail. Check `mmgbsa.csv` manually and see if it has a "DELTA Energy Terms" section with a "TOTAL" column.

**"Correlations are all None"** — this means either the Boltz-2 or MM-GBSA values are missing for most designs. Check `bash ../03_mmgbsa/check_status.sh $SCRATCH/mmgbsa` to confirm all designs completed the full pipeline.

**"One design is missing from the table"** — the analysis script prints warnings when it can't find data for a design. Scroll up in the terminal output to see which design was skipped and why.

For other problems, see `docs/troubleshooting.md`.

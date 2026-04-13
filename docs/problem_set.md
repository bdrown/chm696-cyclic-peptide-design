# Problem Set: Evaluating Cyclic Peptide Binders to MDM2

## Background

MDM2 is an E3 ubiquitin ligase and negative regulator of the tumor suppressor p53. It binds p53 via a shallow hydrophobic cleft on its N-terminal domain and targets p53 for proteasomal degradation. In many cancers, MDM2 is overexpressed or amplified, which effectively silences p53 even when the *TP53* gene itself is intact. Pharmacologically disrupting the p53-MDM2 interaction is therefore a long-standing target for anticancer drug discovery.

Two classes of MDM2 binders have been extensively characterized:

- **The p53 transactivation helix** (residues 17-29 of p53) binds MDM2 via three hydrophobic anchor residues — **Phe19, Trp23, and Leu26** — that dock into subpockets within the MDM2 cleft. This interaction is characterized in the cocrystal structure PDB **1YCR**.
- **Nutlin-3a** is a small molecule cis-imidazoline that was developed by Hoffmann-La Roche as a p53-MDM2 inhibitor. It mimics the three hydrophobic anchors using halogenated aromatic rings. The cocrystal structure is PDB **1RV1** (nutlin-2) or **4HG7** (nutlin-3a-related).

In this problem set, you will evaluate five computationally designed cyclic peptide binders against MDM2 using three qualitatively different methods (generative ranking, ML structure prediction, and physics-based MD) and compare your designs' binding modes to those of p53 and nutlin.

## Deliverables

Submit to Gradescope as a **single PDF** containing:

1. Your answers to the questions below (Parts 1-5)
2. All four figures from Stage 4 (`bar_delta_total.png`, `scatter_iptm_vs_ddg.png`, `components.png`, `per_frame_delta.png`)
3. At least two PyMOL screenshots supporting your answers regarding the binding mode to MDM2 (see Part 4)

Before you start writing, make sure you have completed all four stages of the pipeline and that `bash 03_mmgbsa/check_status.sh $SCRATCH/mmgbsa` shows `done` for all sub-stages on all five designs. Then run `04_analysis/analyze_results.py` to produce `results.md` and supporting figures.

---

## Part 1: Understanding the methods (5 points)

**1.1 (3 pts)** In your own words, describe what question each of the three evaluation methods (ODesign rank, Boltz-2 ipTM, MM-GBSA ΔGbind) is actually answering. Why might they not agree?

**1.2 (2 pts)** The MM-GBSA calculation uses the single-trajectory approach rather than running separate MD simulations of the complex, receptor, and ligand alone. What does single-trajectory gain us in terms of statistical precision, and what does it cost us in terms of physical realism?

---

## Part 2: Interpreting your results (5 points)

Refer to your own `results.md` summary table and `components.png` figure.

**2.1 (2 pts)** Rank your five designs by MM-GBSA ΔG_bind (strongest to weakest). Which design is the strongest binder by this metric, and what is its ΔG_bind value (with uncertainty)?

**2.2 (3 pts)** Examine the MM-GBSA component decomposition (vdW, EEL, EGB, ESURF) in `components.png` for your top two binders. Are they binding through the same type of interactions, or are there qualitative differences in their binding modes? Specifically compare the magnitudes of the vdW and electrostatic (EEL) contributions. Propose a sequence-level explanation for any differences you see.

---

## Part 3: Comparison with p53 and nutlin binding modes (5 points)

Access the PDB structures **1YCR** (p53-MDM2) from the RCSB Protein Data Bank. Open your best-binding design's minimized structure (`$SCRATCH/mmgbsa/<best_design>/complex_min.pdb`) alongside them in PyMOL.

Align all three structures on MDM2 using PyMOL's `align` or `cealign` command so the binding clefts overlap. Take screenshots as needed.

**4 (5 pts)** Examine your best-binding cyclic peptide design. Which residues in your peptide occupy the Phe19, Trp23, and Leu26 subpockets? Does your design successfully recapitulate all three hydrophobic anchors, only some of them, or does it bind in a different mode entirely? Include a PyMOL screenshot showing your design's anchor residues alongside p53's.

---

## Part 5: Synthesis and critical thinking (15 points)

**5.1 (3 pts)** Suppose a classmate claims that Boltz-2's ipTM score is a sufficient metric for prioritizing binder candidates for experimental validation, arguing that "if the model is confident in the structure, the binding must be real." Based on your results, write a short rebuttal that identifies the specific weakness in this reasoning.

**5.3 (2 pts)** The absolute ΔG_bind values from MM-GBSA are typically 10 kcal/mol too negative compared to experimental Kd measurements. List at least two physical effects that are missing from our calculation and explain qualitatively how including them would change the computed binding energy.

---

## Grading notes

- Clarity and correctness of reasoning are weighted more heavily than getting "the right answer" — for several questions there is no single correct answer and a well-defended response will receive full credit.
- You are encouraged to use PyMOL, ChimeraX, or any other structural visualization tool to support your answers in Part 4.
- Cite any external sources you consult (papers, PDB entries, tutorial materials). You do not need to cite the course pipeline itself.

## Helpful commands

Check pipeline status:
```bash
bash 03_mmgbsa/check_status.sh $SCRATCH/mmgbsa
```

Regenerate the analysis outputs if you change anything:
```bash
ml conda
source activate $SCRATCH/conda/envs/openmm
python 04_analysis/analyze_results.py \
    --manifest $SCRATCH/boltz/inputs/odesign_eval_manifest.tsv \
    --mmgbsa-root $SCRATCH/mmgbsa \
    --boltz-root $SCRATCH/boltz/outputs/odesign_eval/boltz_results_odesign_eval/predictions \
    --output-dir $SCRATCH/results
```

Load your best design in PyMOL on a login node:
```bash
/class/bsdrown/apps/pymol/pymol $SCRATCH/mmgbsa/design_XX/complex_min.pdb
```

Fetch reference structures from the PDB within PyMOL:
```
fetch 1YCR, async=0
fetch 4HG7, async=0
align design_XX, 1YCR
```
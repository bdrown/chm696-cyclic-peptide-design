# Pipeline overview

This document explains what each stage of the pipeline does, why it does it that way, and how the stages fit together. Read this once to understand the project as a whole; come back to the individual stage READMEs when you need task-specific details.

## The scientific question

Modern drug discovery increasingly relies on computational methods to propose and prioritize candidate binders before any experimental work begins. Different methods answer different questions about a candidate, and it is easy to treat them as interchangeable when they are not. This pipeline demonstrates three qualitatively different approaches to evaluating whether a candidate cyclic peptide will bind MDM2:

1. **Generative ranking** — trust the order in which a model proposes candidates, on the assumption that earlier-generated designs are more confident.
2. **Structural confidence** — use an ML structure predictor to build the candidate-target complex and trust its internal confidence metrics (ipTM, pTM) as proxies for binding quality.
3. **Physics-based binding energy** — run a molecular dynamics simulation and compute a binding free energy from the ensemble.

These three methods answer slightly different questions:

- **"Can I build a plausible sequence that occupies the binding site?"** — generative ranking.
- **"Given this sequence, can I construct a geometrically reasonable bound complex?"** — Boltz-2 confidence.
- **"Will this complex be thermodynamically stable?"** — MM-GBSA.

When all three agree, you probably have a real binder. When they disagree, the disagreement itself tells you something about the nature of the problem. This pipeline is structured to expose those disagreements as the central pedagogical experience.

## The four stages in detail

### Stage 1: ODesign generates candidates

ODesign is a flexible-receptor generative diffusion model. You specify the target protein, a set of hotspot residues defining the binding site, and constraints like peptide length and whether to cyclize. ODesign returns a collection of candidate peptide sequences each paired with a predicted bound conformation.

For this pipeline we use MDM2 (residues 17-125 of UniProt Q00987, matching the 1YCR crystal structure) with 12 hotspot residues surrounding the p53 binding cleft, and we request cyclic peptides of 8-12 residues. ODesign runs with 3 random seeds and 5 samples per seed to give 15 total candidates with some diversity.

**Critical point**: ODesign is a backbone-only model. Its output CIF files contain Cα, N, C, and O atoms but no side chain atoms. This is a deliberate design choice — placing side chains conditional on a backbone is a different problem than generating backbones, and ODesign focuses on the latter. The missing side chains are added in Stage 2.

**Another critical point**: ODesign does not provide per-candidate confidence scores. The "ranking" we use in the downstream analysis is simply the order in which candidates were generated (sample 1, sample 2, ...). This is not a meaningful confidence ranking, and the near-zero correlation between ODesign rank and the other two methods in our results is therefore expected, not surprising.

### Stage 2: Boltz-2 refines to full atom

We take the top 5 ODesign candidates that pass a minimum-length filter (skipping any 3-4 residue designs that would be synthetically inaccessible) and re-predict each one's complex with MDM2 using Boltz-2. Boltz-2 is an open-source AlphaFold3-style co-folding model that jointly predicts backbone geometry, side chains, and complex structure.

We deliberately enable the `cyclic: true` flag on the peptide chain in the Boltz input YAML, which tells the model to treat the peptide as head-to-tail cyclic via a modified positional encoding. Boltz-2's handling of cyclic peptides is known to be weaker than its handling of linear peptides — its training data includes relatively few cyclic peptides — but it is the best we have.

For each design, Boltz-2 runs 5 diffusion samples and writes one CIF per sample along with per-model confidence JSONs. We use `model_0` (the top-ranked sample) for downstream analysis and extract `iptm` as the primary confidence metric. ipTM is specifically designed to score the interface between chains in multi-chain predictions, which is exactly what we want for a binder-target complex.

**What we do not use**: Boltz-2's binding affinity prediction head. Its authors have stated explicitly that this head is trained on drug-like small molecules with fewer than ~50 heavy atoms and should not be used for larger peptidic molecules. A 10-residue cyclic peptide has roughly 100 heavy atoms, so using the affinity head would be a methodological error. We rely on MM-GBSA in Stage 3 for binding energy estimates instead.

### Stage 3: MM-GBSA via molecular dynamics

This is the longest and most technically involved stage. The high-level workflow for each design is:

1. Convert the Boltz-2 CIF to Amber topology/coordinate files with explicit cyclization bonds
2. Energy-minimize the complex to resolve any strain from tLeap's automated hydrogen placement
3. Heat to 300 K and run 10 ns of production MD in implicit solvent
4. Run MMPBSA.py on the trajectory to compute ΔG_bind and its component decomposition

The three hardest technical problems in this stage are:

**Handling head-to-tail cyclization in tLeap.** AmberTools does not automatically detect that a peptide is cyclic. By default, tLeap assigns N-terminal and C-terminal residue patches to the first and last residues of any chain it loads, which gives them the wrong atom types (a protonated N3 nitrogen instead of a regular peptide bond nitrogen, and an OXT-containing C-terminus instead of a regular peptide bond carbon). When we then try to add the closing bond, ff14SB has no parameters for a bond between an N3 and a regular C and tLeap reports a cascade of missing-parameter errors.

The workaround is to build the peptide from sequence first using tLeap's `sequence { GLY ALA ... }` command, which creates a chain where every residue uses the internal (non-terminal) template. Then we snap coordinates from the PDB onto the pre-built chain with `loadpdbusingseq`, which transfers atom positions without re-assigning residue names. Finally we add an explicit `bond pep.N.C pep.1.N` to close the ring. The closing bond is then between atom types `N` and `C` — a regular amide bond that ff14SB has parameters for.

**Consistent GB model across stages.** MMPBSA.py computes binding energies using a GB implicit solvent model. For the numbers to be meaningful, the trajectory must have been generated with the same GB model and the same atomic radii. We use `igb=5` (GB-OBC2) with `mbondi2` radii throughout: tLeap sets the radii when generating the prmtop, OpenMM uses the matching `OBC2` model during MD, and MMPBSA.py uses the same `igb=5` for scoring. Any mismatch between these would silently produce wrong binding energies.

**Single-trajectory MM-GBSA.** There are two ways to compute MM-GBSA binding energies. The dual-trajectory approach runs separate MD simulations of the complex, the receptor alone, and the ligand alone, then computes ΔG_bind = G_complex - G_receptor - G_ligand by averaging over each separately. This captures any conformational reorganization upon binding but is noisy because the internal energies of receptor and ligand do not exactly cancel between bound and unbound states.

The single-trajectory approach runs MD only on the complex and extracts receptor-only and ligand-only "snapshots" from each frame of the complex trajectory. Internal energies cancel exactly by construction, which dramatically reduces noise. The cost is that you ignore conformational reorganization — the ligand is assumed to adopt the same conformation bound and free. For a relative comparison among closely related binders, this is the correct tradeoff and is what most published MM-GBSA studies use. We adopt it here.

### Stage 4: Cross-method comparison

The analysis script reads outputs from all three previous stages and joins them into a single comparison. It computes Spearman rank correlations between each pair of methods (using rank-based statistics because n=5 is too small for meaningful Pearson correlations), produces a summary table, and makes four plots:

- **Bar chart** of ΔG_bind by design with error bars
- **Scatter plot** of Boltz-2 ipTM vs MM-GBSA ΔG_bind (the key method-comparison figure)
- **Grouped bar chart** of MM-GBSA energy components (vdW, EEL, EGB, ESURF)
- **Time traces** of per-frame ΔG for each design, which diagnoses trajectory stability

The per-frame time trace is the single most pedagogically valuable figure. It lets you see, in real data, when a peptide dissociates from its target during MD — a physical failure mode that pure ML-based methods are blind to.

## The idempotence and dependency-chain design of Stage 3

Stage 3 has 4 sub-stages × 5 designs = 20 jobs. Managing these by hand would be painful. The stage is organized so that:

1. **Each sub-stage checks its own outputs before running.** If the expected files already exist and look valid, the job exits successfully without doing anything. This means re-running after a partial completion is cheap.

2. **The submission wrapper checks outputs at submission time** and only submits jobs that actually need to run. If stages 1-2 are done for a design but stages 3-4 failed, only stages 3-4 get submitted for the retry.

3. **SLURM dependencies chain the sub-stages with `afterok`.** Downstream jobs only start if their upstream predecessor succeeded with exit 0. A failure at any sub-stage cancels everything downstream for that design without wasting compute.

4. **A status reporter queries SLURM and the filesystem** to show students the current state of all 5 designs × 4 sub-stages in a single table. This is the primary tool for understanding "where am I" during a multi-hour pipeline run.

The combination of per-stage idempotence and dependency chaining makes failure recovery trivial in the common case: fix whatever went wrong, re-run `submit_all_designs.sh`, and everything that was already done is skipped while everything that was broken gets retried.

## A note on what this pipeline is for

This pipeline is a teaching tool, not a production drug discovery workflow. Specifically:

- **10 ns MD is short.** Real MM-GBSA analysis of protein-peptide binding typically uses 100 ns to 1 μs of sampling. Ten nanoseconds is enough to get a rough estimate and to distinguish obvious non-binders from plausible binders, but you should not trust small differences between designs.
- **Implicit solvent is an approximation.** GB solvation models are calibrated to reproduce bulk water effects on average but miss specific water-mediated interactions (bridging waters, etc.) that can matter at protein-peptide interfaces. An explicit-solvent protocol would be more rigorous but much slower.
- **Entropy is ignored.** We do not compute the normal-mode or quasi-harmonic entropy contribution to ΔG_bind. For closely related compounds entropies tend to be similar and partially cancel, so ignoring them is tolerable for relative ranking. For absolute affinity prediction, ignoring entropy introduces errors of 10-30 kcal/mol.
- **MDM2 histidine protonation is default.** pdb4amber assigns HIE to all histidines without checking local environment. For MDM2 this is fine because no key binding-site histidines have unusual protonation, but for other targets you might need to handle this explicitly.
- **No water or salt is modeled explicitly.** Whole classes of interactions (bridging waters, salt-bridge screening dynamics) are therefore missing from the simulation.

These limitations are part of the pedagogical value. Students should understand what their numbers do and do not represent, and the report questions in `docs/problem_set.md` are designed to surface this understanding.

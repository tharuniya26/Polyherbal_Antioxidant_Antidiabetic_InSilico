# 3L4W/MIG diagnostic native-ligand redocking

**Experiment:** 12 October 2026 UTC (13 October in Sri Lanka)  
**Scope:** pose-recovery test for deposited miglitol (MIG) in human N-terminal maltase-glucoamylase, [PDB 3L4W](https://www.rcsb.org/structure/3L4W). This is a local-pocket diagnostic, not production phytochemical docking or biochemical validation. The [3L4W chemistry review](vina_3L4W_chemistry_review.md) defines the input microstate and its limitations.

## Inputs and pre-run checks

The primary receptor was `3L4W_noNAG_D443m_D542m_H600HID_pH6p5_dry_sensitivity_Meeko0.8.0.pdbqt`: a **dry, no-NAG sensitivity model** with Asp443⁻, Asp542⁻ and neutral His600-HID. The ligand was the unchanged, neutral [MIG crystal-reference PDBQT](../../docking/validation/vina/3L4W_MIG_crystal_pose_reference.pdbqt); the unchanged [native SDF](../../docking/validation/vina/3L4W_MIG_native_pose.sdf) was the RMSD reference. The script verified these files and the deposited PDB and both atom maps against the recorded SHA256 hashes before docking. The receptor PDBQT SHA256 was `7a76e0fa3b7bf2197714a57a88c60d2f4eefedccb18636951f3f643473afadf3`; MIG SDF was `ea673a59043bea05abbe7d66d8bd29c3786021f564c394027a8b0201385503f7`; MIG PDBQT was `c7417f6bcb661e0f2217c6ba3619b6fd6c03ffb029fe5840d4b0d13c1bf16169`. The [manifest](../../docking/validation/vina/3L4W_diagnostic_redocking/experiment_manifest.json) records all audited hashes, including the executable.

RDKit and Meeko readback checked MIG as neutral C₈H₁₇NO₅ with **14 heavy atoms**, the complete expected molecular graph, four defined carbon stereocentres, valid 3D coordinates and the deposited heavy-atom naming/coordinate map. The four reference carbon configurations were 9S, 10R, 11R and 12R in the audited SDF indexing. Graph enumeration found one allowed heavy-atom automorphism, so no alternative symmetry mapping altered the reported RMSDs. The selected microstate remains a computational hypothesis; neither the CCD nor diffraction fixes the amine's bound-state protonation.

The center was the centroid of the **deposited MIG heavy atoms**, **(45.38664286, 92.35157143, 34.82200000) Å**. Box dimensions were **24 × 18 × 18 Å**. Each axis extends at least 5 Å beyond the deposited MIG atoms and catalytic Asp443/Asp542/His600 side-chain heavy atoms, with half-extents rounded up to whole Å. MIG heavy-atom extrema were (42.501, 89.475, 32.147) to (48.486, 94.054, 36.873) Å. The site and ligand geometry, rather than an optimized score, defined this box.

## Search design and reproducibility

The crystal pose was never used as a starting pose or overwritten. For each seed, NumPy PCG64 generated a deterministic random rigid quaternion rotation and translation within ±3 Å of the native centroid; atom ordering, bonding, stereochemistry and protonation were unchanged. Each accepted start lay within the docking box and was at least 5 Å **unfitted heavy-atom RMSD** from the crystal. The five measured start RMSDs were 5.2422, 5.0564, 5.2932, 5.4181 and 5.3926 Å in seed order. These files provide a visible anti-bias check. The script also supplied the stated seed to Vina's global search, whose internal search randomizes ligand configurations independently of the input placement. The installed Vina 1.2.7 `--randomize_only` path generates a fresh internal seed rather than honoring the supplied seed, so it was not used to construct the reproducible starts ([versioned Vina source](https://github.com/ccsb-scripps/AutoDock-Vina/blob/v1.2.7/src/lib/vina.cpp); [Vina FAQ](https://autodock-vina.readthedocs.io/en/stable/faq.html)).

Five sequential Vina 1.2.7 runs used seeds **1101, 2202, 3303, 4404, 5505**; `--scoring vina --exhaustiveness 32 --num_modes 20 --energy_range 5.0 --min_rmsd 1.0 --cpu 1`. There were no setting deviations. The requested 20 modes are an upper bound: seed 4404 returned 18. The [results CSV](../../data/metadata/vina_3L4W_redocking_results.csv) contains each exact command, runtime, input/output hashes, seed, exit code and status. Logs, randomized inputs, pose outputs, per-seed metrics and [all 98 pose metrics](../../docking/validation/vina/3L4W_diagnostic_redocking/all_pose_metrics.csv) are preserved in the [diagnostic output directory](../../docking/validation/vina/3L4W_diagnostic_redocking/).

Reproduce or independently reanalyse with the project Python and [script](../../scripts/run_vina_3L4W_redocking.py):

```powershell
.\.venv\Scripts\python.exe scripts\run_vina_3L4W_redocking.py
```

The script verifies existing artifacts and analyses them without overwriting an earlier docking output. Software versions were Python 3.12.10, RDKit 2026.03.6, Meeko 0.8.0 and AutoDock Vina 1.2.7; the Vina executable SHA256 is in the manifest. An initial post-run analysis expected separate Meeko molecules for separate PDBQT modes. Meeko 0.8.0 instead read them as conformers of one molecule; the parser was corrected and **the already completed first docking run was not repeated**. This was an analysis error, not a Vina run failure.

## RMSD and ligand quality control

Each PDBQT output was read through Meeko's PDBQT parser and reconstructed with RDKit. For every returned mode, the script checked the 14-heavy-atom identity/map, graph, stereocentres, bond-length sanity and finite coordinates. Coordinates were compared to the preserved SDF in the deposited receptor frame, using the audited SDF-to-PDBQT atom map and any chemically allowed graph automorphism. **No independent pose fitting or superposition** was applied. RMSD is therefore a crystal-frame, heavy-atom pose-recovery measure. Vina's own inter-mode RMSD column was not substituted for this calculation. A ≤2.0 Å top-pose threshold was specified as a conventional initial screen, not a proof of scoring accuracy.

## Results

**Five completed runs; zero docking failures; 98 returned poses.** Energies below are Vina scores in kcal/mol, not measured binding free energies or inhibition constants.

| Seed | Modes | Rank-1 score | Rank-1 RMSD (Å) | Best RMSD (Å) | Best-RMSD rank | Its score |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1101 | 20 | −5.548 | 1.8775 | 0.6863 | 2 | −5.449 |
| 2202 | 20 | −5.535 | 1.8691 | **0.6806** | 2 | −5.431 |
| 3303 | 20 | −5.544 | **1.8629** | 1.8629 | 1 | −5.544 |
| 4404 | 18 | −5.540 | 1.8818 | 0.7227 | 2 | −5.363 |
| 5505 | 20 | −5.531 | 1.8849 | 1.8849 | 1 | −5.531 |

All five rank-1 poses crossed the prespecified ≤2.0 Å screen. Their RMSDs ranged from 1.8629 to 1.8849 Å (mean 1.8752 Å); pairwise rank-1 heavy-atom RMSDs were approximately 0.01–0.05 Å in the common receptor frame. This supports **search convergence to one rank-1 binding mode** under the tested settings. The best overall pose was seed 2202, **rank 2 at 0.6806 Å**; the best rank-1 pose was seed 3303 at **1.8629 Å**. A more crystal-like rank-2 pose occurred in three of five runs, yet Vina ranked it behind a roughly 1.87 Å pose by 0.099–0.177 kcal/mol. This small score gap does not establish physical superiority of either mode.

Only 8 of 98 returned poses were below 2 Å; 90 were at or above 2 Å, including many distant alternative modes. The successful rank-1 convergence therefore coexists with substantial alternative-pose sampling. In the seed-1101 comparison, both the native-like and rank-1 poses remain near the crystallographic catalytic pocket: nearest heavy-atom distances to Asp443 were 2.758 Å in the crystal, 3.140 Å for rank 1 and 3.108 Å for rank 2; for Asp542, 2.629, 2.700 and 2.599 Å, respectively. These geometric contacts do not establish hydrogen-bond direction, catalysis or affinity.

## Interpretation and next validation gate

**Diagnostic native-pose recovery was achieved for this dry no-NAG, neutral-MIG receptor/ligand model** under the initial ≤2 Å rank-1 criterion, with five reproducible runs. It is not full validation of the crystallographic receptor, Vina's energy ranking, glycan treatment, biochemical inhibition, or the planned phytochemical library. [The structural study](https://doi.org/10.1021/bi9016457) and [RCSB entry](https://www.rcsb.org/structure/3L4W) establish the experimental complex; they do not verify this computational microstate.

The model omits six covalently linked NAG units, including NAG A2005 nearest the catalytic region, and all active-site waters; it uses provisional Asp542⁻, His600-HID and neutral MIG states at a preparation pH of 6.5. Gln837 and other deposited segments remain unobserved; no missing coordinates were invented. Diffraction does not resolve the relevant hydrogens. Static-water retention, glycan omission and ligand protonation can shift the binding-site model. The second-ranked near-native poses in three runs also indicate that search recovery and score ordering should be assessed separately.

The receptor may advance **only to a controlled validation sensitivity stage**: repeat multiple seeds at higher exhaustiveness to check mode ranking and convergence, then assess predeclared MIG-amine/Asp542/His600 and selected-water variants, and determine whether a chemically valid glycan-retaining receptor can be tested. Apply the same atom-map RMSD method to every variant. Record any failure or changed ranking. Do not infer potency from raw scores; do not advance this receptor to production phytochemical docking on this experiment alone. No phytochemical docking or ADMET was performed.

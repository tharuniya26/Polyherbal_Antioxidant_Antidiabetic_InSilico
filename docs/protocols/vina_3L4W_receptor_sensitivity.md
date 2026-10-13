# 3L4W receptor-protonation sensitivity of native MIG pose recovery

**Scope:** human N-terminal maltase-glucoamylase (ntMGAM), deposited [3L4W/MIG complex](https://www.rcsb.org/structure/3L4W); conditional pH-6.5 receptor microstates, using AutoDock Vina 1.2.7 and Meeko 0.8.0. This is diagnostic crystallographic pose recovery, not a measurement of inhibitor affinity. The [previous validation assessment](vina_3L4W_validation_assessment.md) established five-seed recovery for the dry, no-NAG Asp443⁻/Asp542⁻/His600-HID receptor with neutral and modelled amine-protonated miglitol (MIG). All previous receptors, ligands, crystal coordinates, atom maps and docking outputs remain unchanged.

## Mechanistic basis and protonation evidence

[Primary structural work on human ntMGAM and ntSI](https://pmc.ncbi.nlm.nih.gov/articles/PMC2878540/) identifies **Asp443 as the catalytic nucleophile** and **Asp542 as the general acid/base**; His600 is a conserved substrate-binding residue. The [3L4W structural study](https://doi.org/10.1021/bi9016457) establishes the native MIG complex, not the positions of catalytic hydrogens. Thus Asp443⁻ is retained in every tested model, while Asp542 protonation is a plausible sensitivity variable. No diffraction-derived hydrogen or experimental site-specific pKa assignment is claimed.

The earlier [chemistry review](vina_3L4W_chemistry_review.md) ran PROPKA 3.5.1 on the same deposited heavy coordinates. Predicted Asp542 pKa changed from **7.76 without MIG** to **5.59 with neutral MIG**; at the operational pH **6.5**, these imply opposite dominant-state tendencies. Predicted Asp443 values were **5.19/3.41** and His600 values **0.83/0.01** (apo/neutral-MIG input). These are model-dependent predictions on a no-NAG protein with an unobserved Gln837 gap. They justify testing an Asp542 acid form but do not establish its population in the actual bound enzyme. The pH is the deposited crystallization value and does not describe every intestinal condition.

| Receptor hypothesis | Assessment | Decision |
| --- | --- | --- |
| Asp443⁻/Asp542⁻/His600-HID | Previously prepared dry no-NAG baseline; five-seed E64 pose recovery already recorded. | Preserve as comparator; do not alter. |
| **Asp443⁻/Asp542H(OD2)/His600-HID** | Standard Meeko ASH residue template, with the proton on deposited **OD2**. Ligand-free OpenMM 8.6.1 placed its direction; the O–H vector was normalized to 0.98 Å. Same deposited protein heavy coordinates and all other hydrogens were retained. | **Prepared conditional sensitivity variant**; diagnostic redocking performed only after checks below. |
| Asp542H on deposited OD1 | OD1 lies 2.629 Å from MIG OAC. A directly donor-oriented 0.98 Å OD1–H would lie approximately **0.734 Å from the existing modeled MIG OAC hydrogen**, an unacceptable H–H overlap. A different joint ligand/receptor hydroxyl orientation and validated OD1 acid template would be needed. | No variant or docking; this rotamer is unresolved, not chemically ruled out. |
| His600-HIE | A ligand-free OpenMM HIE placement, normalized to 1.0 Å N–H, came **1.127 Å from the modeled MIG OAD hydrogen** in either supplied MIG SDF. The earlier review found a 0.940 Å overlap in another trial. Ligand hydroxyl H is not crystallographically fixed, so joint hydrogen reorientation might resolve this. | No receptor-only HIE comparison with the unchanged ligand hydrogen model. |
| His600-HIP | Present in an earlier automatic preparation, but doubly protonated His600 is disfavored by both PROPKA estimates and lacks direct crystallographic support. | No HIP redocking. |

The prepared Asp542H(OD2) geometry is **one acid rotamer**, not all possible Asp542 protonation arrangements. Its OH faces away from MIG. This avoids a modeled donor–donor clash and avoids choosing H orientation to optimize the native ligand pose, but it may differ from a catalytic proton-transfer configuration.

## New receptor quality control

The [reproducibility script](../../scripts/run_vina_3L4W_receptor_sensitivity.py) hashes the audited baseline and produces two newly named files: [Asp542H working PDB](../../structures/proteins/prepared/3L4W_noNAG_D443m_D542H_OD2_H600HID_pH6p5_dry_sensitivity.pdb) and [Meeko PDBQT](../../structures/proteins/vina_prepared/3L4W_noNAG_D443m_D542H_OD2_H600HID_pH6p5_dry_sensitivity_Meeko0.8.0.pdbqt). Meeko was explicitly instructed to use its packaged `ASH` template at A542. That template defines protonated OD2 connectivity and donor hydrogen rather than interpreting an untyped extra H. No residue was deleted or silently renamed.

- All **6,912 protein heavy coordinates and atom names** matched the baseline to PDB precision; every pre-existing hydrogen was unchanged. One HD2 atom was added at **(47.208, 88.638, 33.437) Å**. Its OD2–H length is 0.98 Å and CG–OD2–HD2 angle is 112.16°.
- The new H is at least **3.413 Å from another protein residue** and **2.555 Å from either modeled MIG ligand**. This is a clash screen, not a hydrogen-bond-energy calculation.
- Meeko PDBQT readback retained the 6,912 heavy coordinates to 0.001 Å and contained 8,429 receptor atoms. The Asp542 residue charge sums to approximately **0**, Asp443 to **−1**, and His600-HID to **0**; OD2 is an acceptor oxygen with an attached donor H (`HD`). These partial-charge checks establish template consistency, not physical electrostatic accuracy in standard Vina.
- Compared atom by atom with the unchanged 8,428-atom baseline PDBQT, the only added identity is Asp542 HD2; **all atom types, partial charges and coordinates outside Asp542 match** within 0.001 Å/charge unit. This verifies that the controlled comparison changes only the intended receptor residue in Meeko's exported model.
- Vina 1.2.7 parsed the new receptor with both native-MIG PDBQTs in a **non-scoring zero-grid guard probe**. Both stopped at the expected positive-grid requirement, after receptor and ligand parsing. Actual input usability is further tested by the diagnostic runs below.

Input hashes, ASH template hash, H geometry, software versions, commands and the exact deposited-ligand-derived box are stored in the [experiment manifest](../../docking/validation/vina/3L4W_receptor_sensitivity/experiment_manifest.json). The PDB SHA256 is `8e0872b236cb6975dc81ef131ffd3aacddbd8edd3bbc178e45eb1a18710192fd`; the PDBQT SHA256 is `220e3585558e4bc261e6e76ddb036738dea827543e1b135399506e95ea9f93fd`. To reproduce preparation and quality checks without docking:

```powershell
.\.venv\Scripts\python.exe scripts\run_vina_3L4W_receptor_sensitivity.py --prepare-only
```

## Glycans, water and other structural limitations

The six deposited NAG units are **protein-linked glycans**, not free sugars. The nearest, A2005 linked to Asn209, is **16.409 Å from MIG** and **10.316 Å from the Asp443/Asp542/His600 catalytic group** in the deposited geometry. The dry no-NAG model preserves the local *heavy coordinates by construction*; that does not test glycan-dependent conformational or electrostatic effects. The installed Meeko 0.8.0 packaged residue templates contain no NAG entry, and the available project templates contain no validated four-linkage, stereochemistry- and charge-checked glycan model for these six linked units. Treating isolated NAG CCD units as complete or deleting the link is chemically unsound. A glycan-retaining PDBQT therefore remains **blocked pending a validated polymer-linkage template and charges**; the present receptor is still a local sensitivity model.

HOH A4002/A4004 contact the deposited MIG site, but their hydrogen directions, residence and exchange during ligand binding are not established. Standard rigid-receptor Vina fixes waters and does not provide the [specialized hydrated-docking treatment](https://autodock-vina.readthedocs.io/en/stable/docking_hydrated.html). A static two-water PDBQT exists from prior preparation; running it as if it established solvent sensitivity would overstate the model. Both receptor protonation cases compared here remain **dry**. Gln837 is unobserved far from MIG and was not reconstructed.

## Controlled redocking comparison

The diagnostic test uses the **same neutral and protonated MIG PDBQTs**, deposited heavy-atom RMSD reference, symmetry-aware atom map, and box center **(45.38664286, 92.35157143, 34.82200000) Å**, size **24 × 18 × 18 Å** as the earlier E64 runs. Seeds are **1101, 2202, 3303, 4404, 5505**, with Vina scoring, exhaustiveness 64, `num_modes 20`, `energy_range 5.0`, `min_rmsd 1.0` and one CPU. Each ligand begins from its independently randomized, ≥5 Å-from-crystal start. The previous script's Meeko readback, explicit charged-N geometry checks and unfitted crystal-frame 14-heavy-atom RMSD calculation are reused. All returned modes, logs and hashes are retained under [3L4W receptor sensitivity outputs](../../docking/validation/vina/3L4W_receptor_sensitivity/).

The **neutral MIG** group completed all five seeds. Each top-ranked pose passed the initial ≤2 Å screen. Values in the comparator columns come from the immutable [earlier E64 results](../../data/metadata/vina_3L4W_sensitivity_results.csv), with the same ligand, seeds, box and settings.

| Seed | Asp542⁻ rank-1 RMSD (Å) | Asp542H rank-1 RMSD (Å) | Asp542⁻ best RMSD (Å) | Asp542H best RMSD (Å) | Asp542H best rank | Asp542H rank-1 Vina score |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1101 | 1.8775 | 1.8814 | 0.6997 | 0.6989 | 2 | −5.547 |
| 2202 | 1.8691 | 1.8916 | 0.6806 | **0.6813** | 2 | −5.531 |
| 3303 | 1.8629 | 1.8813 | 1.8629 | 1.8813 | 1 | −5.533 |
| 4404 | 1.8782 | 1.8634 | 0.7227 | 1.8634 | 1 | −5.546 |
| 5505 | 1.8849 | 1.8632 | 0.6926 | 0.6816 | 2 | −5.548 |

Neutral-MIG **rank-1 recovery was 5/5 for both receptor states**. The top-pose RMSD difference by seed was at most 0.0225 Å. The closer approximately 0.7 Å secondary pose occurred in 4/5 baseline seeds and 3/5 Asp542H seeds; its absence in seed 4404 of the latter is a sampling/ranking difference that should not be hidden. Vina score changes were small, but scores are not experimental affinities and cannot establish a protonation-state preference.

The **modelled amine-protonated MIG** group also completed all five seeds. The charged reference uses the same 14 deposited heavy-atom coordinates and unchanged stereodefined carbon framework; its protonated-N geometry was checked in every output, as described in the [prior ligand sensitivity assessment](vina_3L4W_validation_assessment.md).

| Seed | Asp542⁻ rank-1 RMSD (Å) | Asp542H rank-1 RMSD (Å) | Asp542⁻ best RMSD (Å) | Asp542H best RMSD (Å) | Asp542H best rank | Asp542H rank-1 Vina score |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1101 | 1.8759 | 1.8759 | 0.6183 | 0.6183 | 2 | −5.770 |
| 2202 | 1.8360 | 1.8360 | 0.6223 | 0.6223 | 2 | −5.796 |
| 3303 | 1.8302 | 1.8302 | 0.6166 | 0.6166 | 2 | −5.799 |
| 4404 | 1.8669 | 1.8669 | 0.6136 | **0.6136** | 2 | −5.793 |
| 5505 | 1.8628 | 1.8628 | 0.6149 | 0.6149 | 2 | −5.793 |

Charged-MIG **rank-1 recovery was 5/5 under both receptor states**, and its RMSDs matched to the displayed 0.0001 Å precision for every seed. In all five new charged runs the more crystal-like 0.61–0.62 Å pose ranked **second**, 0.043–0.102 Vina score units behind rank 1. This is a persistent score-ordering issue, not a failed search. The OD2 proton points away from MIG; its small effect in this scoring model must not be generalized to a ligand-facing Asp542 acid rotamer or to physical proton-coupled binding energetics. Standard Vina scoring does not use supplied partial charges as a full electrostatic model ([Vina FAQ](https://autodock-vina.readthedocs.io/en/stable/faq.html)).

**Ten of ten new docking runs completed, zero failed**, returning **197 poses** (98 neutral and 99 charged). Rank-1 pairwise RMSD across seeds was 0.0152–0.0510 Å for neutral MIG and 0.0319–0.0938 Å for charged MIG in the common receptor frame: both converge tightly to their favored mode. The more native-like secondary mode appeared in **three of five neutral runs** and **all five charged runs**. Every returned mode was retained in the [pose-level table](../../docking/validation/vina/3L4W_receptor_sensitivity/all_pose_metrics.csv); the [10-row run table](../../data/metadata/vina_3L4W_receptor_sensitivity_results.csv) records its top-ranked and best-of-all-poses RMSD separately, command, seed, hashes, status and comparator. `num_modes 20` was requested, not guaranteed: neutral seeds 1101/3303 and charged seed 5505 returned 19 modes. Output, log and start hashes were checked against the results CSV, and the script reanalysed the files on rerun without launching new Vina jobs.

Elapsed wall times varied substantially because the workspace had extended idle or suspended intervals during some runs. They are preserved in logs for audit but **should not be used as comparative search-efficiency evidence**. No E128 test was added: the E64 rank-1 mode was already convergent, while the remaining concerns are score ordering and receptor chemistry.

## Interpretation and decision

**Native MIG rank-1 pose recovery is robust to the two tested receptor microstates**—Asp542⁻ and Asp542H(OD2)—for both neutral MIG and one modelled protonated-MIG geometry. This establishes a limited *local pose-recovery* result. It does not establish which catalytic microstate is populated, verify Vina's score ordering, or show biochemical inhibition. The ligand-charge comparison and the receptor-charge comparison are separate experiments; neither raw score difference is a measured affinity or bound-state pKa.

The Asp542H test is limited to the standard ASH **OD2** rotamer, whose added proton points away from the ligand. The ligand-facing OD1 acid arrangement, His600-HIE/HIP, glycan-retaining receptor, structurally uncertain waters, and missing Gln837 remain unresolved as described above. The reported success does not validate those alternatives or make this a complete biological-receptor model.

**Qualified exploratory phytochemical pose analysis may begin only as hypothesis generation** using these explicitly labelled dry/no-NAG Asp542⁻ and Asp542H(OD2) models, with the same target-specific receptor audit and controls. Report both receptor states, failed poses and sensitivity of pocket placement; avoid potency ranking or mechanistic conclusions from Vina scores. Work intended to compare phytochemical affinities, propose inhibition constants, or claim formulation effects still needs a chemically validated glycan-aware model or a defensible assessment of glycan omission, additional catalytic-state and water sensitivity, and independent biochemical controls. No phytochemical docking or ADMET was performed in this step.

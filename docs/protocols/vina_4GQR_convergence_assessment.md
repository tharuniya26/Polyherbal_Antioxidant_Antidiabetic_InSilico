# 4GQR/MYC exhaustiveness-64 convergence assessment

This is a controlled follow-up to the [five-seed exhaustiveness-32 (E32) diagnostic experiment](vina_4GQR_redocking_report.md) using the deposited myricetin (CCD **MYC**) pose in human pancreatic α-amylase [4GQR](https://www.rcsb.org/structure/4GQR). The prespecified screen is **rank-1, symmetry-aware, crystal-frame 23-heavy-atom RMSD ≤2.0 Å**. It tests pose recovery for an explicitly dry, no-NAG local model; it does not measure binding affinity, enzyme inhibition or antidiabetic activity.

## Frozen inputs and baseline verification

The [E64 script](../../scripts/run_vina_4GQR_convergence.py) checks all five E32 [CSV rows](../../data/metadata/vina_4GQR_redocking_results.csv), saved commands, logs, SHA256 hashes, PDBQT outputs and per-pose RMSD tables before any new docking. It reparses each E32 output through Meeko and verifies its rank-1 and best-pose values. The original E32 inputs and results remain unchanged.

| Constant across E32 and E64 | Value |
| --- | --- |
| Receptor | [Dry/no-NAG 4GQR local PDBQT](../../structures/proteins/vina_prepared/4GQR_dry_noNAG_CA_CL_PCA_pH7p5_diagnostic_Meeko0.8.0.pdbqt), SHA256 `71bee218c2cfede049ad8254bc39b33cb7c916055c0355e17236a2227a3e52b4` |
| Retained structural components | Calcium A502, chloride A503, terminal PCA A1; protein heavy coordinates and provisional pH-7.5 hydrogens unchanged |
| Ligand | Same neutral MYC PDBQT, unchanged [deposited-pose SDF](../../docking/validation/vina/4GQR_MYC_native_pose.sdf), 23-heavy-atom [bridge](../../docking/validation/vina/4GQR_MYC_sdf_to_pdbqt_atom_map.csv) |
| Box center and dimensions | **(13.44317391, 15.13530435, 39.43613043) Å; 26 × 24 × 28 Å** |
| Scoring and reporting | Vina 1.2.7 standard score; `num_modes 20`, `energy_range 5.0`, `min_rmsd 1.0`, `cpu 1` |
| Seeds and starts | 1101, 2202, 3303, 4404, 5505; the E32 PCG64/quaternion/translation procedure regenerates starts **byte-for-byte identical** to each frozen E32 start |
| Only experimental change | Exhaustiveness **32 → 64** |

The [manifest](../../docking/validation/vina/4GQR_convergence_E64/experiment_manifest.json) records software versions, executable/input hashes, settings, receptor quality checks and the E32 CSV hash. Each new command, log, randomized start, docking output and pose table is saved under [4GQR_convergence_E64](../../docking/validation/vina/4GQR_convergence_E64/). The reproducibility script refuses to overwrite changed completed artifacts and reanalyses saved outputs on a later run.

## Pose calculation and decision rule

Each returned PDBQT mode is reconstructed with Meeko/RDKit and checked for the neutral MYC graph, 23 mapped heavy atoms, atom identity/type, finite coordinates and original bond geometry. MYC has no stereogenic atom assignment in this reference. The SDF/PDBQT atom bridge and both allowed chemical graph automorphisms give the **minimum** heavy-atom RMSD against unchanged deposited coordinates in the common receptor frame, **without independently fitting each ligand pose**. Rank-1 RMSD and best RMSD among *all* saved modes are different outcomes. Vina's inter-mode RMSD is not substituted for crystal-reference RMSD. Every returned mode is retained, including alternative binding modes and failed-screen poses.

## Structural limits carried into E64

The unchanged local model omits one covalently linked NAG about 50.93 Å from MYC and all eight waters retained in an earlier wet PDBQT. The wet hydrogen model had an MYC–water H–H clash of about 0.907 Å; the dry choice is explicit and **does not establish the bound-solvent state**. Calcium and chloride are structurally retained, but standard Vina does not provide validated ion-coordination/electrostatic physics simply because it accepts their PDBQT atom types. Asp197/Glu233/Asp300/His201 and neutral MYC protonation remain computational choices. The [official MYC validation record](https://www.rcsb.org/ligand-validation/4gqr/MYC) reports 0.40 occupancy, real-space correlation 0.662 and geometry outliers, limiting certainty in the crystal reference. Those issues are unchanged by a larger search effort.

## E32 versus E64 outcomes

**All five E64 Vina runs and pose reconstructions completed; none failed.** E64 returned **99 poses** (20 for seeds 1101–4404 and 19 for 5505). `num_modes 20` is a request, not a guaranteed count. The [five-row result table](../../data/metadata/vina_4GQR_convergence_results.csv) contains exact commands, hashes, elapsed times, exit codes and E32 comparators; individual logs, starts, output PDBQTs and per-pose metrics are preserved in the new E64 directory. All recorded input/output/log hashes were independently checked after completion. Elapsed wall times varied substantially, especially seed 1101, and may include host idle periods; they are audit metadata, not evidence of relative search efficiency.

| Seed | E32 rank-1 RMSD (Å) | E64 rank-1 RMSD (Å) | E32 best RMSD (Å), rank | E64 best RMSD (Å), rank | E64 rank-1 score | E64 best-pose score |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1101 | 2.4926 | 2.4926 | 2.1490, 3 | 2.1490, 4 | −7.738 | −7.642 |
| 2202 | 2.4957 | 2.4957 | 2.1012, 3 | 2.1367, 3 | −7.744 | −7.686 |
| 3303 | 2.5407 | 2.5081 | 2.4679, 3 | 2.1312, 4 | −7.768 | −7.666 |
| 4404 | 2.5546 | 2.4913 | 2.1227, 5 | 2.1228, 5 | −7.776 | −7.633 |
| 5505 | 2.5601 | 2.5650 | **1.5328, 17** | **2.1103, 3** | −7.750 | −7.675 |

**Rank-1 ≤2.0 Å recovery stayed at 0/5 at both E32 and E64.** The E64 rank-1 RMSDs span **2.4913–2.5650 Å**; their mean decreased slightly from **2.5287 to 2.5105 Å**. That small numeric change does not cross the predefined threshold. E64 rank-1 poses form a tight geometric cluster: pairwise symmetry-aware, same-frame 23-heavy-atom RMSDs span **0.0070–0.1517 Å**. This supports repeatability of the *favored non-native mode* at these two search efforts, not successful native-pose recovery.

At E32, **one of five** seeds returned *any* pose below 2 Å: seed 5505, rank 17, 1.5328 Å. At E64, **zero of five** did. The best E64 pose overall was **2.1103 Å at rank 3**, also seed 5505. It is not the earlier rank-17 near-native geometry merely promoted in ranking: the nearest E64 pose to that E32 rank-17 pose in each seed remains approximately **1.81–1.83 Å away from it**, after allowing the two MYC graph automorphisms. Thus the specific 1.5328 Å alternative was **not reproduced or made more frequent** by E64. The best E64 poses in all five seeds were scored **0.058–0.143 Vina units less favorably** than their respective rank-1 poses. This is a score-ordering observation under this model, not a measured energetic difference.

## Scientific interpretation and next action

**The 4GQR/MYC local-model validation conclusion remains failed.** Doubling exhaustiveness produced some changes in lower-ranked poses and small changes in two rank-1 RMSDs, but it did not improve the **0/5** rank-1 pose-recovery fraction. The repeated approximately 2.5 Å rank-1 mode and absence of sub-2 Å E64 poses make insufficient E32 search effort alone an inadequate explanation for the failed screen. Yet these outputs cannot uniquely separate a scoring preference from missing waters, provisional receptor or ligand protonation, the low-occupancy MYC reference, or other receptor-model error. No causal assignment is justified from these two exhaustiveness settings alone. Vina scores are neither binding affinities nor phytochemical inhibitory potencies.

The next justified investigation is a **separate receptor/reference chemistry audit** before further search escalation: inspect the deposited MYC electron-density fit and geometry; determine whether a chemically sound, nonclashing binding-site water model can be made; review catalytic and MYC phenol microstates at an explicit pH; and assess an ion-aware or independently validated scoring method for Ca/Cl-containing α-amylase. Any alternative receptor or scoring model must be labelled and validated independently against its own crystallographic reference and controls. An E128 run is **not** automatically warranted because E64 did not solve the pose-recovery failure and the rank-1 mode is already geometrically repeatable. This 4GQR Vina model should **not** advance to production phytochemical docking or potency ranking on the present evidence. The completed 3L4W experiment and all original structures were preserved.

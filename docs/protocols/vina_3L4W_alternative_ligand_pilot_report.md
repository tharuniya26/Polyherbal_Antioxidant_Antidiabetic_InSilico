# Step 73 — frozen alternative-ligand 3L4W pose-only pilot

## Scope and frozen inputs

This pilot tests neutral rutin (LIG02; *Aegle marmelos* association) and neutral luteolin (LIG08; *Vitex negundo* association) against the two existing dry/no-NAG human ntMGAM 3L4W receptor hypotheses: Asp443−/Asp542−/His600-HID and Asp443−/Asp542H(OD2)/His600-HID. The frozen [Step 72 input manifest](../../data/metadata/project4_alternative_ligand_input_manifest.csv) supplies the exact source structures, one-record ligand PDBQTs, rigid-core maps, randomized starts, receptor paths and SHA256 hashes. The original gallic-acid, agnuside and quercetin docking outputs are unchanged; the Boerhavia association of quercetin is not a new experiment.

The selected rutin conformer set is one-based SDF records **1, 3, 4, 5 and 9**. The actual frozen seed allocation in the manifest is **1101→1, 2202→4, 3303→5, 4404→3, 5505→9**. The task request lists the selected records in ascending order beside the seed list; it also requires using the exact frozen starts and forbids reallocating conformers. This run preserves the manifest allocation. Luteolin uses source record 1 with five independently randomized, frozen starts. The same start per seed is paired with both receptor states.

## Mandatory preflight

The [Step 72 validator](../../scripts/validate_project4_alternative_ligand_inputs.py) passed with 30 manifest rows, five selected rutin records, ten new start files, ten previous enlarged-box native-MIG controls and all 40 earlier quercetin run records verified. Its audit checks all source SDF records against the chosen reference graphs, formulas, charges and 3D stereochemistry; Meeko-to-RDKit PDBQT reconstruction; complete heavy-atom maps; all ten rutin stereocentres; selected conformer coordinates and energies; ligand/receptor/start hashes; and the prior control inputs, outputs and logs. The [Step 73 runner](../../scripts/run_vina_3L4W_alternative_ligand_pilot.py) rechecked the 20 exact PASS receptor/start rows, fixed search settings, frozen files and hashes, Vina 1.2.7 executable/version, native MIG reference and the 10/10 prior positive-control outcomes. It did not repeat native-MIG docking.

All **20/20** exact receptor/start pairs passed the Vina `--score_only` input parse and grid-compatibility check with exit code 0. The [preflight record](../../docking/pilot/vina/3L4W_alternative_ligands_pose_only/vina_pair_parse_preflight.json) contains commands, stdout/stderr, exit codes, durations and hashes. These Vina scores are technical parsing evidence only and do not enter the scientific decision.

## Prespecified search and geometric screen

AutoDock Vina **1.2.7**, `vina` scoring, exhaustiveness **64**, 20 modes, energy range **5.0**, minimum RMSD **1.0 Å**, CPU **1** per job; at most **two** concurrent processes. Five seeds per compound and receptor yield **20 planned searches**. The center is **(45.38664286, 92.35157143, 34.82200000) Å** and the box is **30 × 24 × 24 Å**. Inputs and settings were not revised after viewing results.

For every returned mode, Meeko/RDKit reconstruction checks molecular graph, formal state, complete heavy-atom mapping, stereochemistry and geometry. The individually frozen three-ring flavonoid rigid-core map is used for each compound; rutin sugars remain in all-heavy pocket, clash and box checks. Pocket membership requires all heavy atoms inside the box, rigid-core centroid within **6.0 Å** of deposited MIG and at least half the core atoms within **4.0 Å** of a MIG atom. Severe clash means a ligand/receptor pair below **1.8 Å**, van-der-Waals overlap above **0.8 Å**, or a nonbonded internal heavy-atom pair below **1.6 Å**. A receptor-level geometric **PASS** requires at least **4/5 rank-1 poses** to pass pocket and clash checks and be within **2.0 Å unfitted symmetry-aware rigid-core RMSD** of that compound/receptor state's medoid. Face margin below **1.0 Å** is a warning. Deposited crystallographic waters are a separate distance diagnostic and do not change this frozen rule. Vina scores are retained as output metadata, never used to choose favorable runs.

## Run outcomes and decision

Results pending completion of the frozen searches. The [run table](../../data/metadata/vina_3L4W_alternative_ligand_pilot_results.csv) and [all-pose table](../../docking/pilot/vina/3L4W_alternative_ligands_pose_only/all_pose_metrics.csv) will provide the per-seed and per-mode audit trail, with raw PDBQT outputs, process logs and hashes under the same output directory.

## Scientific boundary

A geometric PASS can mean only reproducible computed placement in this prepared pocket. It does not establish affinity, enzyme inhibition, clinical antidiabetic effectiveness, formulation synergy or detection of either compound in this project's original aqueous leaf extracts. Dry/no-NAG receptor treatment, omitted active-site waters and glycans, Asp542 protonation, neutral ligand microstate assumptions and especially the five-start/16-torsion rutin conformer sample limit interpretation. The source literature associations are separate from extract identification.

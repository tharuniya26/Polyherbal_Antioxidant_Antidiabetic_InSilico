# Vina native-ligand redocking plan and current hold

**Status (12 October 2026): redocking has not begun and is not yet cleared for either receptor.** The [receptor results](../../data/metadata/vina_receptor_preparation_results.csv) show one preliminary NAG-free 2QMJ sensitivity PDBQT and no valid 1XD0 PDBQT. The crystallographic reference SDF files and their original atom maps remain unchanged. The new PDBQTs in [`docking/validation/vina/`](../../docking/validation/vina/) are **crystal-pose references for coordinate and chemistry checks**, never initial search poses for a meaningful redocking test.

## Native ligand identity and conversion

| Receptor | Exact deposited inhibitor | Original reference | New reference PDBQT | Heavy atoms / Vina torsions | Conversion check |
|---|---|---|---|---:|---|
| 2QMJ | Complete linked **alpha-acarbose, PRD_900007** | [SDF](../../docking/validation/2QMJ_PRD_900007_native_pose.sdf) | [crystal-pose PDBQT](../../docking/validation/vina/2QMJ_PRD_900007_crystal_pose_reference.pdbqt) | 44 / 22 | Meeko parse and RDKit molecular graph/stereochemistry roundtrip pass; all mapped heavy coordinates preserved |
| 1XD0 | Deposited **ARE acarbose-derived pentasaccharide**, not intact acarbose | [SDF](../../docking/validation/1XD0_ARE_native_pose.sdf) | [crystal-pose PDBQT](../../docking/validation/vina/1XD0_ARE_crystal_pose_reference.pdbqt) | 55 / 28 | Same checks pass |

Each native SDF contains one explicit-H, three-dimensional record. Meeko 0.8.0 used `MoleculePreparation(add_index_map=True, charge_model='gasteiger')` and `PDBQTWriterLegacy` with an index map. The outputs have one parseable torsion tree and preserve the reference isomeric molecular graph. The deposited SDF heavy-atom index and PDB atom name are joined to the PDBQT serial in the [2QMJ bridge map](../../docking/validation/vina/2QMJ_PRD_900007_sdf_to_pdbqt_atom_map.csv) and [1XD0 bridge map](../../docking/validation/vina/1XD0_ARE_sdf_to_pdbqt_atom_map.csv). All **44/55** mapped heavy-atom coordinates agree within 0.0011 Å. The crystal-pose reference hydrogens, neutral microstate and Meeko charges remain **computational assumptions**; their chemical suitability in the catalytic sites is unreviewed. No alternate inhibitor, including intact acarbose for 1XD0, was substituted.

## Site center and provisional box

The centers below are heavy-atom centroids of the original deposited reference ligands, recorded in [receptor preparation results](../../data/metadata/receptor_preparation_results.csv). The proposed rectangular boxes fit the full ligand heavy-atom extent with approximately 6 Å clearance on each face, rounded to whole Ångströms. They require visual pocket/neighbor inspection once receptor chemistry is ready.

| Receptor | Center x, y, z (Å) | Ligand heavy-atom extent x, y, z (Å) | Proposed Vina size x, y, z (Å) |
|---|---|---|---|
| 2QMJ | −20.808, −6.586, −5.074 | x −25.849 to −17.604; y −9.368 to −3.278; z −13.380 to 2.819 | **23 × 19 × 29** |
| 1XD0 | 10.441, 15.880, 40.721 | x 0.531 to 17.346; y 10.523 to 20.246; z 34.040 to 45.574 | **32 × 23 × 26** |

The 1XD0 x dimension exceeds 30 Å because ARE is extended. The [Vina FAQ](https://autodock-vina.readthedocs.io/en/latest/faq.html) cautions that large boxes require greater search effort. Both ligands are flexible carbohydrates (**22 and 28 Meeko torsions**), so Vina's finite search may miss the crystallographic pose even with chemically correct inputs. Treat box size, exhaustiveness and independent seeds as prospective settings, chosen before assessing RMSD; do not infer validity from a favorable score. Review whether the entire distal ligand is resolved and whether crystal packing contributes to the 1XD0 pose.

## Redocking procedure after the receptor gate

1. Resolve the chemical blockers in [Vina receptor preparation](vina_receptor_preparation.md): linked NAG/Asn in 2QMJ; peptide-linked PCA in 1XD0; Ca/Cl and coordinating-water typing; catalytic protonation and short generated-H contacts. Regenerate and independently validate each receptor PDBQT. The 2QMJ NAG-free variant may serve only as a documented sensitivity arm.
2. Reassess reference-ligand protonation and tautomer state for the selected docking pH without changing the deposited heavy-atom reference. Generate independent starting conformers away from the crystal pose, while maintaining the exact PRD_900007 or ARE graph and stereochemistry. Do not use the saved crystal-pose PDBQTs as starting conformers for pose-recovery claims.
3. Define Vina box coordinates explicitly; visually confirm site and ligand coverage. Predeclare scoring function, receptor and ligand hashes, box, exhaustiveness, seed set, number of modes, and any fixed-water/ion choices. Use multiple independent runs with greater-than-default search effort for these high-torsion ligands. Preserve all outputs and failed runs.
4. Map each redocked **heavy atom** through the stored SDF-to-PDBQT bridge, confirm the output molecular graph and stereochemistry, and calculate full-ligand positional RMSD in the unchanged receptor frame. State any symmetry correction separately. Inspect per-sugar deviations, catalytic contacts, clashes, water/ion interactions and ring conformations. Do not align the ligand to itself before positional RMSD.

No PDBQT reference conversion or future RMSD alone verifies receptor chemical correctness, inhibition, botanical occurrence in the original aqueous extracts, or DPPH/ABTS antioxidant behavior. See the [prior native-ligand plan](native_ligand_redocking_plan.md) for the deposited atom-name mapping and crystallographic caveats.

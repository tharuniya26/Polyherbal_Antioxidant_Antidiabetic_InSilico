# Meeko receptor-template blocker resolution: 2QMJ and 1XD0

**Assessment date:** 12 October 2026. **Scope:** Meeko 0.8.0 receptor chemistry and PDBQT quality control before AutoDock Vina 1.2.7 native-ligand redocking. **No redocking, scoring, production docking, or ADMET was performed.** The [results table](../../data/metadata/vina_receptor_blocker_resolution_results.csv) and [reproduction script](../../scripts/resolve_vina_receptor_templates.py) provide source hashes, variant decisions, settings and output hashes. The previous [Vina preparation assessment](vina_receptor_preparation.md) and [native-ligand plan](vina_native_redocking_plan.md) remain historical records.

## Meeko template mechanism and provenance

The installed `mk_prepare_receptor` CLI accepts additional residue-template JSON through `--add_templates`; a linked residue needs chemically defined atom names, a stereodefined graph, and a **link label** that matches the inter-residue bond. Automatic CCD lookup in Meeko 0.8.0 cannot resolve every modified backbone or glycan linker. The documented [Meeko template construction guide](https://github.com/forlilab/Meeko/blob/develop/docs/source/py_build_temp.rst) and [CLI options](https://github.com/forlilab/Meeko/blob/develop/docs/source/rec_cli_options.rst) describe this path. No `--delete_residues`, `--delete_bad_res`, residue substitution, or unreviewed SDF-to-template shortcut was used. The installed Meeko default templates contain CA, CL and HOH, but no PCA or NAG.

The script verifies the SHA256 of the original 1XD0 PDB, all five previous GOLD-stage receptor inputs and the earlier Vina PDBQT before processing. Its PCA source is the official [RCSB PCA chemical-component definition](https://files.rcsb.org/ligands/download/PCA.cif), saved byte-for-byte as [`PCA_RCSB_CCD_2024-09-27.cif`](../../docking/configs/vina_templates/PCA_RCSB_CCD_2024-09-27.cif) (SHA256 `2c2aa7b8ee9b3acc6ad5ec66bb88f4701c7260f3ec60681799a1b8ec68366238`). A checksum change stops the script for manual review. The custom [`PCA_C_terminal_Meeko_0.8.0.json`](../../docking/configs/vina_templates/PCA_C_terminal_Meeko_0.8.0.json) is generated from that component; it contains no invented bond, stereocentre, partial charge or atom type. Meeko assigns its usual receptor Gasteiger/template partial charges during preparation. The Vina executable reported **1.2.7** with `--version`; it was not used to score or dock.

## 2QMJ: linked NAG glycans remain blocked

The [deposited 2QMJ PDB](https://files.rcsb.org/download/2QMJ.pdb) `LINK` records identify **all four** NAG units as members of protein-linked glycans:

| NAG unit | Deposited covalent link | Site relevance |
|---|---|---|
| A2005 | C1–Asn209 ND2, 1.46 Å | About **3.895 Å** from the complete alpha-acarbose inhibitor; relevant to pocket shape |
| X1 | C1–Asn393 ND2, 1.44 Å | First unit of an Asn-linked disaccharide |
| X2 | C1–X1 O4, 1.44 Å | Second unit of the same protein-linked glycan |
| A2001 | C1–Asn741 ND2, approximately 1.45–1.46 Å in the two deposited alternates | Glycan attached to alternate Asn741; prior selection kept alternate A |

Thus the chain-X sugars are **not free crystallographic NAG ligands**. Meeko's default template set has no NAG, and the CLI rejected the retained NAG variants as unknown linking fragments. A valid replacement would need the RCSB NAG stereochemistry and site-specific donor/acceptor states, explicit Asn–NAG and NAG–NAG linker definitions, hydrogen placement, atom typing and charge validation, especially at A2005. Those linked glycan states were not adequately established here, so **no NAG template or NAG-retaining PDBQT was created**. The two prior glycan-retaining PDBs remain coordinate controls. The existing [NAG-free, site-water 2QMJ PDBQT](../../structures/proteins/vina_prepared/2QMJ_sitewaters_noNAG_preliminary_sensitivity.pdbqt) was verified unchanged and remains **only a sensitivity model**; glycan absence cannot be treated as the deposited biological receptor.

The prior pH 7 OpenMM model gives Asp443/Asp542 deprotonated carboxylates and His600 ND1-H; this is a provisional protein-only microstate, not an established catalytic proton-transfer state. The missing Gln837 and termini are distant from the inhibitor but remain model-completeness caveats. Ten deposited waters were retained in the sensitivity variant; HOH4627's placed H1 lies **1.553 Å** from a native-ligand heavy atom and needs orientation/protonation review. No water or glycan was removed in this resolution stage.

## 1XD0: peptide-linked PCA template constructed and checked

The [deposited 1XD0 PDB](https://files.rcsb.org/download/1XD0.pdb) records PCA A1 **C–Tyr2 N** as a 1.36 Å peptide link, and its PCA residue has the eight heavy atoms `N, CA, CB, CG, CD, OE, C, O`. The [RCSB PCA component](https://files.rcsb.org/ligands/download/PCA.cif) defines the N-terminal pyroglutamate ring, **S** configuration at CA and carboxyl-terminal OXT/HXT as leaving atoms for peptide formation. The template removes only those CCD terminal leaving atoms, retains the deposited eight-atom ring and carbonyl graph, and labels its carbonyl **C** as a Meeko `C-term` link. It does **not** replace PCA with glutamate, glutamine or a detached ligand. Independent checks find CA **S** in the CCD ideal geometry, the custom template and the 1XD0 deposited geometry. Meeko loads the JSON and recognizes the covalent terminal-residue context.

Both existing 1XD0 GOLD-stage variants converted without deleting PCA:

| New PDBQT | Waters | Heavy atoms retained | Parsed PDBQT atoms | Status |
|---|---:|---:|---:|---|
| [All-water PCA/CA/CL variant](../../structures/proteins/vina_prepared/1XD0_allwaters_CA_CL_PCA_template_Meeko0.8.0.pdbqt) | 142 | 4,090 / 4,090 | 5,266 | Chemically prepared input for controlled review; not docking validated |
| [Site-water PCA/CA/CL variant](../../structures/proteins/vina_prepared/1XD0_sitewaters_CA_CL_PCA_template_Meeko0.8.0.pdbqt) | 6 | 3,954 / 3,954 | 4,858 | Same |

Meeko PDBQT read-back found finite coordinates/charges, expected atom types, the eight PCA heavy atoms, and a PCA partial-charge sum of approximately **0.000**. All input heavy-atom **positions** were preserved within 0.0011 Å. Meeko exchanges the names of the two equivalent **Leu496 terminal carboxylate O/OXT** atoms while retaining both positions; the QC script handles this explicitly by coordinate assignment and reports it in the CSV. It is remote from the ARE site. The custom template has S stereochemistry; no PCA stereocentre was silently inverted. Nonpolar hydrogens are merged by Meeko, while polar hydrogens are represented in PDBQT. The input hydrogen placements came from the earlier provisional OpenMM pH 7 step.

**Calcium and chloride:** both outputs retain CA A498 as atom type `Ca` with charge **+2.000** and CL A500 as `Cl` with charge **−1.000**. Ca-coordinating HOH A536 and A546 remain in both variants; the all-water model retains all 142 deposited waters and the site-water model retains six. The deposited calcium shell includes Asn100 OD1 (2.43 Å), Arg158 backbone O (2.42 Å), Asp167 OD1/OD2 (2.46/2.51 Å), His201 backbone O (2.48 Å), and HOH536/546 (2.53/2.62 Å). These coordinates were preserved. PDBQT type/charge support does **not** prove the chosen Vina scoring function models metal coordination or water thermodynamics faithfully. Remote NAG A499 was removed only in the previous GOLD-stage variants, at about 47.4 Å from ARE; the earlier working receptor still preserves it.

## Active-site quality and limits

The nearest input-receptor heavy atom to deposited ARE is **2.551 Å** and nearest water oxygen **2.859 Å** in each 1XD0 variant. The nearest input polar H is **1.944 Å** in the all-water and **1.960 Å** in the site-water variant; these distances warrant contact inspection but do not by themselves prove an invalid clash. Asp197, Glu233 and Asp300 were assigned deprotonated carboxylates and His201 ND1-H by the prior ligand-free pH 7 OpenMM preparation. Their catalytic protonation and neighboring water orientations remain uncertain, and no protonation variant was generated here. The earlier crystal-symmetry audit found ARE O3A/C3A contacts to a crystal-mate Asp433/Trp434 within 4 Å. A single-chain receptor cannot represent that contact, so full ARE pose recovery will have an interpretation limit.

## Redocking gate

The **PCA template blocker is resolved for 1XD0** at the level of a Meeko-prepared, coordinate-checked receptor PDBQT. The **NAG blocker remains unresolved for deposited 2QMJ**. Neither receptor is fully docking validated. Before a controlled 1XD0 redocking test, predeclare the catalytic acid/His microstate, water variant, Ca/Cl treatment, box and search settings; inspect the native ARE contacts and check that Vina 1.2.7 imports the selected PDBQT without atom-type warnings. The existing 2QMJ NAG-free receptor may only be used for an explicitly labelled sensitivity test after resolving its water-H and catalytic-state concerns. No receptor is marked unconditionally redocking-ready in the results table.

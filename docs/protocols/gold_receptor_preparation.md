# Staged receptor preparation for GOLD native-ligand redocking

**Status (11 October 2026):** Five derived receptor variants were generated and structurally checked. **None has passed a GOLD/Hermes import test or native-ligand redocking.** The NAG-retaining 2QMJ variants are intentionally unprotonated. The other variants have provisional OpenMM pH 7 hydrogens; their catalytic microstates and GOLD atom types remain unverified. See [variant results](../../data/metadata/gold_receptor_preparation_results.csv) and the reproducible [staging script](../../scripts/stage_gold_receptors.py).

## Input integrity and tools

Before any editing, both ligand-free [working PDBs](../../structures/proteins/prepared/) and both native-pose [SDFs and atom maps](../../docking/validation/) existed and matched every corresponding SHA-256 value in [receptor_preparation_results.csv](../../data/metadata/receptor_preparation_results.csv). The script rechecks these values and checks that each reference SDF heavy atom agrees with its mapping coordinates. It never writes to those inputs or the original RCSB files. The original reference-pose centroids remain **2QMJ (−20.808, −6.586, −5.074) Å** and **1XD0 (10.441, 15.880, 40.721) Å**.

Tools found in the project virtual environment were **PDBFixer 1.12.0** for missing-atom audit and **OpenMM 8.6.1** for provisional hydrogen addition. The script uses OpenMM's Reference platform and seed `20261011`; rerunning reproduced the checked protonated output byte for byte. RDKit checks the reference SDF/map and Gemmi remains available from the previous structural audit. No `gold`, `gold_auto`, or `hermes` executable was on the command path or in the inspected local program directories, and the CCDC Python API was absent. Consequently **no actual GOLD import or SYBYL atom-typing test was possible**. CCDC states that GOLD needs explicit hydrogens and correct protein/ligand atom and bond types; a readable PDB file is insufficient ([CCDC setup guidance](https://support.ccdc.cam.ac.uk/support/solutions/articles/103000306440/)).

## Structural inspection

PDBFixer found **no missing standard side-chain heavy atoms** in either working receptor. For 2QMJ it detected missing construct residues **1–6 and Gln837** and a missing terminal **OXT** on His870; the termini and Gln837 gap are more than 50 Å from the native inhibitor. For 1XD0 it found no missing protein residues or standard side-chain heavy atoms. PDBFixer identifies glycosylated Asn as nonstandard due to covalent NAG linkage and would propose replacing the modified N-terminal **PCA A1** in 1XD0 with Gln. **Neither change was accepted**, because both would erase deposited chemistry. PDB `LINK`/`CONECT` records and OpenMM topology were inspected for peptide, glycan, disulfide and ion-related connectivity.

2QMJ has six equal-occupancy alternate-conformation residues: Arg68, Cys479, Glu661, Asn741, Lys762 and NAG A2001. The deterministic rule is **highest mean occupancy for the entire residue; choose alphabetical label on a tie**. All six choose **A** at occupancy 0.50. B atoms and B-specific NAG `LINK` records are omitted; A atom coordinates and original occupancy values remain unchanged. No alternate locations were present in 1XD0. In 2QMJ an explicit chain break is inserted between residues 836 and 838: without it OpenMM inferred an impossible **4.13 Å C836–N838 peptide bond** across missing Gln837. The break creates artificial segment ends, so it is a local docking approximation, not a repaired full protein. The deposited Cys26–Cys44 `SSBOND` has a **2.916 Å SG–SG separation**; this unusual original geometry was retained and warrants density/structural review. No retained non-hydrogen topology bond exceeded 2.2 Å in 1XD0; in 2QMJ the deposited Cys26–Cys44 disulfide is the remaining exception. The working files preserve all retained heavy-atom names and coordinates exactly.

## Heteroatom and water decisions

Distances below are minimum heavy-atom distances to the **unchanged native reference ligand**. A 4.0 Å cutoff selects waters for the site-water variants; it is a screening rule for sensitivity models, not proof that other waters are irrelevant. No waters or heteroatoms were removed from the original structures.

| Receptor | Category | Decision and reason |
| --- | --- | --- |
| 2QMJ | Ten site waters: HOH A4627, A4081, A4623, A4625, A4626, A4624, A4135, A4473, A4450, A4287 | Retain in site-water variants; each is within 4 Å of alpha-acarbose. Retain **all 626** in the all-water coordinate control; remove the other 616 only in the site-water variants. Individual mediation/displacement remains untested. |
| 2QMJ | NAG glycosylation, including NAG A2005 at **3.895 Å** | Retain all deposited NAG in two coordinate controls. Create one **NAG-free sensitivity** variant because OpenMM cannot correctly protonate linked NAG/Asn. NAG removal may change the local pocket and is **not** the default biological claim. |
| 2QMJ | Nine GOL; sulfate SO4 A4001 | Remove from derived variants: nearest GOL is **13.864 Å**, sulfate **49.365 Å** from the inhibitor. These remote crystallization species are not site cofactors in this structure. |
| 1XD0 | Four site waters: HOH A588, A766, A511, A634 | Retain in site-water variant because each is within 4 Å of ARE. Retain all **142** in the all-water variant. |
| 1XD0 | HOH A536 and A546 | Retain in **both** variants despite being outside the 4 Å ARE cutoff: both coordinate Ca A498. The site-water variant therefore has **six** waters, and the other 136 are omitted. |
| 1XD0 | Ca A498 and Cl A500 | Retain both. Calcium's coordination shell and the deposited chloride are part of the amylase structural environment; GOLD metal/anion atom typing remains to be verified ([CCDC metal guidance](https://support.ccdc.cam.ac.uk/support/solutions/articles/103000305962/)). |
| 1XD0 | PCA A1; NAG A499 | Retain peptide-linked PCA and its original heavy coordinates. Remove remote NAG A499 from derived sensitivity models; its nearest ARE distance is **47.411 Å**. The parent receptor still preserves it. |

## Hydrogen and protonation status

| Variant | Water count | H status | Interpretation |
| --- | ---: | --- | --- |
| `2QMJ_gold_stage_allwaters_keepNAG_unprotonated.pdb` | 626 | None added | Coordinate control preserving deposited glycans and all waters; **not chemically prepared**. |
| `2QMJ_gold_stage_sitewaters_keepNAG_unprotonated.pdb` | 10 | None added | Site-water coordinate control with glycans; **not chemically prepared**. |
| `2QMJ_gold_stage_sitewaters_noNAG_pH7.pdb` | 10 | 6,603 protein/water H added | NAG-free sensitivity model; provisionally protonated, **not GOLD validated**. |
| `1XD0_gold_stage_allwaters_CA_CL_pH7.pdb` | 142 | 4,002 protein/water H added | Provisional pH 7 model with Ca, Cl, PCA and all waters; **not GOLD validated**. |
| `1XD0_gold_stage_sitewaters_CA_CL_pH7.pdb` | 6 | 3,730 protein/water H added | Provisional pH 7 site-water model with Ca, Cl and PCA; **not GOLD validated**. |

OpenMM used its standard pH 7 assignment for the NAG-free variants: **2QMJ Asp443 and Asp542** and **1XD0 Asp197, Glu233 and Asp300** have no side-chain proton; catalytic proton-transfer states were **not established** by this choice. **2QMJ His600** and **1XD0 His201** were assigned ND1-H (`HD1` present, `HE2` absent). These are protein-environment estimates made without the native inhibitor present during hydrogen placement. The original ligand's crystallographic heavy-atom coordinates were never touched. The 2QMJ NAG-retaining variants have no hydrogens because OpenMM adds none to NAG and adds chemically excessive amide hydrogens to NAG-linked Asn; writing that result would falsely imply a valid glycoprotein. Glycan protonation requires a glycoprotein-aware tool or manual chemically checked correction.

For the provisional NAG-free models, a fixed-reference ligand comparison found the closest generated receptor H to a native ligand heavy atom at **1.553 Å** for 2QMJ (HOH A4627 H1) and **1.960 Å** for 1XD0 (Gln63 HE22). The short 2QMJ water-H approach needs orientation review or GOLD water optimization before redocking; it is not evidence of a physical clash in the original X-ray model, which lacks observed hydrogen positions. No deposited heavy-atom pose was minimized. The 1XD0 crystal contact to ARE O3A/C3A described in the [receptor assessment](receptor_structure_assessment.md) still limits isolated-chain pose-recovery interpretation.

## Remaining gates

1. Inspect the deposited validation/electron-density information around the 2QMJ NAG A2005, Cys26–Cys44 region, distal acarbose sugars, catalytic residues and selected waters; review the analogous 1XD0 ARE and ion sites. The present geometric audit cannot assess local density quality.
2. In a glycoprotein-aware preparation tool, add chemically valid H to the retained 2QMJ NAG units and their linked Asn, or document why a NAG-free receptor is adequate after a sensitivity test. Consider the remote 2QMJ sequence gap/terminal treatment without changing the catalytic cleft.
3. Review catalytic Asp/Glu and histidine microstates with the native inhibitor and water network present. Check the 2QMJ HOH4627 H orientation and any other short generated-H contacts. Verify PCA, Ca and Cl typing and metal coordination in 1XD0.
4. Import each candidate receptor and a **separate non-native starting conformer** of the correct inhibitor into the installed GOLD/Hermes version. Inspect the import log, atom/bond typing, explicit H, NAG/PCA/ion recognition, water behavior, site definition and search settings. Record failures. No variant may be called GOLD-compatible until this succeeds.

The [GOLD setup plan](../../docking/configs/gold_redocking_setup.md) records prospective settings. **Native-ligand redocking cannot yet begin as a validated experiment in this environment.**

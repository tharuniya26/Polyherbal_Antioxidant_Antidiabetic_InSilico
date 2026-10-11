# Proposed GOLD setup for native-ligand pose-recovery tests

**Configuration status:** Planning record only. No GOLD configuration file has been executed, no receptor has passed a GOLD import test, and no docking/RMSD result exists. The actual installed GOLD version and supported option names must be captured at import time. Follow the [receptor preparation assessment](../../docs/protocols/gold_receptor_preparation.md) and the existing [redocking validation plan](../../docs/protocols/native_ligand_redocking_plan.md).

## Identity, pose and site anchors

| Receptor | Correct native ligand | Immutable reference pose and map | Heavy-atom centroid (Å) | Furthest deposited heavy atom from centroid | Provisional sphere radius |
| --- | --- | --- | --- | ---: | ---: |
| 2QMJ | Whole linked alpha-acarbose `PRD_900007`, 44 heavy atoms | [SDF](../validation/2QMJ_PRD_900007_native_pose.sdf), [atom map](../validation/2QMJ_PRD_900007_atom_map.csv) | **(−20.808, −6.586, −5.074)** | 9.15 Å | **12 Å**, subject to visual site review |
| 1XD0 | ARE acarbose-derived pentasaccharide, 55 heavy atoms; **not intact acarbose** | [SDF](../validation/1XD0_ARE_native_pose.sdf), [atom map](../validation/1XD0_ARE_atom_map.csv) | **(10.441, 15.880, 40.721)** | 10.10 Å | **13 Å**, subject to visual site review |

The provisional radii provide about 3 Å beyond the observed full-ligand extent. In GOLD/Hermes, display the entire native ligand and catalytic site before accepting these numbers. CCDC allows site definition from a point, atom, reference ligand, or residue set ([official site guidance](https://support.ccdc.cam.ac.uk/support/solutions/articles/103000306183-how-do-i-define-my-active-site-)). If a single sphere captures excess space around 2QMJ's distal sugars or misses relevant catalytic contacts, revise and document the site **before** looking at redocking results. Record final center, radius, atom/residue coverage and a screenshot or exported site definition.

## Candidate receptor inputs and water sensitivity

Use the new [staged receptor files](../../structures/proteins/prepared/) only after the remaining chemistry and import checks. For 2QMJ, the **NAG-retaining** all-water and 10-site-water models are coordinate-only and need glycoprotein protonation. The 10-water **NAG-free** pH 7 model is a sensitivity candidate; its omission of nearby NAG A2005 (3.895 Å from the inhibitor) must be considered when interpreting pose recovery. For 1XD0, compare the **142-water** and **six-water** pH 7 models, keeping Ca A498, Cl A500 and calcium-coordinating HOH A536/A546 in both. The selected four other site waters are HOH A588/A766/A511/A634. The parent structures and native SDFs remain unchanged. Water retention and optimization should follow explicit GOLD options and be logged; CCDC notes that removing mediating waters can prevent crystal-pose recovery ([CCDC guidance](https://support.ccdc.cam.ac.uk/support/solutions/articles/103000306441-why-is-it-not-possible-to-return-the-crystallographic-binding-pose-of-the-ligand-)).

## Proposed trial settings, pending GOLD availability

- **Input ligand:** Generate an independent 3D starting conformer from the verified *same* connected graph and stereochemistry. Record seed and conformer method. Do not use the immutable native-pose SDF coordinates as the initial pose. Preserve glycosidic linkages, ring stereochemistry and appropriate protonation. Use multiple independent starting conformers for these flexible carbohydrates.
- **Search:** Begin with a full ligand-flexibility search, approximately **30 independent runs/poses** per receptor–variant–starting-conformer combination, and GOLD's long/high-effort search option where available. Verify the installed version's actual efficiency controls before writing a runnable configuration. Do not impose native-pose-derived distance constraints or covalent docking.
- **Scoring:** Propose **ChemPLP** for an initial pose-recovery test and **GoldScore** as a sensitivity analysis, subject to the installed license/version. Scores are optimization outputs, not measured affinity or inhibition. Record all parameters and seeds before viewing poses.
- **Waters:** Compare all-water and reduced-water receptors where chemically prepared; review each selected site water and specify whether fixed, rotatable or toggleable in GOLD. Keep the 1XD0 calcium-coordination waters with calcium. A NAG-free 2QMJ outcome cannot alone validate the NAG-containing pocket.
- **Protein/ligand preparation:** Check all H, histidine tautomers, catalytic carboxylates, ion/metal coordination, glycosylation/PCA chemistry and atom/bond types in GOLD/Hermes. Do not proceed if the import log reports an untyped atom, broken bond, missing H or unsupported heteroatom.
- **Acceptance assessment:** Save every output pose and log. Before runs, set a pose-recovery criterion in heavy-atom RMSD terms, using the existing atom maps and reporting both full-ligand and diagnostic regional RMSD. Inspect the catalytic contact network, stereochemistry, ring conformation, water interactions and 1XD0 ARE crystal-contacted region. Report failures and any protocol changes transparently; no criterion has yet been tested.

## Stop conditions

Do not launch GOLD while the required version/import test is unavailable, while 2QMJ retained NAG remains unprotonated, or while receptor/ligand atom typing and site coverage are unverified. The current reference SDFs support future RMSD calculation, but they do **not** constitute docking validation. Production phytochemical docking and ADMET remain outside this stage.

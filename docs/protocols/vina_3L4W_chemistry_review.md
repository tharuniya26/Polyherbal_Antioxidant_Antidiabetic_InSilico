# 3L4W/MIG chemistry review for diagnostic Vina redocking

**Assessment date:** 12 October 2026. **Scope:** human N-terminal maltase-glucoamylase (ntMGAM), [PDB 3L4W](https://www.rcsb.org/structure/3L4W), with its deposited miglitol (MIG) pose. This is an input-chemistry assessment, not a redocking result or a claim of biological validation. No docking search, score, RMSD, ADMET calculation, or phytochemical prediction was performed.

## Source integrity and methods

The audited official `3L4W.pdb` and `3L4W.cif`, [MIG CCD definition](https://files.rcsb.org/ligands/download/MIG.cif), existing native SDF/PDBQT, atom maps, and existing local receptor are preserved. Their source SHA256 values match the [prior preparation audit](../../data/metadata/vina_alternative_receptor_preparation_results.csv). The [new script](../../scripts/review_3l4w_chemistry.py) refuses changed source hashes and writes each new variant under a new filename. It uses RDKit 2026.03.6, Meeko 0.8.0, and the installed Vina 1.2.7 executable. Structural distances below were measured from deposited heavy atoms; Gemmi 0.7.3 was used for symmetry-neighbor inspection. Existing receptor hydrogens were added by OpenMM 8.6.1 in the prior preparation and are *models*, not crystallographic assignments.

The preparation pH is **6.5**, the [deposited crystallization pH](../../structures/proteins/original/3L4W.pdb) (`REMARK 280`). It is an operational preparation value, not a measured pH for intestinal fluid or the later docking environment. The PDB model is at 2.00 Å and does not resolve hydrogen positions or site-specific protonation. The [primary structural study](https://doi.org/10.1021/bi9016457) and [RCSB's binding description](https://pdb101.rcsb.org/global-health/diabetes-mellitus/drugs/alpha-glucosidase-inhibitors/drug/miglitol/miglitol) identify Asp443 as nucleophile, Asp542 as general acid/base, and contacts involving His600, Asp327 and Arg526.

## Missing protein segment and glycans

`REMARK 465` identifies **Gln837 as unobserved**. Its flanking modeled residues Ser836 and Thr838 are, respectively, **54.145 and 52.722 Å from MIG**, and at least **48.720 and 46.983 Å from Asp443/Asp542/His600**. The Ser836 C–Thr838 N gap is **4.343 Å**, so it is not an ordinary peptide bond. The working receptor keeps a chain break at this position. No Gln837 coordinates or side-chain atoms were invented. The artificial Ser836/Thr838 termini can alter global electrostatics, even though the gap is far from the site. Residues 1–6 and 871–875 are also unobserved; the reported D772N sequence variation remains as deposited. PDBFixer found no missing standard side-chain atoms in the prior local model; Ser836 OXT was absent at the artificial break.

All **six omitted NAG units are covalently linked glycans**, not independent free sugars. Deposited `LINK` records connect Asn209–NAG A2005, Asn372–NAG A2006, Asn393–NAG C1–C2, and Asn741–NAG B1–B2. Four protein glycosylation sites are also identified in the [RCSB entry](https://www.rcsb.org/structure/3L4W). Each NAG has 14 deposited heavy atoms rather than the isolated CCD's 15; the linkage chemistry must therefore be assessed as a glycoside, not inferred by treating six isolated NAGs as complete monosaccharides. The closest deposited heavy-atom distances are:

| NAG | Protein linkage | To MIG (Å) | To Asp443/Asp542/His600 (Å) |
| --- | --- | ---: | ---: |
| A2005 | Asn209 | 16.409 | 10.316 |
| A2006 | Asn372 | 19.306 | 21.458 |
| B1 | Asn741 | 29.112 | 23.568 |
| B2 | B1 O4–C1 | 31.744 | 26.420 |
| C1 | Asn393 | 25.994 | 19.960 |
| C2 | C1 O4–C1 | 30.701 | 24.590 |

The closest NAG does not directly contact MIG, but A2005 is close enough to the catalytic region that a conformational or electrostatic effect cannot be excluded by distance alone. No validated Meeko 0.8.0 template for all four protein-linked glycans was available in this work. The **no-NAG receptor is an explicitly local sensitivity model**, with Asn chemistry different from the glycosylated crystal. The four deposited glycerols are also absent from this model; the nearest glycerol atom is 13.809 Å from MIG. A full glycan-retaining receptor remains a separate preparation and validation task.

## Protonation evidence and decisions

The existing OpenMM pH-6.5 model had Asp443⁻, Asp542⁻ and **doubly protonated His600** (HD1 plus HE2). The histidine assignment arose from an automatic hydrogen-addition algorithm and is not established by X-ray diffraction. His600 NE2 is **2.851 Å from MIG OAD**; Asp443 OD2 is **2.758 Å from MIG NAN**, and Asp542 OD1 is **2.629 Å from MIG OAC**. These are heavy-atom contacts, not proof of hydrogen-bond direction or proton transfer.

PROPKA 3.5.1 was installed in the project virtual environment for this review. Its [documented method](https://propka.readthedocs.io/en/stable/command.html) predicts pKa values from a 3D structure; it does not experimentally determine microstates. Two runs used the same deposited protein heavy coordinates, with alternate conformations already resolved by the prior occupancy rule, no NAG or waters, and either absent or neutral CCD-reference MIG. To reproduce the inputs, retain the `ATOM` lines whose element is not H from `3L4W_validation_local_noNAG_pH6p5.pdb`; for the second input append the 14 original `HETATM` lines for MIG A1001 from `3L4W.pdb`. End each file with `END`. Run `.venv/Scripts/propka3.exe -o 6.5 -q <input.pdb>` from a temporary directory. No resulting pKa file was written into this repository. Both runs warned about the artificial Thr838 terminus; this is remote from MIG. Values are:

| Group | Protein only, predicted pKa | Protein + neutral deposited MIG, predicted pKa | Interpretation at pH 6.5 |
| --- | ---: | ---: | --- |
| Asp443 | 5.19 | 3.41 | Deprotonated is the stronger working hypothesis. |
| Asp542 | 7.76 | 5.59 | Ligand inclusion reverses the pH-6.5 implication; **conditional** deprotonated working state, with protonated Asp542 reserved for later sensitivity work. |
| His600 | 0.83 | 0.01 | Both runs oppose the doubly protonated default; a **neutral** histidine is the working hypothesis. These unusually low estimates are model dependent. |

PROPKA's protein-plus-MIG run also emitted an apparent MIG NAN pKa of **13.84** from a PDB input without CCD bond-order annotation. That value was **rejected as a reliable miglitol solution pKa**; ligand typing and the binding environment are insufficient to use it to fix amine charge. [PubChem CID 441314](https://pubchem.ncbi.nlm.nih.gov/compound/441314) reports a tabulated miglitol pKa of **5.9**, sourced to DrugBank and *The Merck Index*; this was not independently confirmed here as a primary potentiometric measurement and may not describe the bound-site pKa. **If** it represents the amine conjugate-acid pKa, pH 6.5 would favor the **neutral** species in bulk solution, while a protonated fraction remains plausible and binding can shift the balance. The [MIG CCD](https://www.rcsb.org/ligand/MIG) and deposited pose encode the neutral formula **C8H17NO5**, formal charge 0, but do not experimentally determine the amine's bound-state protonation.

The locally favorable neutral His600 tautomer is **HID** (HD1 retained, NE2 unprotonated). In the existing neutral MIG reference, the modeled MIG OAD hydrogen is **1.889 Å from His600 NE2**, compatible with a donor-to-acceptor hypothesis. An initial HIE test gave a **0.940 Å H–H overlap** between His HE2 and the modeled MIG OAD hydrogen; that test failed the local hydrogen check and its new files were removed. Ligand hydroxyl hydrogen positions are computational, so HID is an explicit *working hypothesis*, not a measured tautomer. Asp443⁻/Asp542⁻/His600-HID was used for the new variants; Asp542 protonation and alternate His tautomer remain sensitivity questions. No original heavy atom moved.

## MIG charge and stereochemical variants

The original native MIG SDF and PDBQT remain unchanged as the **neutral CCD-reference pose** (C8H17NO5, charge 0, 14 heavy atoms, four defined carbon stereocentres, eight Meeko torsions). A separate [protonated-amine SDF](../../docking/validation/vina/3L4W_MIG_NHplus_crystal_pose_sensitivity.sdf) and [PDBQT](../../docking/validation/vina/3L4W_MIG_NHplus_crystal_pose_sensitivity.pdbqt) were generated from that SDF: C8H18NO5⁺, formal charge +1, 14 heavy atoms, eight Meeko torsions. RDKit and Meeko readback confirmed the original four carbon stereocentres, molecular connectivity, charge and unchanged deposited heavy coordinates to 0.001 Å. Protonating the tertiary amine creates an **unassigned nitrogen pyramidal stereocentre** in RDKit; no nitrogen configuration was invented. Added hydrogens and their orientations are preparation hypotheses; the charged SDF is **not** the deposited CCD microstate. The original SDF-to-PDBQT atom mapping remains the reference for its neutral pose; the charged PDBQT contains its own Meeko `INDEX MAP` records and must be mapped explicitly if later used for RMSD.

## Waters, crystal contacts and local quality checks

Of 434 deposited waters, HOH A4002 and A4004 lie within 4 Å of MIG (**2.719 and 2.907 Å**, respectively). A4002 is **2.482 Å** from Asp443 OD2; A4004 is **2.238 Å** from Asp571 OD2. Their deposited occupancies are 1.00 and B factors are 29.00 and 21.96 Å². These positions support plausible site-water roles but not fixed orientation or persistence in solution. The prior OpenMM-added water O–H lengths are **1.174–1.200 Å**, so the two-water model is a sensitivity case. Ordinary rigid-receptor Vina does not model water exchange; the [Vina hydrated-docking protocol](https://autodock-vina.readthedocs.io/en/stable/docking_hydrated.html) requires separate treatment. A dry local model is the primary diagnostic baseline and a two-water model is retained solely to measure dependence on this choice. Neither is a validated solvation model.

Gemmi symmetry-neighbor search used the deposited **P 21 21 21** cell and all symmetry mates. No crystal-mate protein heavy atom was within 12 Å of MIG; the nearest found across the full ligand was **16.793 Å** (MIG OAA to symmetry-related Thr804 O). Thus no direct crystal-mate MIG contact was detected, although lattice packing may still affect the receptor conformation. The crystallographic MIG pose and deposited cell were not edited.

The new HID dry and two-water receptor PDBQTs passed Meeko readback with finite coordinates/charges and retained all **6,912 protein heavy-atom coordinates** to 0.001 Å; the two-water model additionally retains two water oxygens. His600 is typed as ND1 donor (`N`, HD1 retained) and NE2 acceptor (`NA`); both catalytic Asp side chains remain deprotonated in these variants. Summed Meeko PDBQT partial charges are **−1.000** for each Asp443/Asp542 residue and approximately **0** for His600; the neutral and protonated MIG PDBQTs sum to approximately **0** and **+1**, respectively. These bookkeeping checks do not establish Vina's physical treatment of electrostatics. The closest MIG–protein heavy pair is **2.538 Å** (MIG OAE–Asp327 OD2). For the neutral and charged ligand SDFs, the shortest computed receptor–ligand pair involving a hydrogen is **1.669 Å** in the dry model and **1.614 Å** in the two-water model; no sub-1.2 Å overlap was found. This is a distance screen, not an energy minimization or a complete hydrogen-bond validation.

## New inputs and Vina parse preflight

| Variant | Receptor PDB | Meeko PDBQT | Intended use |
| --- | --- | --- | --- |
| D443⁻/D542⁻/H600-HID, dry, no NAG | [PDB](../../structures/proteins/prepared/3L4W_noNAG_D443m_D542m_H600HID_pH6p5_dry_sensitivity.pdb) | [PDBQT](../../structures/proteins/vina_prepared/3L4W_noNAG_D443m_D542m_H600HID_pH6p5_dry_sensitivity_Meeko0.8.0.pdbqt) | Primary **local-pocket diagnostic** input; paired first with the unchanged neutral MIG reference. |
| Same catalytic states, HOH4002/4004 retained | [PDB](../../structures/proteins/prepared/3L4W_noNAG_D443m_D542m_H600HID_pH6p5_sitewaters_sensitivity.pdb) | [PDBQT](../../structures/proteins/vina_prepared/3L4W_noNAG_D443m_D542m_H600HID_pH6p5_sitewaters_sensitivity_Meeko0.8.0.pdbqt) | Static-water sensitivity only. |

Vina 1.2.7 has no ordinary parse-only CLI flag. For **each of the two receptors × neutral/protonated MIG PDBQT** combinations, the script passed the real files to Vina with an intentionally zero-size box. Its [1.2.7 CLI source](https://raw.githubusercontent.com/ccsb-scripps/AutoDock-Vina/v1.2.7/src/main/main.cpp) calls `set_receptor` and `set_ligand_from_file` before `compute_vina_maps`, where the zero-size box raises `Grid box dimensions must be greater than 0 Angstrom`. All four probes reached exactly that guard, returned exit code 1, and created no pose, map, or score output. This is **parser evidence only**; no valid grid was calculated and no docking or scoring path was entered. For a future real calculation, restore a scientifically defined positive-size box centered on the deposited ligand, not this guard box.

New-file SHA256 hashes are recorded in the [results CSV](../../data/metadata/vina_3L4W_chemistry_review_results.csv). The audit remains source-pinned and rerunnable with:

```powershell
.\.venv\Scripts\python.exe scripts/review_3l4w_chemistry.py
```

## Decision and remaining gates

**A controlled diagnostic pose-recovery test can responsibly be attempted with the dry HID receptor and the unchanged neutral MIG reference**, with the separately prepared MIG-NH⁺ and two-water models predeclared as sensitivity cases. The parse preflight and input-chemistry checks passed. This eligibility is limited to a **diagnostic on a no-NAG local model**; it does not satisfy the full glycan-retaining native-complex validation gate in the [earlier strategy](vina_validation_strategy.md). No redocking or RMSD has yet been performed.

Before interpreting any scores or ranking phytochemicals, assess the full glycan-retaining receptor or quantify the effect of its omission, test Asp542 protonation and His600 tautomer sensitivity, examine hydroxyl/water hydrogen orientations, and evaluate whether Vina's scoring function is adequate for this highly polar inhibitor. [Vina documents](https://autodock-vina.readthedocs.io/en/stable/faq.html) that its standard scoring ignores supplied partial charges, so charged-versus-neutral MIG comparisons must not be read as rigorous protonation free energies. Use randomized ligand starts, multiple seeds, the preserved deposited pose and atom mapping for heavy-atom RMSD, as specified in the [validation strategy](vina_validation_strategy.md). Report failed pose recovery or large state dependence explicitly. The unaffected 4GQR, 1XD0 and 2QMJ files and prior results remain unchanged.

# RDKit 3D generation of ten phytochemical reference ligands

**Run date:** 11 October 2026  
**Outcome:** 10/10 **neutral PubChem reference structures** prepared with automated SDF quality checks; 384/384 generated conformers converged and passed export read-back. These are reference models, **not** confirmed constituents of the original aqueous leaf extracts and **not** GOLD-ready docking results.

## Inputs and scope

The run used the 10 distinct `structure_id` values from [final_preliminary_ligand_library.csv](../../data/compounds/final_preliminary_ligand_library.csv), checked against [ligand_structure_verification.csv](../../data/metadata/ligand_structure_verification.csv). The input `isomeric_SMILES` and PubChem CIDs were left unchanged. One file was produced for each unique chemical; the source library's **12 plant–compound associations** remain intact. Eleven associations point to assigned reference structures, including two links to the **same rutin file** (AM02 and SA01). The twelfth, AM_KR01, is an identity-unresolved *Aegle* kaempferol-diglycoside hypothesis with no assigned structure and was not prepared.

The preparation script is [prepare_ligands.py](../../scripts/prepare_ligands.py). Its per-ligand settings, counts, file paths, SHA-256 hashes, evidence fields, and quality-check results are in [ligand_preparation_results.csv](../../data/metadata/ligand_preparation_results.csv). The ten SDFs are in [structures/ligands/prepared](../../structures/ligands/prepared/). Each **individual ligand SDF contains multiple conformer records of the same chemical**; the SDF record count is not a count of distinct compounds or botanical associations.

## Reproducible method

The run used **Python 3.12**, **RDKit 2026.3.6**, `ETKDGv3` embedding with one thread and fixed seeds `20261011 + ligand number`, and **MMFF94s** optimization with a 1,000-iteration cap. The installed RDKit environment was the repository's `.venv`; RDKit was installed there with `python -m pip install rdkit==2026.3.6`. The repository's original `requirements.txt` was not changed. Reproduce with:

```powershell
.\.venv\Scripts\python.exe scripts\prepare_ligands.py
```

The code requests **24 conformers** for each smaller molecule and **60** for each flexible glycoside (rutin, mangiferin, agnuside, and nicotiflorin). ETKDG pruning was disabled so attempted and embedded counts are directly auditable. MMFF94s minimized every conformer. Results are sorted by MMFF energy **within each chemical and microstate**, but energy is neither a binding score nor comparable biological activity. Because no RMS diversity filter was applied, some exported conformers may be geometrically similar; later docking design should select or cluster them explicitly. No pH-specific deprotonated acids, phenolate states, or alternative tautomers were generated. All records retain the neutral PubChem reference graph and explicit hydrogens.

## Per-ligand result

| ID | PubChem reference | Association(s) | Seed | Embedded / requested | MMFF94s converged | SDF records passing QC | Status |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| LIG01 | Gallic acid, CID 370 | AM01 | 20261012 | 24 / 24 | 24 | 24 | Pass |
| LIG02 | Rutin, CID 5280805 | AM02; SA01 | 20261013 | 60 / 60 | 60 | 60 | Pass |
| LIG03 | **trans-Caffeic acid**, CID 689043 | BD01 | 20261014 | 24 / 24 | 24 | 24 | **Pass as E-reference only** |
| LIG04 | Mangiferin, CID 5281647 | SR01 | 20261015 | 60 / 60 | 60 | 60 | Pass |
| LIG05 | **(−)-Epicatechin**, CID 72276 | SR02 | 20261016 | 24 / 24 | 24 | 24 | **Pass as (2R,3R) reference only** |
| LIG06 | Quercetin, CID 5280343 | SA02 | 20261017 | 24 / 24 | 24 | 24 | Pass |
| LIG07 | Kaempferol, CID 5280863 | SA_K01 | 20261018 | 24 / 24 | 24 | 24 | Pass |
| LIG08 | Luteolin, CID 5280445 | SA_L01 | 20261019 | 24 / 24 | 24 | 24 | Pass |
| LIG09 | Agnuside, CID 442416 | VN03 | 20261020 | 60 / 60 | 60 | 60 | **Pass as reference only** |
| LIG10 | Nicotiflorin, CID 5318767 | SA_KR01 | 20261021 | 60 / 60 | 60 | 60 | Pass as *Senna* reserve reference |

**Totals:** 384 requested; 384 embedded; 384 MMFF94s converged; **0 embedding failures, 0 unsupported MMFF parameter sets, 0 nonconverged or failed minimizations, 0 rejected stereo/geometry records**, and 384 records passing SDF read-back. Each row has a SHA-256 digest for its completed file.

## Checks and interpretation

For each ligand, the script confirmed that its CID and source isomeric SMILES agree across the library and verification table and that RDKit calculates the expected neutral formula and formal charge. Before export it used the canonical **isomeric** molecular graph and 3D-derived stereochemistry to reject altered connectivity, inverted centers, or changed alkene geometry. It checked finite 3D coordinates and nontrivial spatial extent. It then reread every SDF record with explicit hydrogens and checked record count, sanitization, atom count, formula, charge, CID property, graph, stereochemistry, and 3D geometry. Every check passed. This is **automated chemical-file QC**; GOLD/Hermes import, atom typing, binding-site protonation, and ligand microstate suitability have not been assessed.

The following distinctions are essential:

- **LIG03** is the underivatized PubChem **E/trans** caffeic-acid reference. The *Boerhavia diffusa* aqueous report used mixed aerial parts, the methanolic-leaf GC-MS report detected a **tri-TMS derivative**, and source E/Z geometry in this project's extracts remains unresolved. Preparing an E reference does not resolve those gaps.
- **LIG05** preserves CID 72276's **(2R,3R) (−)-epicatechin** configuration. The published hot-water-leaf HPLC association does not independently establish the optical isomer by a chiral method.
- **LIG09** preserves the PubChem agnuside stereocentres. The water-only leaf chromatogram and source-level stereochemical detail remain insufficiently accessible; a methanol-first aqueous partition is not direct aqueous extraction.
- **LIG02, LIG04, and LIG10** preserve their respective PubChem glycoside stereocentres and linkages. Their source papers do not necessarily prove every stereo element of those exact reference structures. Nicotiflorin is linked to the *Senna* association only, not to unresolved AM_KR01.
- **LIG01 and LIG03** are neutral acid references. Their physiologically relevant carboxylate states and any other receptor-dependent microstates require a separate, documented decision before GOLD use. The original aqueous extracts were not chemically profiled for any of these ten molecules.

Botanical evidence remains in the source library and in the results CSV's `botanical_evidence_tiers` and `botanical_evidence_status` fields. Preparation success establishes that the **chosen reference graph can be embedded, minimized, and exported without detected structural errors**. It does not upgrade a Tier B/C botanical association, resolve solvent mismatch, establish compound concentration, or imply antioxidant, antidiabetic, or docking activity. No GOLD docking, ADMET, or biological prediction was performed.

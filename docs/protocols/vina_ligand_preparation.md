# Preliminary AutoDock Vina ligand preparation with Meeko 0.8.0

**Status (12 October 2026):** Ten neutral PubChem **reference structures** were converted to preliminary PDBQT inputs and passed the checks below. They are **not released for production docking**. No receptor docking, redocking, score, affinity, or ADMET result was produced. A literature-reported compound is not thereby confirmed in the project's original aqueous leaf extracts.

## Inputs and evidence boundary

The ten original multi-record SDF files remain unchanged at [`structures/ligands/prepared/LIG01.sdf` through `LIG10.sdf`](../../structures/ligands/prepared/). Before and after conversion, their SHA-256 hashes matched [ligand_preparation_results.csv](../../data/metadata/ligand_preparation_results.csv). Each record's `structure_id`, PubChem CID, formula, charge and isomeric molecular graph were checked against the [previously verified PubChem reference register](../../data/metadata/ligand_structure_verification.csv) and all corresponding associations in the [final preliminary library](../../data/compounds/final_preliminary_ligand_library.csv). This was a **local comparison to PubChem metadata already checked on 11 October 2026**, not a new live PubChem retrieval. Rutin is prepared once as `LIG02` despite its two plant associations; the unresolved *Aegle* kaempferol-rutinoside-like association has no structure ID and was not converted.

The exact repeatable workflow is [`scripts/prepare_vina_ligands.py`](../../scripts/prepare_vina_ligands.py), run with the project's `.venv` Python. The script requires **Meeko 0.8.0** and records RDKit **2026.03.6**. It writes the [results table](../../data/metadata/vina_ligand_preparation_results.csv) and ten individually named [PDBQT files](../../structures/ligands/vina_prepared/). Source SDFs, botanical evidence datasets, and reference stereoisomer assignments are not edited.

## Record inspection and conformer selection

RDKit read **every SDF record**, with explicit hydrogens and sanitization enabled: 24 records each for LIG01, LIG03, LIG05–LIG08; 60 each for LIG02, LIG04, LIG09 and LIG10 (**384 records total**). All 384 parsed successfully. For every record, the script checked the reference CID/SMILES property, exact isomeric molecular graph, formula, formal charge, complete explicit hydrogens, finite nondegenerate 3D coordinates and plausible covalent distances. It independently reassigned stereochemistry from 3D and compared the resulting isomeric graph to the stored reference. **All 384 passed** these source checks; this includes the PubChem *E* geometry of trans-caffeic acid and the assigned glycoside and epicatechin stereocentres.

For this **initial conversion test**, the representative was the record with the lowest **stored MMFF94s energy for that same neutral molecule**, breaking ties by earliest SDF record. The source SDFs happen to be energy-sorted, so the selected record is index 1 in each file, but the script calculates the minimum rather than assuming that ordering. The `conformer_number` is the original source property and differs from the SDF record index. These isolated-molecule energies are **not comparable binding energies**, and the lowest-energy isolated conformer is **not presumed to be the bioactive bound conformation**.

| ID | Reference ligand | SDF records | Selected record / source `conformer_number` | MMFF94s energy (kcal/mol) | PDBQT atoms | Torsions (`TORSDOF`) |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| LIG01 | Gallic acid | 24 | 1 / 19 | −12.22774638 | 16 | 5 |
| LIG02 | Rutin | 60 | 1 / 17 | 218.31880727 | 53 | 16 |
| LIG03 | *E*-caffeic acid reference | 24 | 1 / 6 | −24.83531901 | 16 | 5 |
| LIG04 | Mangiferin | 60 | 1 / 21 | 120.47779257 | 38 | 10 |
| LIG05 | PubChem (2R,3R) (−)-epicatechin reference | 24 | 1 / 4 | 42.36955595 | 26 | 6 |
| LIG06 | Quercetin | 24 | 1 / 23 | 43.38272614 | 27 | 6 |
| LIG07 | Kaempferol | 24 | 1 / 16 | 46.36634832 | 25 | 5 |
| LIG08 | Luteolin | 24 | 1 / 21 | 18.01877740 | 25 | 5 |
| LIG09 | Agnuside | 60 | 1 / 29 | 118.69383853 | 39 | 13 |
| LIG10 | Nicotiflorin | 60 | 1 / 45 | 222.56519023 | 51 | 15 |

## Meeko conversion and quality checks

Each selected RDKit molecule was passed to `MoleculePreparation` with **Gasteiger partial charges**, default nonpolar-H merging (`merge_these_atom_types=('H',)`), index-map remarks enabled, no hydration, rigid macrocycles disabled and flexible amides disabled. `PDBQTWriterLegacy.write_string()` produced one PDBQT per structure. These settings follow Meeko's [ligand preparation interface](https://github.com/forlilab/Meeko/blob/develop/docs/source/lig_prep_basic.rst) and its requirement for a 3D, explicitly hydrogenated RDKit input ([Meeko overview](https://github.com/forlilab/Meeko/blob/develop/docs/source/lig_overview.rst)).

For each output, the script required a nonempty file and one valid `ROOT`/`TORSDOF` torsion tree with balanced branches. Meeko's `PDBQTMolecule` parsed the written file, and `RDKitMolCreate` reconstructed **exactly one molecule with the same isomeric graph** as the selected SDF. The Meeko index map linked every source heavy atom to an output PDBQT atom; the maximum coordinate difference across all outputs was **0.0008 Å**, consistent with PDBQT's three-decimal coordinate storage. All Gasteiger charge sums were within **0.003 e** of the neutral formal charge after output rounding. PDBQT atom counts are lower than source SDF atom counts because Meeko merged nonpolar hydrogens; **no heavy atom was lost**. Each file's SHA-256, source record index, original conformer number, energy, input seed, atom counts, PDBQT atom types, charge sum, torsion count and read-back result are in the [results CSV](../../data/metadata/vina_ligand_preparation_results.csv).

All ten conversions passed these **format and molecular-identity** checks. This does not test a Vina run, receptor compatibility, microstate suitability, torsional search convergence, or experimental activity.

## Required docking-specific decisions

- **Acids and phenols:** LIG01 and LIG03 are **neutral carboxylic-acid references**. Their carboxylate states at the chosen site pH require an explicit choice and separately named PDBQT derivatives if included. Quercetin, kaempferol, luteolin and other phenolics may need documented tautomer/protonation review. No alternative states were silently generated.
- **Source stereochemistry:** LIG03 is specifically PubChem *E*/trans caffeic acid, while the relevant *Boerhavia* source does not resolve E/Z in isolated leaves and a reported leaf GC-MS signal is derivatized. LIG05 retains PubChem CID 72276 (2R,3R) stereochemistry, but non-chiral leaf HPLC does not prove that optical isomer. LIG09 is the PubChem agnuside reference while complete water-only leaf evidence remains unavailable. These three are **test files under existing source/identity holds**, not approved production ligands.
- **Flexible glycosides:** LIG02, LIG04, LIG09 and LIG10 have many torsions; one selected conformer and a successful file conversion do not establish adequate bound-state sampling. Preserve their reference sugar stereochemistry and consider multiple independent starting conformers after a receptor-specific protocol is chosen. LIG10 is the *Senna* nicotiflorin reference, not the unresolved *Aegle* diglycoside.
- **Receptor and grid:** Select the actual receptor microstate, binding-site waters/ions, charge treatment and Vina search space before docking. Review Vina and Meeko parsing in the intended run environment. The current PDBQT files are **preliminary technical inputs**, not docking-ready biological claims.

## Reproduction

From the repository root, run:

```powershell
.venv\Scripts\python.exe scripts\prepare_vina_ligands.py
```

The script fails an individual ligand if any source SDF record, CID/reference comparison, Meeko conversion or read-back check fails, and records the exception in `error_message`; it never writes a success row for a failed conversion. Its outputs are fixed, uniquely named `LIG##_CID######_neutral_reference.pdbqt` files. Re-run with the same unchanged SDFs and pinned software environment to regenerate and compare the recorded hashes.

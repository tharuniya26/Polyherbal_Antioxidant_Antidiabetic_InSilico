# AutoDock Vina receptor preparation assessment

**Assessment date:** 12 October 2026. **Status:** no receptor has cleared native-ligand redocking review. The only saved receptor PDBQT is a **2QMJ NAG-free sensitivity model**, not a validated representation of the deposited glycoprotein. The [machine-readable results](../../data/metadata/vina_receptor_preparation_results.csv) and [reproducible audit/export script](../../scripts/prepare_vina_validation_inputs.py) record the files and checks. No Vina docking was run.

## Source integrity and software

The script checked SHA256 values in [receptor preparation results](../../data/metadata/receptor_preparation_results.csv) and [GOLD-stage preparation results](../../data/metadata/gold_receptor_preparation_results.csv) before reading the five existing receptor variants and both deposited reference-ligand SDFs/maps. It did not change the original PDB/mmCIF, working PDB, reference SDF, or GOLD-stage files. The Vina executable at `E:\docking tools\vina_1.2.7_win.exe` reports version **1.2.7**; it was not used to dock. The project environment supplies **Meeko 0.8.0**, RDKit and Gemmi. Meeko's `mk_prepare_receptor` CLI was examined and then run against each variant in a temporary directory with `--read_pdb` and `--write_pdbqt`; no residue-deletion or atom-filtering option was used. Its installed default residue templates include standard protein residues, HOH, CA and CL, but no linked NAG or PCA template. Meeko [receptor CLI documentation](https://github.com/forlilab/Meeko/blob/develop/docs/source/cli_rec_prep.rst) describes custom templates for nonstandard residues; their bonded linkers and charges require chemical validation.

## Prior structural preparation carried forward

The existing [GOLD-stage assessment](gold_receptor_preparation.md) resolved the six 2QMJ alternate residues by highest whole-residue occupancy, then alphabetical label on a tie (A for all six). It reports no missing standard side-chain heavy atoms in either receptor. 2QMJ still lacks construct residues 1–6 and Gln837, has an explicit break across 837, lacks terminal OXT, and has a deposited Cys26–Cys44 sulfur separation of 2.916 Å. These issues are distant from the active site but remain structural caveats. 1XD0 has no detected missing protein residue or alternate conformer. Its N-terminal pyroglutamate, **PCA A1**, remains peptide-linked; replacing it with glutamine would change the deposited chemical structure.

| Receptor variant | Input chemistry | Meeko CLI result | Saved receptor PDBQT |
|---|---|---|---|
| 2QMJ, all waters, NAG retained | 626 waters; NAG X1, X2, A2001, A2005; no added H | Blocked: linked NAG lacks a validated residue template | None |
| 2QMJ, site waters, NAG retained | 10 waters; same NAG; no added H | Blocked: linked NAG lacks a validated residue template | None |
| 2QMJ, site waters, NAG removed | 10 waters; pH 7 OpenMM H estimate | Converted; **preliminary sensitivity model only** | [`2QMJ_sitewaters_noNAG_preliminary_sensitivity.pdbqt`](../../structures/proteins/vina_prepared/2QMJ_sitewaters_noNAG_preliminary_sensitivity.pdbqt) |
| 1XD0, all waters, Ca/Cl retained | 142 waters; PCA A1; Ca A498; Cl A500 | Blocked: peptide-linked PCA lacks a validated residue template | None |
| 1XD0, site waters, Ca/Cl retained | 6 waters including Ca-coordinating HOH536/546; PCA, Ca, Cl | Blocked: peptide-linked PCA lacks a validated residue template | None |

The Meeko error identifies NAG/PCA as unknown **linking fragments**. Its suggestion to use `--delete_residues` was intentionally not used. A separate disposable diagnostic that deleted 1XD0 PCA showed Meeko could emit `Ca` at +2.000 and `Cl` at −1.000, but that chemical deletion invalidates the intact receptor. It proves only template availability for the isolated ions, not correct Ca coordination or chloride treatment in a valid whole-receptor PDBQT. No PCA-deleted file was saved.

## The 2QMJ preliminary PDBQT

The accepted Meeko conversion contains **8,436 parsed PDBQT atoms**, including **6,922 heavy atoms**, equal to all 6,922 heavy atoms in the NAG-free input PDB. A coordinate-multiset comparison retained every heavy-atom position to 0.001 Å; Meeko read-back found finite coordinates and partial charges. The 10 retained waters are explicitly present. This is format and coordinate quality control, **not** validation of charges, protonation, water orientation, glycosylation choice, or native-pose recovery.

Removing NAG A2005 is consequential because it lies about **3.895 Å** from the deposited alpha-acarbose. The NAG-retaining variants were left unprotonated in the earlier preparation because the available OpenMM treatment did not add chemically valid hydrogens across linked NAG/Asn. A glycoprotein-aware, chemically checked NAG/Asn model and a Meeko template with correct linkage, atom types and charge would be needed to represent the deposited receptor. The NAG-free model is only a sensitivity case. Its Asp443/Asp542 carboxylates and His600 ND1-H state are pH 7 estimates made without native ligand present; proton-transfer microstates need review. OpenMM-placed HOH4627 H1 approaches a native-ligand heavy atom by **1.553 Å**, so this water orientation must be resolved before redocking. The nine glycerols and remote sulfate were already omitted in the derived variants because they lie at least 13.864 Å and 49.365 Å from the ligand, respectively; no active-site heteroatom was removed merely to make Meeko succeed.

## The 1XD0 block

Both derived 1XD0 variants retained **Ca A498**, **Cl A500**, and peptide-linked **PCA A1**. The selected site-water variant retains four ligand-near waters plus calcium-coordinating **HOH A536/A546**; the all-water variant retains all 142. A remote NAG A499 was omitted at the earlier GOLD-staging step because its nearest ARE distance is about 47.4 Å. Meeko cannot infer a reliable template for peptide-linked PCA A1; therefore **no 1XD0 receptor PDBQT was saved**. Ca/Cl must be retained and their coordination checked after a chemically validated PCA template is supplied. The pH 7 states for Asp197, Glu233, Asp300 and His201 are provisional. ARE also contacts a crystal mate near O3A/C3A, which limits full-pose interpretation in a single-chain receptor.

## Gate before redocking

Prepare and validate the NAG/Asn linkages for a glycan-retaining 2QMJ receptor or explicitly justify a NAG-free experimental sensitivity comparison; prepare peptide-linked PCA A1 for 1XD0 without substituting or deleting it. Review both catalytic microstates and water H orientations with the native ligand present. Verify full receptor atom mapping, charge and ion/water typing after any new template; keep the original heavy coordinates and reference poses. Only then perform a Vina import/energy check and native-ligand redocking. A successful PDBQT parse alone does not clear this gate.

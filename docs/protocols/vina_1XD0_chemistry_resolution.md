# 1XD0 receptor chemistry assessment for AutoDock Vina

**Assessment date:** 12 October 2026. **Decision:** the deposited ion identities, local contacts, existing PDBQT atom types, selected water identities, and an ARE crystal-mate contact are now documented. The two existing receptor PDBQTs remain **unvalidated for physically interpretable standard-Vina redocking**. No structures, atom mappings, protonation states, scores, or docking poses were changed or generated. The accompanying [issue register](../../data/metadata/vina_1XD0_chemistry_resolution.csv) separates coordinate observations from modelling proposals.

## Inputs and method

The unchanged official [1XD0 entry](https://www.rcsb.org/structure/1XD0) is human pancreatic α-amylase (AMY2A) with the deposited **ARE acarbose-derived pentasaccharide**, calcium and chloride, determined by X-ray diffraction at 2.00 Å. The [primary structural report](https://doi.org/10.1021/bi048334e) describes the acarbose-rearrangement chemistry; ARE must not be replaced with intact acarbose. Local coordinate measurements below use `structures/proteins/original/1XD0.pdb` (SHA256 `029c9fc8e1570649e02cd7fa088b27d44cfc4eda6b78183062ed45316e41c147`) and confirm the ions and 55 ARE atoms in `structures/proteins/original/1XD0.cif` (SHA256 `1ca3fb8a26fa29ac1d71b33033b80c62f41b72b6edaad1331fc224f35701dc4d`). Distances are Euclidean heavy-atom distances rounded to 0.001 Å; contact distance alone does not establish a bond or protonation state.

The existing [preflight results](../../data/metadata/vina_1xd0_input_preflight_results.csv), [preflight report](vina_1xd0_input_preflight.md), [blocker-resolution audit](../../data/metadata/vina_receptor_blocker_resolution_results.csv), [receptor-preparation audit](../../data/metadata/vina_receptor_preparation_results.csv), and [native-redocking plan](vina_native_redocking_plan.md) were read. The requested `vina_1XD0_preflight.md` does not exist; the lowercase [existing report](vina_1xd0_input_preflight.md) is its actual path. The two existing PDBQTs, unchanged deposited [ARE SDF](../../docking/validation/1XD0_ARE_native_pose.sdf), [PDB atom map](../../docking/validation/1XD0_ARE_atom_map.csv), [PDBQT bridge map](../../docking/validation/vina/1XD0_ARE_sdf_to_pdbqt_atom_map.csv), and [ARE crystal-pose PDBQT](../../docking/validation/vina/1XD0_ARE_crystal_pose_reference.pdbqt) were inspected. The preparation and preflight records already verify the mapped ARE graph, stereochemistry and 55 heavy coordinates. The original crystal pose remains the RMSD reference only.

## Calcium and chloride: deposited chemistry versus Vina typing

| Ion | Deposited position and local contacts in 1XD0 | Existing PDBQT field | Interpretation |
| --- | --- | --- | --- |
| **Ca A498**, (10.212, 6.265, 48.613) Å | Asn100 OD1 **2.428** Å; Arg158 backbone O **2.417** Å; Asp167 OD1/OD2 **2.464/2.506** Å; His201 backbone O **2.479** Å; HOH536/HOH546 oxygen **2.531/2.616** Å. | `Ca`, `+2.000` in both variants; both coordinating waters retained. | This is a plausible seven-oxygen coordination shell, consistent with the primary [human pancreatic α-amylase structural description](https://pubmed.ncbi.nlm.nih.gov/8528071/). His201's **backbone carbonyl**, not its imidazole, contacts calcium. Ca is about **10.077 Å** from the nearest ARE heavy atom. |
| **Cl A500**, (7.567, 23.217, 48.802) Å | Arg195 NE/NH2 **3.118/3.304** Å; Arg337 NH2/NH1 **3.135/3.374** Å; Asn298 ND2 **3.391** Å; HOH519 oxygen **3.229** Å. | `Cl`, `−1.000` in both variants; **HOH519 present only in the 142-water variant**. | The arginine/asparagine environment matches the chloride-activation literature. The closest ARE heavy atom is about **5.541 Å** away. Omission of HOH519 changes its deposited hydration environment. |

The [chloride-activation structural study](https://pmc.ncbi.nlm.nih.gov/articles/PMC2279296/) identifies chloride as a catalytic activator and calcium as a structural stabilizer; it reports catalytic-site rearrangement on chloride loss. Retaining both ions is therefore the defensible structural starting point. The distances above come from 1XD0, not from an idealized coordination template.

The existing PDBQT *charge fields and type labels* encode the intended formal ion identities. They do **not** establish correct scoring: the [Vina FAQ](https://autodock-vina.readthedocs.io/en/stable/faq.html) says standard Vina ignores user-supplied partial charges, while [Vina 1.2.7 atom typing](https://github.com/ccsb-scripps/AutoDock-Vina/blob/v1.2.7/src/lib/model.cpp) maps chlorine to a hydrophobic type and metals to a generic type. Thus the parsed Cl atom is not an explicit chloride-electrostatics term and Ca coordination is not calibrated by its `+2` field. Both `--score_only` parse tests passed previously; this is **format compatibility**, not chemical validation. Do not delete ions for an easier run. A separately validated electrostatics/ion-aware scoring protocol, with documented Ca/Cl parameters and native-pose controls, is needed before interpreting ion-sensitive scores. AutoDock4-style maps or another engine are possibilities, not verified solutions in this project.

## Catalytic groups and ARE microstates

The original PDB/mmCIF has no hydrogens that experimentally assign side-chain protonation. The existing OpenMM-prepared pH-7 receptor uses Asp197, Glu233 and Asp300 as deprotonated carboxylates, and His201 as neutral **HD1** (ND1-protonated/NE2-accepting). These are **computational assignments**. OpenMM 8.6.1 is installed, but its pH-based hydrogen placement is not a site-specific pKa calculation. Neither PROPKA nor PDB2PQR is installed in the project environment or found as a command; no residue pKa was calculated. RDKit, Gemmi and Meeko likewise do not resolve the protein catalytic pKa question here.

The original PDB lists crystallization at **pH 7.5**. Use pH **7.5 as a proposed reference preparation condition** for future explicitly named variants, and retain the existing pH-7 structures as historical variants. The biochemical assay pH and the microenvironment of a bound rearranged inhibitor need independent review; pH 7.5 is not a measured residue protonation assignment. A later PROPKA/PDB2PQR or equivalent run should record software/version, receptor and ligand state, ions/waters/glycan inclusion, pH, predicted pKas, warnings and alternative states. Predictions should be checked against catalytic mechanism and hydrogen-bond geometry, then represented in **new** files if justified.

| Group | Observed local environment in deposited coordinates | Defensible microstate assessment |
| --- | --- | --- |
| **Asp197** | OD2–ARE O7H **2.922 Å**; OD1–Arg195 NH1 **2.925 Å**; OD2–HOH513 O **2.916 Å**. | The established catalytic nucleophile favors a carboxylate in a conventional resting/complex model, but the coordinates alone do not prove its charge. Examine the Asp197–Arg195/ARE/HOH513 network before any alternative. |
| **Glu233** | OE1–ARE O3C **2.593 Å**; OE2–ARE N4C **2.816 Å**; OE2–HOH511 O **2.814 Å**. | [Mutagenesis and pH-activity work](https://doi.org/10.1021/bi011821z) supports Glu233 as the general acid/base. Its protonated versus deprotonated state depends on catalytic step and ligand microstate; the existing deprotonated form is a hypothesis. Predefine any protonated-Glu233 sensitivity variant rather than changing this model silently. |
| **Asp300** | OD2–ARE O2H **2.592 Å**; OD1–ARE O3H **2.777 Å**; OD2–HOH511 O **2.797 Å**. | Binding/water network is compatible with a deprotonated acceptor but does not prove it. Evaluate with Glu233 and ARE together; local ion/water changes could alter preferences. |
| **His201** | Side-chain NE2–ARE O2C **2.551 Å**; side-chain ND1–Val157 O **2.955 Å**; **backbone O–Ca 2.479 Å**. | Existing neutral HD1 leaves NE2 as an acceptor, plausible if ARE O2C is a donor. Neutral HE2 or doubly protonated histidine are alternatives requiring explicit assessment. Do not mistake the backbone Ca contact for evidence of imidazole protonation. |
| **ARE N4C** | The reference SDF/PDBQT represents a neutral secondary amine with one H. N4C–Glu233 OE2 **2.816 Å** and N4C–Asp300 OD2 **3.076 Å**. | Neutral versus protonated (+1) amine states require a pKa/microstate assessment at the stated pH. The deposited X-ray heavy atoms cannot place its H or establish its charge. Preserve the deposited **55-heavy-atom pose** and atom map; any protonated ARE is a separate, stereochemistry-preserving preparation. |

The mechanistic literature establishes roles for Asp197, Glu233 and Asp300; it does not give the bound 1XD0 microstate directly. A contact involving O/N atoms is a geometric observation, not proof of hydrogen-bond direction.

## Water inventory and controlled models

Both existing receptors contain HOH536 and HOH546, the two measured Ca-shell waters. The all-water receptor has **142** retained waters and the site-water receptor has **six**: **511, 536, 546, 588, 634, 766**. The latter omits **HOH519** from the Cl shell and **HOH513** near Asp197. All values in the table use deposited water oxygen coordinates; occupancies are 1.00 in the entry, while B factors and contact patterns remain relevant to uncertainty.

| HOH A | Key distances (Å) | Deposited B (Å²) | Assessment |
| --- | --- | ---: | --- |
| 536, 546 | Ca–O **2.531, 2.616**; nearest ARE **10.774, 11.976** | 11.89, 14.60 | Plausible structural Ca-shell waters; retain in any ion-preserving model. |
| 519 | Cl–O **3.229**; nearest ARE **5.743** | 12.12 | Plausible Cl-shell water; absent from six-water receptor. Retain for a chemistry-preserving ion-shell model, with its orientation reviewed. |
| 511 | Asp300 OD2–O **2.797**; Glu233 OE2–O **2.814**; nearest ARE O2H **3.749** | 10.89 | Plausible catalytic-pocket network water, but its role and H orientation are uncertain. |
| 513 | Asp197 OD2–O **2.916**; nearest ARE **5.103** | 12.77 | Plausible Asp197 network water; include as a predefined expanded-network sensitivity case. |
| 588 | ARE O6A–O **2.859** | 29.97 | Direct ARE contact; higher B and distal position make retention uncertain. |
| 634 | ARE O2G–O **3.865** | 34.76 | Marginal/weak direct contact at distal sugar; uncertain. |
| 766 | ARE O3C–O **2.971**; Glu233 **backbone O**–water O **2.920** (side-chain OE1–water O **3.486**) | 32.66 | Potential water-mediated contact; distinguish backbone from catalytic side chain and assess orientation. |

For future **new** receptor variants, predefine: (1) a *fully dry stress-test* with all crystal waters omitted but Ca and Cl retained; this leaves their coordination incomplete and is **not** a primary physical model; (2) an *ion-shell baseline* retaining **519, 536, 546**; (3) a seven-water site sensitivity model adding **511, 588, 634, 766**; and (4) an eight-water expanded-network model also adding **513**. The existing six-water variant is an historical sensitivity model because it lacks HOH519. The 142-water variant is an upper-retention comparison, not a default native-pose receptor. Select among these on structural grounds before reading any docking output. No new variant is prepared here.

Existing waters are `OA` oxygen plus two `HD` hydrogens in PDBQT, with charge fields around `−0.412/+0.206/+0.206`. The preflight measured water O–H lengths around **1.17–1.20 Å**, longer than an ordinary prepared water model and requiring geometry review. Standard rigid-receptor Vina fixes their positions and cannot exchange or displace them; its ordinary scoring does not model water reorientation/entropy or true ionic hydration. Static waters are sensitivity assumptions, **not validated dynamic solvent**. Even a water set that looks geometrically plausible cannot by itself clear the Ca/Cl scoring problem.

## Crystal packing and pose-recovery limit

The deposited cell is **P 21 21 21**, `a=52.9`, `b=68.9`, `c=132.0 Å`. Using the original PDB `SMTRY` matrices on chain-A protein heavy atoms, plus neighboring unit-cell translations from −1 to +1, gives five ARE–symmetry-mate protein heavy-atom pairs within **4.0 Å**. All involve symmetry operation 3 with cell shift `(0, −1, 0)`, equivalent to `−x, y−1/2, −z+1/2` in fractional coordinates:

| Crystal-pose ARE atom | Symmetry-mate protein atom | Distance (Å) |
| --- | --- | ---: |
| O3A | Asp433 backbone O | **2.721** |
| O3A | Asp433 C | **3.429** |
| O3A | Asp433 CB | **3.449** |
| O3A | Trp434 CE3 | **3.866** |
| C3A | Asp433 backbone O | **3.899** |

This is a **crystal-coordinate inference**; it does not establish how much lattice packing determines the solution pose. It does show that the distal A sugar has contacts absent in a single-chain receptor. If future redocking recovers the catalytic core but displaces that sugar, report both full-ligand and per-sugar heavy-atom RMSD, and inspect the missing crystal contact before calling the search a failure or success. Keep the deposited ARE coordinates unchanged as the reference.

## Decision and remaining work

The current **site-water** variant has a chemically incomplete chloride hydration shell; the **all-water** variant has the Cl-shell water but imposes 142 fixed waters. Neither is a defensible *primary* standard-Vina validation receptor in its current state. The best future **diagnostic** 1XD0 candidate is a newly prepared, ion-retaining **seven-water model** (Cl/Ca-shell waters plus the four selected ARE/pocket waters), with a prespecified eight-water Asp197-network comparison and a dry stress-test clearly labeled artificial. It requires a protonation/ARE-microstate assessment, water geometry checks, verified PDBQT readback and a fresh Vina preflight. Even then, standard Vina's Ca/Cl terms limit any score interpretation; native-pose RMSD with independent seeds would be a diagnostic only. **Controlled diagnostic redocking should remain on hold until those inputs and choices are documented.**

The [4GQR myricetin complex](https://www.rcsb.org/structure/4GQR) offers a smaller, experimentally observed active-site ligand for a separate α-amylase pose-recovery test, as described in the [validation strategy](vina_validation_strategy.md), but it also retains Ca/Cl and does not solve Vina's ion model. A validated ion-aware alternative scoring method may ultimately be needed. 1XD0 itself is **not fully validated**; no redocking, RMSD assessment, production ligand docking or ADMET analysis has been performed.

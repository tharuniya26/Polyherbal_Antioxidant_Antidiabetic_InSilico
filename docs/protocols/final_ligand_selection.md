# Final preliminary five-plant ligand selection

**Selection date:** 11 October 2026  
**Status:** Literature-supported exploratory library for later structure assessment. No ligand is declared docking-ready by this document.

## Scope and source hierarchy

The [consolidated CSV](../../data/compounds/final_preliminary_ligand_library.csv) preserves one row per **plant–compound association** and uses `structure_id` to identify one unique reference chemical across plants. It combines the original [candidate register](../../data/compounds/phytochemical_candidates.csv), the [expanded register](../../data/compounds/expanded_compound_candidates.csv), the [five-plant shortlist](../../data/compounds/five_plant_preliminary_shortlist.csv), the [five-plant evidence review](../../literature/five_plant_compound_evidence_review.md), the [focused *Aegle*–*Senna* review](../../literature/aegle_senna_compound_evidence_review.md), and the [evidence corrections](../compound_evidence_corrections.md). The newer source checks govern classifications where the older original register differs; no source file was overwritten.

The library contains **12 plant–compound associations** and **10 unique assigned PubChem reference structures**. One row includes an identity-unresolved *Aegle* kaempferol-diglycoside hypothesis with **no assigned structure ID, PubChem CID, or SMILES**. The same rutin reference (`LIG02`, PubChem CID 5280805) is linked to both *Aegle marmelos* and *Senna auriculata*. The *Senna* kaempferol-3-O-rutinoside reference (`LIG10`, CID 5318767) is **not** copied to the unresolved *Aegle* peak. Rutin and kaempferol-3-O-rutinoside are chemically different: their PubChem reference formulae are C₂₇H₃₀O₁₆ and C₂₇H₃₀O₁₅, respectively. [PubChem rutin](https://pubchem.ncbi.nlm.nih.gov/compound/5280805); [PubChem kaempferol-3-O-rutinoside](https://pubchem.ncbi.nlm.nih.gov/compound/5318767).

The CSV's `canonical_SMILES` column contains PubChem's current **ConnectivitySMILES** (the legacy canonical field, without stereochemical specification); `isomeric_SMILES` contains PubChem's current **SMILES** (the stereochemically specified field where applicable). Values and CIDs were checked directly through [PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest) on 11 October 2026. [PubChem's terminology note](https://pubchem.ncbi.nlm.nih.gov/docs/glossary) explains the field rename. This verifies database reference identities and strings. It does **not** prove that every source paper established a compound's exact stereochemistry, linkage, tautomer, charge state, or occurrence in the original experimental extracts. No structure files or coordinates were downloaded.

## Evidence and preparation decisions

| Structure / association | Plant-specific evidence | Priority | Remaining gate |
| --- | --- | --- | --- |
| `LIG01` gallic acid / AM01 | **Tier B:** standard-matched HPLC-PDA peak from **ethanol leaves**, [Ahmad et al.](https://doi.org/10.3390/plants10122573). | **Primary** | Verify relevance to water-only leaves; later assess carboxylate state. |
| `LIG02` rutin / AM02 and SA01 | **Tier B:** *Aegle* ethanol-leaf HPLC and independent methanol/hydroalcoholic leaf reports; *Senna* isolation from **80% aqueous-acetone leaves** under synonym *Cassia auriculata*, [Ahmad et al.](https://doi.org/10.3390/plants10122573), [Ibrahim et al.](https://doi.org/10.1016/j.toxrep.2018.09.005), [Nakamura et al.](https://doi.org/10.1248/cpb.c14-00420). | **Reserve, one molecule** | No water-only leaf proof for either plant; retain sugar stereochemistry and address size/flexibility. |
| `LIG03` caffeic acid / BD01 | **Tier C:** underivatized acid in a **boiled-water infusion of mixed aerial parts**, [Sinan et al.](https://doi.org/10.3390/antiox10122003); a separate methanolic-leaf [GC-MS study](https://doi.org/10.1016/j.jtcme.2019.02.002) detected a **tri-TMS derivative**, not the modelled parent. | **Primary, conditional** | Leaf-only parent occurrence and source E/Z geometry unverified. PubChem CID 689043 encodes the *E* reference; do not equate it with the derivatized leaf signal. Hold preparation until explicit source/isomer review. |
| `LIG04` mangiferin / SR01 | **Tier A:** compound-specific HPLC determination in **water leaves**, [Yoshino et al.](https://doi.org/10.1271/bbb.80854). | **Primary** | PubChem sugar stereochemistry is a reference; source-level stereochemical detail and original-extract occurrence remain unverified. |
| `LIG05` epicatechin / SR02 | **Tier A at epicatechin level:** HPLC report from **hot-water leaves**, [Sekiguchi et al.](https://doi.org/10.1186/1756-0500-5-198). | **Primary, stereoisomer conditional** | [PubChem CID 72276](https://pubchem.ncbi.nlm.nih.gov/compound/72276) is **(2R,3R) (−)-epicatechin**. The leaf paper does not verify optical configuration by chiral chromatography or rotation. The CID is a provisional reference pending stereochemical review. |
| `LIG06` quercetin / SA02 | **Tier B:** isolated from the **ethyl-acetate fraction** of 80%-aqueous-acetone *Cassia auriculata* leaves, [Nakamura et al.](https://doi.org/10.1248/cpb.c14-00420). | **Primary** | No water-only leaf recovery demonstrated. |
| `LIG07` kaempferol / SA_K01 | **Tier B:** isolated from the same *Senna* leaf extract and fraction, [Nakamura et al.](https://doi.org/10.1248/cpb.c14-00420). | **Primary** | No water-only leaf recovery demonstrated. |
| `LIG08` luteolin / SA_L01 | **Tier B:** isolated from the same *Senna* leaf extract and fraction, [Nakamura et al.](https://doi.org/10.1248/cpb.c14-00420). | **Primary** | No water-only leaf recovery demonstrated. |
| `LIG09` agnuside / VN03 | **Publisher-abstract-only aqueous-leaf comparison:** [Katare et al.](https://doi.org/10.1007/s42250-022-00383-8) name agnuside and compare water with organic solvents in *Vitex negundo* leaves; an independent [leaf HPLC paper](https://doi.org/10.1016/j.jpha.2013.09.008) starts with methanol before water partitioning. | **Primary, conditional** | Full water-specific chromatogram and values unavailable. Optimized 20%-ethanol-water yield is not a water-only result; source stereochemistry and original-extract occurrence remain unverified. Hold preparation until full aqueous-specific data review. |
| `LIG10` kaempferol-3-O-rutinoside / SA_KR01 | **Tier B:** isolated from the ***n*-butanol fraction** of 80%-aqueous-acetone *Senna* leaves, [Nakamura et al.](https://doi.org/10.1248/cpb.c14-00420). | **Reserve** | Check known-compound structural comparison and glycoside stereochemistry; no water-only leaf recovery demonstrated. |
| No assigned structure / AM_KR01 | **Tier B for a tentative leaf assignment:** *Aegle* methanol-leaf LC-MS/MS reports a kaempferol-rutinoside-like peak with [M+H]⁺ 595 and base peak 287, [Ibrahim et al.](https://doi.org/10.1016/j.toxrep.2018.09.005). | **Identity hold; excluded from ligand counts** | Formula/aglycone fragments do not establish the exact 3-O sugar linkage and stereoisomer. No PubChem structure is assigned to this *Aegle* association. |

The [*Senna* taxonomy check](https://powo.science.kew.org/taxon/urn%3Alsid%3Aipni.org%3Anames%3A484161-1) supports treating *Cassia auriculata* as a synonym. The [Food Chemistry aerial-part study](https://doi.org/10.1016/j.foodchem.2010.08.065) corroborates *Senna* species-level chemistry but does not by itself prove leaf-only occurrence. The [Pharmacognosy Journal study](https://doi.org/10.5530/pj.2018.3.75) examines flowers and is not used to upgrade leaf evidence. TPC/TFC standards, aqueous chromatographic mobile phases, post-extraction water fractions, and crude-extract bioassays are not compound-specific aqueous-leaf identifications.

**Primary** means preferred for the *planned* small exploratory library, **reserve** means optional comparison, and **identity hold** means no stereodefined ligand assignment. These are preparation priorities, not docking results or statements of readiness. There are **8 primary reference structures**, **2 reserve reference structures**, and **1 unresolved association without an assigned structure**. The CSV's separate `docking_readiness` field records pending reviews; no row claims a ligand is ready for docking.

## Formulation mapping

The formulation columns reflect [config/formulations.yaml](../../config/formulations.yaml). Values are **aqueous leaf-extract proportions**, not concentrations of any listed compound. Every association inherits its plant's formulation membership; the mapping does not establish that the molecule was in that formulation.

| Plant code and species | C1 extract ratio | C2 extract ratio | C3 extract ratio |
| --- | ---: | ---: | ---: |
| AM — *Aegle marmelos* | 2.0 | 2.0 | Absent |
| BD — *Boerhavia diffusa* | 1.5 | 1.0 | Absent |
| SA — *Senna auriculata* | 2.0 | 2.5 | 3.5 |
| SR — *Salacia reticulata* | 2.0 | 2.0 | 3.5 |
| VN — *Vitex negundo* | 2.5 | 2.5 | 3.0 |

Thus **all five plants are represented in C1 and C2**, while **C3 contains only *Senna*, *Salacia*, and *Vitex***. The library covers all five plants, but the evidence does not establish a five-plant molecular profile for any actual extract. It also gives no basis to infer docking affinity, molecular synergy, antioxidant activity, or antidiabetic effects for individual compounds.

## Gate before any later modelling

Recheck the selected source paper's exact compound assignment and PubChem structure, with special attention to SR02 optical isomer, BD01 alkene geometry and underivatized parent, glycoside stereochemistry/linkage, and VN03 water-specific analytical results. Any later ligand preparation must record the chosen protonation/tautomer and stereochemical assumptions. The *Aegle* kaempferol-rutinoside peak remains outside the assigned structure set until its exact identity is supported. Docking, protein acquisition, and ADMET analysis are outside this selection stage.

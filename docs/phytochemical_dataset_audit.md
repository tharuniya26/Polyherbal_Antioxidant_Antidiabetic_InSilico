# Scientific audit of the preliminary phytochemical candidate dataset

**Audit date:** 11 October 2026  
**Files reviewed in full:** `data/compounds/phytochemical_candidates.csv`, `literature/compound_evidence.md`, and `docs/protocols/compound_selection.md`.  
**Scope:** Evidence and data-quality review only. No compound was tested in this audit; no docking or ADMET analysis was performed.

## Overall assessment

The CSV is a coherent **literature candidate register**, not a chemical profile of the original aqueous leaf extracts. All ten rows have a traceable primary-paper DOI, a plant-specific leaf report, and a PubChem CID/structure-string pair that resolves to the cited PubChem compound record. The strength of the botanical identifications and their relevance to water-only extraction vary substantially. Only two rows have published water-leaf detection, both for *Salacia reticulata*. **Zero rows have direct identification in the original project's extracts.**

The most consequential issue is in the two *Boerhavia diffusa* entries: the cited GC-MS paper reports **trimethylsilyl (TMS) derivatives** under different PubChem CIDs from the underivatized caffeic and ferulic acids in the CSV. The CSV's parent-acid CIDs are chemically valid reference structures, but they are **not the CIDs printed in the experimental source**. This distinction must be preserved when interpreting the analytical evidence and preparing any eventual ligand.

## Dataset integrity and summary counts

| Check | Audit result |
| --- | --- |
| Total records | 10 plant–compound associations |
| Records by species | *A. marmelos* 2; *B. diffusa* 2; *S. reticulata* 2; *S. auriculata* 2; *V. negundo* 2 |
| Distinct names / PubChem CIDs / connectivity strings | 9 / 9 / 9 |
| Intentional cross-plant duplicate | Rutin: AM02 and SA01 share the same name, CID 5280805, and connectivity SMILES; each has separate botanical evidence |
| Duplicate `compound_id` or duplicate within a species | None |
| CSV format | 14 expected headers; all 10 records parse to those 14 fields; no empty cells or malformed rows detected |
| DOI checks | All six unique DOI links resolve to journal-article records with matching titles/journals in Crossref; Kew's taxonomy entry supports the *Cassia auriculata* synonym |
| PubChem checks | All 10 CID values resolve; each stored `canonical_SMILES` exactly matches PubChem's **ConnectivitySMILES** for that CID |
| Published water-only leaf detection | 2 rows (SR01, SR02); 8 rows have leaf evidence from organic or aqueous-organic extraction |
| Confirmed in original aqueous extracts | 0 rows; original extract chromatograms and compound-level measurements are unavailable |

The header `canonical_SMILES` is potentially misleading: its contents are PubChem **connectivity** SMILES, which discard stereochemistry and double-bond geometry. Seven unique CID records have a different isomeric SMILES from the stored string. This matters particularly for (−)-epicatechin, the glycosides, and hydroxycinnamic acids. A matching connectivity string verifies the stored text against PubChem; it does **not** verify that the paper established the exact stereoisomer or that the stored string is a ready-to-use docking ligand.

## Candidate-by-candidate assessment

**AM01 — Gallic acid, *Aegle marmelos*.** [E1] used authenticated leaves, extracted them with ethanol, and compared an HPLC peak with a gallic-acid standard. This supports a reported species/leaf association, with identification limited by reliance on chromatographic matching rather than an independently deposited spectrum. It is **not** water-leaf evidence. The phenolic hydroxyl groups make antioxidant chemistry plausible; the paper's glucose-related and radical-scavenging results concern the extract, not isolated gallic acid. PubChem [CID 370](https://pubchem.ncbi.nlm.nih.gov/compound/370) and the stored connectivity string match. **Docking disposition:** chemically tractable small molecule, but its acid ionization state, chosen target, and direct water-leaf relevance need resolving; retain on hold.

**AM02 — Rutin, *Aegle marmelos*.** [E1] reports HPLC matching in the same ethanol leaf extract. The source does not establish recovery in water alone or a causal role for rutin in the extract assays. PubChem [CID 5280805](https://pubchem.ncbi.nlm.nih.gov/compound/5280805) and the stored string match; this is the **same chemical** as SA01, so a future chemical ligand library should contain one verified rutin structure linked to two botanical records. **Docking disposition:** large, flexible, highly hydroxylated glycoside; a protocol would need careful stereochemistry, conformer sampling, and score interpretation. Hold pending aqueous evidence.

**BD01 — Caffeic acid, *Boerhavia diffusa*.** [E2] analyzed methanolic leaves by derivatized GC-MS and reported a caffeic-acid-related peak. Its Table 2 gives **CID 5376254 for the tri-TMS derivative**, while the CSV gives [CID 689043](https://pubchem.ncbi.nlm.nih.gov/compound/689043) for underivatized trans-caffeic acid. Both CIDs resolve, but they denote different chemical entities; the parent acid is an inferred pre-derivatization candidate, not the directly measured GC-MS species. The reported GC-MS/library assignment should be treated as **putative** without orthogonal confirmation, and no water-only leaf detection is shown. Phenolic structure supports an antioxidant hypothesis; the cited study addresses wound healing and does not establish an antidiabetic mechanism for this compound. **Docking disposition:** small and chemically feasible after verifying the parent identity and acid protonation; keep on hold.

**BD02 — Ferulic acid, *Boerhavia diffusa*.** The same [E2] methanolic-leaf GC-MS study reports ferulic acid as a **di-TMS derivative, CID 5379186**; the CSV's [CID 445858](https://pubchem.ncbi.nlm.nih.gov/compound/445858) is underivatized trans-ferulic acid. This is the same parent-versus-derivative issue as BD01, with putative botanical identification and no water-only evidence. Phenolic antioxidant relevance is plausible; compound-specific antidiabetic activity is not established by [E2]. **Docking disposition:** small molecule, but confirm the parent analyte and relevant ionization/geometry first; keep on hold.

**SR01 — Mangiferin, *Salacia reticulata*.** [E3] prepared a water extract from leaves and determined mangiferin by HPLC, supporting **species + leaf + water** association. Its mouse glucose findings apply to the complete extract; they cannot be assigned to mangiferin. PubChem [CID 5281647](https://pubchem.ncbi.nlm.nih.gov/compound/5281647) and the stored connectivity string match. The polyphenolic xanthone scaffold provides an antioxidant rationale; separate root-derived enzyme literature is mechanistic context, not proof of the leaf extract's active principle. **Docking disposition:** reasonable preliminary candidate once stereochemical structure, protonation, and a biologically justified target are specified. It is still unconfirmed in the project's extracts.

**SR02 — (−)-Epicatechin, *Salacia reticulata*.** [E4] reports an epicatechin HPLC signal in a hot-water leaf extract; [E3] also describes analysis of leaf-water catechins. This supports a published aqueous-leaf association, though chromatographic peak identity and stereochemistry deserve review. PubChem [CID 72276](https://pubchem.ncbi.nlm.nih.gov/compound/72276) resolves to a stereodefined epicatechin record, but the stored connectivity SMILES **cannot encode the “(−)” assignment**. The [E4] biological experiment concerns an inflammatory-cell model, not isolated-compound antidiabetic action. Antioxidant relevance is plausible from its flavanol chemistry. **Docking disposition:** feasible only with the isomeric structure and documented stereochemical choice; preliminarily retain.

**SA01 — Rutin, *Senna auriculata*.** [E5] isolated rutin from *Cassia auriculata* leaves; [T1] confirms this as a synonym of *Senna auriculata*. The source used **80% aqueous acetone** and further fractionation, and also recovered rutin from an n-butanol fraction. Neither finding demonstrates a water-only leaf extract. CID 5280805 and its connectivity string match PubChem. The cited bioassay concerns hepatoprotection, not an antidiabetic target or the original formulation; antioxidant relevance is a plausible chemical hypothesis. **Docking disposition:** same structure and glycoside-related caveats as AM02; hold for water relevance and use only one chemical ligand record if later advanced.

**SA02 — Quercetin, *Senna auriculata*.** [E5] isolated quercetin from the **ethyl-acetate fraction** of an 80% aqueous-acetone leaf extract under the *Cassia* synonym [T1]. Species and leaf identification are strong, but this organic fraction provides weak support for occurrence in a water-only extract. PubChem [CID 5280343](https://pubchem.ncbi.nlm.nih.gov/compound/5280343) and connectivity SMILES match. A flavonol redox rationale is reasonable; the cited paper does not test quercetin in the project's antioxidant or glucose assays. **Docking disposition:** tractable aglycone after target and ionization selection, but deprioritize as a representation of the original aqueous formulation until water-leaf evidence is found.

**VN01 — Chlorogenic acid, *Vitex negundo*.** [E6] used **80% aqueous methanol** on authenticated leaves and reported LC-MS/MS identification and quantification of chlorogenic acid. This is solid species/leaf evidence but not water-only extraction. PubChem [CID 1794427](https://pubchem.ncbi.nlm.nih.gov/compound/1794427) and its connectivity string match. The article relates antioxidant activity to multi-compound fractions, not purified chlorogenic acid, and does not provide a compound-specific antidiabetic mechanism for this dataset. **Docking disposition:** plausible with a positional-isomer check (chlorogenic/neochlorogenic/cryptochlorogenic acid are distinct), stereochemistry, conformers, and carboxylate state documented. Hold for water-only relevance.

**VN02 — Isoorientin, *Vitex negundo*.** [E6] reports LC-MS/MS identification and quantification from the same aqueous-methanol leaves. Its fraction-level antioxidant evidence supports interest but does not isolate isoorientin's contribution, and the cited work does not establish antidiabetic action by this compound. PubChem [CID 114776](https://pubchem.ncbi.nlm.nih.gov/compound/114776) matches the stored connectivity string. **Docking disposition:** chemically definable C-glycoside but stereochemistry, flexibility, and target choice require care; hold for water-only evidence.

## Reference and mechanism appraisal

The six `E1`–`E6` citations in [the evidence register](../literature/compound_evidence.md) correspond to traceable journal articles; the DOI-to-title mappings were checked against [Crossref](https://api.crossref.org/). The [Kew taxon record](https://powo.science.kew.org/taxon/urn%3Alsid%3Aipni.org%3Anames%3A518328-1) supports the *Senna*/*Cassia* synonym. PubChem validates structure identifiers, **not** species occurrence. Source reliability is not uniform at the *compound-identification* level: [E5](https://doi.org/10.1248/cpb.c14-00420) isolated named constituents; [E6](https://doi.org/10.3390/molecules29133133) used LC-MS/MS with quantitative analysis; [E1](https://doi.org/10.3390/plants10122573), [E3](https://doi.org/10.1271/bbb.80854), and [E4](https://doi.org/10.1186/1756-0500-5-198) used targeted HPLC methods; [E2](https://doi.org/10.1016/j.jtcme.2019.02.002) used derivatized GC-MS and warrants the greatest identity caution.

The CSV's antioxidant and antidiabetic columns are mostly **rationales**, not compound-level mechanism evidence. Phenolic hydroxylation makes redox chemistry plausible, and an extract or fraction may score in DPPH/ABTS assays, but that does not identify the operative constituent or prove a protein target. Only a specified, biologically relevant target with a validated binding-site protocol could make docking an interpretable test; a docking score would still not establish antioxidant activity, glucose uptake, or clinical benefit.

## Recommendations before advancing candidates

1. **Preserve the ten botanical records** but treat the two *Boerhavia* identifications as putative parent compounds inferred from TMS derivative peaks. A later dataset revision should annotate both derivative CIDs explicitly and seek independent, underivatized LC-MS/HPLC confirmation.
2. **Seek water-only leaf analytical evidence** for the eight Tier-B rows. Do not promote an aqueous-acetone or aqueous-methanol report to the protocol's aqueous Tier A. Compare extraction temperature and duration as well as solvent, since the original extraction details are not available here.
3. **Re-examine compound identity before structure-based work:** use isomeric PubChem records, confirm epicatechin and sugar stereochemistry, distinguish chlorogenic positional isomers, and specify the protonation state of phenolic acids. One unique molecule should be prepared once even when linked to multiple plants.
4. **Specify a mechanism-led target plan** and suitable reference ligands/controls before assessing docking suitability. Large glycosides merit extra conformational sampling and cautious score comparisons with smaller aglycones.
5. **Keep claims bounded:** the two published aqueous-leaf associations are useful preliminary priorities, but none of the ten compounds was identified in the original extracts. Direct chemical profiling of those extracts is required to make that claim.

## Source links used for this audit

- [E1: *Aegle marmelos* leaf HPLC and ethanol extraction](https://doi.org/10.3390/plants10122573)
- [E2: *Boerhavia diffusa* leaf methanol GC-MS, including Table 2 derivative CIDs](https://doi.org/10.1016/j.jtcme.2019.02.002)
- [E3: *Salacia reticulata* leaf-water extract and mangiferin HPLC](https://doi.org/10.1271/bbb.80854)
- [E4: *Salacia reticulata* hot-water leaves and epicatechin HPLC](https://doi.org/10.1186/1756-0500-5-198)
- [E5: *Cassia auriculata* leaves and isolated rutin/quercetin](https://doi.org/10.1248/cpb.c14-00420)
- [E6: *Vitex negundo* leaves, LC-MS/MS, chlorogenic acid and isoorientin](https://doi.org/10.3390/molecules29133133)
- [Kew: accepted *Senna auriculata* name and *Cassia auriculata* synonym](https://powo.science.kew.org/taxon/urn%3Alsid%3Aipni.org%3Anames%3A518328-1)
- [PubChem PUG REST documentation](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest) and [Crossref API](https://api.crossref.org/) for identifier and metadata checks

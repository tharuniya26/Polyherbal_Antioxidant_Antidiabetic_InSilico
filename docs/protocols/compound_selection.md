# Protocol: preliminary phytochemical candidate selection

## Purpose and scope

This protocol defines a transparent, manageable literature-screening process for a Biomedical Science final-year computational study of five medicinal plants: *Aegle marmelos*, *Boerhavia diffusa*, *Salacia reticulata*, *Senna auriculata*, and *Vitex negundo*. The original experiment used **aqueous leaf extracts**. Literature occurrence in another plant part or solvent is not evidence of occurrence in those extracts.

The output of this stage is [the candidate CSV](../../data/compounds/phytochemical_candidates.csv) and [its evidence register](../../literature/compound_evidence.md). This stage does not select protein targets, prepare ligands, dock molecules, or predict ADMET.

## Search and source appraisal

1. Search peer-reviewed primary studies with the accepted species name and known synonyms, the terms *leaf/leaves*, *water/aqueous/hot water*, and compound-analysis terms such as HPLC, LC-MS/MS, or isolation. Search the synonym *Cassia auriculata* for *Senna auriculata* and verify taxonomy against an authoritative plant database.
2. Read the methods and results of each primary paper. Record the exact species, plant part, extraction solvent, analytical identification method, compound name, and DOI or stable publisher URL. Prefer papers that disclose a voucher/authentication method, reference standards, spectra, or an isolation/structure-elucidation procedure.
3. Distinguish detected, tentatively annotated, and merely discussed compounds. A named calibration standard, a compound in a review, or an assay expressed as gallic acid/quercetin equivalents is **not** evidence that that compound was detected in the sample.
4. Look up each proposed chemical name in PubChem. Record the CID and a machine-readable structure only if the name, stereochemistry, glycosylation, and positional isomer match. Record the query and access date in the evidence register. PubChem does not establish plant occurrence.
5. Cross-check whether biological evidence pertains to the isolated compound, an enriched fraction, or a crude extract. State the actual level of evidence in the CSV and register.

## Inclusion criteria

A preliminary row may be included when **all** of the following hold:

- A traceable primary publication reports the named compound in the correct plant species, with the synonym relationship documented when applicable.
- The plant part and extraction method are stated precisely enough to judge relevance to the original aqueous leaf extraction.
- The compound has a sufficiently defined chemical identity for a PubChem check or is explicitly marked as needing identity verification with CID and SMILES blank.
- There is a defensible antioxidant- or antidiabetic-related reason to consider it, with the strength of that reason stated without assigning crude-extract effects to the molecule.
- The set remains small and balanced across all five plants, initially about two plant–compound rows per species. A molecule occurring in more than one species receives separate rows because its botanical evidence is different.

### Evidence tiers and decisions

| Tier | Evidence | CSV treatment |
| --- | --- | --- |
| A: species + leaf + water | Primary analytical evidence of the named molecule in a water-only leaf extract | `evidence_status=species_leaf_aqueous_verified`; `inclusion_decision=include_preliminary` |
| B: species + leaf + other solvent | Primary analytical evidence from a leaf extract prepared with methanol, ethanol, acetone mixture, or a subsequent organic-solvent fraction | `evidence_status=species_leaf_nonaqueous_water_unverified`; `inclusion_decision=hold_for_aqueous_verification` |
| C: species, other or unspecified part | Compound is reported for the species but leaf occurrence is unverified | Keep out of the preliminary CSV unless there is a compelling reason to track it; if tracked, mark `requires_verification` and do not advance it |
| D: genus-only, review-only, unclear identification, or citation mismatch | No sufficiently traceable species-level primary evidence | Exclude from the candidate CSV until verified |

An aqueous-organic mixture belongs in Tier B. Detection in the project's original extracts would require direct chemical analysis of those extracts; even Tier A is not such confirmation.

## Exclusion and deferral criteria

Exclude or defer a proposed candidate when it appears only as a TPC/TFC assay standard, only in another species, only in a root/stem/flower study without leaf evidence, or only in a review without a checkable primary source. Defer ambiguous isomers, glycosides, salts, and stereoisomers until their chemical identity can be matched to the structure record. Exclude poorly supported GC-MS library matches from the final ligand list unless the analytical evidence can be reappraised. A lack of PubChem CID or a plausible structure is a reason to pause structure-based work, not a reason to guess a CID or SMILES.

## Data recording and review

- Use one CSV row per **plant–compound association** with a stable `compound_id`. Keep the 14 specified columns in their present order; quote CSV fields containing commas.
- Put the primary-paper citation key in `evidence_source`, the exact DOI or article URL in `DOI_or_URL`, and the solvent and plant part as reported, rather than normalizing them to the project's water extraction.
- Use `antioxidant_relevance` and `antidiabetic_relevance` to describe *why a molecule might be studied*. Label an untested relationship as a hypothesis. Keep measurement values out unless they are verified and sourced.
- A blank `PubChem_CID` or `canonical_SMILES` means the identity was not sufficiently checked; a blank DOI/URL means source verification is incomplete. Mark `evidence_status=requires_verification` and `inclusion_decision=defer` in that case.
- Review at least the paper's methods and relevant results table/figure before changing a Tier B row to Tier A. Record the new paper and rationale in the evidence register. Do not replace older evidence silently.
- Before any later docking, recheck PubChem identity and retrieve an **isomeric** structure with appropriate charge/protonation and stereochemistry. The CSV's connectivity SMILES is a reference, not a prepared ligand.

## Current stop point

The preliminary dataset contains ten rows: two for each plant. Two rows have published water-leaf evidence and eight need water-only leaf verification. None is confirmed in the original experimental extracts. Compound selection can be refined after source and extract checks; docking and ADMET remain separate future stages.

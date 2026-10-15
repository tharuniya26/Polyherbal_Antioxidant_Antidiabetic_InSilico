# Project 4 manuscript-readiness assessment

## Decision

**Ready for a carefully limited, reproducible computational methods-and-pose-results draft; not yet ready for a full integrated experimental–computational research article.** The [scientific synthesis](project4_five_plant_scientific_synthesis.md) and [association table](../../data/metadata/project4_five_plant_evidence_summary.csv) provide traceable candidate provenance and preserve both positive and negative pose screens. The original dissertation, raw assay data, detailed aqueous-extraction method, replicate counts and uncertainty estimates are not deposited. Therefore an article cannot yet verify quantitative links to the original in-vitro work or claim an extract-to-ligand mechanism.

**Proposed accurate title:** “Literature-qualified phytochemicals from five medicinal plants: exploratory pose-placement and receptor-sensitivity analysis at human maltase-glucoamylase.” A subtitle may identify the work as a computational extension of an aqueous-leaf polyherbal dissertation. The title should avoid “inhibitors,” “efficacy,” “synergy,” or “validated antidiabetic activity.”

## Claims supported by existing evidence

1. The repository records three aqueous-leaf formulation ratios and an assay inventory for five plants. The README reports qualitative formulation-level DPPH, ABTS, glucose-adsorption and yeast-uptake trends; their source values are unavailable for independent checking.
2. The final preliminary library contains **10 unique structure-assigned references across 11 plant–compound associations** plus one identity-unresolved Aegle association. Botanical occurrence evidence ranges from compound-specific published water-leaf analysis to organic-leaf or mixed-aerial-part evidence. None is a measured analyte in the dissertation's own extracts.
3. The two tested dry/no-NAG 3L4W receptor microstates each recover native neutral MIG at rank 1 in **5/5** seeds in both the original and enlarged box; the enlarged-box control totals **10/10** with RMSD ≤2 Å. This validates a narrow crystal-pose recovery screen only.
4. The original and enlarged mangiferin/quercetin pilots each have **4/8 condition-level geometric PASS calls**, but two particular calls change. Mangiferin and quercetin findings must retain their exact ligand and receptor microstate and box. A PASS is a pocket-localized pose hypothesis, not experimental binding.
5. Gallic acid and agnuside remain **INCONCLUSIVE** under their frozen rank-1 rules in both receptors. Their lower-ranked pocket modes and overlaps with omitted crystallographic waters are reportable structural diagnostics, not rescued PASS outcomes.
6. The current 4GQR/MYC standard-Vina model fails rank-1 native-pose recovery at E32, E64 and the prespecified Glu233H sensitivity test, each **0/5**; it is suspended for phytochemical docking. Boerhavia LIG03 remains on botanical/E–Z-isomer HOLD and was not docked.

All claims must cite the relevant immutable run table and original output provenance, with the date and model version. Exact Vina scores may be reported as *within-run ranking metadata*, without converting them to measured affinity, inhibition constants or potency orders.

## Proposed article structure

| Section | Content that can be written now | Evidence still needed before a fully integrated article |
| --- | --- | --- |
| Introduction | Original aqueous-leaf formulation context; digestive-enzyme rationale; why botanical occurrence, molecular identity and predicted pose are distinct. | Original dissertation text and checked primary background references for final citation. |
| Methods | Five-plant literature screening and evidence tiers; PubChem identity checks; 3D/Meeko preparation; 3L4W/4GQR receptor chemistry; native controls; frozen Vina boxes, seeds, microstates, RMSD mapping, pocket/clash/medoid rule; preserved technical failures. | Complete original extraction and assay methods if an integrated experimental Methods section is claimed. |
| Results | Candidate-evidence matrix; 3L4W/MIG and 4GQR/MYC validation; original-versus-enlarged box pilot; gallic acid and agnuside inconclusive results; off-pocket/water/NAG diagnostics. | Raw experimental observations and statistics for numerical in-vitro results. |
| Discussion | Narrow geometric interpretations, two changed box classifications, receptor/ionization effects, dry-model water conflicts, failed α-amylase model, source-specific botanical gaps. | Experimental profiling or independent inhibition data before proposing compound-level causal mechanisms. |
| Limitations | Missing original source data, no compound detection in project extracts, conditional stereochemistry, receptor chemistry, scoring and sampling limits, technical exit uncertainties. | None needed to state the limitations; they should remain visible after future work. |
| Conclusion | Reproducible native-MIG pose control plus conditional and inconclusive phytochemical placement; no biochemical or clinical claim. | Broader conclusion only after new independent evidence. |

## Required tables and figures from existing data

| Item | Existing source | Required depiction |
| --- | --- | --- |
| Table 1: original formulation and assay context | [`config/formulations.yaml`](../../config/formulations.yaml), [README](../../README.md) | C1/C2/C3 extract ratios and assays; label qualitative findings as unverified source summaries. |
| Table 2: five-plant candidate provenance | [five-plant evidence CSV](../../data/metadata/project4_five_plant_evidence_summary.csv), [ligand library](../../data/compounds/final_preliminary_ligand_library.csv) | All 12 association rows, distinguishing 11 assigned links, ten unique structures and one unresolved identity; exact part/solvent and no original-extract detection. |
| Table 3: native-ligand validation | [3L4W target summary](vina_target_validation_summary.md), [large-box MIG results](../../data/metadata/vina_3L4W_large_box_mig_results.csv), [4GQR run tables](../../data/metadata/vina_4GQR_redocking_results.csv) | Per-complex run counts, rank-1/best-pose RMSD, passes and failures; separate 4GQR E64/Glu233H tables in supplement. |
| Table 4: frozen pose screen | [original pilot CSV](../../data/metadata/vina_3L4W_phytochemical_pilot_results.csv), [enlarged pilot CSV](../../data/metadata/vina_3L4W_box_adequacy_results.csv), [gallic CSV](../../data/metadata/vina_3L4W_gallic_acid_pilot_results.csv), [agnuside CSV](../../data/metadata/vina_3L4W_agnuside_pilot_results.csv) | Separate ligand microstates/receptor states/boxes; show numerator out of five, not a pooled ligand PASS. |
| Figure 1: study-evidence flow | Existing library, audit, and target records | Diagram from source assay context → literature candidate → verified reference structure → native control → conditional pose-only test. No invented compound concentration arrow. |
| Figure 2: paired box outcomes | Original and enlarged pilot CSVs | Eight paired condition calls and box-face margin distributions; preserve two changed calls and exit uncertainties. |
| Figure 3: native controls | 3L4W and 4GQR per-run CSVs | Per-seed crystal-frame RMSD with the prospective 2 Å line; display 4GQR failures. |
| Figure 4: deposited-water diagnostics | [gallic off-pocket table](../../data/metadata/vina_3L4W_gallic_acid_off_pocket_investigation.csv), [agnuside structural table](../../data/metadata/vina_3L4W_agnuside_structural_investigation.csv), deposited 3L4W PDB | Receptor-frame rank-1 and representative lower-ranked pose coordinates beside deposited water oxygen positions; mark omissions, distances, ranks and box faces. |

Figure generation can use existing coordinate and CSV data; **no figure has been created by this Step 69**. A future figure script should preserve the receptor frame, disclose whether a pose is rank 1 or secondary, and use only deposited/reference structures and archived outputs. A static deposited water is not a validated dynamic-solvent representation.

## Claims to exclude

- That gallic acid, caffeic acid, mangiferin, quercetin or agnuside was chemically detected in the original research extracts, or that a formulation ratio is a compound concentration.
- That a Vina PASS or score measures maltase-glucoamylase binding, α-amylase inhibition, glucose uptake, antioxidant activity, antidiabetic efficacy, clinical benefit or polyherbal synergy.
- That 4GQR passed validation, or that its isolated lower-ranked native-like MYC mode rescues the failed 0/5 rank-1 result.
- That one phenolic ionization state, Asp542 state, water occupancy, glycan treatment or glycoside stereoisomer was experimentally confirmed by the docking runs.
- That all five plants are equally docked or equally supported by water-leaf analytical evidence. Boerhavia's LIG03 was held; most reference structures were never docked.
- That DPPH/ABTS activity is reproduced by enzyme docking. No antioxidant quantum-chemical or ADMET analysis has been performed.

## Repository-documentation check and next task

The current [README](../../README.md) predates the computational extension. It says target structures and docking results are merely planned, no analysis scripts or output files exist, and the named directories are empty. Those status claims conflict with the present repository. This assessment **does not edit** the README or prior reports. The first documentation task is to update the README with the actual evidence tiers, prepared files, completed native controls, limited pose-only findings, 4GQR suspension and current limitations, while retaining the original report's authorship and qualitative-experiment caveats.

For a manuscript, next assemble a source-linked draft of Tables 1–4 and the Methods/Results text from existing CSVs, then obtain the original dissertation and raw assay records before writing quantitative experimental Results. Reconcile the original bibliographic details and experimental inconsistencies from those source records; never fill missing values from the docking output. This can proceed without more docking. Any later independent α-amylase complex must pass its own chemically justified native-ligand validation before inclusion.

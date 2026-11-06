# Step 91 revised manuscript-results verification checklist

This checklist separates regenerated computational facts from interpretations that remain experimentally unverified. It applies to `project4_insilico_manuscript_draft.md` and the figures in `results/figures/project4_manuscript/`.

| Check | Status | Evidence |
|---|---|---|
| Candidate-search count is 130 | PASS | Six authoritative source result CSVs; automated final-summary test |
| Candidate-pose count is 2,559 | PASS | Six authoritative all-pose CSVs; automated final-summary test |
| Exit status is 126 zero and four complete-output/unconfirmed | PASS | Archived run rows and attempt metadata; qualification retained in text and Figure 5 |
| Registered native controls total 40 searches | PASS | Eight-row `vina_target_validation_summary.csv`: 25 3L4W and 15 4GQR searches |
| Separate enlarged-box 3L4W controls total 10 searches | PASS | `vina_3L4W_large_box_mig_results.csv`: five seeds in each of two receptor states |
| Combined archived native-control calculations total 50 | PASS | The manuscript and Figures 1–2 retain the 40 + 10 scope distinction |
| Final summary contains 26 unique conditions | PASS | Deterministic final summary and schema test |
| Final classification is 10 PASS / 16 INCONCLUSIVE | PASS | Recalculated from `geometric_classification` in the final summary |
| Four-of-five geometric rule is unchanged | PASS | Authoritative synthesis, manuscript Methods and automated checks agree |
| Gallic acid, rutin and agnuside remain INCONCLUSIVE | PASS | All eight associated final-summary rows remain INCONCLUSIVE |
| Luteolin remains PASS 5/5 in both receptor states | PASS | Both LIG08 final-summary rows |
| Quercetin remains microstate/receptor dependent | PASS | Eight LIG06 rows; only 4′-O−/Asp542H conditions PASS |
| Mangiferin remains ligand-state/receptor/box dependent | PASS | Eight LIG04 rows; six PASS and two INCONCLUSIVE |
| No Vina score is used as an affinity or potency estimate | PASS | Classification uses geometry only; manuscript states the limitation |
| Waters and N-glycans are reported as omitted | PASS | Methods, captions and limitations |
| Structural panels identify PDB, ligand/receptor state, seed and rank | PASS | Figure 2 and Figure 4 captions plus provenance CSV |
| Structural panels preserve archived coordinates | PASS | Figure script projects coordinates without molecular fitting or coordinate writes |
| Hydrogen bonds are not assigned from distance alone | PASS | Figure 4 and its caption explicitly report nearest heavy-atom distances only |
| Figure selection rules are non-cherry-picked | PASS | Figure 2 uses a median rule; Figure 4 uses the frozen medoid; Figures 3/5/6 show all rows or seeds |
| Manuscript body remains within the focused draft range | PASS | Automated test requires 4,000–6,500 words before the References heading |
| References are real and DOI-checked | PASS | Twenty-seven-entry `reference_verification.csv`; primary structural, botanical and method papers prioritized |
| Figure 2 contains no synthetic midpoint marker | PASS | Only observed rank-1 ranges and best returned modes are plotted; automated source test |
| Public availability statement matches baseline commit | PASS | Public archived analysis is distinguished from local untracked manuscript release-candidate files |
| Original-extract compound occurrence | NOT VERIFIED | No new chromatography or mass spectrometry was performed on the original aqueous extracts |
| ntMGAM inhibition, binding affinity and in vivo efficacy | NOT VERIFIED | No new biochemical, biophysical, animal or clinical experiment was performed |
| Experimental repeatability | BLOCKED | Original wet-laboratory experiment was not independently repeated in this computational step |

No docking, rescoring, molecular dynamics or ADMET calculation was executed during manuscript and figure preparation.

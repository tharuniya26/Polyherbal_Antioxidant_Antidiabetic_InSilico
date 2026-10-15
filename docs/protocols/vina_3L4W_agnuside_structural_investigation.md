# Step 68 — Structural investigation of the 3L4W agnuside pilot

## Decision

**B — Retain neutral agnuside as INCONCLUSIVE in both receptor states and proceed to project synthesis.** The completed experiment has **0/5 qualifying rank-1 poses** with Asp542− and **0/5** with Asp542H(OD2), under the unchanged MIG-pocket, severe-clash, core-clustering and 4/5 rules. The search found 13 lower-ranked pocket members, but their existence does not repair rank-1 failure. A defensible water-aware or conformer-specific follow-up is not presently validated: the dry receptors omit waters that overlap many saved poses, and the five starts confound conformer choice with seed. No new docking or score-based model adjustment is recommended from these data.

The [23-row structural table](../../data/metadata/vina_3L4W_agnuside_structural_investigation.csv) records every rank-1 pose and all 13 pocket-localized secondary poses. It links to the preserved [run table](../../data/metadata/vina_3L4W_agnuside_pilot_results.csv), [all-pose table](../../docking/pilot/vina/3L4W_agnuside_pose_only/all_pose_metrics.csv), [pilot report](vina_3L4W_agnuside_pilot_report.md), [frozen plan](vina_3L4W_agnuside_pilot_plan.md) and original deposited [3L4W PDB](../../structures/proteins/original/3L4W.pdb). No earlier result or classification was changed.

## Verification and measurement

All **ten** completed output PDBQT files matched their archived SHA256 values. Their **200** saved modes were re-read through the existing Meeko/RDKit parser; stored scores, pocket flags and box margins were cross-checked. All **70** archived file/hash pairs for receptor, ligand, source SDF, crystal reference, randomized start, log and output matched. The deposited PDB hash also matched the experiment manifest. Each completed search retains its zero exit status. Two earlier seed-4404 attempts were interrupted without a finalized output or confirmed exit; they remain separate from the ten completed searches and were not treated as successful. The later same-input retries are the completed seed-4404 results.

Distances below are minimum **heavy-atom Euclidean distances in the deposited receptor coordinate frame**, without pose fitting. The deposited PDB contributed **6,949 chain-A protein heavy atoms**, **434 water oxygen records** and **84 heavy atoms in six NAG units**. For each pose, protein contacts were grouped by residue; water and NAG distances used deposited coordinates that are absent from the dry docking receptor. A distance is a geometric observation, not an assigned hydrogen bond, favorable interaction, affinity or inhibition. The CSV retains individual output hashes, ranks, coordinates-derived measures, technical attempts and provenance hashes.

## Recurrent rank-1 regions and deposited waters

The ten rank-1 poses occupy related but nonidentical sites adjacent to the MIG reference volume. Their nearest protein residues are primarily **Asp327** (six poses), **Asp203** (two), or **Gln603** (one); the remaining pose is nearest Asp327 as well, giving **Asp327 seven, Asp203 two and Gln603 one**. Other recurrent nearby residues include Thr205, Tyr299, Phe575, His600, Asp542 and Trp406. These residues establish a local receptor-surface neighborhood, not a validated alternative agnuside binding groove. Rank-1 iridoid-core centroids are **5.98–7.63 Å** from the MIG heavy-atom centroid; several have a centroid inside 6 Å yet fail the additional frozen core-overlap condition. Thus “off pocket” here means **outside the prespecified MIG-localization rule**, not necessarily on a remote protein surface.

| Receptor | Seed / SDF record | Nearest protein residue (Å) | Nearest deposited water O (Å) | Distinct water-O overlaps <1.8 Å | Nearest box face / margin (Å) |
| --- | --- | --- | --- | ---: | --- |
| Asp542− | 1101 / 1 | Asp203 (2.708) | A:4113 (1.200) | 3 | Y− / 3.860 |
| Asp542− | 2202 / 5 | Asp203 (2.786) | A:4268 (1.109) | 1 | Y− / 2.514 |
| Asp542− | 3303 / 8 | Gln603 (2.704) | A:4113 (1.298) | 3 | Y− / 4.018 |
| Asp542− | 4404 / 11 | Asp327 (2.708) | A:4447 (1.322) | 3 | Y− / 1.087 |
| Asp542− | 5505 / 17 | Asp327 (2.775) | A:4113 (0.746) | 5 | Y− / 3.137 |
| Asp542H(OD2) | 1101 / 1 | Asp327 (2.636) | A:4412 (0.631) | 4 | Y− / 3.951 |
| Asp542H(OD2) | 2202 / 5 | Asp327 (2.731) | A:4447 (1.133) | 3 | **Y− / 0.859** |
| Asp542H(OD2) | 3303 / 8 | Asp327 (2.656) | A:4412 (0.629) | 4 | Y− / 3.853 |
| Asp542H(OD2) | 4404 / 11 | Asp327 (2.702) | A:4447 (1.222) | 3 | Y− / 1.044 |
| Asp542H(OD2) | 5505 / 17 | Asp327 (2.743) | A:4113 (0.795) | 5 | Y− / 3.140 |

All ten rank-1 poses overlap at least one omitted deposited water oxygen by the prespecified **<1.8 Å** diagnostic; the smallest separation is **0.629 Å**. Across all saved modes, **186/200** have such an overlap. These are severe geometric incompatibilities *if those deposited water positions are retained unchanged*; they are not clashes with the actual dry Vina receptor. The pattern could reflect omitted-water sensitivity, possible water displacement, or limitations of the static crystallographic solvent model. It does not establish which waters would persist or be displaced on binding. The nearest search-box face is Y− for all rank-1 poses, but only **one** has the predefined <1 Å face warning; **15/200** modes carry that warning overall. The data do not support treating the Y− Vina search face as a biochemical pocket boundary or concluding that box truncation alone explains the result.

## All 13 secondary MIG-pocket poses

Six pocket members arose with Asp542− (seeds 2202, 3303 and 4404) and seven with Asp542H(OD2) (seeds 1101, 2202, 3303 and 4404). No mode from seed 5505 met pocket membership in either receptor. The table gives each saved pocket mode, not a selected favorable subset. Δscore is that mode's Vina score minus its own run's rank-1 score; positive values mean Vina ranked it less favorably. Catalytic distances are the smallest heavy-atom separations to Asp443, Asp542 or His600, not assigned catalytic interactions. Exact per-residue distances are in the CSV.

| Receptor | Seed | Rank | Δscore | Nearest protein (Å) | Nearest water O (Å); overlaps <1.8 Å | Nearest catalytic residue (Å) | Box margin (Å) |
| --- | ---: | ---: | ---: | --- | --- | --- | ---: |
| Asp542− | 2202 | 18 | +0.623 | Asp542 (2.727) | A:4304 (0.906); 4 | Asp542 (2.727) | 3.399 |
| Asp542− | 3303 | 12 | +0.468 | Asp542 (3.001) | A:4304 (0.952); 3 | Asp542 (3.001) | 3.715 |
| Asp542− | 3303 | 19 | +0.736 | Asp542 (2.767) | A:4304 (0.391); 2 | Asp542 (2.767) | 1.716 |
| Asp542− | 4404 | 11 | +0.493 | His600 (3.070) | A:4412 (2.723); 0 | His600 (3.070) | 2.922 |
| Asp542− | 4404 | 14 | +0.559 | Asp542 (2.720) | A:4304 (0.283); 3 | Asp542 (2.720) | 4.283 |
| Asp542− | 4404 | 16 | +0.635 | His600 (2.806) | A:4304 (0.768); 2 | His600 (2.806) | 2.791 |
| Asp542H(OD2) | 1101 | 17 | +0.600 | Asp542 (2.875) | A:4304 (1.318); 2 | Asp542 (2.875) | 3.791 |
| Asp542H(OD2) | 2202 | 6 | +0.292 | Asp542 (2.720) | A:4304 (0.563); 3 | Asp542 (2.720) | 4.136 |
| Asp542H(OD2) | 2202 | 15 | +0.640 | Thr205 (2.701) | A:4447 (0.460); 2 | Asp542 (2.777) | 3.390 |
| Asp542H(OD2) | 3303 | 5 | +0.334 | His600 (2.811) | A:4304 (0.892); 3 | His600 (2.811) | 3.717 |
| Asp542H(OD2) | 4404 | 5 | +0.350 | Asp542 (2.763) | A:4304 (0.322); 1 | Asp542 (2.763) | 4.194 |
| Asp542H(OD2) | 4404 | 8 | +0.420 | Asp327 (2.699) | A:4268 (0.425); 1 | His600 (2.757) | 1.816 |
| Asp542H(OD2) | 4404 | 19 | +0.784 | Asp542 (2.829) | A:4447 (1.019); 2 | Asp542 (2.829) | 2.641 |

The alternative poses are often closer to Asp542/His600 than the corresponding rank-1 positions and have **1.72–4.28 Å** minimum face margins. **Twelve of thirteen** overlap at least one deposited water oxygen, usually A:4304; the sole nonoverlapping mode is Asp542− seed 4404 rank 11. Their rank penalties span **+0.292 to +0.784 Vina score units**. This is evidence of a *ranking difference in this dry Vina model*, not evidence that a lower-ranked pose is physically preferable. Neither catalytic proximity nor pocket membership alone validates a pose.

## Omitted NAG units and conformer/seed pattern

The original PDB has six NAG units. Deposited `LINK` records connect **NAG A2005 to Asn209**, **A2006 to Asn372**, **C1–C2 to Asn393**, and **B1–B2 to Asn741**. Thus these are linked glycans rather than six independent free ligands. Across the 23 decision-relevant poses, the closest omitted NAG atom is always in **A2005**. The nearest NAG separations are **7.35–10.60 Å** for rank-1 poses and **7.54–12.03 Å** for pocket alternatives; no NAG atom lies within **4 Å** of any of these poses. The deposited glycan coordinates therefore show no direct local steric overlap with these saved poses. Omission still leaves possible indirect effects on receptor conformations and solvent untested; these distances cannot validate a no-NAG receptor model.

The seed/conformer allocation was fixed as **1101/record 1, 2202/5, 3303/8, 4404/11, 5505/17** and the same start was used across receptor states. Rank-1 cores from seeds 1101 and 3303 form one close pair (**0.19 Å** Asp542−; **0.14 Å** Asp542H); seeds 2202 and 4404 form another (**0.78 Å** and **0.20 Å**). The between-pair separation is approximately **3–4 Å**, and seed 5505 does not join either close pair. Under the Asp542H state, the 1101/3303 positions have MIG-centroid distances **7.60–7.63 Å**, whereas the 2202/4404 positions are **5.98–6.13 Å** but still fail the overlap portion of the frozen pocket rule. This pattern is consistent with multiple accessible placements. It **cannot attribute cause to starting conformer versus random seed**, because each conformer was paired with only one seed. Five starts also do not prove adequate sampling of a 13-torsion glycoside.

## Follow-up feasibility and scientific boundary

A water-aware diagnostic would require a documented rule for which crystallographic waters to retain, chemically credible hydrogen orientations or a method allowing water displacement, atom typing and scoring support, and an independent native-ligand control under that exact receptor model. Simply adding static water oxygens to standard Vina or deleting overlapping waters after seeing poses would not establish solvent sensitivity. A conformer-sensitivity diagnostic would need a crossed conformer-by-seed design fixed before any scores are viewed, enough searches to assess a 13-torsion ligand, and the same independent chemistry and pocket rules. Neither protocol is currently qualified by the saved artifacts. The existing native-MIG enlarged-box control passed **10/10** pose-recovery runs, but a smaller native ligand cannot establish agnuside glycoside sampling or water physics.

Accordingly, **no additional agnuside run is justified now**. Retain **INCONCLUSIVE 0/5 in each receptor state**, the 13 secondary pocket modes, ten rank-1 water overlaps, the seed-4404 interruption/retry history and all source hashes as reported. LIG09 is a stereodefined PubChem reference; source-specific agnuside stereochemistry and identification in the project's original aqueous *Vitex negundo* extracts were not established. **LIG03 trans-caffeic acid from *Boerhavia diffusa* remains on its botanical/E–Z-isomer HOLD.** These poses cannot establish binding affinity, enzyme inhibition, comparative potency, antidiabetic effect or formulation synergy. Dry/no-NAG receptor construction, omitted water, provisional Asp542/His600 protonation and rigid-receptor Vina scoring remain limitations.

**Next concrete task:** synthesize the existing 3L4W pose-only findings and compound-occurrence evidence across the studied plants, preserving per-compound inconclusive outcomes and the LIG03 HOLD. Any later water-aware or crossed-conformer experiment should first produce an independent, prespecified preparation and native-control validation protocol; it is not part of this Step 68 decision.

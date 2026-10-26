# Project 4 final five-plant scientific synthesis

## Scope and decision language

This document consolidates the completed molecular-docking investigations preserved through Git commit `8253265`. It covers six unique literature-associated reference compounds across five plants: gallic acid and rutin for *Aegle marmelos*; the existing quercetin molecular results shared by *Boerhavia diffusa* and *Senna auriculata*; mangiferin for *Salacia reticulata*; and agnuside and luteolin for *Vitex negundo*. The Boerhavia association adds botanical provenance to LIG06 and **does not constitute another docking experiment**.

The authoritative quantitative sources are the archived run CSVs, all-pose tables and experiment manifests listed in the final [26-row machine-readable synthesis](../../data/metadata/project4_final_five_plant_docking_summary.csv). That table is deterministically rebuilt by [`scripts/build_project4_final_five_plant_summary.py`](../../scripts/build_project4_final_five_plant_summary.py). The generator recalculates the frozen rank-1 decision from the saved pocket, clash and medoid-RMSD fields and stops if it disagrees with a stored classification.

In this report, **PASS** means only that at least four of five prespecified rank-1 poses met the frozen geometric rule in one particular ligand state, receptor state and box. It does not mean experimentally demonstrated binding, inhibition, potency, antidiabetic activity, clinical efficacy or extract composition. **INCONCLUSIVE** is not evidence of inactivity. No named reference compound was analytically demonstrated in the original project's aqueous leaf extracts.

## Source verification and frozen rule

The source audit reconciled **130 analyzable candidate searches** and **2,559 valid returned poses**:

- Original-box mangiferin/quercetin pilot: 40 complete outputs, 781 poses, 38 confirmed zero exits and two complete outputs with unconfirmed exits.
- Enlarged-box mangiferin/quercetin comparison: 40 complete outputs, 798 poses, 38 confirmed zero exits and two different complete outputs with unconfirmed exits.
- Gallic-acid pilot: 20/20 confirmed zero exits and 382 poses.
- Agnuside pilot: 10/10 confirmed zero exits and 200 poses; two earlier interrupted attempts remain separate technical records.
- Rutin/luteolin pilot: 20/20 confirmed zero exits and 398 poses.

The recomputation found **zero classification mismatches**. Across the 26 unique molecular conditions, 10 rows are conditional geometric PASS and 16 are INCONCLUSIVE. Quercetin rows are molecular conditions shared between two botanical associations and must not be double-counted.

The fixed candidate rule required every ligand heavy atom inside the box; rigid-core centroid no more than 6.0 Å from deposited miglitol (MIG); at least half of core atoms within 4.0 Å of a MIG heavy atom; no receptor–ligand distance below 1.8 Å, van-der-Waals overlap above 0.8 Å or specified internal nonbonded distance below 1.6 Å; and no more than 2.0 Å unfitted symmetry-aware rigid-core RMSD from the within-condition rank-1 medoid. At least 4/5 rank-1 seeds had to pass all tests. A box-face margin below 1.0 Å and deposited-water overlaps were diagnostics, not post hoc decision criteria. Scores never entered the classification.

The principal human-readable sources are the [target-validation summary](vina_target_validation_summary.md), [receptor-preparation assessment](vina_receptor_preparation.md), [ligand-preparation assessment](vina_ligand_preparation.md), [integrated mangiferin/quercetin assessment](vina_3L4W_integrated_pilot_assessment.md), [gallic-acid pilot](vina_3L4W_gallic_acid_pilot_report.md), [gallic-acid structural investigation](vina_3L4W_gallic_acid_off_pocket_investigation.md), [agnuside pilot](vina_3L4W_agnuside_pilot_report.md), [agnuside structural investigation](vina_3L4W_agnuside_structural_investigation.md), [rutin/luteolin pilot](vina_3L4W_alternative_ligand_pilot_report.md), and [alternative-compound evidence review](project4_alternative_compound_selection.md). Their run-level and per-pose tables remain the numerical authority.

## Target validation and model quality

### 3L4W native-MIG control

Human maltase-glucoamylase 3L4W is the only target advanced for candidate pose placement. Neutral MIG recovered a rank-1 pose within the prespecified 2.0 Å crystal-frame RMSD threshold in **5/5 seeds** with the original Asp542− receptor at exhaustiveness 32 and 64, **5/5** with the Asp542H(OD2) receptor, and **10/10** across both receptors in the enlarged 30 × 24 × 24 Å box. Original-box neutral rank-1 RMSDs were 1.8629–1.8849 Å; enlarged-box values were 1.8671–1.9362 Å. More crystal-like modes near 0.6–0.7 Å sometimes ranked second, warning that native-like geometry was not always scored best.

This is reproducible native-pose placement for MIG in the tested local models. It does not validate binding energies, phytochemical score ranking, protonation, water displacement or inhibition.

### Failed 4GQR control

The human pancreatic α-amylase 4GQR/MYC model failed the same prospective rank-1 screen: **0/5** at exhaustiveness 32, **0/5** at exhaustiveness 64 and **0/5** after the Glu233H(OE2) sensitivity change. The best E32 native-like pose was rank 17. The current 4GQR model therefore remains suspended for candidate docking and supplies negative validation evidence rather than a second target supporting the plant hypotheses.

### Receptor and binding-site assumptions

Candidate studies used two rigid, dry/no-NAG pH-6.5 hypotheses: Asp443−/Asp542−/His600-HID and Asp443−/Asp542H(OD2)/His600-HID. The six deposited linked NAG units and crystallographic waters were omitted. The site and box centre were defined by the deposited 14-heavy-atom MIG centroid `(45.38664286, 92.35157143, 34.82200000) Å`. The original candidate box was 24 × 18 × 18 Å; the preregistered enlarged comparison and later pilots used 30 × 24 × 24 Å.

The original box was not geometry-independent: 36/40 original mangiferin/quercetin rank-1 poses were within 1 Å of a face, versus 1/40 in the enlarged comparison, and two of eight classifications changed. The larger box passed the native-MIG control and reduced edge crowding, but it is only one qualified exploratory search space, not an optimized physical binding model.

## Quantitative synthesis

`Exit 0/U` gives confirmed zero exits and complete outputs with unconfirmed exits. Mean Vina scores are descriptive output values in kcal/mol and must not be interpreted as measured affinities or biological potency. The full CSV also records median score, returned-pose count, source hashes, water provenance and detailed limitations.

| Plant association | Compound and state | Receptor | Box (Å) | Exit 0/U | Pocket rank 1 | Full passes | Mean Vina | Core RMSD max (Å) | Classification | Main condition limit |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | --- |
| *Aegle marmelos* | Gallic acid, neutral acid | Asp542− | 30×24×24 | 5/0 | 0/5 | 0/5 | −6.413 | 0.0794 | INCONCLUSIVE | Tight off-pocket Z− face cluster |
| *Aegle marmelos* | Gallic acid, neutral acid | Asp542H | 30×24×24 | 5/0 | 0/5 | 0/5 | −6.418 | 0.0716 | INCONCLUSIVE | Tight off-pocket Z− face cluster |
| *Aegle marmelos* | Gallic acid, carboxylate | Asp542− | 30×24×24 | 5/0 | 0/5 | 0/5 | −6.457 | 0.0487 | INCONCLUSIVE | Tight off-pocket Z− face cluster |
| *Aegle marmelos* | Gallic acid, carboxylate | Asp542H | 30×24×24 | 5/0 | 0/5 | 0/5 | −6.447 | 0.0660 | INCONCLUSIVE | Tight off-pocket Z− face cluster |
| *Aegle marmelos* | Rutin, neutral reference | Asp542− | 30×24×24 | 5/0 | 1/5 | 0/5 | −8.756 | 10.9953 | INCONCLUSIVE | Flexible glycoside; poor rank-1 localization/convergence |
| *Aegle marmelos* | Rutin, neutral reference | Asp542H | 30×24×24 | 5/0 | 1/5 | 0/5 | −8.944 | 11.1674 | INCONCLUSIVE | Flexible glycoside; poor rank-1 localization/convergence |
| *Boerhavia diffusa*; *Senna auriculata* | Quercetin, neutral | Asp542− | 24×18×18 | 5/0 | 0/5 | 0/5 | −7.158 | 0.0370 | INCONCLUSIVE | Tight but outside frozen MIG-overlap rule |
| Same shared LIG06 experiment | Quercetin, neutral | Asp542− | 30×24×24 | 5/0 | 0/5 | 0/5 | −7.160 | 0.0459 | INCONCLUSIVE | Enlarging box does not rescue localization |
| Same shared LIG06 experiment | Quercetin, neutral | Asp542H | 24×18×18 | 5/0 | 0/5 | 0/5 | −7.212 | 0.0849 | INCONCLUSIVE | Tight but outside frozen MIG-overlap rule |
| Same shared LIG06 experiment | Quercetin, neutral | Asp542H | 30×24×24 | 5/0 | 0/5 | 0/5 | −7.221 | 0.0943 | INCONCLUSIVE | Enlarging box does not rescue localization |
| Same shared LIG06 experiment | Quercetin, 4′-O monoanion | Asp542− | 24×18×18 | 5/0 | 0/5 | 0/5 | −6.865 | 2.0543 | INCONCLUSIVE | Receptor-state-dependent off-rule placement |
| Same shared LIG06 experiment | Quercetin, 4′-O monoanion | Asp542− | 30×24×24 | 5/0 | 0/5 | 0/5 | −6.850 | 2.0851 | INCONCLUSIVE | Receptor-state-dependent off-rule placement |
| Same shared LIG06 experiment | Quercetin, 4′-O monoanion | Asp542H | 24×18×18 | 5/0 | 5/5 | 5/5 | −7.227 | 0.0653 | PASS | Restricted to one receptor/microstate hypothesis |
| Same shared LIG06 experiment | Quercetin, 4′-O monoanion | Asp542H | 30×24×24 | 5/0 | 5/5 | 5/5 | −7.228 | 0.0301 | PASS | Stable across boxes, not across receptor states |
| *Salacia reticulata* | Mangiferin, neutral | Asp542− | 24×18×18 | 3/2 | 5/5 | 5/5 | −6.899 | 0.0382 | PASS | Original-box face crowding |
| *Salacia reticulata* | Mangiferin, neutral | Asp542− | 30×24×24 | 5/0 | 5/5 | 5/5 | −7.494 | 0.0500 | PASS | Site persists but cores shift 1.83–1.87 Å |
| *Salacia reticulata* | Mangiferin, neutral | Asp542H | 24×18×18 | 5/0 | 5/5 | 5/5 | −7.131 | 0.0541 | PASS | Does not survive enlarged-box test |
| *Salacia reticulata* | Mangiferin, neutral | Asp542H | 30×24×24 | 3/2 | 5/5 | 3/5 | −7.686 | 2.3680 | INCONCLUSIVE | Five pocket poses split into 3+2 cores |
| *Salacia reticulata* | Mangiferin, 6-O monoanion | Asp542− | 24×18×18 | 5/0 | 0/5 | 0/5 | −6.423 | 0.0973 | INCONCLUSIVE | Strong box dependence |
| *Salacia reticulata* | Mangiferin, 6-O monoanion | Asp542− | 30×24×24 | 5/0 | 4/5 | 4/5 | −7.230 | 6.4997 | PASS | One off-pocket, near-face seed remains |
| *Salacia reticulata* | Mangiferin, 6-O monoanion | Asp542H | 24×18×18 | 5/0 | 5/5 | 5/5 | −6.791 | 0.0560 | PASS | Detailed pose changes with box |
| *Salacia reticulata* | Mangiferin, 6-O monoanion | Asp542H | 30×24×24 | 5/0 | 5/5 | 5/5 | −7.462 | 0.0442 | PASS | Core shifts 2.43–2.45 Å across boxes |
| *Vitex negundo* | Agnuside, neutral reference | Asp542− | 30×24×24 | 5/0 | 0/5 | 0/5 | −7.886 | 3.0669 | INCONCLUSIVE | Conformer and seed effects are confounded |
| *Vitex negundo* | Agnuside, neutral reference | Asp542H | 30×24×24 | 5/0 | 0/5 | 0/5 | −7.922 | 4.1802 | INCONCLUSIVE | Conformer and seed effects are confounded |
| *Vitex negundo* | Luteolin, neutral reference | Asp542− | 30×24×24 | 5/0 | 5/5 | 5/5 | −7.948 | 0.0592 | PASS | Only neutral ligand state tested |
| *Vitex negundo* | Luteolin, neutral reference | Asp542H | 30×24×24 | 5/0 | 5/5 | 5/5 | −8.174 | 0.0295 | PASS | Only neutral ligand state tested |

## Five-plant interpretation

### 1. *Aegle marmelos*: gallic acid and rutin

Both associations are **Tier B**: the principal Aegle leaf evidence uses ethanol, with additional organic or aqueous-organic support for rutin. Neither compound was measured in the project's aqueous Aegle extract.

Gallic acid was tested as a neutral-acid reference and a carboxylate pH-6.5 hypothesis. All four receptor/charge combinations are **INCONCLUSIVE 0/5**. Their small within-condition RMSDs describe repeatable ranking of the same off-pocket placement, not pocket convergence. Every rank-1 pose lies 0.008–0.067 Å from the enlarged box's Z− face and conflicts with six deposited water positions; every run also contains a rank-2/3 MIG-pocket alternative with no <1.8 Å deposited-water overlap. Selecting those lower-ranked alternatives would violate the prospective rank-1 rule.

Neutral rutin is **INCONCLUSIVE 0/5** in both receptors. Only seed 2202 is pocket-localized in each state, but it is more than 8.5 Å from the medoid; other rank-1 cores are outside the frozen pocket rule. Forty-four of 200 modes are pocket members, 42 at lower ranks. Four of five rank-1 poses in each receptor overlap deposited water positions, and three of five have a box-face warning. Five conformer-linked starts cannot fully characterize a 16-torsion glycoside.

### 2. *Boerhavia diffusa*: qualified quercetin association

The quercetin occurrence tier is **A with analytical qualification**: a separately prepared boiled-water leaf extract received an HPLC-DAD quercetin assignment, with earlier HPLC-PAD-ESI-MS/MS leaf profiling as context. HPLC-DAD is less definitive than purified-compound NMR, and this literature evidence does not establish quercetin in the project's aqueous Boerhavia extract.

No Boerhavia-specific docking was performed. The same LIG06 neutral and 4′-O monoanion results already generated for the molecular reference are reused. Neutral quercetin is **INCONCLUSIVE 0/5** for both receptors and both boxes: the rank-1 clusters are tight but have only 27.27% rigid-core MIG overlap, below the frozen 50% threshold. The monoanion is **INCONCLUSIVE 0/5** with Asp542− but **PASS 5/5** with Asp542H(OD2) in both boxes. That persistent PASS is a restricted receptor/microstate hypothesis; the cross-receptor medoids remain about 5.7 Å apart. Deposited-water overlap was not quantified in the Step 58/59 source tables.

### 3. *Salacia reticulata*: mangiferin

Mangiferin has **Tier A** published compound-specific water-leaf HPLC evidence. The source did not independently resolve every stereocentre, and the compound was not measured in the project's own Salacia extract.

The neutral and modeled xanthone 6-O monoanion states show substantial receptor and box dependence. Neutral/Asp542− is the only mangiferin condition with **PASS 5/5 in both boxes** and persistent site localization, but same-seed cores shift about 1.83–1.87 Å. Neutral/Asp542H changes from **PASS 5/5** to **INCONCLUSIVE 3/5** because the enlarged-box rank-1 poses split into 3+2 core placements. The monoanion/Asp542− changes from **INCONCLUSIVE 0/5** to **PASS 4/5**, with one enlarged-box off-pocket, near-face seed. Monoanion/Asp542H remains **PASS 5/5**, but the core shifts 2.43–2.45 Å and the contact pattern changes. Mangiferin therefore supplies conditional site hypotheses, not one box- and receptor-independent pose. Deposited-water overlap was not quantified in these earlier candidate tables.

### 4. *Senna auriculata*: quercetin

The Senna association is **Tier B**, based on leaf isolation after 80% aqueous-acetone extraction and fractionation. It is not water-only evidence and does not establish quercetin in the project's aqueous Senna extract.

The molecular results are exactly the shared LIG06 results described for Boerhavia. They count once scientifically: neutral quercetin remains INCONCLUSIVE in all tested submodels; the 4′-O monoanion supplies a reproducible pocket placement only with the Asp542H receptor. The computation is not an independent experiment for each plant and cannot discriminate botanical source, extract recovery or concentration.

### 5. *Vitex negundo*: agnuside and luteolin

Agnuside is **Tier B with abstract-only aqueous evidence**: full water-only chromatographic detail and source stereochemistry remain unresolved. Luteolin is **Tier B**, supported by authenticated leaves extracted with organic solvents and compared with an HPLC standard. Neither was detected in the original aqueous Vitex extract.

Neutral agnuside is **INCONCLUSIVE 0/5** in both receptor states. Thirteen of 200 lower-ranked modes are pocket members, but none is rank 1. All ten rank-1 poses and 186/200 modes overlap at least one deposited water oxygen within 1.8 Å; 12/13 secondary pocket modes also overlap. Each of five conformers was paired with only one seed, confounding conformer and stochastic-search effects for a 13-torsion glycoside. The saved poses are too uncertain for a mechanistic binding-site claim.

Neutral luteolin is **PASS 5/5** in both receptor states. All ten rank-1 poses are pocket-localized and clash-free; maximum within-state core RMSD is 0.0592 Å and 0.0295 Å, minimum face clearance exceeds 4.3 Å, and no rank-1 pose overlaps a deposited water position within 1.8 Å. Cross-receptor medoid RMSD is 0.1166 Å, with matched-seed values 0.0836–0.1376 Å. This is the strongest cross-receptor geometric result, but only the neutral reference was tested and the result still depends on the dry rigid receptor and Vina scoring.

## Water, clashes, conformers and scoring

No rank-1 candidate pose triggered the frozen severe-clash rule. Six lower-ranked severe-clash modes occurred among the 781 original-box mangiferin/quercetin poses; the enlarged comparison had none. Gallic acid, agnuside, rutin and luteolin received explicit deposited-water distance audits because their dry receptors omitted those waters. Mangiferin and quercetin did not receive an equivalent per-pose water audit; no overlap count should be inferred for them.

Water overlap is especially limiting for the gallic-acid off-pocket rank-1 site and agnuside. Rutin also shows frequent overlap, whereas luteolin has none at rank 1. These coordinate comparisons cannot decide whether a real ligand would retain, reorganize or displace a water. A valid water-aware study would need prespecified water selection, credible hydrogen orientations or a displacement-capable method, compatible scoring and a fresh native-ligand control.

Flexible-glycoside results remain sampling-limited. Rutin uses five selected conformers across seeds; agnuside assigns five conformers one-to-one with five seeds. The agnuside design cannot separate conformer from random-seed effects. Native recovery of smaller MIG does not validate glycoside conformational coverage.

Ordinary Vina scores are empirical ranking outputs, not binding free energies. Scores cannot select a ligand ionization state, receptor protonation state, botanical candidate, potency or biological mechanism. The summary retains score means for traceability only and does not rank compounds by them.

## Answers to the final scientific questions

1. **Reliably reproduced:** all archived inputs, commands, outputs and classification fields are hash-linked; the frozen decision rule reproduces every stored candidate classification; 3L4W/MIG rank-1 pose recovery is repeatable in both receptor hypotheses and both tested box regimes.
2. **Conditional geometric PASS:** neutral luteolin with both receptors; quercetin 4′-O monoanion with Asp542H in both boxes; mangiferin neutral/Asp542− in both boxes; mangiferin monoanion/Asp542H in both boxes; and two box-specific mangiferin outcomes. These are submodel-specific geometric hypotheses.
3. **Inconclusive:** all gallic-acid states, neutral rutin, neutral agnuside, neutral quercetin, quercetin monoanion/Asp542−, and the stated box-specific mangiferin conditions.
4. **Dependence on modelling choices:** strong. Two of eight mangiferin/quercetin calls change with box size; quercetin monoanion is receptor-state dependent; mangiferin changes with receptor, charge and box; only the tested neutral luteolin placement is geometrically stable across the two receptor hypotheses in the enlarged box.
5. **Enzyme-inhibition conclusion:** none. Docking supplies no measured inhibition, affinity, kinetics, dose response or causal link to glucose assays.
6. **Suitable undergraduate-manuscript content:** the transparent native-control strategy, negative 4GQR result, frozen geometric rule, complete positive and negative candidate outcomes, box-sensitivity experiment, hash/provenance audit and explicit extract/docking boundary are suitable. Claims of potency, synergy, clinical effect or compound detection are not.
7. **Further validation required:** chromatographic and preferably MS/MS/standard confirmation of candidates in the actual aqueous extracts; direct human ntMGAM inhibition assays on identified compounds and extracts with appropriate controls; concentration and replicate uncertainty; then, if computational refinement is justified, prespecified water/glycan-aware receptor preparation, experimentally informed protonation, crossed conformer-by-seed sampling for flexible glycosides, and an independently validated analysis protocol.

## Reproducibility assessment

Recorded execution versions are AutoDock Vina **1.2.7**, Meeko **0.8.0**, RDKit **2026.03.6** and Python **3.12.10**. The archived Vina executable SHA256 is `e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5`. The two principal receptor PDBQT hashes are `7a76e0fa3b7bf2197714a57a88c60d2f4eefedccb18636951f3f643473afadf3` for Asp542− and `220e3585558e4bc261e6e76ddb036738dea827543e1b135399506e95ea9f93fd` for Asp542H(OD2). The deposited MIG-reference SDF hash is `ea673a59043bea05abbe7d66d8bd29c3786021f564c394027a8b0201385503f7`; the original 3L4W PDB hash used for water diagnostics is `57638e45666e2f848e37691a914a978d85c377c13eaed9e5f94af24d2c52c95e`.

The source manifests preserve receptor, ligand, source-SDF, randomized-start, output, log, Vina and analysis-table hashes. Archived receptor structures, ligand preparation tables, atom maps, frozen seeds, commands, box definitions and analysis thresholds are present. The final CSV embeds current SHA256 values for every source results table, all-pose table and manifest.

Independent clean-environment reproduction is **not yet demonstrated**. The repository's `requirements.txt` declares only PyYAML 6.0.3, while the archived workflows also import NumPy, SciPy, RDKit, Meeko and Gemmi. Exact NumPy, SciPy and Gemmi versions are not consistently recorded in the experiment manifests. Vina is an external Windows executable rather than a repository dependency, and some saved commands contain machine-specific absolute paths. Before a clean reproduction, the project needs a tested environment/constraints file, installation instructions for Vina 1.2.7 and scientific Python dependencies, platform-neutral path handling, and a read-only verification command that checks archived hashes and rebuilds derived tables without launching docking.

## Defensible conclusion

The strongest result is reproducible **geometric placement**, not demonstrated inhibition: neutral luteolin localizes consistently in the prespecified MIG pocket across both tested 3L4W receptor hypotheses in the enlarged box. Quercetin monoanion/Asp542H and selected mangiferin submodels also produce repeatable placements, but their receptor, ligand-state or box dependence is scientifically material. Gallic acid, rutin and agnuside remain inconclusive under the predefined rank-1 rule.

The project supports an undergraduate computational manuscript when presented as a validation-aware, hypothesis-generating docking study with openly retained negative and conditional results. It does not support claims that the six references occur in the original extracts, inhibit ntMGAM, explain the original antioxidant or glucose-response assays, act synergistically, or have clinical antidiabetic efficacy.

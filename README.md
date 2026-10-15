# Integrated in-vitro and in-silico investigation of five polyherbal medicinal plants

**Research:** *Integrated In-Vitro and In-Silico Investigation of Antioxidant and Antidiabetic-Related Properties of Polyherbal Formulations Containing Aegle marmelos, Boerhavia diffusa, Salacia reticulata, Senna auriculata, and Vitex negundo*

**Researcher:** Tharuniya Anpalakan · **Institution:** Northumbria University, Department of Applied Sciences

**Original in-vitro report date:** 11 August 2026

## Scope and current status

The original study compared **aqueous leaf extracts** of five medicinal plants in three polyherbal formulations using phytochemical, antioxidant-related and glucose-related in-vitro assays. This repository extends that work with an evidence-tiered phytochemical reference library, prepared molecular structures, receptor assessments, native-ligand pose-recovery tests and limited **AutoDock Vina 1.2.7 exploratory pose-placement experiments**. The computational work is a hypothesis-generating extension, not a measurement of enzyme inhibition or clinical antidiabetic efficacy.

The original dissertation, raw experimental measurements and full extraction/assay records are **not deposited here**. The qualitative experimental comparisons below are transcribed project background and cannot be independently recalculated from this repository. **No named phytochemical was analytically confirmed in the original project's aqueous extracts**, and the formulation ratios are proportions of *extracts*, not concentrations of individual molecules.

## Original experimental context

Plant codes and formulation ratios come from [the configuration transcription](config/formulations.yaml). A dash means that plant is absent; each formulation totals ten unspecified ratio parts.

| Formulation | *Aegle marmelos* (AM) | *Boerhavia diffusa* (BD) | *Senna auriculata* (SA) | *Salacia reticulata* (SR) | *Vitex negundo* (VN) | Recorded basis |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| C1 | 2.0 | 1.5 | 2.0 | 2.0 | 2.5 | Equal TPC/TFC/TAC weighting |
| C2 | 2.0 | 1.0 | 2.5 | 2.0 | 2.5 | Greater TAC emphasis; TPC/TFC/TAC weights 0.25/0.25/0.50 |
| C3 | — | — | 3.5 | 3.5 | 3.0 | Selected plants with higher antioxidant-related measurements |

The original work used aqueous **leaf** extracts and recorded qualitative phytochemical screening; total phenolic content (TPC), total flavonoid content (TFC) and total antioxidant capacity (TAC); DPPH and ABTS radical-scavenging assays; glucose adsorption; and yeast glucose uptake. The existing README transcription reports C1 as strongest for DPPH, C2 generally highest for glucose adsorption, and C3 strongest for ABTS with generally high yeast uptake. Those are **reported qualitative trends**, without deposited raw values, replicate counts, uncertainties or statistical tests. Extraction temperature, duration, leaf-to-water ratio, yield and some reported data inconsistencies also remain unreconciled against the original source records.

## Literature-qualified compound library

The [final preliminary library](data/compounds/final_preliminary_ligand_library.csv) holds **10 assigned reference chemical structures across 11 assigned plant–compound associations**, plus **one unresolved Aegle kaempferol-rutinoside association without an assigned ligand**. Rutin (LIG02) is one structure associated with both Aegle and Senna; it is not two different ligands. The [five-plant evidence summary](data/metadata/project4_five_plant_evidence_summary.csv) preserves plant part, extraction solvent, literature tier, structure uncertainty and docking status for every association.

| Plant | Selected reference(s) and occurrence boundary | 3L4W status |
| --- | --- | --- |
| *Aegle marmelos* | Gallic acid and rutin reported from **ethanol-extracted leaves**; a separate kaempferol-rutinoside assignment remains identity-unresolved. | Gallic acid tested; all four neutral/carboxylate × receptor conditions **INCONCLUSIVE**. Rutin and tentative glycoside untested. |
| *Boerhavia diffusa* | Underivatized caffeic acid reported from a boiled-water **mixed-aerial-part** infusion; methanolic-leaf GC-MS identified a **tri-TMS derivative**. Source E/Z geometry and isolated aqueous-leaf identity are unresolved. | Trans-caffeic-acid PubChem reference **on botanical/isomer HOLD; not docked**. |
| *Salacia reticulata* | Mangiferin has compound-specific published **water-leaf** HPLC evidence. Hot-water leaves support epicatechin at the compound level, without proof of the exact PubChem optical isomer. | Mangiferin neutral and 6-O monoanion tested; results depend on receptor state and box. Epicatechin untested. |
| *Senna auriculata* | Quercetin, kaempferol, luteolin, rutin and nicotiflorin have leaf literature using **80% aqueous acetone and fractions**, not water-only leaf identification. | Quercetin neutral and 4′-O monoanion tested; only the monoanion/Asp542H model met the geometric screen in both boxes. Other references untested. |
| *Vitex negundo* | Agnuside has leaf occurrence evidence, but full water-only chromatographic details and source-specific stereochemistry remain unresolved; methanol-first aqueous partition is a different extraction method. | Neutral PubChem agnuside reference tested; both receptor conditions **INCONCLUSIVE**. |

These are associations in **external publications**. A PubChem CID, prepared reference structure or docked pose does not identify a molecule in the original laboratory extract. The [literature reviews](literature/five_plant_compound_evidence_review.md) and [evidence corrections](docs/compound_evidence_corrections.md) document the source-specific boundaries.

## Completed molecular work

### Ligands and receptors

All ten assigned reference ligands have checked 3D SDF structures under [`structures/ligands/prepared/`](structures/ligands/prepared/) and preliminary Meeko PDBQT files under [`structures/ligands/vina_prepared/`](structures/ligands/vina_prepared/). Preparation scripts and per-ligand metadata record source identifiers, stereochemistry, conformers, atom maps and hashes. Some exact botanical isomers or bound protonation states remain conditional even when a reference structure passes file-quality checks. Additional neutral/ionized hypotheses used in pilots have **separately named** files; original references remain preserved.

Official RCSB structures and working receptor files are under [`structures/proteins/`](structures/proteins/). The controlled pose work uses human N-terminal maltase-glucoamylase **3L4W** with native miglitol (MIG), in two pH-6.5 **dry/no-NAG** receptor hypotheses: Asp443−/Asp542−/His600-HID and Asp443−/Asp542H(OD2)/His600-HID. Human pancreatic α-amylase **4GQR** with native myricetin (MYC) was also assessed. Earlier 1XD0 and 2QMJ preparation records are retained; they are not presented as validated phytochemical targets. [Receptor preparation](docs/protocols/vina_4GQR_3L4W_preparation.md), [3L4W chemistry review](docs/protocols/vina_3L4W_chemistry_review.md) and [target summary](docs/protocols/vina_target_validation_summary.md) explain the structural choices.

### Native-ligand validation

The 3L4W **native-MIG positive control** achieved rank-1, symmetry-aware crystal-frame heavy-atom RMSD ≤2.0 Å in **5/5 seeds per receptor** in the original box and **10/10 total** in the enlarged 30 × 24 × 24 Å box. The enlarged-box range was **1.8671–1.9362 Å**. This is a limited **pose-recovery screen**; more crystal-like secondary MIG poses sometimes ranked lower. It does not establish physical binding energetics, water or glycan treatment, or biochemical potency. See the [large-box validation](docs/protocols/vina_3L4W_large_box_mig_validation.md) and [run CSV](data/metadata/vina_3L4W_large_box_mig_results.csv).

The current 4GQR/MYC Vina model **failed** the prespecified rank-1 ≤2.0 Å native-pose screen: **0/5** at exhaustiveness 32, **0/5** at 64, and **0/5** after the one-variable Glu233H sensitivity test. This model is **suspended for phytochemical docking**. A lower-ranked native-like pose in one run does not reverse that result. See the [failure investigation](docs/protocols/vina_4GQR_failure_investigation.md) and [Glu233H report](docs/protocols/vina_4GQR_glu233_sensitivity_report.md).

### Exploratory 3L4W phytochemical poses

The frozen pilot rule evaluates MIG-pocket localization, severe clashes and **unfitted rigid-core clustering** in at least **4/5 rank-1 seeds per ligand/receptor condition**. A geometric **PASS** indicates a repeatable *computed placement in that model only*. It is not a binding or inhibition result. Neutral and ionized ligand hypotheses are reported separately.

| Ligand hypothesis | Asp542−: original → enlarged box | Asp542H(OD2): original → enlarged box |
| --- | --- | --- |
| Mangiferin neutral | PASS 5/5 → PASS 5/5 | PASS 5/5 → INCONCLUSIVE 3/5 |
| Mangiferin 6-O monoanion | INCONCLUSIVE 0/5 → PASS 4/5 | PASS 5/5 → PASS 5/5 |
| Quercetin neutral | INCONCLUSIVE 0/5 → INCONCLUSIVE 0/5 | INCONCLUSIVE 0/5 → INCONCLUSIVE 0/5 |
| Quercetin 4′-O monoanion | INCONCLUSIVE 0/5 → INCONCLUSIVE 0/5 | PASS 5/5 → PASS 5/5 |

The original box was **24 × 18 × 18 Å**; its face lay within 1 Å of **36/40 rank-1 poses**. A single predefined enlarged-box comparison reduced this to **1/40**, but changed **two of eight** condition classifications. Both boxes had 4/8 PASS calls in total. The larger box passed its native-MIG control and is qualified for **limited pose-only work**, while candidate geometry remains **box and receptor-state sensitive**. The original 40-run pilot and enlarged comparison each have two complete pose outputs with **unconfirmed process exits**; their provenance is preserved rather than relabelled as confirmed successes. Read the [integrated assessment](docs/protocols/vina_3L4W_integrated_pilot_assessment.md), [box audit](docs/protocols/vina_3L4W_box_adequacy_report.md) and corresponding [original](data/metadata/vina_3L4W_phytochemical_pilot_results.csv) and [enlarged](data/metadata/vina_3L4W_box_adequacy_results.csv) run tables.

**Gallic acid:** neutral acid and carboxylate each yielded **INCONCLUSIVE 0/5** on both enlarged-box receptors. All 20 rank-1 poses favored a Z− box-face site outside the frozen MIG pocket; every run also sampled a lower-ranked pocket pose. The off-pocket site overlaps omitted deposited waters. See the [pilot report](docs/protocols/vina_3L4W_gallic_acid_pilot_report.md) and [structural investigation](docs/protocols/vina_3L4W_gallic_acid_off_pocket_investigation.md).

**Agnuside:** neutral PubChem CID 442416 yielded **INCONCLUSIVE 0/5** on both receptors. Thirteen lower-ranked modes occupied the frozen MIG pocket; all ten rank-1 poses and 12/13 pocket alternatives overlapped omitted deposited water oxygen positions. Five starting conformers were paired one-to-one with five seeds, so conformer and seed effects cannot be separated. See the [pilot report](docs/protocols/vina_3L4W_agnuside_pilot_report.md) and [structural investigation](docs/protocols/vina_3L4W_agnuside_structural_investigation.md).

## Reproducing and inspecting the record

The project virtual environment was used for the completed work. [`requirements.txt`](requirements.txt) currently pins **PyYAML only**; it is **not a complete docking-environment lockfile**. The archived reports record Vina 1.2.7, Meeko 0.8.0, RDKit/Gemmi usage, settings, seeds, commands and SHA256 values. Vina is a separate executable; CCDC GOLD and BIOVIA Discovery Studio were considered in planning but **no GOLD production docking, antioxidant quantum-chemical analysis, ADMET or molecular dynamics** is reported here. To inspect the YAML and CSV record without rerunning docking:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -c "import yaml; print(yaml.safe_load(open('config/formulations.yaml', encoding='utf-8'))['formulations'])"
```

This installs only the declared dependency. Re-executing molecular preparation or Vina scripts requires separately documented compatible RDKit, Meeko, Gemmi, NumPy and Vina installations; some scripts refer to the original workstation's `E:\docking tools\vina_1.2.7_win.exe` and need local-path configuration. **Do not run a docking script merely to inspect archived results.**

| Research artifact | Location |
| --- | --- |
| Literature evidence and selection | [`literature/`](literature/), [`data/compounds/`](data/compounds/), [five-plant evidence table](data/metadata/project4_five_plant_evidence_summary.csv) |
| Ligand/receptor preparation, structures and metadata | [`structures/`](structures/), [`data/metadata/`](data/metadata/), [`scripts/`](scripts/) |
| Frozen docking configurations, starts, output poses, process records and per-pose tables | [`docking/validation/vina/`](docking/validation/vina/), [`docking/pilot/vina/`](docking/pilot/vina/) |
| Scientific interpretation and reproducibility reports | [`docs/protocols/`](docs/protocols/), [five-plant synthesis](docs/protocols/project4_five_plant_scientific_synthesis.md), [manuscript-readiness assessment](docs/protocols/project4_manuscript_readiness.md) |

Run-level CSVs record input/output paths and hashes; experiment manifests and saved starts link each run to its exact chemistry and settings. Some Vina `.log` files are currently ignored by `.gitignore`; the [Git checkpoint audit](docs/protocols/project4_git_checkpoint_audit.md) identifies the files needed to back up the full provenance. The README describes files present in this working repository; GitHub completeness must be checked against that audit.

## Scientific and publication limits

- The missing original dissertation and raw in-vitro observations prevent independent numerical reproduction, statistical interpretation and verification of reported experimental inconsistencies.
- The candidate library is drawn from external literature with unequal plant-part, solvent and stereochemical evidence. It does **not** profile the original aqueous leaf extracts or establish constituent concentrations in C1–C3.
- Both tested 3L4W receptors omit six linked NAG units and crystallographic waters and use provisional Asp542/His600 protonation. Static deposited water positions reveal geometric conflicts but do not establish dynamic solvent behavior. A successful native-MIG pose screen cannot validate these omitted components or scoring physics.
- Phenolic acid/phenolate and glycoside reference states may differ from solution or bound microstates. Five starts do not establish complete sampling of agnuside's flexible glycoside conformations.
- Vina scores and contact distances cannot establish biochemical inhibition, relative potency, antioxidant mechanism, antidiabetic efficacy or formulation synergy. DPPH/ABTS are solution assays and are not reproduced by protein docking.

The [scientific synthesis](docs/protocols/project4_five_plant_scientific_synthesis.md) states the defensible findings and unresolved gates for a manuscript. Future work should recover the original experimental records, update the software environment specification, and independently validate any new target/solvent model before broader docking.

## Attribution and source record

The original in-vitro work was conducted by **Tharuniya Anpalakan**, Department of Applied Sciences, **Northumbria University**, under **Mr. Visal Chamikara** (supervisor) and **Ms. Upamitha Gunawardane** (co-supervisor). The original report is identified in project documentation as *Formulation and in vitro Evaluation of Polyherbal Combinations of Aegle marmelos, Boerhavia diffusa, Salacia reticulata, Senna auriculata, and Vitex negundo: A Study of Antioxidant and Antidiabetic Properties*, dated **11 August 2026**; the report itself and its bibliography are not deposited. The [configuration](config/formulations.yaml) transcribes its formulation and assay context. Verified compound-specific literature references and their limitations are linked from the [library](data/compounds/final_preliminary_ligand_library.csv) and [evidence reviews](literature/five_plant_compound_evidence_review.md); no citation or experimental value has been inferred to fill the missing dissertation record.

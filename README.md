# Integrated In-Vitro and In-Silico Investigation of Antioxidant and Antidiabetic-Related Properties of Polyherbal Formulations Containing *Aegle marmelos*, *Boerhavia diffusa*, *Salacia reticulata*, *Senna auriculata*, and *Vitex negundo*

| Project detail | Description |
| --- | --- |
| Researcher | Tharuniya Anpalakan |
| Institution | Northumbria University, Department of Applied Sciences |
| Original in-vitro research report | 11 August 2026 |
| Project status | In-vitro findings have been reported; the in-silico extension is planned. |

## Research overview

This project brings together an existing experimental evaluation of aqueous leaf extracts from five medicinal plants and a proposed computational investigation of their antioxidant and antidiabetic-related constituents. The original study compared three polyherbal formulations using phytochemical, radical-scavenging, and glucose-related assays. The computational extension is intended to examine literature-supported plant compounds, their interactions with justified protein targets, antioxidant-related molecular properties, and predicted pharmacokinetic and safety profiles.

The in-vitro work is the source of the qualitative findings described below. A formulation-level assay response does not identify which constituent caused it, and a future docking or ADMET prediction would provide a hypothesis rather than proof of clinical efficacy.

## Objectives

1. Document the five selected plants, the three extract ratios, and the assays used in the original study.
2. Preserve the reported experimental comparisons while identifying data that still need verification against the original report and raw records.
3. Select candidate phytochemicals from traceable literature and structure records, then prepare reproducible ligand and target datasets.
4. Evaluate plausible antidiabetic-related target interactions by molecular docking and explore antioxidant-related molecular properties with an explicitly defined computational method.
5. Assess predicted absorption, distribution, metabolism, excretion, and toxicity (ADMET), and interpret computational results alongside the experimental evidence and its limitations.

## Plants and formulations

The plant codes and ratios below are transcribed from [`config/formulations.yaml`](config/formulations.yaml). They describe **proportions of aqueous leaf extracts**, not concentrations of purified phytochemicals. The config does not specify a mass, volume, or concentration unit for each ratio part.

| Code | Medicinal plant |
| --- | --- |
| AM | *Aegle marmelos* |
| BD | *Boerhavia diffusa* |
| SA | *Senna auriculata* |
| SR | *Salacia reticulata* |
| VN | *Vitex negundo* |

| Formulation | AM | BD | SA | SR | VN | Basis recorded in config |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| Combination 1 | 2.0 | 1.5 | 2.0 | 2.0 | 2.5 | Equal weighting of total phenolic content (TPC), total flavonoid content (TFC), and total antioxidant capacity (TAC) |
| Combination 2 | 2.0 | 1.0 | 2.5 | 2.0 | 2.5 | Greater emphasis on TAC; TPC/TFC/TAC weights of 0.25/0.25/0.50 |
| Combination 3 | — | — | 3.5 | 3.5 | 3.0 | Selected plants with higher antioxidant-related measurements |

An em dash means the plant is absent from that formulation. Each listed ratio totals 10 parts.

## Original in-vitro study

### Aqueous leaf extraction and formulation

The original research used leaf material from the five plants to prepare **aqueous leaf extracts**. The extracts were then combined in the three proportions above for the comparative assays. The materials currently in this repository do not provide the extraction temperature, duration, leaf-to-water ratio, drying method, yield, or complete preparation protocol. These parameters must be recovered from the original laboratory record or report before the extraction can be reproduced precisely.

### Assays

| Assay | Experimental purpose |
| --- | --- |
| Qualitative phytochemical screening | Screen extracts for classes of phytochemicals |
| TPC, TFC, and TAC | Characterize phenolic content, flavonoid content, and total antioxidant capacity; inform formulation selection |
| DPPH and ABTS radical-scavenging assays | Compare antioxidant-related responses under two assay conditions |
| Glucose adsorption | Assess glucose binding/adsorption under the reported in-vitro conditions |
| *Saccharomyces cerevisiae* (yeast) glucose uptake | Assess glucose uptake in a yeast-based assay |

### Reported findings and data quality

| Comparison | Finding reported for the original experiment |
| --- | --- |
| DPPH | Combination 1 performed strongest in the reported comparison. |
| Glucose adsorption | Combination 2 generally showed the highest response. |
| ABTS | Combination 3 performed strongest in the reported comparison. |
| Yeast glucose uptake | Combination 3 generally showed high uptake. |

These are qualitative summaries of the existing experimental research, not values recalculated from deposited data. The original report's underlying tables and raw measurements are not present here, so numerical values, uncertainty, replicate counts, units, and statistical significance cannot be checked. Reported experimental-data inconsistencies remain unresolved; any conflicting values, labels, or trends in the source material should be reconciled against laboratory records before quantitative comparison, pooling, or publication. No missing values have been inferred to resolve those inconsistencies.

## Planned in-silico investigation

The computational work has **not yet been performed**. The proposed workflow is:

1. **Compound selection:** identify plant-associated candidate compounds through citable literature, record the plant-to-compound evidence and chemical identifiers, and avoid treating an extract ratio as a compound concentration.
2. **Ligand acquisition and preparation:** obtain chemical structures and identifiers from [PubChem](https://pubchem.ncbi.nlm.nih.gov/); inspect and prepare three-dimensional ligand structures using Avogadro. Record protonation/tautomer decisions, stereochemistry, and preparation settings.
3. **Target selection and preparation:** choose biologically justified antidiabetic-related and, where suitable, antioxidant-related protein targets; obtain experimental structures from the [RCSB Protein Data Bank](https://www.rcsb.org/). Document PDB IDs, chain selection, missing residues, cofactors, waters, bound ligands, and binding-site definitions. No target or PDB entry has been selected in this repository.
4. **Docking:** use CCDC GOLD to dock prepared ligands into defined binding sites. Record software version, scoring function, search settings, reference-ligand validation where possible, poses, and scores. Inspect contacts and figures in BIOVIA Discovery Studio Visualizer. Docking scores must be interpreted within a consistent protocol, not as direct measures of activity.
5. **Antioxidant-related computational analysis:** define and justify molecular descriptors or an appropriate mechanistic calculation relevant to antioxidant behavior, including the chemical state and solvent assumptions where applicable. Compare these results cautiously with DPPH, ABTS, and TAC observations; the exact calculation protocol remains to be selected.
6. **ADMET prediction:** evaluate candidate compounds with a documented prediction tool and version for properties such as absorption, distribution, metabolism, excretion, and toxicity. Report model applicability, confidence or uncertainty when available, and any conflicting predictions. No ADMET predictions are available yet.
7. **Analysis and reporting:** use Python to organize compound provenance, formulation metadata, docking outputs, antioxidant-related calculations, and ADMET predictions. Keep raw inputs, settings, and processed results linked so analyses can be repeated.

GOLD and Discovery Studio are separate licensed applications; they are not installed by the Python requirements file. PubChem and RCSB PDB are source databases, while Avogadro supports molecular preparation and inspection. None of these tools implies that a docking or ADMET result already exists.

## Software and Python environment

The repository currently lists only `PyYAML==6.0.3` in [`requirements.txt`](requirements.txt). It is sufficient to read the YAML configuration with Python. No analysis scripts, notebooks, docking inputs, or computational output files have been deposited, and no Python version is pinned. A future computational workflow will need its own documented software versions, dependencies, and run commands.

### Installation

Clone or download the repository, then create a Python virtual environment in its root directory. For example, in PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

On macOS or Linux, activate the environment with `source .venv/bin/activate` before installing requirements. To inspect the currently reproducible formulation configuration:

```powershell
python -c "import yaml; from pathlib import Path; data = yaml.safe_load(Path('config/formulations.yaml').read_text(encoding='utf-8')); print(data['formulations'])"
```

This command reads and prints the recorded ratios; it does not run an assay or computational analysis.

## Repository structure

```text
.
├── README.md                 Project documentation
├── requirements.txt          Python dependency (PyYAML)
├── config/
│   └── formulations.yaml     Plants, extract ratios, and assay inventory
├── admet/                    Reserved for future ADMET work
├── data/                     Reserved for research data
├── docking/                  Reserved for docking inputs and outputs
├── docs/                     Reserved for supporting documentation
├── literature/               Reserved for literature records
├── manuscript/               Reserved for manuscript materials
├── notebooks/                Reserved for analysis notebooks
├── results/                  Reserved for result tables and figures
├── scripts/                  Reserved for analysis scripts
├── structures/               Reserved for ligand and protein structures
└── tests/                    Reserved for tests
```

At present, the named research directories contain no deposited files; empty directories may not appear in a Git clone. The tracked project content is the README, formulation YAML, requirements file, and `.gitignore`.

## Reproducibility and results status

| Component | Current status |
| --- | --- |
| Plant identities, formulation ratios, assay inventory | Recorded in `config/formulations.yaml`; readable with the installation steps above |
| Original experimental comparisons | Qualitatively reported; raw data, detailed methods, and source tables are not deposited |
| Ligand and target selection | Planned; no curated compound list or PDB structures deposited |
| Molecular docking and antioxidant-related calculations | Planned; no protocols, scripts, poses, scores, or calculation results deposited |
| ADMET analysis | Planned; no predictions deposited |

Consequently, the repository currently supports reproduction of the **configuration record**, not independent reproduction of the laboratory experiments or execution of the proposed computational workflow.

## Scientific limitations

- Aqueous extract assays measure mixtures; the contribution and concentration of any individual compound are unknown.
- The available repository material lacks complete extraction and assay protocols, raw observations, replicate information, uncertainty estimates, and resolved experimental inconsistencies.
- DPPH, ABTS, glucose adsorption, and yeast uptake are screening assays. They do not by themselves establish mechanism, efficacy in humans, or clinical safety.
- Plant-to-compound assignments from literature require source-specific verification; a compound reported in a plant is not automatically confirmed in these aqueous leaf extracts.
- Docking is sensitive to target structure, ligand state, binding-site setup, and scoring choices. Antioxidant-related calculations and ADMET predictions also depend on their models and assumptions. Any future predictions will need experimental follow-up.

## Source and references

1. Anpalakan, Tharuniya. *Formulation and in vitro Evaluation of Polyherbal Combinations of Aegle marmelos, Boerhavia diffusa, Salacia reticulata, Senna auriculata, and Vitex negundo: A Study of Antioxidant and Antidiabetic Properties.* Final-year Biomedical Science research report, Northumbria University, Department of Applied Sciences, 11 August 2026. The report is identified in the existing project material but is not included in this repository; its bibliography cannot be verified here.
2. [`config/formulations.yaml`](config/formulations.yaml). Project transcription of formulation ratios, plant codes, and assay categories, attributed in its header to Section 3.0 and Appendix 3 of the original report.

The original report's individual literature references are not available in this repository, so they are not reproduced or supplemented with unverified citations here. Future compound and target choices should carry their own traceable references.

## Authorship and supervision

The original in-vitro research was conducted by **Tharuniya Anpalakan** at the **Department of Applied Sciences, Northumbria University**. **Mr. Visal Chamikara** supervised the in-vitro research, with **Ms. Upamitha Gunawardane** as co-supervisor. The report date above refers to that original experimental work; it is not a completion date for the planned computational extension.

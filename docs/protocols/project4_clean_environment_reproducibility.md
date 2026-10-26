# Project 4 clean-environment reproducibility audit

## Audit scope and status

This Step 86 audit was performed on 26 October 2026 without launching AutoDock Vina, rescoring poses, changing thresholds, or modifying archived scientific outputs. The source working tree was on `main` at commit `82532653d17273d18265eadf555d766ed1b4008b`. The configured origin was `https://github.com/tharuniya26/Polyherbal_Antioxidant_Antidiabetic_InSilico.git`.

The analysis-only reproduction of the **current working-tree synthesis** is **PASS** in an isolated Windows environment after applying the accompanying `.gitattributes` policy. A clone of the currently published commit alone is **BLOCKED** from reproducing the final synthesis because the Step 85 generator, CSV and synthesis report, and the Step 86 environment/test files are not present in commit `8253265`. They must be committed before a GitHub clean-clone PASS can be claimed.

## Tested platform and environment

- Operating system: Windows 10 64-bit, build 19045.
- Python: CPython 3.12.10, 64-bit.
- Package installer: pip 25.0.1.
- Tested package specification: `requirements-project4-analysis.txt` constrained by `constraints-project4-py312-windows.txt`.
- `pip check`: PASS (`No broken requirements found`).
- Scientific imports: PASS for Gemmi, Meeko, NumPy, OpenMM, PDBFixer, PROPKA, PyYAML, RDKit and SciPy.

The exact tested package set was:

| Package | Version |
|---|---:|
| colorama | 0.4.6 |
| Gemmi | 0.7.3 |
| Meeko | 0.8.0 |
| NumPy | 2.5.4 |
| OpenMM | 8.6.1 |
| PDBFixer | 1.12.0 |
| Pillow | 12.3.0 |
| PROPKA | 3.5.1 |
| PyYAML | 6.0.3 |
| RDKit | 2026.3.6 |
| SciPy | 1.18.1 |
| tqdm | 4.70.1 |

The final synthesis generator uses only the Python standard library. The scientific package set is retained because the archived preparation, input-validation and pose-analysis scripts import these packages.

## Installation

From a checkout containing the Step 86 files:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-project4-analysis.txt
.\.venv\Scripts\python.exe -m pip check
```

The constraints file is a tested Windows/CPython 3.12 environment record, not a claim that every version combination has been tested on Linux or macOS. A different Python minor version may not have matching wheels for RDKit, OpenMM or Gemmi.

## External AutoDock Vina dependency

The frozen runs record AutoDock Vina 1.2.7 with executable SHA256:

```text
e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5
```

The audited workstation executable was external to the repository at `E:\docking tools\vina_1.2.7_win.exe`; `vina_1.2.7_win.exe --version` returned `AutoDock Vina v1.2.7`, and its recalculated SHA256 matched the recorded value. The executable is not bundled. Obtain version 1.2.7 from the official AutoDock Vina release distribution and verify the binary before any separately authorized docking reproduction. Vina is **not required** to regenerate the final CSV or run the Step 86 tests.

## Required analysis inputs

The generator reads the following tracked result, all-pose and experiment-manifest triplets:

1. `data/metadata/vina_3L4W_gallic_acid_pilot_results.csv`, `docking/pilot/vina/3L4W_gallic_acid_pose_only/all_pose_metrics.csv`, and its `experiment_manifest.json`.
2. `data/metadata/vina_3L4W_alternative_ligand_pilot_results.csv`, `docking/pilot/vina/3L4W_alternative_ligands_pose_only/all_pose_metrics.csv`, and its `experiment_manifest.json`.
3. `data/metadata/vina_3L4W_phytochemical_pilot_results.csv`, `docking/pilot/vina/3L4W_phytochemical_pose_only/all_pose_metrics.csv`, and its `experiment_manifest.json`.
4. `data/metadata/vina_3L4W_box_adequacy_results.csv`, `docking/pilot/vina/3L4W_box_adequacy_sensitivity/all_pose_metrics.csv`, and its `experiment_manifest.json`.
5. `data/metadata/vina_3L4W_agnuside_pilot_results.csv`, `docking/pilot/vina/3L4W_agnuside_pose_only/all_pose_metrics.csv`, and its `experiment_manifest.json`.

It also reads `data/metadata/vina_3L4W_gallic_acid_off_pocket_investigation.csv`, `data/metadata/vina_target_validation_summary.csv`, botanical evidence metadata, and the archived receptor, ligand, start, output and attempt-record paths cited by the result tables. These primary inputs and archived outputs are tracked at the audited baseline. The local `.venv/` and runner diagnostic directory are ignored and are not scientific inputs. AutoDock Vina is the only required external executable for rerunning docking; it is not an analysis-only dependency.

## Analysis-only reproduction commands

These commands do not invoke a docking runner:

```powershell
.\.venv\Scripts\python.exe scripts\build_project4_final_five_plant_summary.py
Get-FileHash -Algorithm SHA256 data\metadata\project4_final_five_plant_docking_summary.csv
.\.venv\Scripts\python.exe -m unittest discover -s tests -p "test_*.py" -v
```

Expected output:

- CSV: `data/metadata/project4_final_five_plant_docking_summary.csv`
- Rows: 26 unique molecular-condition rows.
- SHA256: `5e3168dd3e5ed289e53eb8b48922a4c809a942f6b80a8dd7f294f184b20da3a1`
- Test result: 19 tests passed, 0 failed (nine Step 86 integrity tests and ten non-scientific Step 79 runner-hardening tests).

## Independent isolated result

The test directory was `C:\Users\strangers\AppData\Local\Temp\project4_step86_20261026`. A fresh Python 3.12 environment was created there and packages were installed from the new specification. A local-only simulated post-checkpoint clone was used because the required Step 85/86 files are not yet committed.

The generator was run twice. Both files were byte-identical and both hashes were:

```text
5e3168dd3e5ed289e53eb8b48922a4c809a942f6b80a8dd7f294f184b20da3a1
```

The tests independently confirmed:

- 130 candidate searches and 2,559 valid candidate poses;
- 126 zero-exit searches and four complete-output searches with unconfirmed exits;
- all referenced candidate receptor, ligand, source, start, output and attempt-record hashes;
- complete `MODEL`/`ENDMDL` counts for the candidate outputs;
- 40 native-control runs represented by eight control conditions;
- unchanged frozen classifications and the four-of-five geometric rule;
- gallic acid, rutin and agnuside remain INCONCLUSIVE;
- luteolin remains PASS 5/5 in both tested receptor states;
- mangiferin remains box-, ligand-state- and receptor-state-dependent;
- quercetin remains microstate- and receptor-state-dependent;
- the shared quercetin molecular experiment is represented once, with two botanical associations rather than two experiments;
- no unsupported binding, inhibition, synergy, antidiabetic-effect or clinical-efficacy conclusion is introduced.

## Line-ending provenance and platform limitation

The first ordinary Windows clone converted provenance-critical LF files through the global `core.autocrlf=true` setting. It preserved the numerical records and classifications, but changed archived input/manifest byte hashes and produced final CSV hash `3842d9d0eac598f6c11198a3242cf14ff71304b16d57b1803447c95f2a0b986f`. This was a genuine byte-reproducibility failure.

The added `.gitattributes` records the mixed historical line-ending conventions: result CSVs, prepared receptor PDBQT files, archived docked outputs and prepared source SDFs use CRLF; JSON records, prepared ligand PDBQT files, frozen start files and crystallographic reference SDFs preserve LF bytes. In a simulated clone made after that attribute file was committed, all archived hashes and the expected final CSV hash passed. The policy itself must be committed before it can protect a real GitHub clone. Linux and macOS checkout behavior was not independently tested.

## Analysis reproducibility versus docking reproducibility

This audit establishes **analysis reproducibility** from the archived files: deterministic aggregation, provenance checks, frozen-rule recomputation and scientific-claim boundaries. It does not repeat the docking searches and therefore does not independently establish **docking execution reproducibility**. A future docking reproduction would require the verified Vina 1.2.7 executable, the frozen inputs and parameters, comparable CPU/platform behavior, and explicit authorization to run the 170 searches. Small floating-point or runtime differences can occur across hardware or builds even when seeds and parameters are held fixed.

## Publication blocker

The current GitHub commit does not contain the Step 85 synthesis or Step 86 reproducibility additions. Until the intended files are reviewed, staged, committed and pushed, GitHub clean-clone reproducibility remains **BLOCKED**. Unrelated untracked root entries (`pace'` and the file whose displayed name begins `taged files`) were neither opened, altered, deleted nor included in the audit file set.

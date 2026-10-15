# Step 60 — native miglitol redocking in the enlarged 3L4W box

**Decision: supported for limited pose-only follow-up in both tested receptor states.** This decision applies to the prespecified 30 × 24 × 24 Å box and the two dry, no-NAG human ntMGAM receptor hypotheses. It is a native-ligand *pose-recovery* result, not validation of binding affinity, solution protonation, glycan omission, solvent treatment, or phytochemical inhibition.

## Why this test was run

The [Step 59 box audit](vina_3L4W_box_adequacy_report.md) found that 36/40 rank-1 mangiferin/quercetin pilot poses approached a face of the original 24 × 18 × 18 Å search box within 1 Å. One predefined enlargement to 30 × 24 × 24 Å reduced this to 1/40 but changed two of eight geometric pilot classifications. Before interpreting those enlarged-box candidate poses, the deposited native inhibitor had to recover its pose in the *same* enlarged box for both previously tested receptor microstates. No new box or chemistry was selected from the results below.

## Frozen inputs and pre-run verification

The [reproducibility script](../../scripts/run_vina_3L4W_large_box_mig_validation.py) checked both receptor PDBQT files against the [frozen Step 57 manifest](../../data/metadata/vina_3L4W_pilot_input_manifest.csv):

| Dry, no-NAG receptor microstate | PDBQT SHA256 | Prior native-MIG E64 control |
| --- | --- | --- |
| Asp443− / Asp542− / His600-HID | `7a76e0fa3b7bf2197714a57a88c60d2f4eefedccb18636951f3f643473afadf3` | 5/5 rank-1 RMSD ≤2 Å in 24 × 18 × 18 Å box |
| Asp443− / Asp542H(OD2) / His600-HID | `220e3585558e4bc261e6e76ddb036738dea827543e1b135399506e95ea9f93fd` | 5/5 rank-1 RMSD ≤2 Å in 24 × 18 × 18 Å box |

The unchanged neutral-MIG crystallographic reference PDBQT has SHA256 `c7417f6bcb661e0f2217c6ba3619b6fd6c03ffb029fe5840d4b0d13c1bf16169`; the preserved deposited crystal SDF has SHA256 `ea673a59043bea05abbe7d66d8bd29c3786021f564c394027a8b0201385503f7`. The script also checked the original PDB, deposited atom map, SDF-to-PDBQT heavy-atom map, receptor source PDBs, Vina executable/version, and **all ten previous E64 neutral-MIG start/output/log hashes**. The SDF contains the audited neutral C8H17NO5 graph, 14 heavy atoms, and four stereodefined carbons. Its 14 heavy coordinates match the deposited MIG atoms to PDB precision. The previously published native-MIG parser found the chemically allowed graph automorphism set; the crystal pose and atom maps were not altered.

For each seed the script regenerated the audited NumPy PCG64 quaternion rotation and translation with the original starting-pose rule, required ≥5 Å unfitted crystal-frame start RMSD, and required **byte-for-byte identity** with the previous E64 randomized start. The same copied start was used for both receptor states; all starts also lie within the enlarged box. Vina's seeded global search was then run from these non-native inputs. The deposited reference SDF was used only for subsequent RMSD calculation.

## Prespecified search and analysis

- AutoDock Vina 1.2.7, `vina` scoring, exhaustiveness **64**, `num_modes 20`, `energy_range 5.0`, `min_rmsd 1.0`, one CPU per job; two jobs ran concurrently.
- Seeds **1101, 2202, 3303, 4404, 5505** for each receptor. Box center **(45.38664286, 92.35157143, 34.82200000) Å**, the deposited MIG heavy-atom centroid; size **30 × 24 × 24 Å**.
- Every returned PDBQT mode was reconstructed through Meeko/RDKit, checked against the original ligand graph, heavy-atom identity/map, stereochemistry, finite coordinates and bond lengths, then compared to the unchanged crystal pose using symmetry-aware **14-heavy-atom RMSD in the original receptor coordinate frame**. Docked poses were **not** independently fitted to the reference.
- The prospective screen was **rank-1 RMSD ≤2.0 Å in each of the five seeds per receptor**. Best RMSD among *all* returned modes was separately recorded. A Vina score was used only to describe ranking within its own run.

The [10-row results CSV](../../data/metadata/vina_3L4W_large_box_mig_results.csv) records commands, settings, process exits, runtimes, hashes, scores, rank-1 and best-of-all-pose RMSDs, old-box comparators and errors. The [raw output directory](../../docking/validation/vina/3L4W_large_box_MIG/) contains each randomized start, Vina output, log, all **200** per-pose metrics and the preregistered reproducibility manifest. All source structures and earlier validation/pilot results remain unchanged.

## Results and direct old-box comparison

**Ten of ten planned runs completed with confirmed Vina exit code zero; zero failed or timed out.** Vina returned 20 modes per run, and all 200 passed the specified Meeko/RDKit readback and geometry checks. All ten enlarged-box rank-1 poses passed the original ≤2.0 Å screen.

| Receptor | Seed | Old 24 × 18 × 18 rank-1 RMSD (Å) | New 30 × 24 × 24 rank-1 RMSD (Å) | Old best RMSD (Å), rank | New best RMSD (Å), rank | New rank-1 Vina score |
| --- | ---: | ---: | ---: | --- | --- | ---: |
| Asp542− | 1101 | 1.8775 | 1.8865 | 0.6997, 2 | 1.8865, 1 | −5.530 |
| Asp542− | 2202 | 1.8691 | 1.9338 | 0.6806, 2 | 1.9338, 1 | −5.511 |
| Asp542− | 3303 | 1.8629 | 1.8859 | 1.8629, 1 | **0.4658, 2** | −5.536 |
| Asp542− | 4404 | 1.8782 | 1.8755 | 0.7227, 2 | 1.8755, 1 | −5.528 |
| Asp542− | 5505 | 1.8849 | 1.9362 | 0.6926, 2 | 1.9362, 1 | −5.512 |
| Asp542H(OD2) | 1101 | 1.8814 | 1.8865 | 0.6989, 2 | 1.8865, 1 | −5.542 |
| Asp542H(OD2) | 2202 | 1.8916 | 1.8943 | 0.6813, 2 | 1.8943, 1 | −5.538 |
| Asp542H(OD2) | 3303 | 1.8813 | 1.8693 | 1.8813, 1 | 1.8693, 1 | −5.524 |
| Asp542H(OD2) | 4404 | 1.8634 | 1.8675 | 1.8634, 1 | 1.8675, 1 | −5.537 |
| Asp542H(OD2) | 5505 | 1.8632 | 1.8671 | 0.6816, 2 | **0.4897, 2** | −5.550 |

Asp542− rank-1 RMSDs span **1.8755–1.9362 Å** in the enlarged box, compared with **1.8629–1.8849 Å** in the original box. Asp542H(OD2) spans **1.8671–1.8943 Å**, compared with **1.8632–1.8916 Å** before. Thus both states maintain **5/5** rank-1 recovery, with modest seed-dependent changes in pose geometry. Rank-1 Vina scores changed by −0.007 to +0.024 score units relative to the old box; this is neither a measured free-energy difference nor a basis to select a receptor state.

**Secondary-pose sampling changed.** A more crystal-like rank-2 pose appeared in **7/10** original-box comparisons but only **2/10** enlarged-box runs. The latter two modes have 0.4658 and 0.4897 Å RMSD but score **0.144** and **0.141** Vina units worse than rank 1, respectively. The original box returned 17 poses at ≤2 Å among 196 modes for these two neutral-MIG receptor cases; the enlarged box returned **12 among 200**. This does not invalidate the predefined rank-1 screen, but it shows that the wider search space changed discovery or ranking of secondary modes. The data alone do not distinguish search coverage from scoring preference or receptor-chemistry effects. Increasing exhaustiveness or tuning the box after seeing these outputs was not part of this test.

## Scientific decision and next task

**Both Asp542 receptor models maintain reproducible native-MIG rank-1 pose recovery in the single preregistered enlarged box.** The **30 × 24 × 24 Å box is supported for a limited, explicitly pose-only follow-up** with these same dry/no-NAG receptor variants. It does not erase Step 59's **BOX-SENSITIVE** classification of the mangiferin/quercetin pilot: two of eight candidate geometric classifications changed with box size, and several candidate clusters depend on receptor or ligand microstate. Passing MIG redocking establishes a positive local pose-control check for this box, not candidate binding or potency.

The receptors still omit six covalently linked NAG units and all active-site waters; Asp542− versus Asp542H(OD2), His600-HID, and neutral MIG are computational microstate hypotheses. Ligand-facing Asp542H rotamers, alternative histidine tautomer, glycan-retaining preparation, and defensible water treatment remain unresolved. Fixed-water Vina and raw Vina scores cannot validate solvent dynamics, enzyme inhibition, binding affinity, formulation synergy or clinical effects. No new phytochemical docking, ADMET, molecular dynamics or potency ranking was performed.

**Exact next task:** perform a read-only integrated review of the existing Step 58 and Step 59 mangiferin/quercetin pose-only outputs under the now checked 30 × 24 × 24 Å box, applying the *already frozen* pocket, clash, rigid-core and 4/5 convergence rules to each ligand/receptor microstate, with special attention to the split mangiferin-neutral/Asp542H cluster and the off-pocket mangiferin-monoanion/Asp542− seed. Decide whether any limited ligand-specific follow-up is needed before considering another chemical. Expansion to all five medicinal plants is not yet justified by this native-MIG control alone.

Reanalyse immutable outputs without repeating docking:

```powershell
.\.venv\Scripts\python.exe scripts\run_vina_3L4W_large_box_mig_validation.py
```

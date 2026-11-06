# Project 4 final validation summary

**Validation date:** 7 November 2026
**Repository baseline:** `cbfe0f0039d4a89f915dc46db41c4aae48a8f8ac` (`main`, matching `origin/main` before local release-candidate changes)
**Publication state:** author-approved public release; published through a normal commit after the validation recorded here

## Final disposition

The repository is ready for scientific supervisor review and author approval as a carefully limited computational pose-geometry study. It is not evidence of compound occurrence in the original aqueous extracts, measured binding affinity, enzyme inhibition, antioxidant activity, antidiabetic efficacy, clinical benefit or polyherbal synergy. The original dissertation, raw assay observations, replicate-level statistics and direct chemical profiling remain unavailable.

No new docking, rescoring, molecular dynamics, ADMET or wet-laboratory analysis was performed during finalization. Existing archived calculations were inspected, deterministically synthesized, tested, interpreted and documented without changing their frozen decision rule. Publication uses a normal commit on `main`; Git history is not rewritten.

## Verified computational record

- The candidate archive contains **130 analyzable searches and 2,559 returned poses**.
- **126 searches** have recorded process exit zero. Four mangiferin searches have complete, chemistry-checked outputs with unconfirmed historical process exits; they were not relabelled as confirmed exit-zero runs.
- The deterministic final table contains **26 unique ligand-state/receptor-state/box conditions: 10 PASS and 16 INCONCLUSIVE**.
- The summary builder recalculated the frozen rank-1 rule from archived fields and reported **zero classification mismatches**.
- The summary SHA256 is `5e3168dd3e5ed289e53eb8b48922a4c809a942f6b80a8dd7f294f184b20da3a1`.

### Native controls and target scope

- The registered native-control record contains 40 searches: 25 for 3L4W and 15 for 4GQR. Ten additional 3L4W MIG searches assess the separate enlarged-box condition, for 50 archived native-control calculations in total.
- Human ntMGAM 3L4W recovered rank-1 native MIG within the prospective 2.0 Å criterion in every tested condition, including 10/10 enlarged-box runs. This validates only limited crystal-pose recovery for the tested model.
- Human pancreatic α-amylase 4GQR failed rank-1 native MYC recovery at exhaustiveness 32, exhaustiveness 64 and the Glu233H sensitivity condition (0/5 each). It remains suspended for candidate screening.

### Candidate outcomes

- Neutral luteolin produced geometric PASS results in both 3L4W receptor states (5/5 each).
- Quercetin was conditional: only the modelled 4′-O monoanion/Asp542H(OD2) combination passed, in both boxes. Neutral quercetin and Asp542− monoanion conditions were inconclusive.
- Mangiferin produced six PASS and two INCONCLUSIVE conditions and remained sensitive to ligand state, receptor state and box size.
- Gallic acid, rutin and agnuside remained INCONCLUSIVE under the frozen rank-1 rule. Lower-ranked pocket poses did not override the prespecified classification.
- Quercetin is one molecular experiment shared by the *Boerhavia diffusa* and *Senna auriculata* botanical associations; it was not double-counted.

PASS means only that at least four of five prespecified rank-1 poses met the frozen geometric criteria in one defined model condition. It is not a biochemical or therapeutic result. INCONCLUSIVE is not evidence of inactivity.

## Manuscript and figure package

- The revised Markdown manuscript contains **5,613 words before References** and 27 DOI-registered references.
- Six figures were regenerated as both 300-dpi PNG and editable SVG files from archived repository data.
- `figure_provenance.csv` records source paths, SHA256 values and selection rules for all six figures.
- The final PDF is 13 A4 pages with a Unicode serif font, numbered pages, figures, tables, live DOI links and the footer **“Supervisor-review draft • Not peer-reviewed.”**
- PDF path: `manuscript/Tharuniya_Anpalakan_Polyherbal_Computational_Manuscript.pdf`
- PDF size: 2,211,784 bytes
- PDF SHA256: `d32448fa57e2042208da345d10e554d92a26d6029c4e6c6fa9add577298c9f64`
- Markdown SHA256: `c6186b738b382825bf9e90935c57d8247046edc3bbe0d5aca0ac45547c8eb29b`

Every rendered PDF page was visually inspected after the final layout revision. No clipped text, overlapping captions, broken tables or footer collisions were observed. The first build's footer collision and avoidable final-page whitespace were corrected before these hashes were recorded.

## Validation executed

| Check | Result |
|---|---|
| Deterministic final-summary regeneration | PASS — 26 verified rows |
| Manuscript figure regeneration | PASS — six PNG and six SVG files plus provenance |
| Complete repository test suite | PASS — 29/29 in CPython 3.12.10 with cache disabled |
| Local Markdown links | PASS — README, manuscript and figure-caption targets resolve |
| PDF structure | PASS — 13 A4 pages, expected text and figure captions extractable |
| PDF visual inspection | PASS — all 13 pages reviewed |
| Git whitespace check | PASS — no patch errors; Git only reported expected Windows line-ending notices |
| Sensitive-value scan | PASS — no credentials, passwords, API keys or private keys detected |
| Unrelated tracked scientific files | PASS — no tracked docking, structure or archived-result diff |

The documented CPython 3.12.10 `.venv` contains the archived analysis stack and ReportLab 5.0.1 installed through `requirements-project4-manuscript.txt` and its tested constraints. The final publication PDF was rebuilt in that environment.

## README and publication preparation

The README now reports the final six-compound screen, including rutin and luteolin; the 130-search/2,559-pose archive; the 10 PASS/16 INCONCLUSIVE synthesis; the failed 4GQR gate; manuscript scope; figure provenance; PDF link; rebuild commands; and unresolved scientific limitations.

Two verified accidental top-level terminal captures were removed: a 17,182-byte `less` help dump and a 1,082,074-byte coloured Git whitespace dump. No research input, output, script, structure or documentation file was deleted.

Files included in the author-approved publication commit are:

- `.gitattributes`
- `README.md`
- `requirements-project4-manuscript.txt`
- `constraints-project4-manuscript-py312-windows.txt`
- `manuscript/project4_insilico_manuscript_draft.md`
- `manuscript/figure_captions.md`
- `manuscript/reference_verification.csv`
- `manuscript/results_verification_checklist.md`
- `manuscript/Tharuniya_Anpalakan_Polyherbal_Computational_Manuscript.pdf`
- `results/figures/project4_manuscript/`
- `scripts/build_project4_manuscript_figures.py`
- `scripts/build_project4_manuscript_pdf.py`
- `tests/test_project4_manuscript_outputs.py`
- this validation summary

The publication manuscript uses the author-approved professional correspondence address `tharuniya08@gmail.com`. The preserved Step 89 draft retains the former institutional address locally and is deliberately excluded from the public commit, together with the internal Step 90 review and Step 91 response. No institutional correspondence address is included in the staged publication set.

## Remaining scientific limitations

1. The six named reference compounds were not analytically demonstrated in the original project's aqueous leaf extracts.
2. The original experimental measurements, replicate counts, uncertainties, statistics, detailed extraction method and dissertation bibliography are not deposited.
3. The 3L4W work uses rigid dry/no-NAG receptor hypotheses, omits crystallographic waters and linked glycans, and samples limited receptor and ligand microstates.
4. Native MIG recovery does not validate Vina affinity rankings across different phytochemical scaffolds.
5. Rutin and agnuside conformers are confounded with seed, and flexible glycoside sampling may be incomplete.
6. Four complete candidate outputs retain unconfirmed historical process exits.
7. The failed 4GQR model cannot support α-amylase candidate claims.
8. No biochemical, biophysical, cell, animal or clinical validation was performed.

The next defensible experimental priority is targeted LC-MS/HPLC profiling of the original extracts, followed by concentration-matched ntMGAM inhibition assays for analytically confirmed constituents.

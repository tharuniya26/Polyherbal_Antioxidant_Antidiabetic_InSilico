# Project 4 Step 79 — runner hardening validation

Date: 2026-10-20  
Frozen checkpoint: `a55389009774f8ebe466130454f1a44ed812f4fa`  
Scope: operational hardening and non-docking validation only

## Outcome

The Step 73 runner was completed as a fail-closed, write-once recovery runner. No production docking, Vina `score_only` calculation, completed-job rerun, or scientific-output modification was performed during Step 79.

The required read-only inventory passed:

- `REUSE: 8`
- `MISSING: 12`
- `BLOCKED: 0`

All 16 archived run artifacts (eight PDBQT outputs and eight attempt JSON records) are byte-identical to checkpoint `a553890`.

## Review of the partial Step 78 change

`step78_partial_changes.txt` records the partial work inherited by Step 79. It added imports and constants for UTC timestamps, faulthandler, atomic writes, thread-safe session logging, a fixed Vina SHA256, a 7200-second timeout, two workers, and strict attempt-name matching. It also made input verification capable of suppressing external checks.

Those changes were useful but incomplete. The active execution path still used `subprocess.run`, trusted an exit-code/output pair without validating the recorded output hash or model structure, had no `--dry-run` inventory, did not create launch/PID evidence, did not use the new atomic writer for attempts, and could proceed to aggregate analysis after incomplete work. The partial additions did not change a scientific input, search parameter, output, or analysis rule.

## Step 79 operational changes

The completed runner now:

- builds job descriptions without creating directories;
- verifies exact commands, frozen manifest/input hashes, exit code, recorded output SHA256, and balanced consecutive `MODEL`/`ENDMDL` records before reuse;
- treats an output without a valid attempt JSON as `BLOCKED`;
- treats a launch record without a final attempt record as `BLOCKED`;
- records a PID and exact command in a write-once launch sidecar immediately after process creation;
- writes final attempt JSON atomically without replacement;
- preserves legacy completed records that predate launch sidecars, but subjects them to all other integrity checks;
- records persistent UTC JSONL session events and a unique faulthandler diagnostic file;
- kills a timed-out child before writing `timeout_killed` status;
- captures worker exceptions and continues collecting the statuses of other submitted jobs;
- withholds aggregate analysis unless all 20 jobs are verified complete;
- uses at most two worker processes;
- checks the Windows process table for an already running `vina_1.2.7_win.exe` before any retry;
- supplies a read-only `--dry-run` path that does not invoke Vina, create output directories, or write control evidence.

The initial process-table implementation used `tasklist`, which returned `Access denied` in this environment. It was replaced with read-only PowerShell `Get-Process`; the subsequent required dry-run passed. The failed inventory attempt stopped before classification and wrote no repository artifact.

## Scientific parameter integrity

| Setting | Frozen baseline | Hardened runner | Status |
|---|---:|---:|---|
| AutoDock Vina | 1.2.7; executable SHA256 `e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5` | same | unchanged |
| Scoring | `vina` | `vina` | unchanged |
| Exhaustiveness | 64 | 64 | unchanged |
| `num_modes` | 20 | 20 | unchanged |
| `energy_range` | 5.0 | 5.0 | unchanged |
| `min_rmsd` | 1.0 Å | 1.0 Å | unchanged |
| CPU per search | 1 | 1 | unchanged |
| Maximum concurrent searches | 2 | 2 (`MAX_WORKERS`) | unchanged |
| Box center | (45.38664286, 92.35157143, 34.82200000) Å | same | unchanged |
| Box dimensions | 30 × 24 × 24 Å | same | unchanged |
| Seeds | 1101, 2202, 3303, 4404, 5505 | same | unchanged |
| Receptors | two frozen human 3L4W states | manifest paths and SHA256 values revalidated | unchanged |
| Ligands/starts | frozen qualified neutral LIG02/LIG08 conformers and starts | manifest paths and SHA256 values revalidated | unchanged |
| Native-MIG controls | frozen 10-run control table and native pose | hashes and 10/10 status criteria revalidated | unchanged |
| Pose analysis/classification | Step 73 implementation at `a553890` | analysis function unchanged | unchanged |

The dry-run validated the frozen input manifest against the saved 20/20 preflight record and validated all current receptor, qualified-ligand, randomized-start, MIG-reference, and MIG-control hashes. It did not invoke the validator subprocess or Vina.

## Isolated tests

Command:

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -p 'test_step79_runner_hardening.py' -v
```

Result: `Ran 9 tests ... OK`.

Covered cases:

1. successful mock subprocess and verified reuse;
2. finalized nonzero subprocess exit;
3. timeout, child kill, and finalized timeout evidence;
4. launch sidecar without a final attempt record;
5. rejection after output-hash alteration;
6. invalid/incomplete attempt JSON;
7. output without an attempt record and proof that no process starts or output changes;
8. atomic create-without-overwrite behavior and temporary-file cleanup;
9. correct `REUSE`/`MISSING`/`BLOCKED` counting.

Every generated mock PDBQT contains `REMARK MOCK NON-SCIENTIFIC` and exists only in a temporary test directory. The tests do not reference or execute the real Vina binary.

## Required real-project dry-run

Command:

```powershell
.\.venv\Scripts\python.exe scripts\run_vina_3L4W_alternative_ligand_pilot.py --dry-run
```

Result:

- `REUSE`: LIG02 seeds 1101, 2202, 3303, and 4404 for both `Asp542-` and `Asp542H(OD2)` (8).
- `MISSING`: LIG02 seed 5505 for both states (2), plus LIG08 seeds 1101, 2202, 3303, 4404, and 5505 for both states (10).
- `BLOCKED`: none.

The dry-run exited 0. `_runner_control` did not exist afterward, confirming that this inventory path wrote no control, launch, diagnostic, attempt, or output file.

## Archived artifact verification

`git ls-tree` found exactly 16 archived run files at `a553890`. A path-limited `git diff --name-only a553890` returned zero changed run files. The runner additionally matched each PDBQT SHA256 to its attempt record and verified exit code 0 plus 20 balanced, consecutively ranked models.

| Receptor state | Seed | Attempt JSON SHA256 | PDBQT SHA256 |
|---|---:|---|---|
| Asp542H(OD2) | 1101 | `abd95d401015a4255e07a76c9d4ab5fa51f9161aaf5ffadf1183473f8164fa05` | `04feea544c0f8dce8cc01a474524a81b6845dbe082cdb325308ed5ad0dbe9538` |
| Asp542H(OD2) | 2202 | `d610a7026e7cf398437d6f74d59fb4b3522a1a9503f3b5d728e3a762a62e8095` | `39e284646305fac0bc17170332225413adfe71a6db0bddc9ec1eb6793ed58d57` |
| Asp542H(OD2) | 3303 | `8fa8c519dad14b3dc7cd0861345a845897ec2fe1d93b590e18e249ee6116948a` | `909eab4df80d2f998cd62e64e2aeff634f3acac06a80a15577206f9bb81659ec` |
| Asp542H(OD2) | 4404 | `8fcb2d98e90422be4799fe752c1e4faf8ba4a681e5ab2cb40d54541e3a4e6bb6` | `90b876e6f3c00f8014a44d64939c350042a77ff2f86dadd655cefd67d76a5638` |
| Asp542- | 1101 | `aa3f4b633292c39118ca59363c8662094594e0d7f245a3db30ed8c4f66081f67` | `cf37dae4b5c4a7dd0fc187106937a39c505b257f3f8e75d5640c893fd98a3a27` |
| Asp542- | 2202 | `082de6ba66a5e3894d064569ad2beb953355dfe491358b44a1efa00d4970cf1f` | `af60ccb062ca6c0c347b7438d4103de86132d1d75a2e383f36202b7d286be151` |
| Asp542- | 3303 | `52888c347340517a7e60a97a7e649fcee3131b29767427c2b3b7ba9faa29a5ac` | `699bb55d6151445e07b11ba2ba45bc85bd2c2ca61f19a8869dc7d0f0a46a4178` |
| Asp542- | 4404 | `0abc50469ff82ab5f560dc03eebfe67639ed6ad70dd4adda9c5c8e2523df3a77` | `3179e5fe07dc2868527c53c486df2737e44afb5744bf4d0a66a3f8eae29b7826` |

## Remaining risks and limitations

- The underlying cause of the prior `python312.dll` access violation remains unknown. Faulthandler improves Python-side diagnostics but cannot guarantee a useful traceback for every native crash.
- A very small unavoidable interval remains between Windows returning a child PID and durable creation of its launch sidecar. A parent failure in that interval is mitigated by the independent process-table check and output-without-log blocking, but cannot be completely eliminated by the runner.
- Orphan detection is deliberately conservative and image-name based. Any running `vina_1.2.7_win.exe` blocks all missing jobs; it does not infer ownership from inaccessible command-line details.
- There is an unavoidable race between the process-table check and later process creation. Only one runner instance should be started.
- A hard interruption after launch leaves evidence and blocks automatic retry. That is intentional: stale launch evidence must be reconciled manually after confirming the recorded PID is gone and inspecting any output; it must not be deleted merely to make the runner proceed.
- The production path was not exercised with real Vina in Step 79. Only mock subprocess behavior and the read-only real-project inventory were tested.
- Faulthandler files use `.log` and are ignored by the repository-wide `*.log` rule. They persist locally, but a later evidence checkpoint must explicitly stage a selected diagnostic with `git add -f` if it is scientifically relevant. Session JSONL and launch JSON records are not ignored.

## Controlled Step 73 resumption procedure

1. Review and checkpoint only the hardened runner, isolated test, this protocol, and (if desired as provenance) `step78_partial_changes.txt`. Do not include the unrelated root files `pace'` or `taged files…`.
2. Re-run the isolated test command above and require all tests to pass.
3. Run `--dry-run` again immediately before execution and require exactly `REUSE=8 MISSING=12 BLOCKED=0`. Any other count is a stop condition.
4. Confirm that no other terminal, scheduler, Python runner, or `vina_1.2.7_win.exe` process is active. Start only one runner instance from a stable, non-disposable terminal/session.
5. Run the normal Step 73 command once:

   ```powershell
   .\.venv\Scripts\python.exe scripts\run_vina_3L4W_alternative_ligand_pilot.py
   ```

   The runner will revalidate the frozen inputs and saved preflight, reuse the eight completed searches, and submit only the 12 missing searches with at most two Vina children.
6. Preserve `_runner_control/sessions`, `_runner_control/launch_records`, and `_runner_control/diagnostics` together with all attempt records. Do not edit or remove them during execution.
7. If the runner terminates unexpectedly, do not restart it immediately. First run `--dry-run`. If any job is `BLOCKED`, confirm the recorded PID is absent, inspect the launch/session/faulthandler evidence and any PDBQT, and create a separate audit decision before permitting a retry. Never overwrite or delete the evidence or a partial output.
8. Only when all 20 jobs validate does the runner perform the unchanged aggregate pose analysis and write its write-once result CSV/manifest files.

At the end of Step 79, the implementation is ready for controlled resumption subject to an explicit code/documentation checkpoint and the pre-run `8/12/0` gate. It is not evidence that the unknown native crash has been eliminated.

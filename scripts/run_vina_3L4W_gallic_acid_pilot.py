"""Step 63: frozen gallic-acid 3L4W pose-only pilot (20 searches).

Run from the repository root with .venv/Scripts/python.exe. --preflight-only
checks inputs and Vina parsing without docking. Existing raw attempts are never
overwritten; completed attempts are analysed on rerun rather than repeated.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path
import re
import subprocess
import sys
import time

import meeko
import numpy as np
from rdkit import Chem, rdBase

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import run_vina_3L4W_phytochemical_pilot as pilot  # noqa: E402

INPUT_MANIFEST = ROOT / "data/metadata/vina_3L4W_gallic_acid_input_manifest.csv"
INPUT_MANIFEST_SHA = "148c510888ead5dabb7d7245e7648f60a2cd1737d29217e160529573852c035c"
INPUT_VALIDATOR = ROOT / "scripts/validate_vina_3L4W_gallic_acid_inputs.py"
OUTPUT_DIR = ROOT / "docking/pilot/vina/3L4W_gallic_acid_pose_only"
RESULTS_CSV = ROOT / "data/metadata/vina_3L4W_gallic_acid_pilot_results.csv"
POSE_CSV = OUTPUT_DIR / "all_pose_metrics.csv"
RUN_MANIFEST = OUTPUT_DIR / "experiment_manifest.json"
PREFLIGHT = OUTPUT_DIR / "vina_pair_parse_preflight.json"
INTERRUPTED = OUTPUT_DIR / "interrupted_attempts.json"
SEEDS = (1101, 2202, 3303, 4404, 5505)
CENTER = np.array([45.38664286, 92.35157143, 34.82200000])
SIZE = np.array([30.0, 24.0, 24.0])
pilot.SIZE = SIZE  # The reused frozen metrics and start routine refer to this constant.
pilot.CENTER = CENTER
RECEPTOR_LABEL = pilot.RECEPTOR_LABEL
RUN_FIELDS = [*pilot.RUN_FIELDS, "source_state_sdf_path", "source_state_sdf_sha256", "meeko_input_parse_status", "attempt_count", "attempt_log_paths", "attempt_log_hashes", "interrupted_attempt_record_path", "interrupted_attempt_record_sha256", "technical_interrupted_attempts_before_completion", "state_formal_charge", "box_size_A", "all_modes_charge_qc"]
POSE_FIELDS = [*pilot.POSE_FIELDS, "pdbqt_charge_qc", "nearest_box_face", "nearest_face_sdf_atom_1based", "nearest_face_element"]


def write_once(path: Path, content: bytes) -> None:
    pilot.write_once(path, content)


def rel(path: Path) -> str:
    return pilot.rel(path)


def verified_rows() -> tuple[list[dict[str, str]], str, Path]:
    check = subprocess.run([sys.executable, str(INPUT_VALIDATOR)], cwd=ROOT, capture_output=True, text=True, timeout=120)
    if check.returncode:
        raise RuntimeError(f"Step 62 validator failed:\n{check.stdout}\n{check.stderr}")
    if pilot.digest(INPUT_MANIFEST) != INPUT_MANIFEST_SHA:
        raise RuntimeError("Frozen Step 62 manifest SHA256 changed")
    rows = pilot.rows(INPUT_MANIFEST)
    expected = {(state, receptor) for state in ("neutral_acid_reference", "carboxylate_pH6p5_hypothesis") for receptor in RECEPTOR_LABEL}
    if len(rows) != 4 or {(r["microstate"], r["receptor_state"]) for r in rows} != expected:
        raise RuntimeError("Expected exactly four frozen LIG01 input pairs")
    executable = Path(r"E:\docking tools\vina_1.2.7_win.exe")
    if pilot.digest(executable) != rows[0]["vina_executable_sha256"]:
        raise RuntimeError("Vina executable SHA256 changed")
    version = subprocess.run([str(executable), "--version"], capture_output=True, text=True, timeout=15)
    if version.returncode or "1.2.7" not in version.stdout + version.stderr:
        raise RuntimeError("Vina is not 1.2.7")
    for r in rows:
        if (r["structure_id"], r["decision"], r["heavy_atoms"], r["stereocentres"], r["box_center_A"], r["box_size_A"], r["scoring"], r["exhaustiveness"], r["num_modes"], r["energy_range_kcal_mol"], r["min_rmsd_A"], r["cpu"], r["seeds"]) != (
            "LIG01", "PASS", "12", "0", "45.38664286;92.35157143;34.82200000", "30;24;24", "vina", "64", "20", "5.0", "1.0", "1", "1101;2202;3303;4404;5505"):
            raise RuntimeError("Frozen chemistry, box or search settings changed")
        for path_key, sha_key in (
            ("source_sdf_path", "source_sdf_sha256"), ("source_pdbqt_path", "source_pdbqt_sha256"),
            ("prepared_sdf_path", "prepared_sdf_sha256"), ("prepared_pdbqt_path", "prepared_pdbqt_sha256"),
            ("receptor_source_pdb_path", "receptor_source_pdb_sha256"), ("receptor_pdbqt_path", "receptor_pdbqt_sha256"),
            ("native_mig_sdf_path", "native_mig_sdf_sha256"), ("native_mig_pdbqt_path", "native_mig_pdbqt_sha256"),
            ("large_box_mig_results_path", "large_box_mig_results_sha256"),
        ):
            path = ROOT / r[path_key]
            if not path.is_file() or pilot.digest(path) != r[sha_key]:
                raise RuntimeError(f"Frozen SHA256 mismatch: {path}")
        # Adapt only field names required by the archived pose-analysis module.
        r["prepared_variant_sdf_path"] = "" if r["microstate"] == "neutral_acid_reference" else r["prepared_sdf_path"]
        r["prepared_variant_pdbqt_path"] = r["prepared_pdbqt_path"]
        r["prepared_variant_pdbqt_sha256"] = r["prepared_pdbqt_sha256"]
    controls = pilot.rows(ROOT / rows[0]["large_box_mig_results_path"])
    if len(controls) != 10 or any(c["status"] != "completed" or c["exit_code"] != "0" or float(c["rank1_RMSD_A"]) > 2.0 for c in controls):
        raise RuntimeError("Enlarged-box native MIG controls no longer pass 10/10")
    return rows, INPUT_MANIFEST_SHA, executable


def parse_check(r: dict[str, str], start: Path, executable: Path) -> dict:
    command = [str(executable), "--receptor", str(ROOT / r["receptor_pdbqt_path"]), "--ligand", str(start),
               "--center_x", "45.38664286", "--center_y", "92.35157143", "--center_z", "34.82200000",
               "--size_x", "30", "--size_y", "24", "--size_z", "24", "--scoring", "vina", "--score_only", "--cpu", "1"]
    t0 = time.perf_counter()
    completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=180)
    return {"microstate": r["microstate"], "receptor_state": r["receptor_state"], "command": subprocess.list2cmdline(command),
            "start_path": rel(start), "start_sha256": pilot.digest(start), "receptor_sha256": r["receptor_pdbqt_sha256"],
            "ligand_source_sha256": r["prepared_pdbqt_sha256"], "exit_code": completed.returncode,
            "runtime_seconds": round(time.perf_counter() - t0, 3), "stdout": completed.stdout, "stderr": completed.stderr,
            "interpretation": "Vina input/grid parse only; score not used or compared"}


def preflight() -> tuple[list[dict[str, str]], str, Path, dict, dict, np.ndarray, list[dict]]:
    frozen, manifest_sha, executable = verified_rows()
    mig = next(iter(Chem.SDMolSupplier(str(ROOT / frozen[0]["native_mig_sdf_path"]), removeHs=False)))
    if mig is None:
        raise RuntimeError("MIG crystal reference cannot be parsed")
    mig_xyz = np.asarray([mig.GetConformer().GetAtomPosition(a.GetIdx()) for a in mig.GetAtoms() if a.GetAtomicNum() > 1])
    if len(mig_xyz) != 14 or np.max(np.abs(mig_xyz.mean(axis=0) - CENTER)) > 1e-7:
        raise RuntimeError("MIG crystal centroid changed")
    contexts = {r["microstate"]: pilot.ligand_context(r) for r in frozen}
    receptor_contexts = {r["receptor_state"]: pilot.receptor_context(ROOT / r["receptor_pdbqt_path"]) for r in frozen}
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if PREFLIGHT.exists():
        prior = json.loads(PREFLIGHT.read_text(encoding="utf-8"))
        if prior.get("frozen_manifest_sha256") != manifest_sha or len(prior.get("parse_checks", [])) != 4:
            raise RuntimeError("Previous Vina parse-preflight manifest is incomplete or changed")
        parse_results = prior["parse_checks"]
        for item in parse_results:
            if pilot.digest(ROOT / item["start_path"]) != item["start_sha256"]:
                raise RuntimeError("Saved Vina parse-preflight start SHA256 changed")
    else:
        parse_results = []
        for r in frozen:
            ctx = contexts[r["microstate"]]
            start_data, _ = pilot.random_start(r, SEEDS[0], ctx)
            state_label = "neutral" if r["formal_charge"] == "0" else "carboxylate"
            preflight_start = OUTPUT_DIR / f"parse_{state_label}_seed_1101_start.pdbqt"
            write_once(preflight_start, start_data)
            parse_results.append(parse_check(r, preflight_start, executable))
        write_once(PREFLIGHT, (json.dumps({"frozen_manifest_sha256": manifest_sha, "parse_checks": parse_results}, indent=2) + "\n").encode("utf-8"))
    if {(x["microstate"], x["receptor_state"]) for x in parse_results} != {(r["microstate"], r["receptor_state"]) for r in frozen}:
        raise RuntimeError("Vina parse-preflight pair set changed")
    if any(x["exit_code"] != 0 for x in parse_results):
        raise RuntimeError("At least one receptor/ligand pair failed Vina input parsing; see preflight JSON")
    return frozen, manifest_sha, executable, contexts, receptor_contexts, mig_xyz, parse_results


def run_one(job: dict) -> dict:
    output, log = job["output"], job["log"]
    if log.exists() or output.exists():
        if not log.is_file() or not output.is_file():
            return {"status": "failed_prior_incomplete_artifacts", "exit_code": "", "runtime": "", "error": "Existing attempt has only one of log/output; preserved", "attempt_logs": [log] if log.exists() else []}
        previous = json.loads(log.read_text(encoding="utf-8"))
        if previous.get("return_code") == 0:
            return {"status": "completed", "exit_code": 0, "runtime": previous.get("elapsed_seconds", ""), "error": "", "attempt_logs": [log]}
        return {"status": "failed_prior_attempt", "exit_code": previous.get("return_code", ""), "runtime": previous.get("elapsed_seconds", ""), "error": previous.get("stderr", ""), "attempt_logs": [log]}
    start_time = time.perf_counter()
    try:
        completed = subprocess.run(job["cmd"], cwd=ROOT, capture_output=True, text=True, timeout=5400)
        elapsed = round(time.perf_counter() - start_time, 3)
        record = {"command": job["cmd"], "command_display": subprocess.list2cmdline(job["cmd"]), "start_sha256": pilot.digest(job["start"]),
                  "ligand_source_sha256": job["row"]["prepared_pdbqt_sha256"], "receptor_sha256": job["row"]["receptor_pdbqt_sha256"],
                  "return_code": completed.returncode, "elapsed_seconds": elapsed, "stdout": completed.stdout, "stderr": completed.stderr}
        write_once(log, (json.dumps(record, indent=2) + "\n").encode("utf-8"))
        status = "completed" if completed.returncode == 0 and output.is_file() else "vina_failed"
        return {"status": status, "exit_code": completed.returncode, "runtime": elapsed, "error": completed.stderr.strip(), "attempt_logs": [log]}
    except subprocess.TimeoutExpired as exc:
        elapsed = round(time.perf_counter() - start_time, 3)
        write_once(log, (json.dumps({"command": job["cmd"], "timeout_seconds": 5400, "elapsed_seconds": elapsed,
                                     "stdout": str(exc.stdout), "stderr": str(exc.stderr)}, indent=2) + "\n").encode("utf-8"))
        return {"status": "timeout_with_saved_output" if output.is_file() else "vina_timeout", "exit_code": "unconfirmed", "runtime": elapsed,
                "error": "Timeout; process exit unconfirmed", "attempt_logs": [log]}


def pose_charge_and_face(segment: str, ctx: dict, xyz: np.ndarray, formal_charge: int) -> dict:
    atoms = pilot.parse_pdbqt_atoms(segment)
    original = ctx["atoms"]
    if set(atoms) != set(original) or any(abs(atoms[s]["charge"] - original[s]["charge"]) > 0.0011 for s in original):
        raise ValueError("Pose PDBQT atom set or partial charges differ from input")
    if abs(sum(a["charge"] for a in atoms.values()) - formal_charge) > 0.02:
        raise ValueError("Pose PDBQT partial-charge sum differs from formal-charge hypothesis")
    margins = np.column_stack((
        xyz[:, 0] - (CENTER[0] - SIZE[0] / 2), (CENTER[0] + SIZE[0] / 2) - xyz[:, 0],
        xyz[:, 1] - (CENTER[1] - SIZE[1] / 2), (CENTER[1] + SIZE[1] / 2) - xyz[:, 1],
        xyz[:, 2] - (CENTER[2] - SIZE[2] / 2), (CENTER[2] + SIZE[2] / 2) - xyz[:, 2],
    ))
    i, j = np.unravel_index(np.argmin(margins), margins.shape)
    i, j = int(i), int(j)
    return {"pdbqt_charge_qc": "pass", "nearest_box_face": ("X-", "X+", "Y-", "Y+", "Z-", "Z+")[j],
            "nearest_face_sdf_atom_1based": i + 1, "nearest_face_element": ctx["heavy"].GetAtomWithIdx(i).GetSymbol()}


def analyse(frozen: list[dict], manifest_sha: str, contexts: dict, receptors: dict, mig_xyz: np.ndarray, jobs: list[dict], executable: Path) -> dict:
    interruption_record = json.loads(INTERRUPTED.read_text(encoding="utf-8")) if INTERRUPTED.exists() else {"attempts": []}
    interruptions = {(x["microstate"], x["receptor_state"], int(x["seed"])): x for x in interruption_record["attempts"]}
    if len(interruptions) != len(interruption_record["attempts"]):
        raise RuntimeError("Duplicate technical interruption record")
    run_rows: list[dict] = []
    pose_rows: list[dict] = []
    rank1_coords: dict[tuple[str, str], dict[int, np.ndarray]] = {}
    pose_coords: dict[tuple[str, str, int, int], np.ndarray] = {}
    for job in jobs:
        r, ctx, info = job["row"], job["ctx"], job["run_info"]
        state, receptor, seed = r["microstate"], r["receptor_state"], job["seed"]
        interrupted = (state, receptor, seed) in interruptions
        row = {"structure_id": "LIG01", "compound_name": "Gallic acid", "microstate": state, "receptor_state": receptor, "seed": seed,
               "run_status": info["status"], "receptor_path": r["receptor_pdbqt_path"], "receptor_sha256": r["receptor_pdbqt_sha256"],
               "ligand_path": r["prepared_pdbqt_path"], "ligand_sha256": r["prepared_pdbqt_sha256"],
               "source_sdf_path": r["source_sdf_path"], "source_sdf_sha256": r["source_sdf_sha256"],
               "source_state_sdf_path": r["prepared_sdf_path"], "source_state_sdf_sha256": r["prepared_sdf_sha256"],
               "crystal_reference_path": r["native_mig_sdf_path"], "crystal_reference_sha256": r["native_mig_sdf_sha256"],
               "frozen_manifest_sha256": manifest_sha, "vina_executable_sha256": r["vina_executable_sha256"],
               "start_path": rel(job["start"]), "start_sha256": pilot.digest(job["start"]),
               "start_attempt": job["random_meta"]["attempt"],
               "start_translation_A": ";".join(map(str, job["random_meta"]["translation_A"])),
               "start_rotation_quaternion": ";".join(map(str, job["random_meta"]["rotation_quaternion"])),
               "command": subprocess.list2cmdline(job["cmd"]), "log_path": rel(job["log"]) if job["log"].exists() else "",
               "log_sha256": pilot.digest(job["log"]) if job["log"].exists() else "",
               "output_path": rel(job["output"]) if job["output"].exists() else "",
               "output_sha256": pilot.digest(job["output"]) if job["output"].exists() else "",
               "exit_code": info["exit_code"], "runtime_seconds": info["runtime"], "error": info["error"],
               "meeko_input_parse_status": "pass", "attempt_count": len(info["attempt_logs"]) + int(interrupted),
               "attempt_log_paths": ";".join(rel(p) for p in info["attempt_logs"]),
               "attempt_log_hashes": ";".join(pilot.digest(p) for p in info["attempt_logs"]),
               "interrupted_attempt_record_path": rel(INTERRUPTED) if interrupted else "",
               "interrupted_attempt_record_sha256": pilot.digest(INTERRUPTED) if interrupted else "",
               "technical_interrupted_attempts_before_completion": int(interrupted),
               "state_formal_charge": r["formal_charge"], "box_size_A": "30;24;24"}
        if job["output"].is_file() and job["log"].is_file():
            try:
                parsed, coords = pilot.parse_output(job["output"], ctx, receptors[receptor], mig_xyz, row)
                row["n_modes_returned"] = len(parsed)
                segments = re.findall(r"(?ms)^MODEL\s+(\d+)\s*\n(.*?)^ENDMDL\s*$", job["output"].read_text(encoding="utf-8"))
                if len(segments) != len(parsed):
                    raise ValueError("Mode count changed during charge/face QC")
                for pose, (_, segment) in zip(parsed, segments):
                    if pose["pose_status"] == "valid":
                        try:
                            pose.update(pose_charge_and_face(segment, ctx, coords[int(pose["pose_rank"])], int(r["formal_charge"])))
                        except Exception as exc:
                            pose["pose_status"] = "invalid"
                            pose["error"] = str(exc)
                pose_rows.extend(parsed)
                row["all_modes_charge_qc"] = "pass" if all(p["pdbqt_charge_qc"] == "pass" for p in parsed if p["pose_status"] == "valid") and all(p["pose_status"] == "valid" for p in parsed) else "failed"
                for pose in parsed:
                    rank = int(pose["pose_rank"])
                    if pose["pose_status"] == "valid":
                        pose_coords[(state, receptor, seed, rank)] = coords[rank]
                top = parsed[0]
                if top["pose_status"] == "valid":
                    rank1_coords.setdefault((state, receptor), {})[seed] = coords[1]
                    row.update({"rank1_score": top["score"], "rank1_pocket_member": top["pocket_member"],
                                "rank1_severe_clash": top["severe_clash"],
                                "rank1_core_centroid_MIG_distance_A": top["core_centroid_MIG_distance_A"],
                                "rank1_core_MIG_overlap_fraction": top["core_MIG_overlap_fraction"],
                                "rank1_min_box_margin_A": top["min_box_margin_A"],
                                "rank1_catalytic_contacts": top["catalytic_contacts_within_4A"]})
                if any(p["pose_status"] != "valid" for p in parsed):
                    row["run_status"] = "pose_QC_failed"
                    row["error"] = "; ".join(p.get("error", "") for p in parsed if p["pose_status"] != "valid")
                elif info["status"] == "timeout_with_saved_output":
                    row["run_status"] = "complete_pose_output_exit_unconfirmed"
                elif info["status"] == "completed":
                    row["run_status"] = "completed"
            except Exception as exc:
                row["run_status"] = "pose_parse_failed"
                row["error"] = str(exc)
        run_rows.append(row)
    decisions: dict[tuple[str, str], dict] = {}
    for r in frozen:
        key = (r["microstate"], r["receptor_state"])
        if key in decisions:
            continue
        ctx = contexts[key[0]]
        available = rank1_coords.get(key, {})
        medoid = min(sorted(available), key=lambda s: (sum(pilot.symmetry_core_rmsd(available[s], x, ctx) for x in available.values()), s)) if available else None
        passing = 0
        for row in (x for x in run_rows if (x["microstate"], x["receptor_state"]) == key):
            if medoid is not None and int(row["seed"]) in available:
                rmsd = pilot.symmetry_core_rmsd(available[int(row["seed"])], available[medoid], ctx)
                row["rank1_core_RMSD_to_medoid_A"] = f"{rmsd:.4f}"
                row["medoid_seed"] = medoid
                if row["run_status"] in ("completed", "complete_pose_output_exit_unconfirmed") and row["rank1_pocket_member"] is True and row["rank1_severe_clash"] is False and rmsd <= 2.0:
                    passing += 1
        group = [x for x in run_rows if (x["microstate"], x["receptor_state"]) == key]
        if len(available) < 4 or any(x["run_status"] not in ("completed", "complete_pose_output_exit_unconfirmed") for x in group):
            status = "BLOCKED"
        else:
            status = "PASS" if passing >= 4 else "INCONCLUSIVE"
        decisions[key] = {"status": status, "passing_rank1_of_5": passing, "medoid_seed": medoid}
        for row in group:
            row["combination_status"] = status
            row["passing_rank1_of_5"] = passing
        for pose in (x for x in pose_rows if (x["microstate"], x["receptor_state"]) == key):
            pose["medoid_seed"] = medoid or ""
            pose_key = (key[0], key[1], int(pose["seed"]), int(pose["pose_rank"]))
            if pose["pose_status"] == "valid" and medoid is not None and pose_key in pose_coords:
                pose["core_RMSD_to_medoid_A"] = f"{pilot.symmetry_core_rmsd(pose_coords[pose_key], available[medoid], ctx):.4f}"
    write_once(RESULTS_CSV, pilot.csv_bytes(run_rows, RUN_FIELDS))
    write_once(POSE_CSV, pilot.csv_bytes(pose_rows, POSE_FIELDS))
    metadata = {
        "experiment": "3L4W gallic acid pose-only pilot", "script_path": rel(Path(__file__)), "script_sha256": pilot.digest(Path(__file__)),
        "frozen_input_manifest_path": rel(INPUT_MANIFEST), "frozen_input_manifest_sha256": manifest_sha,
        "validator_path": rel(INPUT_VALIDATOR), "validator_sha256": pilot.digest(INPUT_VALIDATOR),
        "vina_executable": str(executable), "vina_executable_sha256": pilot.digest(executable), "vina_version": "1.2.7",
        "python_version": sys.version.split()[0], "rdkit_version": rdBase.rdkitVersion, "meeko_version": meeko.__version__,
        "box_center_A": CENTER.tolist(), "box_size_A": SIZE.tolist(), "scoring": "vina", "exhaustiveness": 64,
        "num_modes": 20, "energy_range": 5.0, "min_rmsd": 1.0, "cpu": 1, "concurrent_processes": 1,
        "seeds": list(SEEDS), "start_method": "NumPy PCG64 SeedSequence(seed, PubChem CID 370, charge code); seeded random quaternion and translation [-3,3] A; same seed-specific start for both receptors of each ligand state; >=0.5 A box clearance at start",
        "pose_rules": {"pocket": "all heavy atoms inside box; core centroid <=6 A from MIG centroid; >=50% core atoms within 4 A of MIG heavy atoms",
                       "clash": "receptor-ligand heavy distance <1.8 A OR VdW overlap >0.8 A; internal nonbonded >=4-bond heavy distance <1.6 A",
                       "convergence": ">=4/5 rank1 poses per combination pocket+no severe clash+<=2 A symmetry-aware unfitted core RMSD to medoid"},
        "MIG_control_results_path": frozen[0]["large_box_mig_results_path"], "MIG_control_results_sha256": frozen[0]["large_box_mig_results_sha256"],
        "pair_parse_preflight_path": rel(PREFLIGHT), "pair_parse_preflight_sha256": pilot.digest(PREFLIGHT),
        "interrupted_attempt_record_path": rel(INTERRUPTED) if INTERRUPTED.exists() else "",
        "interrupted_attempt_record_sha256": pilot.digest(INTERRUPTED) if INTERRUPTED.exists() else "",
        "technical_interrupted_attempts": len(interruptions),
        "source_hashes": {"/".join((r["microstate"], r["receptor_state"])): {
            "source_sdf": r["source_sdf_sha256"], "state_sdf": r["prepared_sdf_sha256"],
            "ligand_pdbqt": r["prepared_pdbqt_sha256"], "receptor_pdbqt": r["receptor_pdbqt_sha256"],
            "MIG_reference_sdf": r["native_mig_sdf_sha256"]} for r in frozen},
        "results_csv_path": rel(RESULTS_CSV), "results_csv_sha256": pilot.digest(RESULTS_CSV),
        "all_pose_csv_path": rel(POSE_CSV), "all_pose_csv_sha256": pilot.digest(POSE_CSV),
        "confirmed_exit_zero_runs": sum(r["run_status"] == "completed" for r in run_rows),
        "exit_unconfirmed_runs": sum(r["run_status"] == "complete_pose_output_exit_unconfirmed" for r in run_rows),
        "failed_runs": sum(r["run_status"] not in ("completed", "complete_pose_output_exit_unconfirmed") for r in run_rows),
        "returned_modes": len(pose_rows),
        "decisions": {"/".join(k): v for k, v in decisions.items()},
    }
    write_once(RUN_MANIFEST, (json.dumps(metadata, indent=2, sort_keys=True) + "\n").encode("utf-8"))
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--prefill-tail", action="store_true", help="Run only the last nine frozen jobs; main runner later reads their logs")
    args = parser.parse_args()
    frozen, manifest_sha, executable, contexts, receptors, mig_xyz, checks = preflight()
    print(f"Preflight PASS: {len(checks)}/4 Vina pair parses; Step 62 and MIG hashes verified", flush=True)
    if args.preflight_only:
        return 0
    jobs = []
    for r in frozen:
        ctx = contexts[r["microstate"]]
        state_label = "neutral" if r["formal_charge"] == "0" else "carboxylate"
        pair_dir = OUTPUT_DIR / f"LIG01_{state_label}_{RECEPTOR_LABEL[r['receptor_state']]}"
        pair_dir.mkdir(parents=True, exist_ok=True)
        for seed in SEEDS:
            prefix = pair_dir / f"seed_{seed}"
            start = prefix.with_name(prefix.name + "_start.pdbqt")
            output = prefix.with_name(prefix.name + "_docked.pdbqt")
            log = prefix.with_name(prefix.name + "_vina.json")
            start_data, meta = pilot.random_start(r, seed, ctx)
            write_once(start, start_data)
            cmd = [str(executable), "--receptor", str(ROOT / r["receptor_pdbqt_path"]), "--ligand", str(start),
                   "--center_x", "45.38664286", "--center_y", "92.35157143", "--center_z", "34.82200000",
                   "--size_x", "30", "--size_y", "24", "--size_z", "24", "--scoring", "vina",
                   "--exhaustiveness", "64", "--num_modes", "20", "--energy_range", "5.0", "--min_rmsd", "1.0",
                   "--seed", str(seed), "--cpu", "1", "--out", str(output), "--verbosity", "1"]
            jobs.append({"row": r, "ctx": ctx, "seed": seed, "start": start, "output": output, "log": log,
                         "cmd": cmd, "random_meta": meta})
    if len(jobs) != 20:
        raise RuntimeError("Expected exactly twenty planned searches")
    if args.prefill_tail:
        for n, job in enumerate(jobs[11:], 12):
            job["run_info"] = run_one(job)
            print(f"Tail {n}/20 {job['row']['microstate']} {job['row']['receptor_state']} seed {job['seed']}: {job['run_info']['status']}", flush=True)
        return 0
    for n, job in enumerate(jobs, 1):
        job["run_info"] = run_one(job)
        print(f"{n}/20 {job['row']['microstate']} {job['row']['receptor_state']} seed {job['seed']}: {job['run_info']['status']}", flush=True)
    result = analyse(frozen, manifest_sha, contexts, receptors, mig_xyz, jobs, executable)
    print(f"Analysed {result['confirmed_exit_zero_runs']}/20 confirmed runs; {result['failed_runs']} failed; {result['exit_unconfirmed_runs']} exit-unconfirmed; {result['returned_modes']} modes", flush=True)
    return 0 if not result["failed_runs"] else 2


if __name__ == "__main__":
    raise SystemExit(main())

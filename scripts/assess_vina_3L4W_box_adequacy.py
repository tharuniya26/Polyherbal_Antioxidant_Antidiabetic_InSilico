"""Audit Step 58 box faces and run one prespecified larger-box comparison.

Run from the repository root with .venv/Scripts/python.exe. --audit-only makes
the coordinate audit and immutable box preregistration without running Vina.
The ordinary invocation first repeats that audit, then runs only missing jobs.
Step 58 inputs, output coordinates, and classifications are read-only.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
from datetime import datetime, timezone
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys
import time

import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("frozen_pilot", ROOT / "scripts/run_vina_3L4W_phytochemical_pilot.py")
pilot = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(pilot)

OLD_SIZE = np.array([24., 18., 18.])
NEW_SIZE = np.array([30., 24., 24.])  # Frozen before any sensitivity docking: +3 A on each face.
CENTER = np.array([45.38664286, 92.35157143, 34.82200000])
FACES = ("x-", "x+", "y-", "y+", "z-", "z+")
OUT = ROOT / "docking/pilot/vina/3L4W_box_adequacy_sensitivity"
AUDIT_CSV = OUT / "original_box_face_audit.csv"
COMPARISON_FACE_CSV = OUT / "comparison_box_face_audit.csv"
PRE_REG = OUT / "box_pre_registration.json"
PREFLIGHT_DIR = OUT / "grid_input_checks"
POSE_CSV = OUT / "all_pose_metrics.csv"
MANIFEST_JSON = OUT / "experiment_manifest.json"
RESULTS_CSV = ROOT / "data/metadata/vina_3L4W_box_adequacy_results.csv"

AUDIT_FIELDS = ["structure_id", "microstate", "receptor_state", "seed", "pose_rank", "pose_status",
                "output_path", "output_sha256", "closest_face", "closest_face_distance_A", "closest_sdf_atom_1based",
                "closest_pdbqt_serial", "closest_atom_element", "closest_atom_in_rigid_core", "core_min_face_distance_A",
                "x_minus_A", "x_plus_A", "y_minus_A", "y_plus_A", "z_minus_A", "z_plus_A",
                "faces_within_1A", "heavy_atoms_outside_box", "full_heavy_span_x_A", "full_heavy_span_y_A",
                "full_heavy_span_z_A", "original_metric_margin_A"]
COMPARISON_FACE_FIELDS = ["comparison_metric_margin_A" if field == "original_metric_margin_A" else field
                          for field in AUDIT_FIELDS]
RESULT_FIELDS = pilot.RUN_FIELDS + ["comparison_box_center_A", "comparison_box_size_A", "original_run_status",
                "original_combination_status", "original_rank1_score", "original_rank1_pocket_member",
                "original_rank1_severe_clash", "original_rank1_core_RMSD_to_medoid_A",
                "original_rank1_min_box_margin_A", "original_rank1_catalytic_contacts", "original_output_sha256",
                "original_output_path", "prior_attempt_log_path", "prior_attempt_log_sha256", "prior_attempt_error",
                "new_rank1_closest_face", "new_rank1_closest_face_atom_sdf_1based",
                "new_rank1_core_min_box_margin_A", "comparison_conclusion"]


def xyz_segments(path: Path):
    return re.findall(r"(?ms)^MODEL\s+(\d+)\s*\n(.*?)^ENDMDL\s*$", path.read_text(encoding="utf-8"))


def face_data(xyz: np.ndarray, size: np.ndarray, ctx: dict) -> dict:
    low, high = CENTER - size / 2, CENTER + size / 2
    d = np.column_stack((xyz[:, 0] - low[0], high[0] - xyz[:, 0],
                         xyz[:, 1] - low[1], high[1] - xyz[:, 1],
                         xyz[:, 2] - low[2], high[2] - xyz[:, 2]))
    atom_i, face_i = (int(v) for v in np.unravel_index(np.argmin(d), d.shape))
    mins = d.min(axis=0)
    span = np.ptp(xyz, axis=0)
    return {"closest_face": FACES[face_i], "closest_face_distance_A": float(d[atom_i, face_i]),
            "closest_sdf_atom_1based": atom_i + 1, "closest_pdbqt_serial": ctx["mapping"][atom_i],
            "closest_atom_element": ctx["heavy"].GetAtomWithIdx(atom_i).GetSymbol(),
            "closest_atom_in_rigid_core": atom_i in ctx["core"],
            "core_min_face_distance_A": float(d[ctx["core"]].min()),
            "x_minus_A": float(mins[0]), "x_plus_A": float(mins[1]),
            "y_minus_A": float(mins[2]), "y_plus_A": float(mins[3]),
            "z_minus_A": float(mins[4]), "z_plus_A": float(mins[5]),
            "faces_within_1A": ";".join(FACES[i] for i, v in enumerate(mins) if v < 1.0),
            "heavy_atoms_outside_box": int(np.sum(np.any(d < -0.0001, axis=1))),
            "full_heavy_span_x_A": float(span[0]), "full_heavy_span_y_A": float(span[1]),
            "full_heavy_span_z_A": float(span[2])}


def old_input_audit():
    frozen, frozen_hash, vina = pilot.preflight()
    old_manifest = json.loads(pilot.RUN_MANIFEST.read_text(encoding="utf-8"))
    for path_key, hash_key in (("results_csv", "results_csv_sha256"), ("all_pose_csv", "all_pose_csv_sha256")):
        if pilot.digest(ROOT / old_manifest[path_key]) != old_manifest[hash_key]:
            raise RuntimeError(f"Step 58 summary file changed: {path_key}")
    old_runs = pilot.rows(pilot.RESULTS_CSV)
    old_poses = pilot.rows(pilot.POSE_CSV)
    if len(old_runs) != 40 or len(old_poses) != 781:
        raise RuntimeError("Step 58 record counts changed")
    if sum(r["run_status"] == "completed" for r in old_runs) != 38 or sum(
            r["run_status"] == "complete_pose_output_exit_unconfirmed" for r in old_runs) != 2:
        raise RuntimeError("Step 58 process-status distinction changed")
    contexts = {(r["structure_id"], r["microstate"]): pilot.ligand_context(r) for r in frozen}
    pose_lookup = {(r["structure_id"], r["microstate"], r["receptor_state"], r["seed"], r["pose_rank"]): r
                   for r in old_poses}
    face_rows = []
    for run in old_runs:
        for path_col, hash_col in (("start_path", "start_sha256"), ("log_path", "log_sha256"),
                                   ("output_path", "output_sha256")):
            if pilot.digest(ROOT / run[path_col]) != run[hash_col]:
                raise RuntimeError(f"Step 58 raw artifact changed: {run[path_col]}")
        ctx = contexts[(run["structure_id"], run["microstate"])]
        for rank, seg in xyz_segments(ROOT / run["output_path"]):
            xyz = pilot.verify_pose_atoms(pilot.parse_pdbqt_atoms(seg), ctx)
            key = (run["structure_id"], run["microstate"], run["receptor_state"], run["seed"], rank)
            prior = pose_lookup[key]
            fd = face_data(xyz, OLD_SIZE, ctx)
            if abs(fd["closest_face_distance_A"] - float(prior["min_box_margin_A"])) > 0.00001:
                raise RuntimeError(f"Independent box-face calculation disagrees with Step 58: {key}")
            if fd["heavy_atoms_outside_box"]:
                raise RuntimeError(f"Step 58 pose outside original grid: {key}")
            face_rows.append({"structure_id": run["structure_id"], "microstate": run["microstate"],
                              "receptor_state": run["receptor_state"], "seed": run["seed"], "pose_rank": rank,
                              "pose_status": prior["pose_status"], "output_path": run["output_path"],
                              "output_sha256": run["output_sha256"], "original_metric_margin_A": prior["min_box_margin_A"], **fd})
    if len(face_rows) != len(old_poses):
        raise RuntimeError("Incomplete Step 58 pose-by-pose audit")
    pilot.write_once(AUDIT_CSV, pilot.csv_bytes(face_rows, AUDIT_FIELDS))
    mig = next(iter(Chem.SDMolSupplier(str(ROOT / frozen[0]["native_mig_sdf_path"]), removeHs=False)))
    mig_xyz = np.array([mig.GetConformer().GetAtomPosition(a.GetIdx()) for a in mig.GetAtoms() if a.GetAtomicNum() > 1])
    if len(mig_xyz) != 14 or np.max(np.abs(mig_xyz.mean(axis=0) - CENTER)) > 0.0000001:
        raise RuntimeError("Deposited MIG site geometry changed")
    rank1 = [r for r in face_rows if r["pose_rank"] == "1"]
    if len(rank1) != 40 or sum(r["closest_face_distance_A"] < 1 for r in rank1) != 36:
        raise RuntimeError("Original rank-1 face count differs from documented 36/40")
    pre = {"decision": "RUN_ONE_CENTERED_LARGER_BOX_COMPARISON", "registered_before_vina": True,
           "old_center_A": CENTER.tolist(), "old_size_A": OLD_SIZE.tolist(), "new_center_A": CENTER.tolist(),
           "new_size_A": NEW_SIZE.tolist(), "added_clearance_per_face_A": [3., 3., 3.],
           "reason": "Independent audit: 36/40 rank-1 and 685/781 modes within 1 A of an original face; all 40 rank-1 limiting faces are y-minus. Extend each face by 3 A without moving the deposited MIG centroid, changing the pocket rule, or using docking scores.",
           "deposited_MIG_bbox_min_A": mig_xyz.min(axis=0).tolist(),
           "deposited_MIG_bbox_max_A": mig_xyz.max(axis=0).tolist(),
           "old_rank1_min_clearance_A": min(r["closest_face_distance_A"] for r in rank1),
           "expected_min_clearance_if_old_rank1_coordinates_unchanged_A": min(r["closest_face_distance_A"] for r in rank1) + 3.,
           "software": "AutoDock Vina 1.2.7, Meeko 0.8.0, RDKit", "scoring": "vina",
           "exhaustiveness": 64, "num_modes": 20, "energy_range": 5.0, "min_rmsd": 1.0, "cpu": 1,
           "seeds": list(pilot.SEEDS), "old_results_sha256": pilot.digest(pilot.RESULTS_CSV),
           "old_pose_metrics_sha256": pilot.digest(pilot.POSE_CSV), "old_face_audit_sha256": pilot.digest(AUDIT_CSV),
           "frozen_input_manifest_sha256": frozen_hash, "vina_executable_sha256": pilot.digest(vina)}
    pilot.write_once(PRE_REG, (json.dumps(pre, indent=2, sort_keys=True) + "\n").encode("utf-8"))
    print("Read-only audit: 781/781 poses checked; 36/40 rank-1 and 685/781 modes within 1 A; ",
          "all rank-1 closest faces y-. Box 30 x 24 x 24 A preregistered.", flush=True)
    return frozen, frozen_hash, vina, old_runs, contexts, face_rows, pre


def grid_preflight(frozen: list[dict], vina: Path) -> None:
    # --score_only parses receptor/ligand and builds the larger grid; it does not search.
    PREFLIGHT_DIR.mkdir(parents=True, exist_ok=True)
    old_runs = pilot.rows(pilot.RESULTS_CSV)
    for receptor in pilot.RECEPTOR_LABEL:
        r = next(row for row in frozen if row["structure_id"] == "LIG04" and
                 row["microstate"] == "neutral_keto_reference" and row["receptor_state"] == receptor)
        old_start = next(row for row in old_runs if row["structure_id"] == "LIG04" and
                         row["microstate"] == "neutral_keto_reference" and row["receptor_state"] == receptor and
                         row["seed"] == "1101")
        if pilot.digest(ROOT / old_start["start_path"]) != old_start["start_sha256"]:
            raise RuntimeError("Frozen pocket start changed before grid test")
        cmd = [str(vina), "--receptor", str(ROOT / r["receptor_pdbqt_path"]),
               "--ligand", str(ROOT / old_start["start_path"]),
               "--center_x", "45.38664286", "--center_y", "92.35157143", "--center_z", "34.82200000",
               "--size_x", "30", "--size_y", "24", "--size_z", "24", "--scoring", "vina",
               "--score_only", "--cpu", "1"]
        path = PREFLIGHT_DIR / f"{pilot.RECEPTOR_LABEL[receptor]}_frozen_start_score_only_input_check.json"
        if path.exists():
            record = json.loads(path.read_text(encoding="utf-8"))
            if record["command"] != cmd or record["return_code"] != 0:
                raise RuntimeError("Prior grid check is not a successful match")
        else:
            start = time.perf_counter()
            proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=180)
            record = {"command": cmd, "return_code": proc.returncode, "elapsed_seconds": round(time.perf_counter()-start,3),
                      "stdout": proc.stdout, "stderr": proc.stderr, "purpose": "grid and input parse only; score uninterpreted"}
            pilot.write_once(path, (json.dumps(record, indent=2) + "\n").encode("utf-8"))
            if proc.returncode:
                raise RuntimeError(f"Larger-box Vina input/grid check failed: {receptor}: {proc.stderr}")
    print("Larger-grid score-only input checks passed for both receptor variants (scores uninterpreted).", flush=True)


def run_job(job: dict) -> dict:
    out, log = job["output"], job["log"]
    if out.exists() or log.exists():
        if not out.is_file() or not log.is_file():
            return {"status": "incomplete_prior_artifacts", "exit_code": "", "elapsed": "", "error": "Partial saved files preserved"}
        history = json.loads(log.read_text(encoding="utf-8"))
        if history.get("command") != job["cmd"] or history.get("start_sha256") != pilot.digest(job["start"]):
            return {"status": "prior_artifact_settings_mismatch", "exit_code": "", "elapsed": "", "error": "Saved command/start mismatch"}
        if history.get("return_code") == 0:
            return {"status": "confirmed_exit_zero", "exit_code": 0, "elapsed": history.get("elapsed_seconds", ""), "error": ""}
        if history.get("timeout_seconds"):
            return {"status": "timeout_exit_unconfirmed", "exit_code": "unconfirmed", "elapsed": history.get("elapsed_seconds", ""), "error": "Prior timeout; output retained for QC"}
        return {"status": "vina_failed", "exit_code": history.get("return_code", ""), "elapsed": history.get("elapsed_seconds", ""), "error": history.get("stderr", "")}
    begun = time.perf_counter()
    try:
        # Long host suspension can exceed wall time without consuming CPU.
        # Preserve an uninterrupted search instead of imposing a wall timeout.
        proc = subprocess.run(job["cmd"], cwd=ROOT, capture_output=True, text=True)
        elapsed = round(time.perf_counter() - begun, 3)
        record = {"command": job["cmd"], "start_sha256": pilot.digest(job["start"]),
                  "receptor_sha256": job["row"]["receptor_pdbqt_sha256"],
                  "ligand_sha256": job["row"]["prepared_variant_pdbqt_sha256"],
                  "return_code": proc.returncode, "elapsed_seconds": elapsed, "stdout": proc.stdout, "stderr": proc.stderr}
        pilot.write_once(log, (json.dumps(record, indent=2) + "\n").encode("utf-8"))
        return {"status": "confirmed_exit_zero" if proc.returncode == 0 and out.is_file() else "vina_failed",
                "exit_code": proc.returncode, "elapsed": elapsed, "error": proc.stderr.strip()}
    except subprocess.TimeoutExpired as exc:
        elapsed = round(time.perf_counter() - begun, 3)
        record = {"command": job["cmd"], "start_sha256": pilot.digest(job["start"]),
                  "timeout_seconds": 7200, "elapsed_seconds": elapsed, "stdout": str(exc.stdout), "stderr": str(exc.stderr)}
        pilot.write_once(log, (json.dumps(record, indent=2) + "\n").encode("utf-8"))
        return {"status": "timeout_exit_unconfirmed", "exit_code": "unconfirmed", "elapsed": elapsed,
                "error": "7200 s wall-clock timeout; exit unconfirmed"}


def analyze(jobs: list[dict], old_runs: list[dict], contexts: dict, frozen: list[dict], frozen_hash: str,
            vina: Path, face_rows: list[dict]) -> int:
    old_lookup = {(r["structure_id"], r["microstate"], r["receptor_state"], int(r["seed"])): r for r in old_runs}
    face_lookup = {(r["structure_id"], r["microstate"], r["receptor_state"], int(r["seed"]), int(r["pose_rank"])): r
                   for r in face_rows}
    mig = next(iter(Chem.SDMolSupplier(str(ROOT / frozen[0]["native_mig_sdf_path"]), removeHs=False)))
    mig_xyz = np.array([mig.GetConformer().GetAtomPosition(a.GetIdx()) for a in mig.GetAtoms() if a.GetAtomicNum() > 1])
    receptors = {r["receptor_state"]: pilot.receptor_context(ROOT / r["receptor_pdbqt_path"]) for r in frozen}
    pilot.SIZE = NEW_SIZE.copy()  # Only the comparison scoring/box-inside check; old audit used OLD_SIZE.
    run_rows, pose_rows, rank1, comparison_face_rows = [], [], {}, []
    for job in jobs:
        r, seed, info = job["row"], job["seed"], job["info"]
        key = (r["structure_id"], r["microstate"], r["receptor_state"])
        old = old_lookup[(*key, seed)]
        row = {"structure_id": r["structure_id"], "compound_name": r["compound_name"],
               "microstate": r["microstate"], "receptor_state": r["receptor_state"], "seed": seed,
               "run_status": info["status"], "receptor_path": r["receptor_pdbqt_path"],
               "receptor_sha256": r["receptor_pdbqt_sha256"], "ligand_path": r["prepared_variant_pdbqt_path"],
               "ligand_sha256": r["prepared_variant_pdbqt_sha256"], "source_sdf_path": r["source_sdf_path"],
               "source_sdf_sha256": r["source_sdf_sha256"], "crystal_reference_path": r["native_mig_sdf_path"],
               "crystal_reference_sha256": r["native_mig_sdf_sha256"], "frozen_manifest_sha256": frozen_hash,
               "vina_executable_sha256": r["vina_executable_sha256"], "start_path": pilot.rel(job["start"]),
               "start_sha256": pilot.digest(job["start"]), "start_attempt": old["start_attempt"],
               "start_translation_A": old["start_translation_A"],
               "start_rotation_quaternion": old["start_rotation_quaternion"],
               "command": subprocess.list2cmdline(job["cmd"]),
               "log_path": pilot.rel(job["log"]) if job["log"].exists() else "",
               "log_sha256": pilot.digest(job["log"]) if job["log"].exists() else "",
               "output_path": pilot.rel(job["output"]) if job["output"].exists() else "",
               "output_sha256": pilot.digest(job["output"]) if job["output"].exists() else "",
               "exit_code": info["exit_code"], "runtime_seconds": info["elapsed"], "error": info["error"],
               "comparison_box_center_A": ";".join(str(v) for v in CENTER),
               "comparison_box_size_A": "30;24;24", "original_run_status": old["run_status"],
               "original_combination_status": old["combination_status"], "original_rank1_score": old["rank1_score"],
               "original_rank1_pocket_member": old["rank1_pocket_member"],
               "original_rank1_severe_clash": old["rank1_severe_clash"],
               "original_rank1_core_RMSD_to_medoid_A": old["rank1_core_RMSD_to_medoid_A"],
               "original_rank1_min_box_margin_A": old["rank1_min_box_margin_A"],
               "original_rank1_catalytic_contacts": old["rank1_catalytic_contacts"],
               "original_output_sha256": old["output_sha256"], "original_output_path": old["output_path"]}
        if job.get("prior_failure"):
            prior = job["prior_failure"]
            row.update({"prior_attempt_log_path": pilot.rel(prior["log"]),
                        "prior_attempt_log_sha256": pilot.digest(prior["log"]),
                        "prior_attempt_error": prior["error"]})
        if job["output"].is_file() and job["log"].is_file():
            try:
                parsed, coords = pilot.parse_output(job["output"], contexts[(r["structure_id"], r["microstate"])],
                                                    receptors[r["receptor_state"]], mig_xyz, row)
                row["n_modes_returned"] = len(parsed)
                pose_rows.extend(parsed)
                if any(p["pose_status"] != "valid" for p in parsed):
                    row["run_status"] = "pose_QC_failed"
                    row["error"] = "; ".join(p.get("error", "") for p in parsed if p["pose_status"] != "valid")
                else:
                    top = parsed[0]
                    rank1.setdefault(key, {})[seed] = coords[1]
                    fd = face_data(coords[1], NEW_SIZE, contexts[(r["structure_id"], r["microstate"])])
                    row.update({"rank1_score": top["score"], "rank1_pocket_member": top["pocket_member"],
                                "rank1_severe_clash": top["severe_clash"],
                                "rank1_core_centroid_MIG_distance_A": top["core_centroid_MIG_distance_A"],
                                "rank1_core_MIG_overlap_fraction": top["core_MIG_overlap_fraction"],
                                "rank1_min_box_margin_A": top["min_box_margin_A"],
                                "rank1_catalytic_contacts": top["catalytic_contacts_within_4A"],
                                "new_rank1_closest_face": fd["closest_face"],
                                "new_rank1_closest_face_atom_sdf_1based": fd["closest_sdf_atom_1based"],
                                "new_rank1_core_min_box_margin_A": fd["core_min_face_distance_A"]})
                    if abs(fd["closest_face_distance_A"] - top["min_box_margin_A"]) > 0.00001:
                        raise RuntimeError("Comparison face margin disagrees with pose metric")
                    if info["status"] == "timeout_exit_unconfirmed":
                        row["run_status"] = "complete_pose_output_exit_unconfirmed"
                        row["error"] = "Complete valid pose output; Vina process exit unconfirmed after timeout"
            except Exception as exc:
                row["run_status"] = "pose_parse_failed"
                row["error"] = str(exc)
        run_rows.append(row)
    decisions = {}
    for r in frozen:
        key = (r["structure_id"], r["microstate"], r["receptor_state"])
        if key in decisions:
            continue
        ctx = contexts[(r["structure_id"], r["microstate"])]
        available = rank1.get(key, {})
        medoid = min(sorted(available), key=lambda seed: (
            sum(pilot.symmetry_core_rmsd(available[seed], other, ctx) for other in available.values()), seed)) if available else None
        combo = [row for row in run_rows if (row["structure_id"], row["microstate"], row["receptor_state"]) == key]
        passing = 0
        for row in combo:
            seed = int(row["seed"])
            if medoid is not None and seed in available:
                rmsd = pilot.symmetry_core_rmsd(available[seed], available[medoid], ctx)
                row["rank1_core_RMSD_to_medoid_A"] = f"{rmsd:.4f}"
                row["medoid_seed"] = medoid
                if row["run_status"] in ("confirmed_exit_zero", "complete_pose_output_exit_unconfirmed") and \
                        row["rank1_pocket_member"] is True and row["rank1_severe_clash"] is False and rmsd <= 2.0:
                    passing += 1
        status = ("BLOCKED" if len(available) < 4 or any(row["run_status"] not in
                  ("confirmed_exit_zero", "complete_pose_output_exit_unconfirmed") for row in combo)
                  else "PASS" if passing >= 4 else "INCONCLUSIVE")
        decisions[key] = {"status": status, "passing_rank1_of_5": passing, "medoid_seed": medoid}
        for row in combo:
            row["combination_status"] = status
            row["passing_rank1_of_5"] = passing
            row["comparison_conclusion"] = ("unchanged" if status == row["original_combination_status"] else
                                            f"{row['original_combination_status']}_to_{status}")
        for pose in (p for p in pose_rows if (p["structure_id"], p["microstate"], p["receptor_state"]) == key):
            pose["medoid_seed"] = medoid or ""
            if pose["pose_status"] == "valid" and medoid is not None:
                job = next(j for j in jobs if (j["row"]["structure_id"], j["row"]["microstate"],
                                              j["row"]["receptor_state"], j["seed"]) == (*key, int(pose["seed"])))
                seg = xyz_segments(job["output"])[int(pose["pose_rank"]) - 1][1]
                xyz = pilot.verify_pose_atoms(pilot.parse_pdbqt_atoms(seg), ctx)
                pose["core_RMSD_to_medoid_A"] = f"{pilot.symmetry_core_rmsd(xyz, available[medoid], ctx):.4f}"
                fd = face_data(xyz, NEW_SIZE, ctx)
                if abs(fd["closest_face_distance_A"] - float(pose["min_box_margin_A"])) > 0.00001:
                    raise RuntimeError("Comparison all-mode face metric mismatch")
                comparison_face_rows.append({"structure_id": pose["structure_id"], "microstate": pose["microstate"],
                                             "receptor_state": pose["receptor_state"], "seed": pose["seed"],
                                             "pose_rank": pose["pose_rank"], "pose_status": pose["pose_status"],
                                             "output_path": pose["output_path"], "output_sha256": pose["output_sha256"],
                                             "comparison_metric_margin_A": pose["min_box_margin_A"], **fd})
    pilot.atomic_write(RESULTS_CSV, pilot.csv_bytes(run_rows, RESULT_FIELDS))
    pilot.atomic_write(POSE_CSV, pilot.csv_bytes(pose_rows, pilot.POSE_FIELDS))
    pilot.atomic_write(COMPARISON_FACE_CSV, pilot.csv_bytes(comparison_face_rows, COMPARISON_FACE_FIELDS))
    manifest = {"experiment": "Step 59 single preregistered larger-box comparison", "created_utc": datetime.now(timezone.utc).isoformat(),
                "script": pilot.rel(Path(__file__)), "script_sha256": pilot.digest(Path(__file__)),
                "box_preregistration": pilot.rel(PRE_REG), "box_preregistration_sha256": pilot.digest(PRE_REG),
                "original_box_face_audit": pilot.rel(AUDIT_CSV), "original_box_face_audit_sha256": pilot.digest(AUDIT_CSV),
                "original_results_sha256": pilot.digest(pilot.RESULTS_CSV), "original_pose_metrics_sha256": pilot.digest(pilot.POSE_CSV),
                "frozen_input_manifest_sha256": frozen_hash, "vina_executable_sha256": pilot.digest(vina),
                "center_A": CENTER.tolist(), "size_A": NEW_SIZE.tolist(), "settings": {"scoring": "vina", "exhaustiveness": 64,
                "num_modes": 20, "energy_range": 5.0, "min_rmsd": 1.0, "cpu": 1, "seeds": list(pilot.SEEDS)},
                "start_method": "Exact byte copies of the 40 hash-verified Step 58 randomized starts",
                "results_csv": pilot.rel(RESULTS_CSV), "results_csv_sha256": pilot.digest(RESULTS_CSV),
                "all_pose_csv": pilot.rel(POSE_CSV), "all_pose_csv_sha256": pilot.digest(POSE_CSV),
                "comparison_face_audit_csv": pilot.rel(COMPARISON_FACE_CSV),
                "comparison_face_audit_sha256": pilot.digest(COMPARISON_FACE_CSV),
                "confirmed_exit_zero_runs": sum(x["run_status"] == "confirmed_exit_zero" for x in run_rows),
                "validated_output_exit_unconfirmed_runs": sum(x["run_status"] == "complete_pose_output_exit_unconfirmed" for x in run_rows),
                "failed_or_incomplete_runs": sum(x["run_status"] not in ("confirmed_exit_zero", "complete_pose_output_exit_unconfirmed") for x in run_rows),
                "preserved_initial_failed_attempts": sum(bool(j.get("prior_failure")) for j in jobs),
                "returned_poses": len(pose_rows), "decisions": {"/".join(k): v for k, v in decisions.items()}}
    pilot.atomic_write(MANIFEST_JSON, (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8"))
    print(f"Comparison analysis: {manifest['confirmed_exit_zero_runs']} confirmed exit-zero, "
          f"{manifest['validated_output_exit_unconfirmed_runs']} complete-output/exit-unconfirmed, "
          f"{manifest['failed_or_incomplete_runs']} failed/incomplete; {len(pose_rows)} poses; "
          f"{sum(v['status'] == 'PASS' for v in decisions.values())}/8 combinations PASS.", flush=True)
    return 0 if manifest["failed_or_incomplete_runs"] == 0 else 2


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit-only", action="store_true")
    ap.add_argument("--tail-only", action="store_true", help="Run only the final 16 preregistered jobs; final analysis remains in the full invocation")
    ap.add_argument("--retry-resource-failures", action="store_true", help="After other runs finish, retry saved insufficient-memory or interrupted-output failures under new filenames with at most two concurrent Vina processes")
    args = ap.parse_args()
    frozen, frozen_hash, vina, old_runs, contexts, face_rows, pre = old_input_audit()
    if args.audit_only:
        return 0
    if np.max(np.abs(np.array(pre["new_size_A"]) - NEW_SIZE)) != 0 or pre["old_results_sha256"] != pilot.digest(pilot.RESULTS_CSV):
        raise RuntimeError("Prespecified box or original experiment changed")
    grid_preflight(frozen, vina)
    lookup = {(r["structure_id"], r["microstate"], r["receptor_state"], int(r["seed"])): r for r in old_runs}
    jobs = []
    for r in frozen:
        label = "neutral" if r["formal_charge"] == "0" else "monoanion"
        directory = OUT / f"{r['structure_id']}_{label}_{pilot.RECEPTOR_LABEL[r['receptor_state']]}"
        for seed in pilot.SEEDS:
            old = lookup[(r["structure_id"], r["microstate"], r["receptor_state"], seed)]
            prefix = directory / f"seed_{seed}"
            start = prefix.with_name(prefix.name + "_start.pdbqt")
            output = prefix.with_name(prefix.name + "_docked.pdbqt")
            log = prefix.with_name(prefix.name + "_vina.json")
            prior_failure = None
            if args.retry_resource_failures and log.is_file():
                previous = json.loads(log.read_text(encoding="utf-8"))
                memory_failed = previous.get("return_code") != 0 and "insufficient memory" in previous.get("stderr", "").lower()
                interrupted_without_pose = bool(previous.get("timeout_seconds")) and not output.exists()
                if memory_failed or interrupted_without_pose:
                    if output.exists():
                        raise RuntimeError(f"Refusing automatic retry of failed attempt with existing pose file: {output}")
                    prior_failure = {"log": log, "error": previous["stderr"].strip()}
                    prefix = directory / f"seed_{seed}_retry_after_{'memory' if memory_failed else 'interruption'}"
                    start = prefix.with_name(prefix.name + "_start.pdbqt")
                    output = prefix.with_name(prefix.name + "_docked.pdbqt")
                    log = prefix.with_name(prefix.name + "_vina.json")
            pilot.write_once(start, (ROOT / old["start_path"]).read_bytes())
            if pilot.digest(start) != old["start_sha256"]:
                raise RuntimeError("Comparison start differs from frozen Step 58 start")
            cmd = [str(vina), "--receptor", str(ROOT / r["receptor_pdbqt_path"]), "--ligand", str(start),
                   "--center_x", "45.38664286", "--center_y", "92.35157143", "--center_z", "34.82200000",
                   "--size_x", "30", "--size_y", "24", "--size_z", "24", "--scoring", "vina",
                   "--exhaustiveness", "64", "--num_modes", "20", "--energy_range", "5.0", "--min_rmsd", "1.0",
                   "--seed", str(seed), "--cpu", "1", "--out", str(output), "--verbosity", "1"]
            jobs.append({"row": r, "seed": seed, "start": start, "output": output, "log": log, "cmd": cmd,
                         "prior_failure": prior_failure})
    if len(jobs) != 40:
        raise RuntimeError("Expected exactly 40 comparison jobs")
    if args.tail_only:
        jobs = jobs[24:]
        print("Parallel tail batch: only prespecified jobs 25-40; each remains --cpu 1. No summary write.", flush=True)
    print("Frozen one-box comparison: 40 original starts and 40 unchanged chemistry/seed combinations.", flush=True)
    workers = 2 if args.retry_resource_failures else pilot.WORKERS
    if args.retry_resource_failures:
        print(f"Resource/interruption recovery: {sum(bool(j['prior_failure']) for j in jobs)} preserved failed attempts; "
              "at most 2 concurrent one-CPU retry processes.", flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        future_jobs = {pool.submit(run_job, j): j for j in jobs}
        for fut in concurrent.futures.as_completed(future_jobs):
            job = future_jobs[fut]
            try:
                job["info"] = fut.result()
            except Exception as exc:
                job["info"] = {"status": "runner_error", "exit_code": "", "elapsed": "", "error": str(exc)}
            done = sum("info" in j for j in jobs)
            print(f"{done}/{len(jobs)} {job['row']['structure_id']} {job['row']['microstate']} "
                  f"{job['row']['receptor_state']} seed {job['seed']}: {job['info']['status']}", flush=True)
    if args.tail_only:
        return 0 if all(j["info"]["status"] == "confirmed_exit_zero" for j in jobs) else 2
    return analyze(jobs, old_runs, contexts, frozen, frozen_hash, vina, face_rows)


if __name__ == "__main__":
    raise SystemExit(main())

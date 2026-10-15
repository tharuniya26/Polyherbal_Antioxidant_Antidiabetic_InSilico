"""Execute the frozen Step 66 agnuside pose-only pilot; never tune inputs.

Run from the project Python environment. --preflight-only performs Vina
score_only parsing checks but does not run a docking search. Existing attempts
are read and never overwritten. All RMSDs are unfitted in the receptor frame.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import re
import subprocess
import sys
import time

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree

import run_vina_3L4W_phytochemical_pilot as base

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "data/metadata/vina_3L4W_agnuside_pilot_manifest.csv"
FROZEN_SHA = "1fbba07b7b68d1ffd3417e97d6c1ee56aa18198b97fa55a89ebd9accc8613ac5"
VALIDATOR = ROOT / "scripts/validate_vina_3L4W_agnuside_pilot_inputs.py"
OUT = ROOT / "docking/pilot/vina/3L4W_agnuside_pose_only"
PREFLIGHT = OUT / "vina_pair_parse_preflight.json"
POSE_CSV = OUT / "all_pose_metrics.csv"
RUN_CSV = ROOT / "data/metadata/vina_3L4W_agnuside_pilot_results.csv"
EXPERIMENT = OUT / "experiment_manifest.json"
INTERRUPTED = OUT / "interrupted_attempts.json"
VINA = Path(r"E:\docking tools\vina_1.2.7_win.exe")
SEEDS = (1101, 2202, 3303, 4404, 5505)
STATES = ("Asp542-", "Asp542H(OD2)")
CENTER = np.array([45.38664286, 92.35157143, 34.82200000])
SIZE = np.array([30., 24., 24.])
base.CENTER = CENTER
base.SIZE = SIZE

RUN_FIELDS = base.RUN_FIELDS + ["source_record_1based", "source_conformer_number", "formal_charge",
    "preflight_status", "attempts", "attempt_log_paths", "attempt_log_hashes", "start_reused_across_receptors",
    "rank1_nearest_face", "rank1_nearest_face_atom", "rank1_nearest_water_id", "rank1_nearest_water_A",
    "rank1_water_under_1p8", "rank1_geometric_pass"]
POSE_FIELDS = base.POSE_FIELDS + ["source_record_1based", "formal_charge", "charge_qc",
    "nearest_box_face", "nearest_face_atom_1based", "nearest_water_id", "nearest_water_A",
    "water_under_1p8", "rank1_geometric_pass"]


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def check_inputs() -> tuple[list[dict], np.ndarray]:
    checked = subprocess.run([sys.executable, str(VALIDATOR)], cwd=ROOT, capture_output=True,
                             text=True, timeout=180)
    if checked.returncode:
        raise RuntimeError(f"Step 66 validator failed: {checked.stdout}\n{checked.stderr}")
    if base.digest(FROZEN) != FROZEN_SHA:
        raise RuntimeError("Frozen Step 66 manifest SHA256 changed")
    rows = base.rows(FROZEN)
    expected = {(seed, state) for seed in SEEDS for state in STATES}
    if len(rows) != 10 or {(int(r["seed"]), r["receptor_state"]) for r in rows} != expected:
        raise RuntimeError("Frozen ten-pair seed/receptor allocation differs")
    if base.digest(VINA) != rows[0]["Vina_executable_SHA256"]:
        raise RuntimeError("Vina executable SHA256 differs")
    version = subprocess.run([str(VINA), "--version"], capture_output=True, text=True, timeout=20)
    if version.returncode or "1.2.7" not in version.stdout + version.stderr:
        raise RuntimeError("Vina executable version differs")
    by_seed = {}
    for r in rows:
        if (r["pilot_decision"], r["PubChem_CID"], r["formal_charge"], r["heavy_atoms"],
            r["defined_stereocentres"], r["pdbqt_torsions"], r["box_center_A"], r["box_size_A"],
            r["scoring"], r["exhaustiveness"], r["num_modes"], r["energy_range_kcal_mol"],
            r["min_rmsd_A"], r["cpu"]) != (
            "READY_CONDITIONAL_REFERENCE_POSE_ONLY", "442416", "0", "33", "9", "13",
            "45.38664286;92.35157143;34.82200000", "30;24;24", "vina", "64", "20", "5.0", "1.0", "1"):
            raise RuntimeError("Frozen chemical or search setting changed")
        for path_key, sha_key in (("source_sdf_path", "source_sdf_SHA256"),
                                 ("ligand_sdf_path", "ligand_sdf_SHA256"),
                                 ("ligand_pdbqt_path", "ligand_pdbqt_SHA256"),
                                 ("original_pdbqt_path", "original_pdbqt_SHA256"),
                                 ("start_path", "start_SHA256"),
                                 ("receptor_pdbqt_path", "receptor_pdbqt_SHA256"),
                                 ("receptor_source_pdb_path", "receptor_source_pdb_SHA256"),
                                 ("native_MIG_sdf_path", "native_MIG_sdf_SHA256"),
                                 ("native_MIG_control_path", "native_MIG_control_SHA256")):
            path = ROOT / r[path_key]
            if not path.is_file() or base.digest(path) != r[sha_key]:
                raise RuntimeError(f"Frozen file/hash mismatch: {path}")
        seed = int(r["seed"])
        pair = (r["source_record_1based"], r["start_SHA256"])
        if seed in by_seed and by_seed[seed] != pair:
            raise RuntimeError("Seed-specific conformer/start differs across receptors")
        by_seed[seed] = pair
        if int(r["source_record_1based"]) != dict(zip(SEEDS, (1, 5, 8, 11, 17)))[seed]:
            raise RuntimeError("Frozen conformer assignment changed")
        # Adapt the previous audited Meeko/RDKit pose reader's field names.
        r["prepared_variant_sdf_path"] = r["ligand_sdf_path"]
        r["prepared_variant_pdbqt_path"] = r["ligand_pdbqt_path"]
        r["heavy_atom_map_sdf_to_pdbqt"] = r["heavy_atom_map_sdf_to_pdbqt"]
    if len(set(x[0] for x in by_seed.values())) != 5:
        raise RuntimeError("Selected source conformers are not distinct")
    controls = base.rows(ROOT / rows[0]["native_MIG_control_path"])
    if len(controls) != 10 or any(c["status"] != "completed" or c["exit_code"] != "0"
                                  or float(c["rank1_RMSD_A"]) > 2 for c in controls):
        raise RuntimeError("Enlarged-box native-MIG control differs from 10/10 pass")
    mig = next(iter(Chem.SDMolSupplier(str(ROOT / rows[0]["native_MIG_sdf_path"]), removeHs=False)))
    if mig is None:
        raise RuntimeError("Deposited MIG SDF does not parse")
    mig_xyz = np.array([list(mig.GetConformer().GetAtomPosition(a.GetIdx()))
                        for a in mig.GetAtoms() if a.GetAtomicNum() > 1])
    if mig_xyz.shape != (14, 3) or np.max(np.abs(mig_xyz.mean(axis=0) - CENTER)) > 1e-6:
        raise RuntimeError("MIG reference or box center differs")
    return rows, mig_xyz


def parse_check(r: dict) -> dict:
    cmd = [str(VINA), "--receptor", str(ROOT / r["receptor_pdbqt_path"]),
           "--ligand", str(ROOT / r["start_path"]), "--center_x", "45.38664286",
           "--center_y", "92.35157143", "--center_z", "34.82200000",
           "--size_x", "30", "--size_y", "24", "--size_z", "24", "--scoring", "vina",
           "--score_only", "--cpu", "1"]
    t0 = time.perf_counter()
    result = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=180)
    return {"seed": int(r["seed"]), "receptor_state": r["receptor_state"], "command": cmd,
            "exit_code": result.returncode, "elapsed_seconds": round(time.perf_counter()-t0, 3),
            "start_sha256": r["start_SHA256"], "receptor_sha256": r["receptor_pdbqt_SHA256"],
            "stdout": result.stdout, "stderr": result.stderr,
            "interpretation": "input/grid parse only; score ignored"}


def preflight(rows: list[dict]) -> list[dict]:
    if PREFLIGHT.exists():
        data = json.loads(PREFLIGHT.read_text(encoding="utf-8"))
        checks = data["checks"]
        if data["frozen_manifest_sha256"] != FROZEN_SHA or len(checks) != 10:
            raise RuntimeError("Prior parse-preflight artifact is incomplete")
    else:
        checks = [parse_check(r) for r in rows]
        base.write_once(PREFLIGHT, (json.dumps({"frozen_manifest_sha256": FROZEN_SHA,
                                               "checks": checks}, indent=2) + "\n").encode("utf-8"))
    if {(int(c["seed"]), c["receptor_state"]) for c in checks} != {(int(r["seed"]), r["receptor_state"]) for r in rows}:
        raise RuntimeError("Parse-preflight pair set differs")
    for c in checks:
        r = next(r for r in rows if int(r["seed"]) == int(c["seed"]) and r["receptor_state"] == c["receptor_state"])
        if c["start_sha256"] != r["start_SHA256"] or c["receptor_sha256"] != r["receptor_pdbqt_SHA256"]:
            raise RuntimeError("Parse-preflight input hashes differ")
    if any(c["exit_code"] != 0 for c in checks):
        raise RuntimeError("One or more Vina input/grid parse checks failed; see preflight JSON")
    return checks


def water_context() -> tuple[np.ndarray, list[str]]:
    path = ROOT / "structures/proteins/original/3L4W.pdb"
    xyz, labels = [], []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("HETATM") and line[17:20] in ("HOH", "WAT") and line[12:16].strip() == "O":
            xyz.append([float(line[i:i+8]) for i in (30, 38, 46)])
            labels.append(f"{line[21:22]}:{line[22:26].strip()}:{line[16:17].strip() or '-'}")
    if not xyz:
        raise RuntimeError("Original deposited water oxygens missing")
    return np.array(xyz), labels


def face_water_diagnostic(xyz: np.ndarray, waters: tuple[np.ndarray, list[str]]) -> dict:
    faces = np.stack((xyz[:,0]-(CENTER[0]-15), (CENTER[0]+15)-xyz[:,0],
                      xyz[:,1]-(CENTER[1]-12), (CENTER[1]+12)-xyz[:,1],
                      xyz[:,2]-(CENTER[2]-12), (CENTER[2]+12)-xyz[:,2]), axis=1)
    atom, face = np.unravel_index(np.argmin(faces), faces.shape)
    d, index = cKDTree(waters[0]).query(xyz)
    near = int(np.argmin(d))
    return {"nearest_box_face": ("X-","X+","Y-","Y+","Z-","Z+")[int(face)],
            "nearest_face_atom_1based": int(atom)+1, "nearest_water_id": waters[1][int(index[near])],
            "nearest_water_A": round(float(d[near]), 4), "water_under_1p8": bool(d[near] < 1.8)}


def run_one(job: dict) -> dict:
    output, log = job["output"], job["log"]
    if log.exists() or output.exists():
        if not log.is_file():
            return {"status": "blocked_existing_output_without_log", "exit_code": "", "elapsed": "", "error": "Output exists without attempt log"}
        prior = json.loads(log.read_text(encoding="utf-8"))
        if prior.get("return_code") == 0 and output.is_file():
            return {"status": "completed", "exit_code": 0, "elapsed": prior.get("elapsed_seconds", ""), "error": ""}
        return {"status": "failed_prior_attempt", "exit_code": prior.get("return_code", ""),
                "elapsed": prior.get("elapsed_seconds", ""), "error": prior.get("stderr", "Prior attempt failed")}
    t0 = time.perf_counter()
    try:
        completed = subprocess.run(job["cmd"], cwd=ROOT, capture_output=True, text=True, timeout=5400)
        elapsed = round(time.perf_counter()-t0, 3)
        record = {"command": job["cmd"], "start_sha256": job["row"]["start_SHA256"],
                  "receptor_sha256": job["row"]["receptor_pdbqt_SHA256"],
                  "ligand_sha256": job["row"]["ligand_pdbqt_SHA256"],
                  "return_code": completed.returncode, "elapsed_seconds": elapsed,
                  "stdout": completed.stdout, "stderr": completed.stderr}
        base.write_once(log, (json.dumps(record, indent=2)+"\n").encode("utf-8"))
        return {"status": "completed" if completed.returncode == 0 and output.is_file() else "vina_failed",
                "exit_code": completed.returncode, "elapsed": elapsed, "error": completed.stderr.strip()}
    except subprocess.TimeoutExpired as exc:
        elapsed = round(time.perf_counter()-t0, 3)
        base.write_once(log, (json.dumps({"command": job["cmd"], "start_sha256": job["row"]["start_SHA256"],
                                         "receptor_sha256": job["row"]["receptor_pdbqt_SHA256"],
                                         "timeout_seconds": 5400, "elapsed_seconds": elapsed,
                                         "stdout": str(exc.stdout), "stderr": str(exc.stderr)}, indent=2)+"\n").encode("utf-8"))
        return {"status": "exit_unconfirmed_output" if output.is_file() else "vina_timeout",
                "exit_code": "unconfirmed", "elapsed": elapsed, "error": "Vina timeout; process exit unconfirmed"}


def analyse(rows: list[dict], mig_xyz: np.ndarray, jobs: list[dict]) -> dict:
    contexts = {int(r["seed"]): base.ligand_context(r) for r in rows if r["receptor_state"] == STATES[0]}
    receptors = {state: base.receptor_context(ROOT / next(r["receptor_pdbqt_path"] for r in rows if r["receptor_state"] == state)) for state in STATES}
    waters = water_context()
    interruptions = json.loads(INTERRUPTED.read_text(encoding="utf-8"))["attempts"] if INTERRUPTED.exists() else []
    interrupted_keys = {(int(x["seed"]), x["receptor_state"]) for x in interruptions}
    if len(interrupted_keys) != len(interruptions):
        raise RuntimeError("Duplicate interrupted-attempt record")
    runs, poses, coords = [], [], {}
    for job in jobs:
        r, seed, state, info = job["row"], int(job["row"]["seed"]), job["row"]["receptor_state"], job["info"]
        run = {"structure_id": "LIG09", "compound_name": "Agnuside", "microstate": "neutral_agnuside_reference",
               "receptor_state": state, "seed": seed, "source_record_1based": r["source_record_1based"],
               "source_conformer_number": r["source_conformer_number"], "formal_charge": "0",
               "run_status": info["status"], "receptor_path": r["receptor_pdbqt_path"],
               "receptor_sha256": r["receptor_pdbqt_SHA256"], "ligand_path": r["ligand_pdbqt_path"],
               "ligand_sha256": r["ligand_pdbqt_SHA256"], "source_sdf_path": r["source_sdf_path"],
               "source_sdf_sha256": r["source_sdf_SHA256"], "crystal_reference_path": r["native_MIG_sdf_path"],
               "crystal_reference_sha256": r["native_MIG_sdf_SHA256"], "frozen_manifest_sha256": FROZEN_SHA,
               "vina_executable_sha256": r["Vina_executable_SHA256"], "start_path": r["start_path"],
               "start_sha256": r["start_SHA256"], "start_attempt": r["start_attempt"],
               "start_translation_A": r["start_translation_A"], "start_rotation_quaternion": r["start_quaternion_wxyz"],
               "start_reused_across_receptors": True, "preflight_status": "PASS", "command": subprocess.list2cmdline(job["cmd"]),
               "log_path": rel(job["log"]) if job["log"].is_file() else "", "log_sha256": base.digest(job["log"]) if job["log"].is_file() else "",
               "output_path": rel(job["output"]) if job["output"].is_file() else "", "output_sha256": base.digest(job["output"]) if job["output"].is_file() else "",
               "exit_code": info["exit_code"], "runtime_seconds": info["elapsed"],
               "attempts": (1 if job["log"].is_file() else 0) + int((seed, state) in interrupted_keys),
               "attempt_log_paths": ";".join(x for x in (rel(INTERRUPTED) if (seed, state) in interrupted_keys else "", rel(job["log"]) if job["log"].is_file() else "") if x),
               "attempt_log_hashes": ";".join(x for x in (base.digest(INTERRUPTED) if (seed, state) in interrupted_keys else "", base.digest(job["log"]) if job["log"].is_file() else "") if x),
               "error": info["error"]}
        if job["output"].is_file():
            try:
                got, xyz_by_rank = base.parse_output(job["output"], contexts[seed], receptors[state], mig_xyz, run)
                segments = re.findall(r"(?ms)^MODEL\s+(\d+)\s*\n(.*?)^ENDMDL\s*$", job["output"].read_text(encoding="utf-8"))
                if len(segments) != len(got):
                    raise ValueError("Output mode segment count differs")
                for pose, (_, segment) in zip(got, segments):
                    rank = int(pose["pose_rank"])
                    pose["source_record_1based"] = r["source_record_1based"]
                    pose["formal_charge"] = 0
                    if pose["pose_status"] == "valid":
                        atoms = base.parse_pdbqt_atoms(segment)
                        source = contexts[seed]["atoms"]
                        if set(atoms) != set(source) or any(abs(atoms[s]["charge"]-source[s]["charge"]) > 0.0011 for s in source):
                            pose["pose_status"], pose["error"] = "invalid", "PDBQT atom/charge identity differs"
                        else:
                            pose["charge_qc"] = "pass"
                            pose.update(face_water_diagnostic(xyz_by_rank[rank], waters))
                            coords[(state, seed, rank)] = xyz_by_rank[rank]
                    poses.append(pose)
                run["n_modes_returned"] = len(got)
                if got and got[0]["pose_status"] == "valid":
                    top = got[0]
                    run.update({"rank1_score": top["score"], "rank1_pocket_member": top["pocket_member"],
                                "rank1_severe_clash": top["severe_clash"],
                                "rank1_core_centroid_MIG_distance_A": top["core_centroid_MIG_distance_A"],
                                "rank1_core_MIG_overlap_fraction": top["core_MIG_overlap_fraction"],
                                "rank1_min_box_margin_A": top["min_box_margin_A"],
                                "rank1_catalytic_contacts": top["catalytic_contacts_within_4A"],
                                "rank1_nearest_face": top["nearest_box_face"],
                                "rank1_nearest_face_atom": top["nearest_face_atom_1based"],
                                "rank1_nearest_water_id": top["nearest_water_id"],
                                "rank1_nearest_water_A": top["nearest_water_A"],
                                "rank1_water_under_1p8": top["water_under_1p8"]})
                if any(p["pose_status"] != "valid" for p in got):
                    run["run_status"] = "pose_QC_failed"
                    run["error"] = "; ".join(str(p.get("error", "")) for p in got if p["pose_status"] != "valid")
            except Exception as exc:
                run["run_status"], run["error"] = "pose_parse_failed", str(exc)
        runs.append(run)
    decisions = {}
    for state in STATES:
        available = {seed: coords[(state, seed, 1)] for seed in SEEDS if (state, seed, 1) in coords}
        if available:
            medoid = min(sorted(available), key=lambda seed: (sum(base.symmetry_core_rmsd(available[seed], x, contexts[seed]) for x in available.values()), seed))
        else:
            medoid = None
        passing = 0
        for run in (x for x in runs if x["receptor_state"] == state):
            seed = int(run["seed"])
            if medoid is not None and seed in available:
                rmsd = base.symmetry_core_rmsd(available[seed], available[medoid], contexts[seed])
                run["rank1_core_RMSD_to_medoid_A"] = round(rmsd, 4)
                run["medoid_seed"] = medoid
                passed = (run["run_status"] == "completed" and run["rank1_pocket_member"] is True
                          and run["rank1_severe_clash"] is False and rmsd <= 2.0)
                run["rank1_geometric_pass"] = passed
                passing += passed
            else:
                run["rank1_geometric_pass"] = False
        group = [x for x in runs if x["receptor_state"] == state]
        decision = "BLOCKED" if any(x["run_status"] != "completed" for x in group) else "PASS" if passing >= 4 else "INCONCLUSIVE"
        decisions[state] = {"status": decision, "passing_rank1_of_5": passing, "medoid_seed": medoid}
        for run in group:
            run["combination_status"] = decision
            run["passing_rank1_of_5"] = passing
        for pose in (x for x in poses if x["receptor_state"] == state and x["pose_status"] == "valid"):
            seed, rank = int(pose["seed"]), int(pose["pose_rank"])
            pose["medoid_seed"] = medoid or ""
            if medoid is not None and (state, seed, rank) in coords:
                pose["core_RMSD_to_medoid_A"] = round(base.symmetry_core_rmsd(coords[(state, seed, rank)], available[medoid], contexts[seed]), 4)
            pose["rank1_geometric_pass"] = next(x["rank1_geometric_pass"] for x in group if int(x["seed"]) == seed) if rank == 1 else ""
    base.write_once(RUN_CSV, base.csv_bytes(runs, RUN_FIELDS))
    base.write_once(POSE_CSV, base.csv_bytes(poses, POSE_FIELDS))
    manifest = {"experiment": "Step 67 LIG09 agnuside reference pose-only", "frozen_manifest": rel(FROZEN),
                "frozen_manifest_sha256": FROZEN_SHA, "script": rel(Path(__file__)),
                "script_sha256": base.digest(Path(__file__)), "validator_sha256": base.digest(VALIDATOR),
                "vina_path": str(VINA), "vina_sha256": base.digest(VINA), "vina_version": "1.2.7",
                "box_center_A": CENTER.tolist(), "box_size_A": SIZE.tolist(), "scoring": "vina",
                "exhaustiveness": 64, "num_modes": 20, "energy_range": 5.0, "min_rmsd": 1.0,
                "cpu_per_process": 1, "max_parallel_processes": 2, "seeds": SEEDS,
                "preflight_path": rel(PREFLIGHT), "preflight_sha256": base.digest(PREFLIGHT),
                "interrupted_attempts_path": rel(INTERRUPTED) if INTERRUPTED.exists() else "",
                "interrupted_attempts_sha256": base.digest(INTERRUPTED) if INTERRUPTED.exists() else "",
                "interrupted_attempt_count": len(interruptions),
                "results_csv": rel(RUN_CSV), "results_sha256": base.digest(RUN_CSV),
                "all_pose_csv": rel(POSE_CSV), "all_pose_sha256": base.digest(POSE_CSV),
                "original_pdb_waters_sha256": base.digest(ROOT / "structures/proteins/original/3L4W.pdb"),
                "native_MIG_control_sha256": rows[0]["native_MIG_control_SHA256"],
                "completed_runs": sum(x["run_status"] == "completed" for x in runs),
                "failed_runs": sum(x["run_status"] != "completed" for x in runs),
                "returned_modes": len(poses), "decisions": decisions}
    base.write_once(EXPERIMENT, (json.dumps(manifest, indent=2, sort_keys=True)+"\n").encode("utf-8"))
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    rows, mig_xyz = check_inputs()
    preflight(rows)
    print("Preflight PASS: 10/10 Vina input/grid parse checks, Step 66 hashes, 10/10 native-MIG controls", flush=True)
    if args.preflight_only:
        return 0
    jobs = []
    for r in rows:
        seed, state = int(r["seed"]), r["receptor_state"]
        directory = OUT / f"LIG09_neutral_{base.RECEPTOR_LABEL[state]}"
        directory.mkdir(parents=True, exist_ok=True)
        output = directory / f"seed_{seed}_docked.pdbqt"
        log = directory / f"seed_{seed}_vina.json"
        cmd = [str(VINA), "--receptor", str(ROOT / r["receptor_pdbqt_path"]), "--ligand", str(ROOT / r["start_path"]),
               "--center_x", "45.38664286", "--center_y", "92.35157143", "--center_z", "34.82200000",
               "--size_x", "30", "--size_y", "24", "--size_z", "24", "--scoring", "vina",
               "--exhaustiveness", "64", "--num_modes", "20", "--energy_range", "5.0", "--min_rmsd", "1.0",
               "--seed", str(seed), "--cpu", "1", "--out", str(output), "--verbosity", "1"]
        jobs.append({"row": r, "output": output, "log": log, "cmd": cmd})
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = {pool.submit(run_one, j): j for j in jobs}
        for done in as_completed(futures):
            job = futures[done]
            job["info"] = done.result()
            print(f"{job['row']['receptor_state']} seed {job['row']['seed']}: {job['info']['status']}", flush=True)
    result = analyse(rows, mig_xyz, jobs)
    print(json.dumps({k: result[k] for k in ("completed_runs", "failed_runs", "returned_modes", "decisions")}, indent=2), flush=True)
    return 0 if result["failed_runs"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())

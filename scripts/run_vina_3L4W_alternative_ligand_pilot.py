"""Step 73 frozen neutral rutin/luteolin pose-only pilot (20 searches).

--preflight-only performs the 20 exact Vina score_only input/grid checks. Scores
from preflight are never interpreted. Existing attempts are preserved; this
runner resumes completed jobs and never overwrites prior outputs.
"""
from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import math
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
INPUT = ROOT / "data/metadata/project4_alternative_ligand_input_manifest.csv"
VALIDATOR = ROOT / "scripts/validate_project4_alternative_ligand_inputs.py"
OUT = ROOT / "docking/pilot/vina/3L4W_alternative_ligands_pose_only"
PREFLIGHT = OUT / "vina_pair_parse_preflight.json"
RUN_CSV = ROOT / "data/metadata/vina_3L4W_alternative_ligand_pilot_results.csv"
POSE_CSV = OUT / "all_pose_metrics.csv"
EXPERIMENT = OUT / "experiment_manifest.json"
VINA = Path(r"E:\docking tools\vina_1.2.7_win.exe")
SEEDS = (1101, 2202, 3303, 4404, 5505)
STATES = ("Asp542-", "Asp542H(OD2)")
SID = ("LIG02", "LIG08")
CENTER = np.array([45.38664286, 92.35157143, 34.82200000])
SIZE = np.array([30., 24., 24.])
base.CENTER = CENTER
base.SIZE = SIZE

RUN_FIELDS = base.RUN_FIELDS + ["botanical_association", "source_record_1based", "formal_charge",
    "preflight_status", "attempt_log_paths", "attempt_log_hashes", "attempt_count",
    "rank1_geometric_pass", "rank1_nearest_face", "rank1_nearest_face_atom_1based",
    "rank1_nearest_water_id", "rank1_nearest_water_A", "rank1_water_under_1p8"]
POSE_FIELDS = base.POSE_FIELDS + ["botanical_association", "source_record_1based", "formal_charge",
    "charge_qc", "nearest_box_face", "nearest_face_atom_1based", "nearest_water_id",
    "nearest_water_A", "water_under_1p8", "rank1_geometric_pass"]


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def verify_inputs() -> tuple[list[dict], np.ndarray, str]:
    check = subprocess.run([sys.executable, str(VALIDATOR)], cwd=ROOT, capture_output=True,
                           text=True, timeout=180)
    if check.returncode:
        raise RuntimeError(f"Step 72 validator failed: {check.stdout}\n{check.stderr}")
    rows = [r for r in base.rows(INPUT) if r["decision"] == "PASS" and r["structure_id"] in SID]
    expected = {(sid, state, seed) for sid in SID for state in STATES for seed in SEEDS}
    if len(rows) != 20 or {(r["structure_id"], r["receptor_state"], int(r["start_seed"])) for r in rows} != expected:
        raise RuntimeError("Frozen new-ligand input pair set differs from 20 prespecified jobs")
    control_hash = rows[0]["MIG_control_results_sha256"]
    executable_sha = base.digest(VINA)
    if executable_sha != "e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5":
        raise RuntimeError("Vina executable differs from previously audited binary")
    version = subprocess.run([str(VINA), "--version"], capture_output=True, text=True, timeout=20)
    if version.returncode or "1.2.7" not in version.stdout + version.stderr:
        raise RuntimeError("Vina executable version differs")
    per_seed = {}
    for r in rows:
        if (r["microstate"], r["formal_charge"], r["box_center_A"], r["box_size_A"],
            r["vina_version"], r["scoring"], r["exhaustiveness"], r["num_modes"],
            r["energy_range"], r["min_rmsd_A"], r["cpu"]) != (
            "neutral_reference", "0", "45.38664286;92.35157143;34.82200000", "30;24;24",
            "1.2.7", "vina", "64", "20", "5.0", "1.0", "1"):
            raise RuntimeError("Frozen chemical state or search setting differs")
        for path, sha in (("source_sdf", "source_sdf_sha256"),
                          ("source_pdbqt", "source_pdbqt_sha256"),
                          ("qualified_sdf", "qualified_sdf_sha256"),
                          ("qualified_pdbqt", "qualified_pdbqt_sha256"),
                          ("start_pdbqt", "start_sha256"),
                          ("receptor_pdbqt", "receptor_sha256")):
            p = ROOT / r[path]
            if not p.is_file() or base.digest(p) != r[sha]:
                raise RuntimeError(f"Input/hash mismatch: {p}")
        if r["MIG_reference_sha256"] != base.digest(ROOT / "docking/validation/vina/3L4W_MIG_native_pose.sdf"):
            raise RuntimeError("MIG reference changed")
        if r["MIG_control_results_sha256"] != control_hash or control_hash != base.digest(ROOT / "data/metadata/vina_3L4W_large_box_mig_results.csv"):
            raise RuntimeError("MIG control table changed")
        key = r["structure_id"], int(r["start_seed"])
        pair = (r["selected_record_1based"], r["start_sha256"])
        if key in per_seed and per_seed[key] != pair:
            raise RuntimeError("Start or conformer differs across receptors")
        per_seed[key] = pair
        if r["structure_id"] == "LIG02" and int(r["selected_record_1based"]) != dict(zip(SEEDS, (1, 4, 5, 3, 9)))[int(r["start_seed"])]:
            raise RuntimeError("Rutin source record allocation changed")
        if r["structure_id"] == "LIG08" and r["selected_record_1based"] != "1":
            raise RuntimeError("Luteolin source record allocation changed")
        # Adapt Step 72 field names to the audited Step 58 pose parser.
        r["prepared_variant_sdf_path"] = r["qualified_sdf"]
        r["prepared_variant_pdbqt_path"] = r["qualified_pdbqt"]
        r["heavy_atom_map_sdf_to_pdbqt"] = r["heavy_atom_map"]
        r["rigid_core_pdbqt_serials_1based"] = r["rigid_core_pdbqt_serials"]
    mig_rows = base.rows(ROOT / "data/metadata/vina_3L4W_large_box_mig_results.csv")
    if len(mig_rows) != 10 or any(r["status"] != "completed" or r["exit_code"] != "0"
                                 or float(r["rank1_RMSD_A"]) > 2 for r in mig_rows):
        raise RuntimeError("Previous enlarged-box MIG validation no longer passes")
    mig = next(iter(Chem.SDMolSupplier(str(ROOT / "docking/validation/vina/3L4W_MIG_native_pose.sdf"), removeHs=False)))
    if mig is None:
        raise RuntimeError("Deposited MIG reference failed SDF readback")
    mig_xyz = np.array([list(mig.GetConformer().GetAtomPosition(a.GetIdx()))
                        for a in mig.GetAtoms() if a.GetAtomicNum() > 1])
    if mig_xyz.shape != (14, 3) or not np.allclose(mig_xyz.mean(axis=0), CENTER, atol=1e-6):
        raise RuntimeError("Deposited MIG centroid differs from frozen box")
    return rows, mig_xyz, executable_sha


def score_only_check(r: dict) -> dict:
    cmd = [str(VINA), "--receptor", str(ROOT / r["receptor_pdbqt"]),
           "--ligand", str(ROOT / r["start_pdbqt"]),
           "--center_x", "45.38664286", "--center_y", "92.35157143", "--center_z", "34.82200000",
           "--size_x", "30", "--size_y", "24", "--size_z", "24", "--scoring", "vina",
           "--score_only", "--cpu", "1"]
    start = time.perf_counter()
    try:
        done = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=180)
        return {"compound": r["structure_id"], "receptor_state": r["receptor_state"],
                "seed": int(r["start_seed"]), "command": cmd, "exit_code": done.returncode,
                "elapsed_seconds": round(time.perf_counter() - start, 3),
                "start_sha256": r["start_sha256"], "receptor_sha256": r["receptor_sha256"],
                "stdout": done.stdout, "stderr": done.stderr,
                "interpretation": "Vina input/grid parsing only; score ignored"}
    except subprocess.TimeoutExpired as exc:
        return {"compound": r["structure_id"], "receptor_state": r["receptor_state"],
                "seed": int(r["start_seed"]), "command": cmd, "exit_code": "timeout",
                "elapsed_seconds": round(time.perf_counter() - start, 3),
                "start_sha256": r["start_sha256"], "receptor_sha256": r["receptor_sha256"],
                "stdout": str(exc.stdout), "stderr": str(exc.stderr),
                "interpretation": "input/grid parse unconfirmed"}


def preflight(rows: list[dict]) -> list[dict]:
    if PREFLIGHT.exists():
        stored = json.loads(PREFLIGHT.read_text(encoding="utf-8"))
        checks = stored["checks"]
        if stored["frozen_manifest_sha256"] != base.digest(INPUT):
            raise RuntimeError("Saved preflight is for a different frozen manifest")
    else:
        checks = [score_only_check(r) for r in rows]
        base.write_once(PREFLIGHT, (json.dumps({"frozen_manifest_sha256": base.digest(INPUT),
                                               "checks": checks}, indent=2) + "\n").encode("utf-8"))
    wanted = {(r["structure_id"], r["receptor_state"], int(r["start_seed"])) for r in rows}
    if len(checks) != 20 or {(r["compound"], r["receptor_state"], int(r["seed"])) for r in checks} != wanted:
        raise RuntimeError("Incomplete Vina score_only preflight pair set")
    for r in rows:
        c = next(c for c in checks if (c["compound"], c["receptor_state"], int(c["seed"])) ==
                 (r["structure_id"], r["receptor_state"], int(r["start_seed"])))
        if c["start_sha256"] != r["start_sha256"] or c["receptor_sha256"] != r["receptor_sha256"]:
            raise RuntimeError("Preflight input hashes differ")
    if any(c["exit_code"] != 0 for c in checks):
        raise RuntimeError("At least one Vina input/grid parse check failed; see preflight JSON")
    return checks


def water_context() -> tuple[np.ndarray, list[str]]:
    xyz, labels = [], []
    for line in (ROOT / "structures/proteins/original/3L4W.pdb").read_text(encoding="utf-8").splitlines():
        if line.startswith("HETATM") and line[17:20] in ("HOH", "WAT") and line[12:16].strip() == "O":
            xyz.append([float(line[i:i+8]) for i in (30, 38, 46)])
            labels.append(f"{line[21:22]}:{line[22:26].strip()}:{line[16:17].strip() or '-'}")
    if not xyz:
        raise RuntimeError("No deposited 3L4W water oxygens read")
    return np.array(xyz), labels


def face_water(xyz: np.ndarray, waters: tuple[np.ndarray, list[str]]) -> dict:
    faces = np.stack((xyz[:, 0] - (CENTER[0] - 15), (CENTER[0] + 15) - xyz[:, 0],
                      xyz[:, 1] - (CENTER[1] - 12), (CENTER[1] + 12) - xyz[:, 1],
                      xyz[:, 2] - (CENTER[2] - 12), (CENTER[2] + 12) - xyz[:, 2]), axis=1)
    atom, face = np.unravel_index(np.argmin(faces), faces.shape)
    d, index = cKDTree(waters[0]).query(xyz)
    near = int(np.argmin(d))
    return {"nearest_box_face": ("X-", "X+", "Y-", "Y+", "Z-", "Z+")[int(face)],
            "nearest_face_atom_1based": int(atom) + 1,
            "nearest_water_id": waters[1][int(index[near])],
            "nearest_water_A": round(float(d[near]), 4),
            "water_under_1p8": bool(d[near] < 1.8)}


def run_one(job: dict) -> dict:
    output, folder = job["output"], job["folder"]
    prior_logs = sorted(folder.glob(f"seed_{job['seed']}_attempt_*.json"))
    if prior_logs:
        last = json.loads(prior_logs[-1].read_text(encoding="utf-8"))
        if last.get("exit_code") == 0 and output.is_file():
            return {"status": "completed", "exit_code": 0, "elapsed_seconds": last["elapsed_seconds"],
                    "error": "", "logs": prior_logs}
        return {"status": "failed_preserved", "exit_code": last.get("exit_code", "unconfirmed"),
                "elapsed_seconds": last.get("elapsed_seconds", ""),
                "error": last.get("stderr", "Prior attempt failed"), "logs": prior_logs}
    if output.exists():
        return {"status": "blocked_output_without_log", "exit_code": "", "elapsed_seconds": "",
                "error": "Existing output without attempt log", "logs": []}
    cmd = job["cmd"]
    log = folder / f"seed_{job['seed']}_attempt_1.json"
    started = time.perf_counter()
    try:
        done = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=7200)
        record = {"command": cmd, "input_manifest_sha256": base.digest(INPUT),
                  "receptor_sha256": job["row"]["receptor_sha256"],
                  "ligand_sha256": job["row"]["qualified_pdbqt_sha256"],
                  "start_sha256": job["row"]["start_sha256"],
                  "exit_code": done.returncode, "elapsed_seconds": round(time.perf_counter()-started, 3),
                  "stdout": done.stdout, "stderr": done.stderr,
                  "output_sha256": base.digest(output) if output.is_file() else ""}
        base.write_once(log, (json.dumps(record, indent=2)+"\n").encode("utf-8"))
        return {"status": "completed" if done.returncode == 0 and output.is_file() else "vina_failed",
                "exit_code": done.returncode, "elapsed_seconds": record["elapsed_seconds"],
                "error": done.stderr.strip(), "logs": [log]}
    except subprocess.TimeoutExpired as exc:
        record = {"command": cmd, "input_manifest_sha256": base.digest(INPUT),
                  "receptor_sha256": job["row"]["receptor_sha256"],
                  "ligand_sha256": job["row"]["qualified_pdbqt_sha256"],
                  "start_sha256": job["row"]["start_sha256"],
                  "exit_code": "unconfirmed_timeout", "elapsed_seconds": round(time.perf_counter()-started, 3),
                  "stdout": str(exc.stdout), "stderr": str(exc.stderr),
                  "output_sha256": base.digest(output) if output.is_file() else ""}
        base.write_once(log, (json.dumps(record, indent=2)+"\n").encode("utf-8"))
        return {"status": "exit_unconfirmed_output" if output.is_file() else "timeout",
                "exit_code": "unconfirmed", "elapsed_seconds": record["elapsed_seconds"],
                "error": "Vina timeout; exit status unconfirmed", "logs": [log]}


def analyse(rows: list[dict], mig_xyz: np.ndarray, jobs: list[dict], exe_sha: str) -> dict:
    contexts = {(r["structure_id"], int(r["start_seed"])): base.ligand_context(r)
                for r in rows if r["receptor_state"] == STATES[0]}
    receptors = {state: base.receptor_context(ROOT / next(r["receptor_pdbqt"] for r in rows
                                                            if r["receptor_state"] == state)) for state in STATES}
    waters = water_context()
    runs, poses, coords = [], [], {}
    for job in jobs:
        r, info = job["row"], job["info"]
        sid, state, seed = r["structure_id"], r["receptor_state"], int(r["start_seed"])
        run = {"structure_id": sid, "compound_name": r["compound_name"],
               "botanical_association": r["plant_species"], "microstate": "neutral_reference",
               "formal_charge": 0, "receptor_state": state, "seed": seed,
               "source_record_1based": int(r["selected_record_1based"]),
               "run_status": info["status"], "preflight_status": "PASS",
               "receptor_path": r["receptor_pdbqt"], "receptor_sha256": r["receptor_sha256"],
               "ligand_path": r["qualified_pdbqt"], "ligand_sha256": r["qualified_pdbqt_sha256"],
               "source_sdf_path": r["source_sdf"], "source_sdf_sha256": r["source_sdf_sha256"],
               "crystal_reference_path": "docking/validation/vina/3L4W_MIG_native_pose.sdf",
               "crystal_reference_sha256": r["MIG_reference_sha256"],
               "frozen_manifest_sha256": base.digest(INPUT), "vina_executable_sha256": exe_sha,
               "start_path": r["start_pdbqt"], "start_sha256": r["start_sha256"],
               "command": subprocess.list2cmdline(job["cmd"]),
               "attempt_count": len(info["logs"]),
               "attempt_log_paths": ";".join(rel(x) for x in info["logs"]),
               "attempt_log_hashes": ";".join(base.digest(x) for x in info["logs"]),
               "log_path": rel(info["logs"][-1]) if info["logs"] else "",
               "log_sha256": base.digest(info["logs"][-1]) if info["logs"] else "",
               "output_path": rel(job["output"]) if job["output"].is_file() else "",
               "output_sha256": base.digest(job["output"]) if job["output"].is_file() else "",
               "exit_code": info["exit_code"], "runtime_seconds": info["elapsed_seconds"],
               "error": info["error"]}
        if job["output"].is_file():
            try:
                got, mapped = base.parse_output(job["output"], contexts[sid, seed], receptors[state], mig_xyz, run)
                segments = re.findall(r"(?ms)^MODEL\s+(\d+)\s*\n(.*?)^ENDMDL\s*$",
                                      job["output"].read_text(encoding="utf-8"))
                if len(segments) != len(got):
                    raise ValueError("Output mode count differs from parsed pose table")
                source_atoms = contexts[sid, seed]["atoms"]
                for pose, (_, segment) in zip(got, segments):
                    pose.update({"botanical_association": r["plant_species"],
                                 "source_record_1based": int(r["selected_record_1based"]), "formal_charge": 0})
                    rank = int(pose["pose_rank"])
                    if pose["pose_status"] == "valid":
                        atoms = base.parse_pdbqt_atoms(segment)
                        if set(atoms) != set(source_atoms) or any(abs(atoms[s]["charge"] - source_atoms[s]["charge"]) > 0.0011 for s in source_atoms):
                            pose["pose_status"], pose["error"] = "invalid", "PDBQT atom identity/charge mismatch"
                        else:
                            pose["charge_qc"] = "pass"
                            pose.update(face_water(mapped[rank], waters))
                            coords[sid, state, seed, rank] = mapped[rank]
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
                                "rank1_nearest_face_atom_1based": top["nearest_face_atom_1based"],
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
    for sid in SID:
        for state in STATES:
            available = {seed: coords[sid, state, seed, 1] for seed in SEEDS if (sid, state, seed, 1) in coords}
            if available:
                medoid = min(sorted(available), key=lambda seed: (
                    sum(base.symmetry_core_rmsd(available[seed], xyz, contexts[sid, seed])
                        for xyz in available.values()), seed))
            else:
                medoid = None
            group = [r for r in runs if r["structure_id"] == sid and r["receptor_state"] == state]
            passed = 0
            for run in group:
                seed = int(run["seed"])
                if medoid is not None and seed in available:
                    rmsd = base.symmetry_core_rmsd(available[seed], available[medoid], contexts[sid, seed])
                    run["rank1_core_RMSD_to_medoid_A"] = round(rmsd, 4)
                    run["medoid_seed"] = medoid
                    yes = (run["run_status"] == "completed" and run["rank1_pocket_member"] is True
                           and run["rank1_severe_clash"] is False and rmsd <= 2.0)
                    passed += yes
                    run["rank1_geometric_pass"] = yes
                else:
                    run["rank1_geometric_pass"] = False
            decision = "BLOCKED" if any(r["run_status"] not in ("completed", "exit_unconfirmed_output") for r in group) else "PASS" if passed >= 4 else "INCONCLUSIVE"
            decisions[f"{sid}/{state}"] = {"status": decision, "passing_rank1_of_5": passed, "medoid_seed": medoid}
            for run in group:
                run["combination_status"], run["passing_rank1_of_5"] = decision, passed
            for pose in (p for p in poses if p["structure_id"] == sid and p["receptor_state"] == state and p["pose_status"] == "valid"):
                seed, rank = int(pose["seed"]), int(pose["pose_rank"])
                pose["medoid_seed"] = medoid or ""
                if medoid is not None and (sid, state, seed, rank) in coords:
                    pose["core_RMSD_to_medoid_A"] = round(base.symmetry_core_rmsd(
                        coords[sid, state, seed, rank], available[medoid], contexts[sid, seed]), 4)
                pose["rank1_geometric_pass"] = next(r["rank1_geometric_pass"] for r in group if int(r["seed"]) == seed) if rank == 1 else ""
    base.write_once(RUN_CSV, base.csv_bytes(runs, RUN_FIELDS))
    base.write_once(POSE_CSV, base.csv_bytes(poses, POSE_FIELDS))
    manifest = {"experiment": "Step 73 neutral rutin/luteolin reference pose-only",
                "frozen_input_manifest": rel(INPUT), "frozen_input_manifest_sha256": base.digest(INPUT),
                "script": rel(Path(__file__)), "script_sha256": base.digest(Path(__file__)),
                "validator_sha256": base.digest(VALIDATOR), "vina_path": str(VINA),
                "vina_sha256": exe_sha, "vina_version": "1.2.7", "box_center_A": CENTER.tolist(),
                "box_size_A": SIZE.tolist(), "scoring": "vina", "exhaustiveness": 64,
                "num_modes": 20, "energy_range": 5.0, "min_rmsd": 1.0, "cpu": 1,
                "max_parallel_processes": 2, "seeds": SEEDS,
                "preflight": rel(PREFLIGHT), "preflight_sha256": base.digest(PREFLIGHT),
                "run_csv": rel(RUN_CSV), "run_csv_sha256": base.digest(RUN_CSV),
                "all_pose_csv": rel(POSE_CSV), "all_pose_csv_sha256": base.digest(POSE_CSV),
                "original_pdb_waters_sha256": base.digest(ROOT / "structures/proteins/original/3L4W.pdb"),
                "MIG_control_sha256": rows[0]["MIG_control_results_sha256"],
                "completed_runs": sum(r["run_status"] == "completed" for r in runs),
                "exit_unconfirmed_runs": sum(r["run_status"] == "exit_unconfirmed_output" for r in runs),
                "failed_runs": sum(r["run_status"] not in ("completed", "exit_unconfirmed_output") for r in runs),
                "returned_modes": len(poses), "decisions": decisions}
    base.write_once(EXPERIMENT, (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8"))
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    rows, mig_xyz, exe_sha = verify_inputs()
    preflight(rows)
    print("Preflight PASS: 20/20 exact Vina score_only parses, source hashes and 10/10 prior MIG controls", flush=True)
    if args.preflight_only:
        return 0
    jobs = []
    for r in rows:
        sid, state, seed = r["structure_id"], r["receptor_state"], int(r["start_seed"])
        folder = OUT / f"{sid}_neutral_{base.RECEPTOR_LABEL[state]}"
        folder.mkdir(parents=True, exist_ok=True)
        output = folder / f"seed_{seed}_docked.pdbqt"
        cmd = [str(VINA), "--receptor", str(ROOT / r["receptor_pdbqt"]),
               "--ligand", str(ROOT / r["start_pdbqt"]),
               "--center_x", "45.38664286", "--center_y", "92.35157143", "--center_z", "34.82200000",
               "--size_x", "30", "--size_y", "24", "--size_z", "24", "--scoring", "vina",
               "--exhaustiveness", "64", "--num_modes", "20", "--energy_range", "5.0",
               "--min_rmsd", "1.0", "--seed", str(seed), "--cpu", "1",
               "--out", str(output), "--verbosity", "1"]
        jobs.append({"row": r, "seed": seed, "folder": folder, "output": output, "cmd": cmd})
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = {pool.submit(run_one, job): job for job in jobs}
        for future in as_completed(futures):
            job = futures[future]
            job["info"] = future.result()
            print(f"{job['row']['structure_id']} {job['row']['receptor_state']} seed {job['seed']}: {job['info']['status']}", flush=True)
    result = analyse(rows, mig_xyz, jobs, exe_sha)
    print(json.dumps({k: result[k] for k in ("completed_runs", "failed_runs", "exit_unconfirmed_runs", "returned_modes", "decisions")}, indent=2), flush=True)
    return 0 if result["failed_runs"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())

"""Step 73 frozen neutral rutin/luteolin pose-only pilot (20 searches).

--preflight-only performs the 20 exact Vina score_only input/grid checks. Scores
from preflight are never interpreted. Existing attempts are preserved; this
runner resumes completed jobs and never overwrites prior outputs.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import faulthandler
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import threading
import time
import traceback
import uuid

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
CONTROL = OUT / "_runner_control"
VINA = Path(r"E:\docking tools\vina_1.2.7_win.exe")
VINA_SHA256 = "e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5"
JOB_TIMEOUT_SECONDS = 7200
MAX_WORKERS = 2
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

ATTEMPT_RE = re.compile(r"^seed_(\d+)_attempt_(\d+)\.json$")
LAUNCH_RE = re.compile(r"^seed_(\d+)_attempt_(\d+)_launch\.json$")


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds").replace("+00:00", "Z")


def atomic_write_new(path: Path, content: bytes) -> None:
    """Durably publish a new file without ever replacing an existing path."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
    try:
        with temporary.open("xb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        # Linking is an atomic create-if-absent operation on the same filesystem.
        os.link(temporary, path)
    except FileExistsError as exc:
        raise FileExistsError(f"Refusing to overwrite existing artifact: {path}") from exc
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_json_new(path: Path, record: dict) -> None:
    atomic_write_new(path, (json.dumps(record, indent=2, sort_keys=True) + "\n").encode("utf-8"))


class SessionLog:
    """Thread-safe append-only JSONL session evidence with an fsync per event."""

    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path
        self._stream = path.open("x", encoding="utf-8", newline="\n", buffering=1)
        self._lock = threading.Lock()

    def emit(self, event: str, **fields: object) -> None:
        record = {"timestamp_utc": utc_now(), "event": event, **fields}
        with self._lock:
            self._stream.write(json.dumps(record, sort_keys=True, default=str) + "\n")
            self._stream.flush()
            os.fsync(self._stream.fileno())

    def close(self) -> None:
        with self._lock:
            if not self._stream.closed:
                self._stream.flush()
                os.fsync(self._stream.fileno())
                self._stream.close()


def verify_inputs(run_external_checks: bool = True) -> tuple[list[dict], np.ndarray, str]:
    if run_external_checks:
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
    if executable_sha != VINA_SHA256:
        raise RuntimeError("Vina executable differs from previously audited binary")
    if run_external_checks:
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


def preflight(rows: list[dict], allow_score_only: bool = True) -> list[dict]:
    if PREFLIGHT.exists():
        stored = json.loads(PREFLIGHT.read_text(encoding="utf-8"))
        checks = stored["checks"]
        if stored["frozen_manifest_sha256"] != base.digest(INPUT):
            raise RuntimeError("Saved preflight is for a different frozen manifest")
    else:
        if not allow_score_only:
            raise RuntimeError("Saved preflight is absent; dry-run will not invoke Vina")
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


def build_jobs(rows: list[dict], output_root: Path = OUT,
               control_root: Path = CONTROL) -> list[dict]:
    """Describe the frozen jobs without creating directories or files."""
    jobs = []
    for r in rows:
        sid, state, seed = r["structure_id"], r["receptor_state"], int(r["start_seed"])
        folder = output_root / f"{sid}_neutral_{base.RECEPTOR_LABEL[state]}"
        output = folder / f"seed_{seed}_docked.pdbqt"
        cmd = [str(VINA), "--receptor", str(ROOT / r["receptor_pdbqt"]),
               "--ligand", str(ROOT / r["start_pdbqt"]),
               "--center_x", "45.38664286", "--center_y", "92.35157143", "--center_z", "34.82200000",
               "--size_x", "30", "--size_y", "24", "--size_z", "24", "--scoring", "vina",
               "--exhaustiveness", "64", "--num_modes", "20", "--energy_range", "5.0",
               "--min_rmsd", "1.0", "--seed", str(seed), "--cpu", "1",
               "--out", str(output), "--verbosity", "1"]
        jobs.append({"row": r, "seed": seed, "folder": folder, "output": output,
                     "cmd": cmd, "control": control_root,
                     "manifest_sha256": base.digest(INPUT)})
    return jobs


def _numbered_records(folder: Path, pattern: re.Pattern, seed: int) -> list[tuple[int, Path]]:
    found = []
    if folder.is_dir():
        for path in folder.iterdir():
            match = pattern.fullmatch(path.name)
            if match and int(match.group(1)) == seed:
                found.append((int(match.group(2)), path))
    return sorted(found)


def launch_folder(job: dict) -> Path:
    return job["control"] / "launch_records" / job["folder"].name


def _validate_record_identity(job: dict, record: dict) -> str | None:
    row = job["row"]
    expected = {
        "input_manifest_sha256": job["manifest_sha256"],
        "receptor_sha256": row["receptor_sha256"],
        "ligand_sha256": row["qualified_pdbqt_sha256"],
        "start_sha256": row["start_sha256"],
    }
    if record.get("command") != job["cmd"]:
        return "recorded command differs from the frozen command"
    for key, value in expected.items():
        if record.get(key) != value:
            return f"recorded {key} differs"
    paths = (("receptor_pdbqt", "receptor_sha256"),
             ("qualified_pdbqt", "qualified_pdbqt_sha256"),
             ("start_pdbqt", "start_sha256"))
    for path_key, hash_key in paths:
        path = ROOT / row[path_key]
        if not path.is_file() or base.digest(path) != row[hash_key]:
            return f"current input/hash mismatch: {path}"
    return None


def output_completeness(path: Path, maximum_modes: int) -> tuple[bool, str, int]:
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        return False, f"output read failed: {exc}", 0
    starts = [int(x) for x in re.findall(r"(?m)^MODEL\s+(\d+)\s*$", text)]
    ends = re.findall(r"(?m)^ENDMDL\s*$", text)
    if not starts or len(starts) != len(ends):
        return False, "MODEL/ENDMDL records are absent or unbalanced", len(starts)
    if starts != list(range(1, len(starts) + 1)):
        return False, "MODEL ranks are not consecutive from 1", len(starts)
    if len(starts) > maximum_modes:
        return False, "output contains more modes than the frozen num_modes", len(starts)
    if not re.search(r"(?m)^ENDMDL\s*$", text.rstrip()):
        return False, "output does not end with ENDMDL", len(starts)
    return True, "complete", len(starts)


def _validate_launch_reference(job: dict, record: dict) -> str | None:
    launch_path = record.get("launch_record_path")
    launch_hash = record.get("launch_record_sha256")
    if bool(launch_path) != bool(launch_hash):
        return "launch-record path/hash pair is incomplete"
    if not launch_path:  # Legacy Step 73 records predate launch sidecars.
        return None
    candidate = Path(launch_path)
    if not candidate.is_absolute():
        candidate = ROOT / candidate
    if not candidate.is_file() or base.digest(candidate) != launch_hash:
        return "launch record is absent or its SHA256 differs"
    try:
        launch_record = json.loads(candidate.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        return f"launch record is unreadable: {exc}"
    mismatch = _validate_record_identity(job, launch_record)
    if mismatch:
        return f"launch record {mismatch}"
    if not isinstance(launch_record.get("pid"), int):
        return "launch record has no integer process ID"
    return None


def verify_completed_job(job: dict, attempt_path: Path) -> tuple[bool, str, dict | None]:
    try:
        record = json.loads(attempt_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        return False, f"attempt JSON is unreadable: {exc}", None
    mismatch = _validate_record_identity(job, record)
    if mismatch:
        return False, mismatch, record
    if record.get("exit_code") != 0:
        return False, f"attempt exit status is {record.get('exit_code')!r}, not 0", record
    if not job["output"].is_file():
        return False, "successful attempt record has no output", record
    actual_hash = base.digest(job["output"])
    if not record.get("output_sha256") or record["output_sha256"] != actual_hash:
        return False, "output SHA256 differs from the attempt record", record
    complete, reason, _ = output_completeness(job["output"], int(job["row"]["num_modes"]))
    if not complete:
        return False, reason, record
    launch_mismatch = _validate_launch_reference(job, record)
    if launch_mismatch:
        return False, launch_mismatch, record
    return True, "verified completed output", record


def classify_job(job: dict, orphan_processes: list[dict] | None = None) -> dict:
    seed = job["seed"]
    attempts = _numbered_records(job["folder"], ATTEMPT_RE, seed)
    launches = _numbered_records(launch_folder(job), LAUNCH_RE, seed)
    if attempts and [number for number, _ in attempts] != list(range(1, attempts[-1][0] + 1)):
        return {"classification": "BLOCKED", "reason": "attempt record numbering has a gap",
                "logs": [p for _, p in attempts], "record": None}
    attempt_numbers = {number for number, _ in attempts}
    incomplete = [(number, path) for number, path in launches if number not in attempt_numbers]
    if incomplete:
        number, path = incomplete[-1]
        pid = "unknown"
        try:
            pid = json.loads(path.read_text(encoding="utf-8")).get("pid", "unknown")
        except (OSError, UnicodeError, json.JSONDecodeError):
            pass
        return {"classification": "BLOCKED", "reason":
                f"attempt {number} has a launch record but no final record (pid {pid})",
                "logs": [p for _, p in attempts], "record": None}
    parsed = []
    launch_numbers = {number for number, _ in launches}
    for number, path in attempts:
        try:
            record = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            return {"classification": "BLOCKED", "reason": f"invalid attempt JSON {path.name}: {exc}",
                    "logs": [p for _, p in attempts], "record": None}
        mismatch = _validate_record_identity(job, record)
        if mismatch:
            return {"classification": "BLOCKED", "reason": f"{path.name}: {mismatch}",
                    "logs": [p for _, p in attempts], "record": record}
        launch_mismatch = _validate_launch_reference(job, record)
        if launch_mismatch:
            return {"classification": "BLOCKED", "reason": f"{path.name}: {launch_mismatch}",
                    "logs": [p for _, p in attempts], "record": record}
        if number in launch_numbers and not record.get("launch_record_path"):
            return {"classification": "BLOCKED",
                    "reason": f"{path.name} does not bind its existing launch record",
                    "logs": [p for _, p in attempts], "record": record}
        parsed.append(record)
    if attempts:
        last_path = attempts[-1][1]
        last = parsed[-1]
        if last.get("exit_code") == 0:
            valid, reason, record = verify_completed_job(job, last_path)
            return {"classification": "REUSE" if valid else "BLOCKED", "reason": reason,
                    "logs": [p for _, p in attempts], "record": record}
        if job["output"].exists():
            return {"classification": "BLOCKED",
                    "reason": "output exists after a non-successful finalized attempt",
                    "logs": [p for _, p in attempts], "record": last}
    elif job["output"].exists():
        return {"classification": "BLOCKED", "reason": "output exists without a valid attempt JSON",
                "logs": [], "record": None}
    if orphan_processes:
        pids = ", ".join(str(p["pid"]) for p in orphan_processes)
        return {"classification": "BLOCKED", "reason": f"Vina process already running (pid {pids})",
                "logs": [p for _, p in attempts], "record": parsed[-1] if parsed else None}
    return {"classification": "MISSING", "reason": "no verified completed output",
            "logs": [p for _, p in attempts], "record": parsed[-1] if parsed else None}


def scan_vina_processes() -> list[dict]:
    """Inspect the process table; this does not invoke the Vina executable."""
    if os.name != "nt":
        return []
    command = (f"@(Get-Process -Name '{VINA.stem}' -ErrorAction SilentlyContinue | "
               "Select-Object Id,ProcessName) | ConvertTo-Json -Compress")
    done = subprocess.run(["powershell.exe", "-NoProfile", "-NonInteractive", "-Command", command],
                          capture_output=True, text=True, timeout=20)
    if done.returncode:
        raise RuntimeError(f"Unable to inspect running Vina processes: {done.stderr.strip()}")
    if not done.stdout.strip():
        return []
    records = json.loads(done.stdout)
    if isinstance(records, dict):
        records = [records]
    return [{"image": f"{record['ProcessName']}.exe", "pid": int(record["Id"])}
            for record in records]


def _record_base(job: dict, started_utc: str, started: float) -> dict:
    return {"command": job["cmd"], "input_manifest_sha256": job["manifest_sha256"],
            "receptor_sha256": job["row"]["receptor_sha256"],
            "ligand_sha256": job["row"]["qualified_pdbqt_sha256"],
            "start_sha256": job["row"]["start_sha256"], "started_utc": started_utc,
            "completed_utc": utc_now(), "elapsed_seconds": round(time.perf_counter() - started, 3)}


def _stop_process(process: subprocess.Popen) -> None:
    if process.poll() is not None:
        return
    process.terminate()
    try:
        process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait(timeout=10)


def run_one(job: dict, session: SessionLog | None = None,
            timeout_seconds: float = JOB_TIMEOUT_SECONDS,
            popen_factory=subprocess.Popen) -> dict:
    before = classify_job(job)
    if before["classification"] == "REUSE":
        record = before["record"] or {}
        return {"status": "completed", "exit_code": 0,
                "elapsed_seconds": record.get("elapsed_seconds", ""), "error": "",
                "logs": before["logs"]}
    if before["classification"] != "MISSING":
        return {"status": "blocked", "exit_code": "", "elapsed_seconds": "",
                "error": before["reason"], "logs": before["logs"]}
    prior_numbers = [number for number, _ in _numbered_records(job["folder"], ATTEMPT_RE, job["seed"])]
    attempt_number = max(prior_numbers, default=0) + 1
    job["folder"].mkdir(parents=True, exist_ok=True)
    log = job["folder"] / f"seed_{job['seed']}_attempt_{attempt_number}.json"
    launch = launch_folder(job) / f"seed_{job['seed']}_attempt_{attempt_number}_launch.json"
    started_utc, started = utc_now(), time.perf_counter()
    process = None
    launch_sha = ""
    try:
        process = popen_factory(job["cmd"], cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                text=True)
        launch_record = {**_record_base(job, started_utc, started), "completed_utc": None,
                         "elapsed_seconds": None, "pid": process.pid,
                         "session_log": str(session.path) if session else ""}
        atomic_json_new(launch, launch_record)
        launch_sha = base.digest(launch)
        if session:
            session.emit("job_launched", structure_id=job["row"]["structure_id"],
                         receptor_state=job["row"]["receptor_state"], seed=job["seed"],
                         attempt=attempt_number, pid=process.pid, launch_record=str(launch))
        try:
            stdout, stderr = process.communicate(timeout=timeout_seconds)
            exit_code = process.returncode
            status = "completed" if exit_code == 0 and job["output"].is_file() else "vina_failed"
            error = stderr.strip() or ("Vina returned success without an output" if exit_code == 0 else "")
        except subprocess.TimeoutExpired:
            process.kill()
            stdout, stderr = process.communicate()
            exit_code, status, error = "timeout_killed", "timeout", "Vina timeout; process killed"
        record = {**_record_base(job, started_utc, started), "pid": process.pid,
                  "exit_code": exit_code, "stdout": stdout or "", "stderr": stderr or "",
                  "output_sha256": base.digest(job["output"]) if job["output"].is_file() else "",
                  "launch_record_path": str(launch), "launch_record_sha256": launch_sha}
        atomic_json_new(log, record)
        if status == "completed":
            valid, reason, _ = verify_completed_job(job, log)
            if not valid:
                status, error = "invalid_output", reason
        if session:
            session.emit("job_finished", structure_id=job["row"]["structure_id"],
                         receptor_state=job["row"]["receptor_state"], seed=job["seed"],
                         attempt=attempt_number, pid=process.pid, exit_code=exit_code,
                         status=status, attempt_record=str(log))
        return {"status": status, "exit_code": exit_code,
                "elapsed_seconds": record["elapsed_seconds"], "error": error, "logs": before["logs"] + [log]}
    except BaseException as exc:
        if process is not None:
            try:
                _stop_process(process)
            except Exception:
                pass
        record = {**_record_base(job, started_utc, started),
                  "pid": getattr(process, "pid", None), "exit_code": "worker_exception",
                  "stdout": "", "stderr": "".join(traceback.format_exception(exc)),
                  "output_sha256": base.digest(job["output"]) if job["output"].is_file() else "",
                  "launch_record_path": str(launch) if launch.exists() else "",
                  "launch_record_sha256": launch_sha}
        try:
            atomic_json_new(log, record)
        except FileExistsError:
            pass
        if session:
            session.emit("job_exception", structure_id=job["row"]["structure_id"],
                         receptor_state=job["row"]["receptor_state"], seed=job["seed"],
                         attempt=attempt_number, pid=getattr(process, "pid", None), error=repr(exc))
        if isinstance(exc, (KeyboardInterrupt, SystemExit)):
            raise
        return {"status": "worker_exception", "exit_code": "worker_exception",
                "elapsed_seconds": record["elapsed_seconds"], "error": repr(exc),
                "logs": before["logs"] + ([log] if log.exists() else [])}


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


def classification_counts(classifications: list[dict]) -> dict[str, int]:
    return {name: sum(item["classification"] == name for item in classifications)
            for name in ("REUSE", "MISSING", "BLOCKED")}


def _show_inventory(jobs: list[dict], classifications: list[dict]) -> dict[str, int]:
    for job, item in zip(jobs, classifications):
        print(f"{item['classification']:7} {job['row']['structure_id']} "
              f"{job['row']['receptor_state']} seed {job['seed']}: {item['reason']}", flush=True)
    counts = classification_counts(classifications)
    print(f"SUMMARY REUSE={counts['REUSE']} MISSING={counts['MISSING']} "
          f"BLOCKED={counts['BLOCKED']}", flush=True)
    return counts


def main() -> int:
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--preflight-only", action="store_true")
    modes.add_argument("--dry-run", action="store_true",
                       help="Read-only job inventory; never invokes Vina or writes files")
    args = parser.parse_args()

    if args.dry_run:
        rows, _, _ = verify_inputs(run_external_checks=False)
        preflight(rows, allow_score_only=False)
        jobs = build_jobs(rows)
        orphan_processes = scan_vina_processes()
        classifications = [classify_job(job, orphan_processes) for job in jobs]
        counts = _show_inventory(jobs, classifications)
        return 0 if counts["BLOCKED"] == 0 else 2

    session_id = f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')}_{os.getpid()}_{uuid.uuid4().hex[:8]}"
    session = SessionLog(CONTROL / "sessions" / f"session_{session_id}.jsonl")
    diagnostic_path = CONTROL / "diagnostics" / f"faulthandler_{session_id}.log"
    diagnostic_path.parent.mkdir(parents=True, exist_ok=True)
    diagnostic = diagnostic_path.open("xb")
    faulthandler.enable(file=diagnostic, all_threads=True)
    session.emit("session_started", argv=sys.argv, pid=os.getpid(), faulthandler=str(diagnostic_path))
    try:
        rows, mig_xyz, exe_sha = verify_inputs()
        preflight(rows)
        session.emit("preflight_passed", jobs=len(rows), vina_sha256=exe_sha)
        print("Preflight PASS: 20/20 exact Vina score_only parses, source hashes and 10/10 prior MIG controls", flush=True)
        if args.preflight_only:
            session.emit("session_finished", status="preflight_only")
            return 0

        jobs = build_jobs(rows)
        orphan_processes = scan_vina_processes()
        classifications = [classify_job(job, orphan_processes) for job in jobs]
        counts = _show_inventory(jobs, classifications)
        session.emit("inventory", counts=counts, orphan_processes=orphan_processes)
        if counts["BLOCKED"]:
            session.emit("session_finished", status="blocked_before_launch")
            return 2

        pending = []
        for job, item in zip(jobs, classifications):
            if item["classification"] == "REUSE":
                record = item["record"] or {}
                job["info"] = {"status": "completed", "exit_code": 0,
                               "elapsed_seconds": record.get("elapsed_seconds", ""),
                               "error": "", "logs": item["logs"]}
            else:
                pending.append(job)
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
            futures = {pool.submit(run_one, job, session): job for job in pending}
            for future in as_completed(futures):
                job = futures[future]
                try:
                    job["info"] = future.result()
                except Exception as exc:
                    job["info"] = {"status": "worker_exception", "exit_code": "worker_exception",
                                   "elapsed_seconds": "", "error": repr(exc), "logs": []}
                    session.emit("future_exception", structure_id=job["row"]["structure_id"],
                                 receptor_state=job["row"]["receptor_state"],
                                 seed=job["seed"], error="".join(traceback.format_exception(exc)))
                print(f"{job['row']['structure_id']} {job['row']['receptor_state']} "
                      f"seed {job['seed']}: {job['info']['status']}", flush=True)
        unfinished = [job for job in jobs if job.get("info", {}).get("status") != "completed"]
        if unfinished:
            session.emit("session_finished", status="incomplete", unfinished=len(unfinished))
            print(f"Analysis withheld: {len(unfinished)} jobs are not verified complete", flush=True)
            return 2
        result = analyse(rows, mig_xyz, jobs, exe_sha)
        summary = {k: result[k] for k in ("completed_runs", "failed_runs",
                                          "exit_unconfirmed_runs", "returned_modes", "decisions")}
        print(json.dumps(summary, indent=2), flush=True)
        session.emit("session_finished", status="completed", result=summary)
        return 0 if result["failed_runs"] == 0 else 2
    except BaseException as exc:
        session.emit("session_fatal", error="".join(traceback.format_exception(exc)))
        if isinstance(exc, (KeyboardInterrupt, SystemExit)):
            raise
        print(f"FATAL: {exc}", file=sys.stderr, flush=True)
        return 2
    finally:
        faulthandler.disable()
        diagnostic.flush()
        os.fsync(diagnostic.fileno())
        diagnostic.close()
        session.close()


if __name__ == "__main__":
    raise SystemExit(main())

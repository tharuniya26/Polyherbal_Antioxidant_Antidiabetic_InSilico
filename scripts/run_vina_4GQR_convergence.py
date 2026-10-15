"""Controlled 4GQR/MYC E64 native-ligand convergence test.

Reuses the audited E32 chemistry, atom mapping, crystal reference, box and
randomized starts. Existing E32 and completed E64 artifacts are immutable.
"""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import subprocess
import sys
import time

import run_vina_4GQR_redocking as base


OUT = base.ROOT / "docking/validation/vina/4GQR_convergence_E64"
RESULTS = base.ROOT / "data/metadata/vina_4GQR_convergence_results.csv"
BASELINE = base.RESULTS
EXHAUSTIVENESS = 64


def verify_baseline(heavy, crystal, mapping, ligand_records, automorphisms):
    if not BASELINE.is_file():
        raise RuntimeError("The E32 results CSV is missing")
    with BASELINE.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 5 or [int(r["seed"]) for r in rows] != list(base.SEEDS):
        raise RuntimeError("E32 baseline does not contain exactly five prespecified seeds")
    required = {
        "receptor_path": base.RECEPTOR,
        "ligand_path": base.LIGAND,
        "crystal_reference_path": base.REFERENCE,
    }
    for row in rows:
        seed = int(row["seed"])
        if (row["status"] != "completed" or row["exit_code"] != "0" or row["error"] or
                row["exhaustiveness"] != "32" or row["num_modes_requested"] != "20" or
                row["energy_range_kcal_mol"] != "5.0" or row["min_rmsd_A"] != "1.0" or
                row["cpu"] != "1" or row["vina_version"] != "1.2.7" or
                row["box_center_xyz_A"] != "13.44317391;15.13530435;39.43613043" or
                row["box_size_xyz_A"] != "26;24;28"):
            raise RuntimeError(f"E32 seed {seed} has unexpected status or settings")
        for key, path in required.items():
            if base.ROOT / row[key] != path:
                raise RuntimeError(f"E32 seed {seed} changed {key}")
        for path_key, hash_key in (
            ("receptor_path", "receptor_SHA256"), ("ligand_path", "ligand_SHA256"),
            ("crystal_reference_path", "crystal_reference_SHA256"),
            ("start_path", "start_SHA256"), ("log_path", "log_SHA256"),
            ("output_path", "output_SHA256"),
        ):
            path = base.ROOT / row[path_key]
            if not path.is_file() or base.digest(path) != row[hash_key]:
                raise RuntimeError(f"E32 seed {seed} missing/altered artifact {path_key}")
        log = (base.ROOT / row["log_path"]).read_text()
        if (f"SEED: {seed}" not in log or "EXIT_CODE: 0" not in log or
                f"Performing docking (random seed: {seed}" not in log or
                row["command"] not in log or f"RECEPTOR_SHA256: {row['receptor_SHA256']}" not in log):
            raise RuntimeError(f"E32 seed {seed} log lacks run/chemistry evidence")
        poses = base.parse_output(base.ROOT / row["output_path"], mapping, ligand_records,
                                  heavy, crystal, automorphisms)
        top, best = poses[0], min(poses, key=lambda p: p["heavy_atom_RMSD_A"])
        if (len(poses) != int(row["n_modes_returned"]) or
                abs(top["heavy_atom_RMSD_A"] - float(row["top_RMSD_A"])) > 0.00006 or
                abs(best["heavy_atom_RMSD_A"] - float(row["best_RMSD_A"])) > 0.00006 or
                best["pose_rank"] != int(row["best_RMSD_rank"]) or
                abs(top["score_kcal_mol"] - float(row["top_score_kcal_mol"])) > 0.0001):
            raise RuntimeError(f"E32 seed {seed} CSV/pose reconstruction mismatch")
        metric_file = base.ROOT / "docking/validation/vina/4GQR_diagnostic_redocking" / f"seed_{seed}_pose_metrics.csv"
        if not metric_file.is_file():
            raise RuntimeError(f"E32 seed {seed} pose metric table missing")
        with metric_file.open(newline="", encoding="utf-8") as handle:
            metric_rows = list(csv.DictReader(handle))
        if len(metric_rows) != len(poses):
            raise RuntimeError(f"E32 seed {seed} pose metric count mismatch")
        for metric, pose in zip(metric_rows, poses):
            if (int(metric["pose_rank"]) != pose["pose_rank"] or
                    abs(float(metric["heavy_atom_RMSD_A"]) - pose["heavy_atom_RMSD_A"]) > 1e-8):
                raise RuntimeError(f"E32 seed {seed} pose metrics mismatch")
    return {int(r["seed"]): r for r in rows}


def main():
    source_hashes = base.audit_inputs()
    heavy, crystal, mapping, ligand_records, automorphisms = base.reference_data()
    receptor_qc = base.prepare_dry_receptor()  # verifies identical existing file; never overwrites it
    center, dims, extent = base.box(crystal)
    baseline = verify_baseline(heavy, crystal, mapping, ligand_records, automorphisms)
    if (tuple(f"{x:.8f}" for x in center) != ("13.44317391", "15.13530435", "39.43613043") or
            dims != (26, 24, 28)):
        raise RuntimeError("Deposited-MYC-derived box differs from E32")
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = {
        "experiment": "4GQR/MYC E64 controlled convergence test",
        "created_UTC": datetime.now(timezone.utc).isoformat(),
        "E32_results_path": str(BASELINE.relative_to(base.ROOT)),
        "E32_results_SHA256": base.digest(BASELINE),
        "source_SHA256": source_hashes,
        "receptor_QC": receptor_qc,
        "receptor_SHA256": base.digest(base.RECEPTOR),
        "ligand_SHA256": base.digest(base.LIGAND),
        "crystal_reference_SHA256": base.digest(base.REFERENCE),
        "Vina_executable_SHA256": base.digest(base.VINA),
        "Vina_version": "1.2.7", "Meeko_version": __import__("meeko").__version__,
        "RDKit_version": base.Chem.rdBase.rdkitVersion,
        "seeds": list(base.SEEDS), "exhaustiveness": EXHAUSTIVENESS,
        "num_modes": base.NUM_MODES, "energy_range_kcal_mol": base.ENERGY_RANGE,
        "min_rmsd_A": base.MIN_RMSD, "cpu": base.CPU,
        "box_center_A": center.tolist(), "box_dimensions_A": list(dims),
        "box_site_extent_A": extent.tolist(),
        "randomized_start_rule": "The E32 PCG64/quaternion/translation method is rerun with the same seed; new start bytes must equal each frozen E32 start exactly.",
        "RMSD_rule": "23 MYC heavy atoms; SDF-to-PDBQT audited bridge; minimum over two graph automorphisms; crystal coordinates fixed; no pose fitting.",
        "model_limits": "Dry/no-NAG local model; retained Ca/Cl/PCA; provisional protonation; standard Vina ion scoring and waters not validated; MYC 0.4 occupancy and geometry outliers.",
    }
    manifest_file = OUT / "experiment_manifest.json"
    if manifest_file.exists():
        manifest["created_UTC"] = json.loads(manifest_file.read_text())["created_UTC"]
    base.immutable(manifest_file, (json.dumps(manifest, indent=2, default=str) + "\n").encode())
    print("Verified five E32 rows, logs, hashes, all pose tables and frozen chemistry", flush=True)
    rows, all_poses = [], []
    source = base.LIGAND.read_text()
    for seed in base.SEEDS:
        e32 = baseline[seed]
        start_text, start_rmsd, shift, attempt = base.randomized_start(
            seed, source, center, dims, crystal, mapping, ligand_records, heavy, automorphisms)
        if start_text.encode("ascii") != (base.ROOT / e32["start_path"]).read_bytes():
            raise RuntimeError(f"E64 seed {seed} start differs from frozen E32 start")
        prefix = f"seed_{seed}"
        start = OUT / f"{prefix}_randomized_start.pdbqt"
        output = OUT / f"{prefix}_docked.pdbqt"
        log = OUT / f"{prefix}_vina.log"
        base.immutable(start, start_text.encode("ascii"))
        cmd = [str(base.VINA), "--receptor", str(base.RECEPTOR), "--ligand", str(start),
               "--center_x", f"{center[0]:.8f}", "--center_y", f"{center[1]:.8f}",
               "--center_z", f"{center[2]:.8f}", "--size_x", str(dims[0]),
               "--size_y", str(dims[1]), "--size_z", str(dims[2]),
               "--scoring", "vina", "--exhaustiveness", "64", "--num_modes", "20",
               "--energy_range", "5.0", "--min_rmsd", "1.0", "--seed", str(seed),
               "--cpu", "1", "--out", str(output), "--verbosity", "1"]
        command = subprocess.list2cmdline(cmd)
        if output.exists() and log.exists():
            stdout = log.read_text()
            code_match = re.search(r"(?m)^EXIT_CODE: (\d+)\s*$", stdout)
            elapsed_match = re.search(r"(?m)^ELAPSED_SECONDS: ([\d.]+)\s*$", stdout)
            if command not in stdout or not code_match or not elapsed_match:
                raise RuntimeError(f"Existing E64 seed {seed} log is incomplete")
            code, elapsed = int(code_match.group(1)), float(elapsed_match.group(1))
            print("Reanalysing E64 seed", seed, flush=True)
        elif output.exists() or log.exists():
            raise RuntimeError(f"Partial E64 seed {seed} artifacts require review")
        else:
            print("Docking E64 seed", seed, "start RMSD", round(start_rmsd, 4), flush=True)
            t0 = time.monotonic()
            try:
                run = subprocess.run(cmd, cwd=base.ROOT, capture_output=True, text=True, timeout=3600)
                code, elapsed = run.returncode, time.monotonic() - t0
                stdout = (f"COMMAND: {command}\nSEED: {seed}\nRECEPTOR_SHA256: {base.digest(base.RECEPTOR)}\n"
                          f"LIGAND_SHA256: {base.digest(base.LIGAND)}\nSTART_SHA256: {base.digest(start)}\n"
                          f"EXIT_CODE: {code}\nELAPSED_SECONDS: {elapsed:.3f}\n"
                          f"--- STDOUT ---\n{run.stdout}\n--- STDERR ---\n{run.stderr}\n")
            except subprocess.TimeoutExpired as exc:
                code, elapsed = 124, time.monotonic() - t0
                stdout = (f"COMMAND: {command}\nSEED: {seed}\nRECEPTOR_SHA256: {base.digest(base.RECEPTOR)}\n"
                          f"LIGAND_SHA256: {base.digest(base.LIGAND)}\nSTART_SHA256: {base.digest(start)}\n"
                          f"EXIT_CODE: {code}\nELAPSED_SECONDS: {elapsed:.3f}\nERROR: {exc}\n")
            base.immutable(log, stdout.encode())
        row = {"PDB_ID": "4GQR", "ligand": "MYC neutral CCD reference", "seed": seed,
               "status": "failed", "vina_version": "1.2.7",
               "receptor_model": e32["receptor_model"],
               "receptor_path": e32["receptor_path"], "receptor_SHA256": base.digest(base.RECEPTOR),
               "ligand_path": e32["ligand_path"], "ligand_SHA256": base.digest(base.LIGAND),
               "crystal_reference_path": e32["crystal_reference_path"],
               "crystal_reference_SHA256": base.digest(base.REFERENCE),
               "start_path": str(start.relative_to(base.ROOT)), "start_SHA256": base.digest(start),
               "start_RMSD_A": f"{start_rmsd:.4f}",
               "box_center_xyz_A": ";".join(f"{x:.8f}" for x in center),
               "box_size_xyz_A": ";".join(map(str, dims)),
               "exhaustiveness": 64, "num_modes_requested": base.NUM_MODES,
               "energy_range_kcal_mol": base.ENERGY_RANGE, "min_rmsd_A": base.MIN_RMSD,
               "cpu": base.CPU, "command": command,
               "log_path": str(log.relative_to(base.ROOT)), "log_SHA256": base.digest(log),
               "output_path": str(output.relative_to(base.ROOT)) if output.exists() else "",
               "output_SHA256": base.digest(output) if output.exists() else "",
               "exit_code": code, "elapsed_seconds": f"{elapsed:.3f}",
               "n_modes_returned": "", "top_score_kcal_mol": "", "top_RMSD_A": "",
               "best_RMSD_A": "", "best_RMSD_rank": "", "best_RMSD_score_kcal_mol": "",
               "top_within_2A": "", "any_within_2A": "",
               "E32_top_RMSD_A": e32["top_RMSD_A"],
               "E32_best_RMSD_A": e32["best_RMSD_A"],
               "E32_best_RMSD_rank": e32["best_RMSD_rank"],
               "E32_top_score_kcal_mol": e32["top_score_kcal_mol"],
               "error": ""}
        try:
            if code or not output.is_file() or not output.stat().st_size:
                raise RuntimeError(f"Vina exit {code}; inspect saved log")
            if f"Performing docking (random seed: {seed}" not in stdout:
                raise RuntimeError("Vina log did not confirm requested seed")
            poses = base.parse_output(output, mapping, ligand_records, heavy, crystal, automorphisms)
            if not poses:
                raise RuntimeError("No chemically validated modes")
            top, best = poses[0], min(poses, key=lambda p: p["heavy_atom_RMSD_A"])
            row.update(status="completed", n_modes_returned=len(poses),
                       top_score_kcal_mol=f"{top['score_kcal_mol']:.3f}",
                       top_RMSD_A=f"{top['heavy_atom_RMSD_A']:.4f}",
                       best_RMSD_A=f"{best['heavy_atom_RMSD_A']:.4f}",
                       best_RMSD_rank=best["pose_rank"],
                       best_RMSD_score_kcal_mol=f"{best['score_kcal_mol']:.3f}",
                       top_within_2A=top["heavy_atom_RMSD_A"] <= 2.0,
                       any_within_2A=best["heavy_atom_RMSD_A"] <= 2.0)
            pose_rows = [{"seed": seed, **p, "within_2A": p["heavy_atom_RMSD_A"] <= 2.0}
                         for p in poses]
            all_poses.extend(pose_rows)
            base.immutable(OUT / f"{prefix}_pose_metrics.csv", base.csv_data(
                pose_rows, ["seed", "pose_rank", "score_kcal_mol", "heavy_atom_RMSD_A", "within_2A"]))
            print(seed, "completed; top RMSD", row["top_RMSD_A"],
                  "best", row["best_RMSD_A"], "rank", row["best_RMSD_rank"], flush=True)
        except Exception as exc:
            row["error"] = str(exc)
            print(seed, "FAILED", exc, flush=True)
        rows.append(row)
    base.immutable(RESULTS, base.csv_data(rows, list(rows[0])))
    if all_poses:
        base.immutable(OUT / "all_pose_metrics.csv", base.csv_data(
            all_poses, ["seed", "pose_rank", "score_kcal_mol", "heavy_atom_RMSD_A", "within_2A"]))
    return 0 if len(rows) == 5 and all(r["status"] == "completed" for r in rows) else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print("FATAL E64 PREFLIGHT/ANALYSIS ERROR:", exc, file=sys.stderr)
        sys.exit(2)

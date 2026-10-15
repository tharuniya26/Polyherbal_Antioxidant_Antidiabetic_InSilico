"""Step 60: fixed large-box native MIG pose-recovery comparison (no phytochemicals).

Run with .venv/Scripts/python.exe scripts/run_vina_3L4W_large_box_mig_validation.py.
All pre-existing source files are read-only. Completed new runs are reanalysed,
never repeated or overwritten. The original E64 randomized starts are reproduced
from the audited seeded algorithm and copied byte-for-byte into this experiment.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys
import time

import numpy as np
from rdkit import Chem
from meeko import PDBQTMolecule, RDKitMolCreate

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import run_vina_3L4W_redocking as base

spec = importlib.util.spec_from_file_location("mig_previous_sensitivity", ROOT / "scripts/run_vina_3L4W_sensitivity.py")
sens = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sens)

OUT = ROOT / "docking/validation/vina/3L4W_large_box_MIG"
RESULTS = ROOT / "data/metadata/vina_3L4W_large_box_mig_results.csv"
FROZEN = ROOT / "data/metadata/vina_3L4W_pilot_input_manifest.csv"
OLD_D_MINUS = ROOT / "data/metadata/vina_3L4W_sensitivity_results.csv"
OLD_D_H = ROOT / "data/metadata/vina_3L4W_receptor_sensitivity_results.csv"
RECEPTORS = {
    "Asp542-": base.RECEPTOR,
    "Asp542H(OD2)": ROOT / "structures/proteins/vina_prepared/3L4W_noNAG_D443m_D542H_OD2_H600HID_pH6p5_dry_sensitivity_Meeko0.8.0.pdbqt",
}
SEEDS = base.SEEDS
CENTER = (45.38664286, 92.35157143, 34.82200000)
SIZE = (30, 24, 24)
SETTINGS = {"scoring": "vina", "exhaustiveness": 64, "num_modes": 20,
            "energy_range": 5.0, "min_rmsd": 1.0, "cpu": 1}
FIELDS = [
    "PDB_ID", "receptor_state", "ligand_state", "seed", "status", "exit_code", "error",
    "vina_version", "vina_executable_SHA256", "receptor_path", "receptor_SHA256",
    "ligand_path", "ligand_SHA256", "reference_path", "reference_SHA256",
    "bridge_map_path", "bridge_map_SHA256", "deposited_map_path", "deposited_map_SHA256",
    "old_start_path", "old_start_SHA256", "start_path", "start_SHA256", "start_RMSD_A",
    "box_center_xyz_A", "box_size_xyz_A", "scoring", "exhaustiveness", "num_modes_requested",
    "energy_range_kcal_mol", "min_rmsd_A", "cpu", "command", "log_path", "log_SHA256",
    "output_path", "output_SHA256", "runtime_seconds", "n_modes_returned",
    "rank1_score_kcal_mol", "rank1_RMSD_A", "best_pose_RMSD_A", "best_pose_rank",
    "best_pose_score_kcal_mol", "rank1_within_2A", "any_pose_within_2A",
    "old_box_rank1_RMSD_A", "old_box_best_RMSD_A", "old_box_best_rank",
    "old_box_rank1_score_kcal_mol", "old_box_output_SHA256", "old_box_status",
]
POSE_FIELDS = ["receptor_state", "ligand_state", "seed", "pose_rank", "score_kcal_mol",
               "heavy_atom_RMSD_A", "within_2A", "output_path", "output_SHA256"]


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def audited(path: Path, expected: str) -> None:
    if not path.is_file():
        raise RuntimeError(f"Missing audited input: {path}")
    observed = base.digest(path)
    if observed.lower() != expected.lower():
        raise RuntimeError(f"SHA256 mismatch: {path}: {observed} != {expected}")


def audited_row_file(row: dict, path_col: str, hash_col: str) -> Path:
    path = ROOT / row[path_col]
    audited(path, row[hash_col])
    return path


def preflight():
    frozen = read_csv(FROZEN)
    if len(frozen) != 8 or set(r["receptor_state"] for r in frozen) != set(RECEPTORS):
        raise RuntimeError("Step 57 frozen receptor manifest is incomplete")
    if any(r["decision"] != "PASS" for r in frozen):
        raise RuntimeError("Frozen pilot manifest has a held pair")
    exe_hashes = {r["vina_executable_sha256"] for r in frozen}
    if len(exe_hashes) != 1:
        raise RuntimeError("Frozen Vina executable hashes disagree")
    audited(base.VINA, exe_hashes.pop())
    version = subprocess.run([str(base.VINA), "--version"], capture_output=True, text=True, timeout=30)
    if version.returncode or "1.2.7" not in version.stdout + version.stderr:
        raise RuntimeError(f"Unexpected Vina version: {version.stdout} {version.stderr}")
    for state, path in RECEPTORS.items():
        rows = [r for r in frozen if r["receptor_state"] == state]
        if len({r["receptor_pdbqt_sha256"] for r in rows}) != 1:
            raise RuntimeError(f"Frozen {state} hashes disagree")
        if any((ROOT / r["receptor_pdbqt_path"]).resolve() != path.resolve() for r in rows):
            raise RuntimeError(f"Frozen {state} receptor path differs")
        audited(path, rows[0]["receptor_pdbqt_sha256"])
        for r in rows:
            audited_row_file(r, "receptor_source_pdb_path", "receptor_source_pdb_sha256")
    for r in frozen:
        for path_col, hash_col in (("native_mig_sdf_path", "native_mig_sdf_sha256"),
                                   ("native_mig_pdbqt_path", "native_mig_pdbqt_sha256")):
            audited_row_file(r, path_col, hash_col)
            target = base.REFERENCE if "sdf" in path_col else base.LIGAND
            if (ROOT / r[path_col]).resolve() != target.resolve():
                raise RuntimeError(f"Frozen native MIG path differs: {path_col}")
    base.verify_hashes()  # Checks deposited PDB, both maps, SDF and neutral PDBQT against audits.
    molecule, heavy, xyz, mapping, input_records, automorphisms = base.get_reference()
    signs = base.stereo_signs(heavy, xyz)
    if not np.allclose(xyz.mean(axis=0), CENTER, atol=1e-6):
        raise RuntimeError("Enlarged-box center differs from deposited MIG centroid")
    if len(mapping) != 14 or heavy.GetNumAtoms() != 14 or not automorphisms:
        raise RuntimeError("Neutral MIG map, heavy atoms or symmetry equivalence invalid")
    old_minus = {int(r["seed"]): r for r in read_csv(OLD_D_MINUS) if r["case"] == "neutral_MIG_dry_HID_E64"}
    old_h = {int(r["seed"]): r for r in read_csv(OLD_D_H) if r["ligand_microstate"] == "neutral"}
    if set(old_minus) != set(SEEDS) or set(old_h) != set(SEEDS):
        raise RuntimeError("Original E64 neutral comparison is incomplete")
    old = {"Asp542-": old_minus, "Asp542H(OD2)": old_h}
    starts = {}
    for seed in SEEDS:
        for state in RECEPTORS:
            row = old[state][seed]
            if row["status"] != "completed" or row["exit_code"] != "0" or row["top_within_2A"] != "True":
                raise RuntimeError(f"Original E64 neutral run not confirmed: {state}/{seed}")
            for pc, hc in (("receptor_path", "receptor_SHA256"), ("ligand_path", "ligand_SHA256"),
                           ("start_path", "start_SHA256"), ("output_path", "output_SHA256"),
                           ("log_path", "log_SHA256")):
                audited_row_file(row, pc, hc)
            if base.digest(RECEPTORS[state]) != row["receptor_SHA256"]:
                raise RuntimeError(f"E64 receptor hash differs: {state}")
            if base.digest(base.LIGAND) != row["ligand_SHA256"]:
                raise RuntimeError(f"E64 ligand hash differs: {state}")
            if row["box_size_xyz_A"] != "24;18;18" or int(row["exhaustiveness"]) != 64:
                raise RuntimeError(f"E64 comparison settings differ: {state}/{seed}")
            start = ROOT / row["start_path"]
            if seed in starts and base.digest(start) != base.digest(starts[seed]):
                raise RuntimeError(f"Original neutral starts differ between receptors: {seed}")
            starts[seed] = start
    # Recreate every previous start with the published seeded quaternion/translation
    # algorithm and the old box constraint; require byte identity before reusing it.
    source = base.LIGAND.read_text()
    target_smiles = Chem.MolToSmiles(heavy, isomericSmiles=True)
    for seed, path in starts.items():
        new_text, rmsd, _, _ = sens.random_start(seed, source, input_records, mapping, heavy,
            xyz, np.asarray(CENTER), (24, 18, 18), automorphisms, signs, target_smiles, molecule)
        if new_text.encode("ascii") != path.read_bytes() or rmsd < 5.0:
            raise RuntimeError(f"Seeded native-MIG start not reproduced byte-for-byte: {seed}")
        records = base.atom_records(new_text)
        coords = base.verify_pose_records(records, mapping, input_records, heavy, xyz, signs)
        if np.any(np.abs(coords - np.asarray(CENTER)) > np.asarray(SIZE) / 2 - 0.5):
            raise RuntimeError(f"Seeded start outside enlarged box: {seed}")
        restored = [m for m in RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(new_text)) if m]
        if len(restored) != 1 or Chem.MolToSmiles(Chem.RemoveHs(restored[0]), isomericSmiles=True) != target_smiles:
            raise RuntimeError(f"Seeded start chemistry failed Meeko readback: {seed}")
    return heavy, xyz, mapping, input_records, automorphisms, signs, old, starts


def run_one(state, seed, old_row, old_start, heavy, xyz, mapping, input_records, automorphisms, signs):
    label = "D542m" if state == "Asp542-" else "D542H_OD2"
    prefix = f"{label}_neutral_MIG_seed_{seed}"
    start, docked, log = [OUT / f"{prefix}_{suffix}" for suffix in
                          ("randomized_start.pdbqt", "docked.pdbqt", "vina.log")]
    base.write_immutable(start, old_start.read_bytes())
    command = [str(base.VINA), "--receptor", str(RECEPTORS[state]), "--ligand", str(start),
               "--center_x", f"{CENTER[0]:.8f}", "--center_y", f"{CENTER[1]:.8f}",
               "--center_z", f"{CENTER[2]:.8f}", "--size_x", str(SIZE[0]),
               "--size_y", str(SIZE[1]), "--size_z", str(SIZE[2]),
               "--scoring", "vina", "--exhaustiveness", "64", "--num_modes", "20",
               "--energy_range", "5.0", "--min_rmsd", "1.0", "--seed", str(seed),
               "--cpu", "1", "--out", str(docked), "--verbosity", "1"]
    command_text = subprocess.list2cmdline(command)
    status = "failed"
    error = ""
    exit_code = ""
    runtime = ""
    if docked.exists() or log.exists():
        if not (docked.exists() and log.exists()):
            raise RuntimeError(f"Partial existing output requires review: {prefix}")
        saved = log.read_text(encoding="utf-8")
        if f"COMMAND: {command_text}" not in saved or f"SEED: {seed}" not in saved:
            raise RuntimeError(f"Existing log command/seed mismatch: {prefix}")
        m = re.search(r"(?m)^EXIT_CODE: (-?\d+)$", saved)
        t = re.search(r"(?m)^ELAPSED_SECONDS: ([\d.]+)$", saved)
        if m is None or t is None:
            raise RuntimeError(f"Existing log has no process status: {prefix}")
        exit_code, runtime = m.group(1), t.group(1)
        print(f"{prefix}: reanalysing immutable output", flush=True)
    else:
        print(f"{prefix}: Vina E64 enlarged-box search", flush=True)
        then = time.monotonic()
        run = subprocess.run(command, cwd=str(ROOT), capture_output=True, text=True)
        runtime = f"{time.monotonic() - then:.3f}"
        exit_code = str(run.returncode)
        saved = (f"COMMAND: {command_text}\nSEED: {seed}\n"
                 f"RECEPTOR_SHA256: {base.digest(RECEPTORS[state])}\n"
                 f"LIGAND_SHA256: {base.digest(base.LIGAND)}\n"
                 f"RANDOMIZED_START_SHA256: {base.digest(start)}\n"
                 f"EXIT_CODE: {exit_code}\nELAPSED_SECONDS: {runtime}\n"
                 f"--- STDOUT ---\n{run.stdout}\n--- STDERR ---\n{run.stderr}\n")
        base.write_immutable(log, saved.encode("utf-8"))
    result = {k: "" for k in FIELDS}
    result.update({
        "PDB_ID": "3L4W", "receptor_state": state, "ligand_state": "neutral_MIG",
        "seed": seed, "status": status, "exit_code": exit_code, "vina_version": "1.2.7",
        "vina_executable_SHA256": base.digest(base.VINA),
        "receptor_path": relative(RECEPTORS[state]), "receptor_SHA256": base.digest(RECEPTORS[state]),
        "ligand_path": relative(base.LIGAND), "ligand_SHA256": base.digest(base.LIGAND),
        "reference_path": relative(base.REFERENCE), "reference_SHA256": base.digest(base.REFERENCE),
        "bridge_map_path": relative(base.MAP), "bridge_map_SHA256": base.digest(base.MAP),
        "deposited_map_path": relative(base.DEPOSITED_MAP), "deposited_map_SHA256": base.digest(base.DEPOSITED_MAP),
        "old_start_path": relative(old_start), "old_start_SHA256": base.digest(old_start),
        "start_path": relative(start), "start_SHA256": base.digest(start),
        "start_RMSD_A": old_row["start_RMSD_A"],
        "box_center_xyz_A": ";".join(f"{x:.8f}" for x in CENTER),
        "box_size_xyz_A": ";".join(map(str, SIZE)), "command": command_text,
        "log_path": relative(log), "log_SHA256": base.digest(log),
        "output_path": relative(docked) if docked.exists() else "",
        "output_SHA256": base.digest(docked) if docked.exists() else "",
        "runtime_seconds": runtime,
        "old_box_rank1_RMSD_A": old_row["top_RMSD_A"],
        "old_box_best_RMSD_A": old_row["best_RMSD_A"],
        "old_box_best_rank": old_row.get("best_RMSD_rank", old_row.get("best_RMSD_pose_rank", "")),
        "old_box_rank1_score_kcal_mol": old_row["top_score_kcal_mol"],
        "old_box_output_SHA256": old_row["output_SHA256"],
        "old_box_status": old_row["status"],
    })
    result.update({"scoring": SETTINGS["scoring"], "exhaustiveness": SETTINGS["exhaustiveness"],
                   "num_modes_requested": SETTINGS["num_modes"], "energy_range_kcal_mol": SETTINGS["energy_range"],
                   "min_rmsd_A": SETTINGS["min_rmsd"], "cpu": SETTINGS["cpu"]})
    poses = []
    try:
        if exit_code != "0" or not docked.exists() or not docked.stat().st_size:
            raise RuntimeError(f"Vina exit {exit_code}; inspect preserved log")
        if f"Performing docking (random seed: {seed}" not in saved:
            raise RuntimeError("Vina log lacks seed confirmation")
        metrics = base.parse_output(docked, mapping, input_records, heavy, xyz, automorphisms, signs)
        if not metrics:
            raise RuntimeError("No chemically verified modes")
        best = min(metrics, key=lambda p: p["heavy_atom_RMSD_A"])
        top = metrics[0]
        result.update({"status": "completed", "n_modes_returned": len(metrics),
                       "rank1_score_kcal_mol": f"{top['score_kcal_mol']:.3f}",
                       "rank1_RMSD_A": f"{top['heavy_atom_RMSD_A']:.4f}",
                       "best_pose_RMSD_A": f"{best['heavy_atom_RMSD_A']:.4f}",
                       "best_pose_rank": best["pose_rank"],
                       "best_pose_score_kcal_mol": f"{best['score_kcal_mol']:.3f}",
                       "rank1_within_2A": top["heavy_atom_RMSD_A"] <= 2.0,
                       "any_pose_within_2A": best["heavy_atom_RMSD_A"] <= 2.0})
        poses = [{"receptor_state": state, "ligand_state": "neutral_MIG", "seed": seed,
                  **p, "within_2A": p["heavy_atom_RMSD_A"] <= 2.0,
                  "output_path": relative(docked), "output_SHA256": base.digest(docked)} for p in metrics]
        print(f"{prefix}: rank1 {top['heavy_atom_RMSD_A']:.4f} A; best {best['heavy_atom_RMSD_A']:.4f} A rank {best['pose_rank']}", flush=True)
    except Exception as exc:
        result["error"] = str(exc)
        print(f"{prefix}: FAILED: {exc}", flush=True)
    return result, poses


def main() -> int:
    heavy, xyz, mapping, input_records, automorphisms, signs, old, starts = preflight()
    OUT.mkdir(parents=True, exist_ok=True)
    print("Preflight PASS: two frozen receptor hashes, deposited 14-atom MIG, exact E64 seeded starts", flush=True)
    manifest = {
        "step": 60, "created_UTC": datetime.now(timezone.utc).isoformat(),
        "purpose": "neutral native MIG pose recovery in one preregistered enlarged box",
        "box_center_A": CENTER, "box_size_A": SIZE, "seeds": SEEDS, "settings": SETTINGS,
        "vina_version": "1.2.7", "vina_executable_SHA256": base.digest(base.VINA),
        "frozen_step57_manifest_SHA256": base.digest(FROZEN),
        "receptor_SHA256": {k: base.digest(v) for k, v in RECEPTORS.items()},
        "ligand_SHA256": base.digest(base.LIGAND), "reference_SHA256": base.digest(base.REFERENCE),
        "bridge_map_SHA256": base.digest(base.MAP), "deposited_map_SHA256": base.digest(base.DEPOSITED_MAP),
        "old_E64_result_SHA256": {"Asp542-": base.digest(OLD_D_MINUS), "Asp542H(OD2)": base.digest(OLD_D_H)},
        "old_seeded_start_SHA256": {str(k): base.digest(v) for k, v in starts.items()},
        "starting_pose_method": "Recompute audited NumPy PCG64(seed) quaternion rigid rotation plus translation using original 24x18x18 constraint; require >=5 A crystal-frame RMSD and byte identity with archived E64 start. Copy identical randomized PDBQT per receptor into new experiment. Vina then performs seeded global search.",
        "RMSD_method": "Audited 14-heavy-atom map and graph automorphisms, deposited crystal frame, no independent pose fit; Meeko/RDKit graph, stereo and bond-geometry readback per returned mode",
        "model_limitations": "Both receptors dry/no-NAG; provisional Asp542 and His600 protonation; fixed-water and glycan effects unresolved; Vina score ordering is not affinity validation",
    }
    manifest_path = OUT / "experiment_manifest.json"
    if manifest_path.exists():
        manifest["created_UTC"] = json.loads(manifest_path.read_text(encoding="utf-8"))["created_UTC"]
    base.write_immutable(manifest_path, (json.dumps(manifest, indent=2) + "\n").encode("utf-8"))
    jobs = [(state, seed) for state in RECEPTORS for seed in SEEDS]
    collected = {}
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = {pool.submit(run_one, state, seed, old[state][seed], starts[seed], heavy, xyz,
                               mapping, input_records, automorphisms, signs): (state, seed)
                   for state, seed in jobs}
        for future in as_completed(futures):
            state, seed = futures[future]
            try:
                collected[(state, seed)] = future.result()
            except Exception as exc:
                print(f"{state}/{seed}: process/QC failure: {exc}", file=sys.stderr, flush=True)
                raise
    rows, all_poses = [], []
    for state, seed in jobs:
        row, poses = collected[(state, seed)]
        rows.append(row)
        all_poses.extend(poses)
    base.write_immutable(OUT / "all_pose_metrics.csv", base.csv_bytes(all_poses, POSE_FIELDS))
    base.write_immutable(RESULTS, base.csv_bytes(rows, FIELDS))
    print(f"Completed {sum(r['status']=='completed' for r in rows)}/10; rank1 <=2 A {sum(r['rank1_within_2A'] is True for r in rows)}/10", flush=True)
    return 0 if all(r["status"] == "completed" for r in rows) else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"STEP 60 ERROR: {exc}", file=sys.stderr, flush=True)
        sys.exit(2)

"""Frozen 3L4W phytochemical pose-only pilot (40 Vina runs, no potency analysis).

Run with .venv/Scripts/python.exe from the repository root. Completed outputs
are read back on rerun and are never overwritten. --preflight-only performs no
docking. All coordinates and rigid-core RMSDs use the deposited receptor frame.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import time

import meeko
import numpy as np
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem, rdBase
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
INPUT_MANIFEST = ROOT / "data/metadata/vina_3L4W_pilot_input_manifest.csv"
INPUT_VALIDATOR = ROOT / "scripts/validate_vina_3L4W_pilot_inputs.py"
OUTPUT_DIR = ROOT / "docking/pilot/vina/3L4W_phytochemical_pose_only"
RESULTS_CSV = ROOT / "data/metadata/vina_3L4W_phytochemical_pilot_results.csv"
POSE_CSV = OUTPUT_DIR / "all_pose_metrics.csv"
RUN_MANIFEST = OUTPUT_DIR / "experiment_manifest.json"
SEEDS = (1101, 2202, 3303, 4404, 5505)
VERSION = "1.2.7"
CENTER = np.array([45.38664286, 92.35157143, 34.82200000])
SIZE = np.array([24., 18., 18.])
EXHAUSTIVENESS = 64
NUM_MODES = 20
ENERGY_RANGE = 5.0
MIN_RMSD = 1.0
CPU = 1
WORKERS = 4  # Four independent Vina processes, each retains frozen --cpu 1.
RECEPTOR_LABEL = {"Asp542-": "D542minus", "Asp542H(OD2)": "D542H_OD2"}
TYPE_ELEMENT = {"A": "C", "C": "C", "N": "N", "NA": "N", "OA": "O", "O": "O", "S": "S", "SA": "S"}
RUN_FIELDS = [
    "structure_id", "compound_name", "microstate", "receptor_state", "seed", "run_status", "combination_status",
    "receptor_path", "receptor_sha256", "ligand_path", "ligand_sha256", "source_sdf_path", "source_sdf_sha256", "crystal_reference_path", "crystal_reference_sha256", "frozen_manifest_sha256", "vina_executable_sha256",
    "start_path", "start_sha256", "start_attempt", "start_translation_A", "start_rotation_quaternion", "command", "log_path", "log_sha256", "output_path", "output_sha256", "exit_code", "runtime_seconds", "n_modes_returned",
    "rank1_score", "rank1_pocket_member", "rank1_severe_clash", "rank1_core_RMSD_to_medoid_A", "rank1_core_centroid_MIG_distance_A", "rank1_core_MIG_overlap_fraction", "rank1_min_box_margin_A", "rank1_catalytic_contacts", "medoid_seed", "passing_rank1_of_5", "error",
]
POSE_FIELDS = [
    "structure_id", "compound_name", "microstate", "receptor_state", "seed", "pose_rank", "pose_status", "score", "pocket_member", "all_heavy_inside_box", "core_centroid_MIG_distance_A", "core_MIG_overlap_fraction", "min_box_margin_A", "box_edge_warning", "severe_clash", "receptor_min_heavy_distance_A", "receptor_max_vdw_overlap_A", "receptor_clash_pairs", "internal_min_nonbonded_distance_A", "internal_clash_pairs", "core_RMSD_to_medoid_A", "medoid_seed", "Asp443_min_distance_A", "Asp542_min_distance_A", "His600_min_distance_A", "catalytic_contacts_within_4A", "heavy_atoms", "stereo_geometry", "bond_geometry", "output_path", "output_sha256", "error",
]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def csv_bytes(records: list[dict], fields: list[str]) -> bytes:
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=fields, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(records)
    return out.getvalue().encode("utf-8")


def atomic_write(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_bytes(content)
    tmp.replace(path)


def write_once(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != content:
            raise RuntimeError(f"Refusing to overwrite differing prior artifact: {path}")
        return
    path.write_bytes(content)


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def parse_pdbqt_atoms(text: str) -> dict[int, dict]:
    result: dict[int, dict] = {}
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        serial = int(line[6:11])
        if serial in result:
            raise ValueError(f"Duplicate atom serial {serial}")
        result[serial] = {
            "name": line[12:16].strip(), "resname": line[17:20].strip(),
            "chain": line[21:22], "residue": line[22:26].strip(),
            "xyz": np.array([float(line[i:i + 8]) for i in (30, 38, 46)]),
            "charge": float(line[70:76]), "type": line.split()[-1],
        }
    if not result or any(not np.isfinite(a["xyz"]).all() or not math.isfinite(a["charge"]) for a in result.values()):
        raise ValueError("Missing/nonfinite PDBQT atoms")
    return result


def canonical(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)


def preflight() -> tuple[list[dict[str, str]], str, Path]:
    check = subprocess.run([sys.executable, str(INPUT_VALIDATOR)], cwd=ROOT, capture_output=True, text=True, timeout=120)
    if check.returncode:
        raise RuntimeError(f"Frozen Step 57 validator failed: {check.stdout}\n{check.stderr}")
    manifest_hash = digest(INPUT_MANIFEST)
    manifest = rows(INPUT_MANIFEST)
    if len(manifest) != 8 or len({(r["structure_id"], r["microstate"], r["receptor_state"]) for r in manifest}) != 8:
        raise RuntimeError("Expected exactly eight unique frozen input pairs")
    expected = {(lig, state, receptor) for lig in ("LIG04", "LIG06")
                for state in (("neutral_keto_reference", "xanthone_6_OH_monoanion_hypothesis") if lig == "LIG04" else ("neutral_keto_reference", "flavonol_4prime_OH_monoanion_hypothesis"))
                for receptor in RECEPTOR_LABEL}
    if {(r["structure_id"], r["microstate"], r["receptor_state"]) for r in manifest} != expected:
        raise RuntimeError("Manifest pair identities differ from the eight prespecified combinations")
    executable = Path(r"E:\docking tools\vina_1.2.7_win.exe")
    version = subprocess.run([str(executable), "--version"], capture_output=True, text=True, timeout=20)
    if version.returncode or VERSION not in version.stdout:
        raise RuntimeError("Incorrect AutoDock Vina executable version")
    for r in manifest:
        if r["decision"] != "PASS" or r["vina_version"] != VERSION:
            raise RuntimeError("An input pair is not frozen PASS/Vina 1.2.7")
        if (r["box_center_A"], r["box_size_A"], r["exhaustiveness"], r["num_modes"], r["energy_range_kcal_mol"], r["min_rmsd_A"], r["cpu"], r["seeds"]) != (
            "45.38664286;92.35157143;34.82200000", "24;18;18", "64", "20", "5.0", "1.0", "1", "1101;2202;3303;4404;5505"):
            raise RuntimeError("Frozen Vina settings differ from Step 57")
        if digest(executable).lower() != r["vina_executable_sha256"].lower():
            raise RuntimeError("Vina executable hash mismatch")
        for pathfield, hashfield in (
            ("source_sdf_path", "source_sdf_sha256"), ("source_pdbqt_path", "source_pdbqt_sha256"),
            ("prepared_variant_sdf_path", "prepared_variant_sdf_sha256"), ("prepared_variant_pdbqt_path", "prepared_variant_pdbqt_sha256"),
            ("receptor_source_pdb_path", "receptor_source_pdb_sha256"), ("receptor_pdbqt_path", "receptor_pdbqt_sha256"),
            ("native_mig_sdf_path", "native_mig_sdf_sha256"), ("native_mig_pdbqt_path", "native_mig_pdbqt_sha256"),
            ("native_control_summary_path", "native_control_summary_sha256")):
            if r[pathfield] and digest(ROOT / r[pathfield]).lower() != r[hashfield].lower():
                raise RuntimeError(f"Frozen hash changed: {r[pathfield]}")
    return manifest, manifest_hash, executable


def selected_molecule(r: dict[str, str]) -> Chem.Mol:
    source = ROOT / (r["prepared_variant_sdf_path"] or r["source_sdf_path"])
    records = [m for m in Chem.SDMolSupplier(str(source), removeHs=False, sanitize=True)]
    index = 0 if r["prepared_variant_sdf_path"] else int(r["selected_record_1based"]) - 1
    if index >= len(records) or records[index] is None:
        raise RuntimeError(f"Ligand SDF record unavailable: {source}")
    m = records[index]
    if sum(a.GetAtomicNum() > 1 for a in m.GetAtoms()) != int(r["heavy_atoms"]) or Chem.GetFormalCharge(m) != int(r["formal_charge"]):
        raise RuntimeError("Ligand heavy count or formal charge differs from manifest")
    return m


def ligand_context(r: dict[str, str]) -> dict:
    mol = selected_molecule(r)
    text = (ROOT / r["prepared_variant_pdbqt_path"]).read_text(encoding="utf-8")
    atoms = parse_pdbqt_atoms(text)
    mapping = {int(pair.split(":")[0]) - 1: int(pair.split(":")[1]) for pair in r["heavy_atom_map_sdf_to_pdbqt"].split(";")}
    heavy = Chem.RemoveHs(mol)
    if set(mapping) != set(range(heavy.GetNumAtoms())):
        raise RuntimeError("Frozen ligand heavy atom mapping is incomplete")
    source_xyz = np.array([atoms[mapping[i]]["xyz"] for i in range(heavy.GetNumAtoms())])
    source_chiral = chirality_signs(heavy, source_xyz)
    core = [int(v) - 1 for v in r["rigid_core_sdf_indices_1based"].split(";")]
    if [mapping[i] for i in core] != [int(v) for v in r["rigid_core_pdbqt_serials_1based"].split(";")]:
        raise RuntimeError("Frozen core map differs from the ligand PDBQT map")
    rebuilt = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(text))
    if len(rebuilt) != 1 or rebuilt[0] is None or canonical(rebuilt[0]) != canonical(mol):
        raise RuntimeError("Ligand input PDBQT chemistry changed")
    automorphisms = heavy.GetSubstructMatches(heavy, uniquify=False, useChirality=True, maxMatches=1000)
    if not automorphisms:
        raise RuntimeError("No ligand graph automorphism")
    return {"mol": mol, "heavy": heavy, "text": text, "atoms": atoms, "mapping": mapping, "source_xyz": source_xyz,
            "chiral_signs": source_chiral, "core": core, "automorphisms": automorphisms}


def chirality_signs(heavy: Chem.Mol, xyz: np.ndarray) -> dict[int, int]:
    result = {}
    centers = Chem.FindMolChiralCenters(heavy, includeUnassigned=True)
    for center, assignment in centers:
        if assignment == "?":
            raise ValueError("Undefined ligand stereocentre")
        neighbors = sorted(a.GetIdx() for a in heavy.GetAtomWithIdx(center).GetNeighbors())
        if len(neighbors) < 3:
            raise ValueError("Cannot verify stereocentre with fewer than three heavy neighbors")
        v = xyz[neighbors[:3]] - xyz[center]
        det = float(np.linalg.det(v))
        if abs(det) < 0.2:
            raise ValueError(f"Collapsed stereocentre geometry at atom {center + 1}")
        result[center] = 1 if det > 0 else -1
    return result


def verify_pose_atoms(records: dict[int, dict], ctx: dict) -> np.ndarray:
    original = ctx["atoms"]
    if set(records) != set(original):
        raise ValueError("Pose atom serial set differs from frozen ligand")
    for serial in original:
        if (records[serial]["name"], records[serial]["type"]) != (original[serial]["name"], original[serial]["type"]):
            raise ValueError(f"Pose atom identity/type differs at serial {serial}")
    xyz = np.array([records[ctx["mapping"][i]]["xyz"] for i in range(ctx["heavy"].GetNumAtoms())])
    if chirality_signs(ctx["heavy"], xyz) != ctx["chiral_signs"]:
        raise ValueError("Pose stereocentre geometry differs from input")
    for bond in ctx["heavy"].GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if abs(np.linalg.norm(xyz[i] - xyz[j]) - np.linalg.norm(ctx["source_xyz"][i] - ctx["source_xyz"][j])) > 0.08:
            raise ValueError(f"Heavy-atom bond geometry changed at {i + 1}-{j + 1}")
    return xyz


def random_start(r: dict[str, str], seed: int, ctx: dict) -> tuple[bytes, dict]:
    # The ligand has no deposited native pose. Each start therefore differs
    # from the isolated prepared coordinates by a seeded global rigid rotation
    # and pocket-centered translation. The same start is used in both receptor
    # states of a given ligand/microstate/seed.
    code = 0 if r["formal_charge"] == "0" else 1
    rng = np.random.default_rng(np.random.SeedSequence([seed, int(r["PubChem_CID"]), code]))
    original = ctx["atoms"]
    xyz_all = np.array([original[s]["xyz"] for s in sorted(original)])
    serials = sorted(original)
    input_center = ctx["source_xyz"].mean(axis=0)
    for attempt in range(1, 2001):
        q = rng.normal(size=4)
        q /= np.linalg.norm(q)
        w, x, y, z = q
        rot = np.array([
            [1 - 2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w)],
            [2*(x*y+z*w), 1-2*(x*x+z*z), 2*(y*z-x*w)],
            [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x*x+y*y)],
        ])
        shift = rng.uniform(-3.0, 3.0, size=3)
        if np.linalg.norm(shift) < 2.0:
            continue
        transformed = (xyz_all - input_center) @ rot.T + CENTER + shift
        transformed_map = dict(zip(serials, transformed))
        ligand_xyz = np.array([transformed_map[ctx["mapping"][i]] for i in range(ctx["heavy"].GetNumAtoms())])
        if np.any(np.abs(ligand_xyz - CENTER) > SIZE / 2 - 0.5):
            continue
        lines = []
        for line in ctx["text"].splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                p = transformed_map[int(line[6:11])]
                line = line[:30] + "".join(f"{v:8.3f}" for v in p) + line[54:]
            lines.append(line)
        output = ("\n".join(lines) + "\n").encode("utf-8")
        check = output.decode("utf-8")
        verify_pose_atoms(parse_pdbqt_atoms(check), ctx)
        reconstructed = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(check))
        if len(reconstructed) != 1 or reconstructed[0] is None or canonical(reconstructed[0]) != canonical(ctx["mol"]):
            raise RuntimeError("Random-start Meeko graph readback failed")
        return output, {"attempt": attempt, "translation_A": [round(float(v), 4) for v in shift],
                        "rotation_quaternion": [round(float(v), 8) for v in q]}
    raise RuntimeError(f"Cannot place seeded full ligand inside frozen box: {r['structure_id']} {seed}")


def receptor_context(path: Path) -> dict:
    atoms = parse_pdbqt_atoms(path.read_text(encoding="utf-8"))
    heavy = [a for a in atoms.values() if a["type"] != "HD"]
    if len(heavy) != 6912 or any(a["type"] not in TYPE_ELEMENT for a in heavy):
        raise RuntimeError("Receptor heavy atom count/type differs from frozen local model")
    xyz = np.array([a["xyz"] for a in heavy])
    period = Chem.GetPeriodicTable()
    radii = np.array([period.GetRvdw(TYPE_ELEMENT[a["type"]]) for a in heavy])
    catalytic = {res: np.array([a["xyz"] for a in heavy if a["chain"] == "A" and a["residue"] == res])
                 for res in ("443", "542", "600")}
    if any(len(v) == 0 for v in catalytic.values()):
        raise RuntimeError("Catalytic residue atoms missing")
    return {"xyz": xyz, "radii": radii, "tree": cKDTree(xyz), "catalytic": catalytic}


def score_pose(xyz: np.ndarray, ctx: dict, receptor: dict, mig_xyz: np.ndarray) -> dict:
    core_xyz = xyz[ctx["core"]]
    margin = float(np.min(SIZE / 2 - np.abs(xyz - CENTER)))
    centroid_dist = float(np.linalg.norm(core_xyz.mean(axis=0) - mig_xyz.mean(axis=0)))
    near_mig = np.linalg.norm(core_xyz[:, None, :] - mig_xyz[None, :, :], axis=2).min(axis=1)
    overlap = float(np.mean(near_mig <= 4.0))
    pocket = margin >= -0.0001 and centroid_dist <= 6.0 and overlap >= 0.5
    period = Chem.GetPeriodicTable()
    ligand_radii = np.array([period.GetRvdw(a.GetSymbol()) for a in ctx["heavy"].GetAtoms()])
    clash_pairs = 0
    min_receptor = float("inf")
    max_overlap = float("-inf")
    nearby = receptor["tree"].query_ball_point(xyz, r=5.0)
    for i, indices in enumerate(nearby):
        if not indices:
            continue
        distances = np.linalg.norm(receptor["xyz"][indices] - xyz[i], axis=1)
        overlap_values = ligand_radii[i] + receptor["radii"][indices] - distances
        min_receptor = min(min_receptor, float(distances.min()))
        max_overlap = max(max_overlap, float(overlap_values.max()))
        clash_pairs += int(np.sum((distances < 1.8) | (overlap_values > 0.8)))
    graph_dist = Chem.GetDistanceMatrix(ctx["heavy"])
    internal_min = float("inf")
    internal_clashes = 0
    for i in range(len(xyz)):
        for j in range(i + 1, len(xyz)):
            if graph_dist[i, j] >= 4:
                d = float(np.linalg.norm(xyz[i] - xyz[j]))
                internal_min = min(internal_min, d)
                internal_clashes += d < 1.6
    contacts = {}
    for residue, cat_xyz in receptor["catalytic"].items():
        d = np.linalg.norm(xyz[:, None, :] - cat_xyz[None, :, :], axis=2)
        contacts[residue] = float(d.min())
    return {"pocket_member": pocket, "all_heavy_inside_box": margin >= -0.0001,
            "core_centroid_MIG_distance_A": centroid_dist, "core_MIG_overlap_fraction": overlap,
            "min_box_margin_A": margin, "box_edge_warning": margin < 1.0,
            "severe_clash": bool(clash_pairs or internal_clashes),
            "receptor_min_heavy_distance_A": min_receptor,
            "receptor_max_vdw_overlap_A": max_overlap,
            "receptor_clash_pairs": clash_pairs,
            "internal_min_nonbonded_distance_A": internal_min,
            "internal_clash_pairs": internal_clashes,
            "Asp443_min_distance_A": contacts["443"], "Asp542_min_distance_A": contacts["542"],
            "His600_min_distance_A": contacts["600"],
            "catalytic_contacts_within_4A": ";".join(res for res, d in contacts.items() if d <= 4.0)}


def symmetry_core_rmsd(xyz_a: np.ndarray, xyz_b: np.ndarray, ctx: dict) -> float:
    indices = ctx["core"]
    allowed = [mapping for mapping in ctx["automorphisms"] if {mapping[i] for i in indices} == set(indices)]
    if not allowed:
        raise RuntimeError("No symmetry mapping preserves the frozen rigid core")
    return min(float(np.sqrt(np.mean(np.sum((xyz_a[[mapping[i] for i in indices]] - xyz_b[indices]) ** 2, axis=1))))
               for mapping in allowed)


def parse_output(path: Path, ctx: dict, receptor: dict, mig_xyz: np.ndarray, run: dict) -> tuple[list[dict], dict[int, np.ndarray]]:
    text = path.read_text(encoding="utf-8")
    segments = re.findall(r"(?ms)^MODEL\s+(\d+)\s*\n(.*?)^ENDMDL\s*$", text)
    if not segments:
        raise ValueError("Vina output has no complete MODEL/ENDMDL records")
    rebuilt = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(text))
    if len(rebuilt) != 1 or rebuilt[0] is None or rebuilt[0].GetNumConformers() != len(segments) or canonical(rebuilt[0]) != canonical(ctx["mol"]):
        raise ValueError("Meeko multi-pose graph/stereochemistry readback mismatch")
    output_hash = digest(path)
    pose_rows: list[dict] = []
    coordinates: dict[int, np.ndarray] = {}
    previous_score = float("-inf")
    for rank_s, segment in segments:
        rank = int(rank_s)
        if rank != len(pose_rows) + 1:
            raise ValueError("Vina mode ranks are not sequential")
        m = re.search(r"(?m)^REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)\s+", segment)
        if not m:
            raise ValueError(f"Mode {rank} has no Vina score")
        score = float(m.group(1))
        if score < previous_score - 0.011:
            raise ValueError("Vina mode score order is inconsistent")
        previous_score = score
        row: dict = {"structure_id": run["structure_id"], "compound_name": run["compound_name"],
                     "microstate": run["microstate"], "receptor_state": run["receptor_state"],
                     "seed": run["seed"], "pose_rank": rank, "score": score,
                     "output_path": rel(path), "output_sha256": output_hash, "heavy_atoms": ctx["heavy"].GetNumAtoms()}
        try:
            xyz = verify_pose_atoms(parse_pdbqt_atoms(segment), ctx)
            coordinates[rank] = xyz
            row.update(score_pose(xyz, ctx, receptor, mig_xyz))
            row.update({"pose_status": "valid", "stereo_geometry": "pass", "bond_geometry": "pass"})
        except Exception as exc:
            row.update({"pose_status": "invalid", "stereo_geometry": "failed_or_unchecked", "bond_geometry": "failed_or_unchecked", "error": str(exc)})
        pose_rows.append(row)
    return pose_rows, coordinates


def run_one(job: dict) -> dict:
    output = job["output"]
    log = job["log"]
    if output.exists() or log.exists():
        if not output.is_file() or not log.is_file():
            return {"status": "failed_incomplete_prior_artifacts", "exit_code": "", "elapsed": "", "error": "Existing output/log pair incomplete; preserved without rerun"}
        history = json.loads(log.read_text(encoding="utf-8"))
        if history.get("return_code") == 0:
            return {"status": "reused_completed_artifacts", "exit_code": 0,
                    "elapsed": history.get("elapsed_seconds", ""), "error": ""}
        if history.get("timeout_seconds") and output.is_file():
            return {"status": "timeout_with_saved_output", "exit_code": "unconfirmed",
                    "elapsed": history.get("elapsed_seconds", ""),
                    "error": "Wall-clock timeout after output file was saved; exit code unconfirmed"}
        return {"status": "failed_prior_artifacts", "exit_code": history.get("return_code", ""),
                "elapsed": history.get("elapsed_seconds", ""), "error": "Prior Vina log did not confirm successful completion"}
    begun = time.perf_counter()
    try:
        completed = subprocess.run(job["cmd"], cwd=ROOT, capture_output=True, text=True, timeout=5400)
        elapsed = time.perf_counter() - begun
        log_data = {"command": job["cmd"], "command_display": subprocess.list2cmdline(job["cmd"]),
                    "start_sha256": digest(job["start"]), "receptor_sha256": job["receptor_sha256"],
                    "ligand_sha256": job["ligand_sha256"], "return_code": completed.returncode,
                    "elapsed_seconds": round(elapsed, 3), "stdout": completed.stdout, "stderr": completed.stderr}
        write_once(log, (json.dumps(log_data, indent=2) + "\n").encode("utf-8"))
        return {"status": "completed" if completed.returncode == 0 and output.is_file() else "vina_failed",
                "exit_code": completed.returncode, "elapsed": round(elapsed, 3), "error": completed.stderr.strip()}
    except subprocess.TimeoutExpired as exc:
        elapsed = time.perf_counter() - begun
        write_once(log, (json.dumps({"command": job["cmd"], "timeout_seconds": 5400,
                                     "elapsed_seconds": round(elapsed, 3), "stdout": str(exc.stdout),
                                     "stderr": str(exc.stderr)}, indent=2) + "\n").encode("utf-8"))
        return {"status": "vina_timeout", "exit_code": "", "elapsed": round(elapsed, 3), "error": "Vina exceeded 5400 s"}


def main() -> int:
    cli = argparse.ArgumentParser()
    cli.add_argument("--preflight-only", action="store_true")
    args = cli.parse_args()
    manifest, manifest_hash, executable = preflight()
    print(f"Preflight PASS: eight input pairs; manifest SHA256 {manifest_hash}", flush=True)
    if args.preflight_only:
        return 0
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    mig = next(iter(Chem.SDMolSupplier(str(ROOT / manifest[0]["native_mig_sdf_path"]), removeHs=False)))
    mig_xyz = np.array([mig.GetConformer().GetAtomPosition(a.GetIdx()) for a in mig.GetAtoms() if a.GetAtomicNum() > 1], dtype=float)
    if len(mig_xyz) != 14 or np.max(np.abs(mig_xyz.mean(axis=0) - CENTER)) > 1e-7:
        raise RuntimeError("Deposited MIG pocket reference changed")
    contexts = {(r["structure_id"], r["microstate"]): ligand_context(r) for r in manifest}
    receptors = {r["receptor_state"]: receptor_context(ROOT / r["receptor_pdbqt_path"]) for r in manifest}
    jobs = []
    for r in manifest:
        ctx = contexts[(r["structure_id"], r["microstate"])]
        state_name = "neutral" if r["formal_charge"] == "0" else "monoanion"
        basename = f"{r['structure_id']}_{state_name}_{RECEPTOR_LABEL[r['receptor_state']]}"
        pair_dir = OUTPUT_DIR / basename
        pair_dir.mkdir(parents=True, exist_ok=True)
        for seed in SEEDS:
            prefix = pair_dir / f"seed_{seed}"
            start = prefix.with_name(prefix.name + "_start.pdbqt")
            output = prefix.with_name(prefix.name + "_docked.pdbqt")
            log = prefix.with_name(prefix.name + "_vina.json")
            start_data, random_meta = random_start(r, seed, ctx)
            write_once(start, start_data)
            cmd = [str(executable), "--receptor", str(ROOT / r["receptor_pdbqt_path"]),
                   "--ligand", str(start), "--center_x", "45.38664286", "--center_y", "92.35157143",
                   "--center_z", "34.82200000", "--size_x", "24", "--size_y", "18", "--size_z", "18",
                   "--scoring", "vina", "--exhaustiveness", "64", "--num_modes", "20",
                   "--energy_range", "5.0", "--min_rmsd", "1.0", "--seed", str(seed),
                   "--cpu", "1", "--out", str(output), "--verbosity", "1"]
            jobs.append({"row": r, "ctx": ctx, "seed": seed, "start": start, "output": output, "log": log,
                         "receptor_sha256": r["receptor_pdbqt_sha256"], "ligand_sha256": r["prepared_variant_pdbqt_sha256"],
                         "cmd": cmd, "random_meta": random_meta})
    if len(jobs) != 40:
        raise RuntimeError("Expected 40 prespecified jobs")
    print("Prepared 40 seeded starts without changing any frozen input", flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futures = {pool.submit(run_one, job): job for job in jobs}
        for future in concurrent.futures.as_completed(futures):
            job = futures[future]
            try:
                job["run_info"] = future.result()
            except Exception as exc:
                job["run_info"] = {"status": "runner_error", "exit_code": "", "elapsed": "", "error": str(exc)}
            print(f"{len([j for j in jobs if 'run_info' in j])}/40 {job['row']['structure_id']} "
                  f"{job['row']['microstate']} {job['row']['receptor_state']} seed {job['seed']}: "
                  f"{job['run_info']['status']}", flush=True)
    run_rows: list[dict] = []
    pose_rows: list[dict] = []
    rank1_coords: dict[tuple[str, str, str], dict[int, np.ndarray]] = {}
    for job in jobs:
        r, ctx, seed = job["row"], job["ctx"], job["seed"]
        info = job["run_info"]
        key = (r["structure_id"], r["microstate"], r["receptor_state"])
        row: dict = {"structure_id": r["structure_id"], "compound_name": r["compound_name"],
                     "microstate": r["microstate"], "receptor_state": r["receptor_state"], "seed": seed,
                     "run_status": info["status"], "receptor_path": r["receptor_pdbqt_path"],
                     "receptor_sha256": r["receptor_pdbqt_sha256"], "ligand_path": r["prepared_variant_pdbqt_path"],
                     "ligand_sha256": r["prepared_variant_pdbqt_sha256"], "source_sdf_path": r["source_sdf_path"],
                     "source_sdf_sha256": r["source_sdf_sha256"], "crystal_reference_path": r["native_mig_sdf_path"],
                     "crystal_reference_sha256": r["native_mig_sdf_sha256"], "frozen_manifest_sha256": manifest_hash,
                     "vina_executable_sha256": r["vina_executable_sha256"], "start_path": rel(job["start"]),
                     "start_sha256": digest(job["start"]), "start_attempt": job["random_meta"]["attempt"],
                     "start_translation_A": ";".join(map(str, job["random_meta"]["translation_A"])),
                     "start_rotation_quaternion": ";".join(map(str, job["random_meta"]["rotation_quaternion"])),
                     "command": subprocess.list2cmdline(job["cmd"]), "log_path": rel(job["log"]) if job["log"].exists() else "",
                     "log_sha256": digest(job["log"]) if job["log"].exists() else "",
                     "output_path": rel(job["output"]) if job["output"].exists() else "",
                     "output_sha256": digest(job["output"]) if job["output"].exists() else "",
                     "exit_code": info["exit_code"], "runtime_seconds": info["elapsed"], "error": info["error"]}
        if job["output"].is_file() and job["log"].is_file():
            try:
                parsed, coords = parse_output(job["output"], ctx, receptors[r["receptor_state"]], mig_xyz, row)
                row["n_modes_returned"] = len(parsed)
                pose_rows.extend(parsed)
                top = parsed[0]
                if top["pose_status"] == "valid":
                    rank1_coords.setdefault(key, {})[seed] = coords[1]
                    row.update({"rank1_score": top["score"], "rank1_pocket_member": top["pocket_member"],
                                "rank1_severe_clash": top["severe_clash"],
                                "rank1_core_centroid_MIG_distance_A": top["core_centroid_MIG_distance_A"],
                                "rank1_core_MIG_overlap_fraction": top["core_MIG_overlap_fraction"],
                                "rank1_min_box_margin_A": top["min_box_margin_A"],
                                "rank1_catalytic_contacts": top["catalytic_contacts_within_4A"]})
                if any(p["pose_status"] != "valid" for p in parsed):
                    row["run_status"] = "pose_QC_failed"
                    row["error"] = "; ".join(p.get("error", "") for p in parsed if p["pose_status"] != "valid")
                elif info["status"] in ("completed", "reused_completed_artifacts"):
                    row["run_status"] = "completed"
                elif info["status"] in ("vina_timeout", "timeout_with_saved_output"):
                    row["run_status"] = "complete_pose_output_exit_unconfirmed"
                    row["error"] = "Vina wall-clock timeout; complete multi-mode output passed Meeko/RDKit pose QC, but process exit code is unconfirmed"
            except Exception as exc:
                row["run_status"] = "pose_parse_failed"
                row["error"] = str(exc)
        run_rows.append(row)
    # One medoid and one prospective 4/5 decision for each frozen combination.
    decisions = {}
    for r in manifest:
        key = (r["structure_id"], r["microstate"], r["receptor_state"])
        if key in decisions:
            continue
        ctx = contexts[(r["structure_id"], r["microstate"])]
        available = rank1_coords.get(key, {})
        medoid_seed = None
        if available:
            medoid_seed = min(sorted(available), key=lambda s: (sum(symmetry_core_rmsd(available[s], other, ctx) for other in available.values()), s))
        passes = 0
        for row in (x for x in run_rows if (x["structure_id"], x["microstate"], x["receptor_state"]) == key):
            if medoid_seed is not None and int(row["seed"]) in available:
                rmsd = symmetry_core_rmsd(available[int(row["seed"])], available[medoid_seed], ctx)
                row["rank1_core_RMSD_to_medoid_A"] = f"{rmsd:.4f}"
                row["medoid_seed"] = medoid_seed
                if row["run_status"] in ("completed", "complete_pose_output_exit_unconfirmed") and row["rank1_pocket_member"] is True and row["rank1_severe_clash"] is False and rmsd <= 2.0:
                    passes += 1
        combo_rows = [x for x in run_rows if (x["structure_id"], x["microstate"], x["receptor_state"]) == key]
        usable = ("completed", "complete_pose_output_exit_unconfirmed")
        if len(available) < 4 or any(x["run_status"] not in usable for x in combo_rows):
            status = "BLOCKED"
        else:
            status = "PASS" if passes >= 4 else "INCONCLUSIVE"
        decisions[key] = (status, passes, medoid_seed)
        for row in combo_rows:
            row["combination_status"] = status
            row["passing_rank1_of_5"] = passes
        for pose in (x for x in pose_rows if (x["structure_id"], x["microstate"], x["receptor_state"]) == key):
            pose["medoid_seed"] = medoid_seed or ""
            if pose["pose_status"] == "valid" and medoid_seed is not None:
                pose_context = next((j for j in jobs if j["row"]["structure_id"] == key[0] and j["row"]["microstate"] == key[1]
                                     and j["row"]["receptor_state"] == key[2] and j["seed"] == int(pose["seed"])), None)
                # Reparse only the coordinate block needed for the per-pose medoid metric.
                if pose_context:
                    text = pose_context["output"].read_text(encoding="utf-8")
                    seg = re.findall(r"(?ms)^MODEL\s+(\d+)\s*\n(.*?)^ENDMDL\s*$", text)[int(pose["pose_rank"]) - 1][1]
                    xyz = verify_pose_atoms(parse_pdbqt_atoms(seg), ctx)
                    pose["core_RMSD_to_medoid_A"] = f"{symmetry_core_rmsd(xyz, available[medoid_seed], ctx):.4f}"
    atomic_write(RESULTS_CSV, csv_bytes(run_rows, RUN_FIELDS))
    atomic_write(POSE_CSV, csv_bytes(pose_rows, POSE_FIELDS))
    run_manifest = {
        "experiment": "3L4W mangiferin/quercetin pose-only pilot", "created_utc": datetime.now(timezone.utc).isoformat(),
        "script_path": rel(Path(__file__)), "script_sha256": digest(Path(__file__)),
        "frozen_input_manifest_path": rel(INPUT_MANIFEST), "frozen_input_manifest_sha256": manifest_hash,
        "input_validator_path": rel(INPUT_VALIDATOR), "input_validator_sha256": digest(INPUT_VALIDATOR),
        "vina_executable": str(executable), "vina_executable_sha256": digest(executable), "vina_version": VERSION,
        "python_version": sys.version.split()[0], "rdkit_version": rdBase.rdkitVersion, "meeko_version": meeko.__version__,
        "settings": {"scoring": "vina", "box_center_A": CENTER.tolist(), "box_size_A": SIZE.tolist(),
                     "exhaustiveness": EXHAUSTIVENESS, "num_modes": NUM_MODES, "energy_range": ENERGY_RANGE,
                     "min_rmsd": MIN_RMSD, "cpu_per_run": CPU, "concurrent_processes": WORKERS, "seeds": list(SEEDS)},
        "start_method": "NumPy PCG64 SeedSequence(seed, PubChem CID, charge code); independent random unit quaternion; translation in [-3,3] A with norm >=2 A; all heavy atoms inside box with >=0.5 A margin; same ligand start in both receptor states",
        "pose_rules": {"pocket": "all heavy inside box; core centroid <=6 A from MIG centroid; >=50% core atoms within 4 A of any MIG heavy atom",
                       "clash": "receptor-ligand heavy distance <1.8 A OR van der Waals overlap >0.8 A; internal >=4-bond heavy distance <1.6 A",
                       "convergence": ">=4/5 rank1 per combination pocket+no severe clash+<=2 A unfitted core RMSD to within-state medoid"},
        "source_input_sha256": {r["structure_id"] + "/" + r["microstate"] + "/" + r["receptor_state"]:
                                {"ligand": r["prepared_variant_pdbqt_sha256"], "receptor": r["receptor_pdbqt_sha256"],
                                 "source_sdf": r["source_sdf_sha256"], "MIG_reference": r["native_mig_sdf_sha256"]} for r in manifest},
        "results_csv": rel(RESULTS_CSV), "results_csv_sha256": digest(RESULTS_CSV),
        "all_pose_csv": rel(POSE_CSV), "all_pose_csv_sha256": digest(POSE_CSV),
        "confirmed_exit_zero_runs": sum(r["run_status"] == "completed" for r in run_rows),
        "validated_output_exit_unconfirmed_runs": sum(r["run_status"] == "complete_pose_output_exit_unconfirmed" for r in run_rows),
        "completed_runs": sum(r["run_status"] in ("completed", "complete_pose_output_exit_unconfirmed") for r in run_rows),
        "failed_runs": sum(r["run_status"] not in ("completed", "complete_pose_output_exit_unconfirmed") for r in run_rows),
        "returned_poses": len(pose_rows),
        "combination_decisions": {"/".join(key): {"status": value[0], "passing_rank1_of_5": value[1], "medoid_seed": value[2]} for key, value in decisions.items()},
    }
    atomic_write(RUN_MANIFEST, (json.dumps(run_manifest, indent=2, sort_keys=True) + "\n").encode("utf-8"))
    print(f"Analysis complete: {run_manifest['completed_runs']}/40 completed, {len(pose_rows)} returned poses; "
          f"{sum(v[0] == 'PASS' for v in decisions.values())}/8 combinations PASS", flush=True)
    return 0 if run_manifest["completed_runs"] == 40 else 2


if __name__ == "__main__":
    raise SystemExit(main())

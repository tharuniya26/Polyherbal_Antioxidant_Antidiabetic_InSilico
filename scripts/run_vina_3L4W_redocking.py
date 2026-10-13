"""Controlled 3L4W neutral-MIG diagnostic redocking, Vina 1.2.7.

Run: .venv/Scripts/python.exe scripts/run_vina_3L4W_redocking.py
Only the audited dry/no-NAG/HID receptor and neutral crystal MIG are used.
Existing source files and completed seed outputs are never overwritten.
"""

from __future__ import annotations

import csv
from datetime import datetime, timezone
from hashlib import sha256
import io
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import time

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from meeko import PDBQTMolecule, RDKitMolCreate


ROOT = Path(__file__).resolve().parents[1]
REVIEW = ROOT / "data/metadata/vina_3L4W_chemistry_review_results.csv"
PREP_AUDIT = ROOT / "data/metadata/vina_alternative_receptor_preparation_results.csv"
RECEPTOR = ROOT / "structures/proteins/vina_prepared/3L4W_noNAG_D443m_D542m_H600HID_pH6p5_dry_sensitivity_Meeko0.8.0.pdbqt"
RECEPTOR_PDB = ROOT / "structures/proteins/prepared/3L4W_noNAG_D443m_D542m_H600HID_pH6p5_dry_sensitivity.pdb"
LIGAND = ROOT / "docking/validation/vina/3L4W_MIG_crystal_pose_reference.pdbqt"
REFERENCE = ROOT / "docking/validation/vina/3L4W_MIG_native_pose.sdf"
MAP = ROOT / "docking/validation/vina/3L4W_MIG_sdf_to_pdbqt_atom_map.csv"
DEPOSITED_MAP = ROOT / "docking/validation/vina/3L4W_MIG_atom_map.csv"
ORIGINAL_PDB = ROOT / "structures/proteins/original/3L4W.pdb"
VINA = Path(r"E:\docking tools\vina_1.2.7_win.exe")
OUTDIR = ROOT / "docking/validation/vina/3L4W_diagnostic_redocking"
RESULTS = ROOT / "data/metadata/vina_3L4W_redocking_results.csv"
SEEDS = (1101, 2202, 3303, 4404, 5505)
EXHAUSTIVENESS = 32
NUM_MODES = 20
ENERGY_RANGE = 5.0
CPU = 1
MIN_RMSD = 1.0
BENCHMARK = 2.0


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def verify_hashes() -> dict[str, str]:
    with REVIEW.open(newline="", encoding="utf-8") as handle:
        rows = {row["component"]: row for row in csv.DictReader(handle)}
    with PREP_AUDIT.open(newline="", encoding="utf-8") as handle:
        prep = next(row for row in csv.DictReader(handle) if row["PDB_ID"] == "3L4W")
    expected = {}
    for key in ("3L4W HID dry receptor", "MIG neutral"):
        row = rows[key]
        for name, hexdigest in zip(row["artifact_paths"].split(";"), row["artifact_SHA256"].split(";")):
            expected[ROOT / name] = hexdigest
    expected.update({
        MAP: prep["native_bridge_map_SHA256"],
        DEPOSITED_MAP: prep["native_atom_map_SHA256"],
        ORIGINAL_PDB: prep["original_PDB_SHA256"],
    })
    if set((RECEPTOR, RECEPTOR_PDB, LIGAND, REFERENCE, MAP, DEPOSITED_MAP, ORIGINAL_PDB)) != set(expected):
        raise RuntimeError("The audited file set differs from the requested primary inputs")
    for path, hexdigest in expected.items():
        if not path.is_file() or digest(path).lower() != hexdigest.lower():
            raise RuntimeError(f"Missing or changed audited input: {path}")
    if not VINA.is_file():
        raise RuntimeError(f"Vina executable absent: {VINA}")
    version = subprocess.run([str(VINA), "--version"], capture_output=True, text=True, timeout=20)
    if version.returncode or "1.2.7" not in version.stdout:
        raise RuntimeError(f"Unexpected Vina version: {version.stdout} {version.stderr}")
    return {str(path.relative_to(ROOT)): digest(path) for path in expected}


def atom_records(text: str) -> dict[int, tuple[str, str, np.ndarray]]:
    result = {}
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        serial = int(line[6:11])
        if serial in result:
            raise ValueError(f"Repeated atom serial {serial}")
        result[serial] = (line[12:16].strip()[0], line.split()[-1],
                          np.array([float(line[i : i + 8]) for i in (30, 38, 46)]))
    return result


def get_reference():
    mols = [m for m in Chem.SDMolSupplier(str(REFERENCE), removeHs=False) if m is not None]
    if len(mols) != 1:
        raise RuntimeError("Expected exactly one native MIG SDF molecule")
    molecule = mols[0]
    heavy = Chem.RemoveHs(molecule)
    if rdMolDescriptors.CalcMolFormula(molecule) != "C8H17NO5" or Chem.GetFormalCharge(molecule) != 0:
        raise RuntimeError("Native MIG formula or charge mismatch")
    if heavy.GetNumAtoms() != 14:
        raise RuntimeError("Native MIG must have 14 heavy atoms")
    cips = dict(Chem.FindMolChiralCenters(heavy, includeUnassigned=True))
    if cips != {9: "S", 10: "R", 11: "R", 12: "R"}:
        raise RuntimeError(f"Native MIG carbon stereochemistry mismatch: {cips}")
    reference_xyz = np.asarray(molecule.GetConformer().GetPositions()[:14], dtype=float)
    with MAP.open(newline="", encoding="utf-8") as handle:
        mapping = {int(row["pdbqt_atom_serial"]): int(row["sdf_atom_index_1based"]) - 1
                   for row in csv.DictReader(handle)}
    if len(mapping) != 14 or set(mapping.values()) != set(range(14)):
        raise RuntimeError("Invalid native SDF/PDBQT heavy-atom bridge")
    with DEPOSITED_MAP.open(newline="", encoding="utf-8") as handle:
        deposited = list(csv.DictReader(handle))
    if len(deposited) != 14:
        raise RuntimeError("Deposited atom map incomplete")
    original_lines = ORIGINAL_PDB.read_text().splitlines()
    for row in deposited:
        idx = int(row["sdf_atom_index_1based"]) - 1
        serial = int(row["PDB_serial"])
        line = next((line for line in original_lines if line.startswith("HETATM")
                     and int(line[6:11]) == serial), None)
        if line is None or line[17:20] != "MIG" or line[12:16].strip() != row["auth_atom_name"]:
            raise RuntimeError(f"Deposited MIG atom mismatch: {serial}")
        xyz = np.array([float(line[i : i + 8]) for i in (30, 38, 46)])
        if np.max(np.abs(xyz - reference_xyz[idx])) > 0.001:
            raise RuntimeError("Crystal SDF differs from deposited MIG coordinates")
    input_records = atom_records(LIGAND.read_text())
    if len(input_records) != 19 or sum(t != "HD" for _, t, _ in input_records.values()) != 14:
        raise RuntimeError("Neutral ligand PDBQT atom count/typing mismatch")
    for serial, idx in mapping.items():
        symbol, atom_type, xyz = input_records[serial]
        if symbol != heavy.GetAtomWithIdx(idx).GetSymbol() or atom_type == "HD":
            raise RuntimeError("Ligand atom identity/element mapping failed")
        if np.max(np.abs(xyz - reference_xyz[idx])) > 0.001:
            raise RuntimeError("Crystal ligand PDBQT pose differs from SDF")
    rebuilt = [m for m in RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(LIGAND.read_text())) if m]
    if len(rebuilt) != 1 or Chem.MolToSmiles(Chem.RemoveHs(rebuilt[0]), isomericSmiles=True) != Chem.MolToSmiles(heavy, isomericSmiles=True):
        raise RuntimeError("Meeko ligand readback graph/stereochemistry differs from crystal SDF")
    automorphisms = heavy.GetSubstructMatches(heavy, uniquify=False, useChirality=True, maxMatches=1000)
    if not automorphisms:
        raise RuntimeError("No valid native MIG graph automorphism")
    return molecule, heavy, reference_xyz, mapping, input_records, automorphisms


def derive_box(reference_xyz: np.ndarray) -> tuple[np.ndarray, tuple[int, int, int], dict]:
    center = reference_xyz.mean(axis=0)
    side_atoms = []
    for line in ORIGINAL_PDB.read_text().splitlines():
        if not line.startswith("ATOM  ") or line[21] != "A":
            continue
        residue = line[22:26].strip()
        name = line[12:16].strip()
        if (residue in ("443", "542") and name in ("CG", "OD1", "OD2")) or (
            residue == "600" and name in ("CG", "ND1", "CD2", "CE1", "NE2")):
            side_atoms.append(np.array([float(line[i : i + 8]) for i in (30, 38, 46)]))
    if len(side_atoms) != 11:
        raise RuntimeError(f"Unexpected catalytic side-chain coverage: {len(side_atoms)}")
    all_site = np.vstack([reference_xyz, side_atoms])
    max_deviation = np.max(np.abs(all_site - center), axis=0)
    dimensions = tuple(int(2 * np.ceil(v + 5.0)) for v in max_deviation)
    if dimensions != (24, 18, 18):
        raise RuntimeError(f"Unexpected box after source verification: {dimensions}")
    return center, dimensions, {
        "ligand_extrema_A": {"min": reference_xyz.min(axis=0).tolist(), "max": reference_xyz.max(axis=0).tolist()},
        "site_max_axis_deviation_A": max_deviation.tolist(),
        "margin_rule": "centroid of deposited MIG heavy atoms; five Angstrom beyond all MIG heavy atoms and Asp443/Asp542/His600 catalytic side-chain heavy atoms; round each half-extent up to an integer",
    }


def symmetry_rmsd(pose: np.ndarray, reference: np.ndarray, automorphisms) -> float:
    return min(float(np.sqrt(np.mean(np.sum((pose[list(match)] - reference) ** 2, axis=1))))
               for match in automorphisms)


def heavy_coordinates(records, mapping) -> np.ndarray:
    xyz = np.empty((14, 3), dtype=float)
    for serial, idx in mapping.items():
        xyz[idx] = records[serial][2]
    return xyz


def stereo_signs(mol: Chem.Mol, coords: np.ndarray) -> dict[int, float]:
    signs = {}
    for center in (9, 10, 11, 12):
        neighbors = sorted(a.GetIdx() for a in mol.GetAtomWithIdx(center).GetNeighbors())
        if len(neighbors) != 3:
            raise RuntimeError(f"Unexpected MIG stereocentre degree: {center}")
        v = coords[neighbors] - coords[center]
        det = float(np.linalg.det(v))
        if abs(det) < 0.2:
            raise RuntimeError(f"Collapsed MIG stereocentre geometry: {center}")
        signs[center] = math.copysign(1, det)
    return signs


def verify_pose_records(records, mapping, input_records, heavy, reference_xyz, reference_signs) -> np.ndarray:
    if set(records) != set(input_records):
        raise RuntimeError("Docked pose atom serial set differs from neutral MIG input")
    for serial in records:
        if records[serial][:2] != input_records[serial][:2]:
            raise RuntimeError(f"Docked atom name/type changed: {serial}")
        if not np.isfinite(records[serial][2]).all():
            raise RuntimeError("Nonfinite docked coordinates")
    coords = heavy_coordinates(records, mapping)
    if stereo_signs(heavy, coords) != reference_signs:
        raise RuntimeError("Docked MIG heavy-atom stereochemistry inverted")
    for bond in heavy.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        original = np.linalg.norm(reference_xyz[a] - reference_xyz[b])
        docked = np.linalg.norm(coords[a] - coords[b])
        if abs(original - docked) > 0.08:
            raise RuntimeError(f"Docked MIG bond geometry changed: {a}-{b} {original:.3f}/{docked:.3f} A")
    return coords


def deterministic_start(seed, source_text, center, dimensions, reference_xyz, mapping,
                        input_records, heavy, automorphisms, reference_signs):
    rng = np.random.default_rng(seed)
    for attempt in range(100):
        q = rng.normal(size=4)
        q /= np.linalg.norm(q)
        w, x, y, z = q
        rotation = np.array([
            [1 - 2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w)],
            [2*(x*y+z*w), 1-2*(x*x+z*z), 2*(y*z-x*w)],
            [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x*x+y*y)],
        ])
        shift = rng.uniform(-3.0, 3.0, size=3)
        if np.linalg.norm(shift) < 3.0:
            continue
        transformed = {serial: (record[2] - center) @ rotation.T + center + shift
                       for serial, record in input_records.items()}
        trial = {serial: (input_records[serial][0], input_records[serial][1], xyz)
                 for serial, xyz in transformed.items()}
        coords = heavy_coordinates(trial, mapping)
        if symmetry_rmsd(coords, reference_xyz, automorphisms) < 5.0:
            continue
        if np.any(np.abs(coords - center) > np.asarray(dimensions) / 2 - 0.5):
            continue
        if stereo_signs(heavy, coords) != reference_signs:
            raise RuntimeError("Rigid transform inverted MIG stereochemistry")
        lines = []
        for line in source_text.splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                serial = int(line[6:11])
                p = transformed[serial]
                line = line[:30] + "".join(f"{v:8.3f}" for v in p) + line[54:]
            lines.append(line)
        result = "\n".join(lines) + "\n"
        parsed = atom_records(result)
        verify_pose_records(parsed, mapping, input_records, heavy, reference_xyz, reference_signs)
        # The REMARK SMILES / INDEX MAP records remain intact after rigid transform.
        rebuilt = [m for m in RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(result)) if m]
        if len(rebuilt) != 1 or Chem.MolToSmiles(Chem.RemoveHs(rebuilt[0]), isomericSmiles=True) != Chem.MolToSmiles(heavy, isomericSmiles=True):
            raise RuntimeError("Randomized starting PDBQT graph readback failed")
        return result, symmetry_rmsd(heavy_coordinates(parsed, mapping), reference_xyz, automorphisms), shift.tolist(), attempt + 1
    raise RuntimeError(f"Could not generate start away from crystal for seed {seed}")


def write_immutable(path: Path, contents: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != contents:
            raise RuntimeError(f"Refusing to overwrite changed diagnostic artifact: {path}")
        return
    path.write_bytes(contents)


def parse_output(path: Path, mapping, input_records, heavy, reference_xyz,
                 automorphisms, reference_signs):
    text = path.read_text()
    # Meeko reads the complete multi-pose Vina PDBQT and rebuilds its chemistry.
    reconstructed = [m for m in RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(text)) if m]
    segments = re.findall(r"(?ms)^MODEL\s+(\d+)\s*\n(.*?)^ENDMDL\s*$", text)
    if not segments:
        raise RuntimeError("Vina output contains no MODEL/ENDMDL pose records")
    if len(reconstructed) != 1 or reconstructed[0].GetNumConformers() != len(segments):
        raise RuntimeError(f"Meeko readback expected one molecule with {len(segments)} conformers; "
                           f"got {[(m.GetNumConformers()) for m in reconstructed]}")
    target_smiles = Chem.MolToSmiles(heavy, isomericSmiles=True)
    metrics = []
    molecule = reconstructed[0]
    for rank_text, segment in segments:
        rank = int(rank_text)
        if rank != len(metrics) + 1:
            raise RuntimeError("Nonsequential Vina mode ranks")
        if Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True) != target_smiles:
            raise RuntimeError(f"Mode {rank} identity/connectivity/stereochemistry failed Meeko readback")
        match = re.search(r"(?m)^REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)\s+", segment)
        if not match:
            raise RuntimeError(f"Mode {rank} lacks a Vina affinity record")
        records = atom_records(segment)
        xyz = verify_pose_records(records, mapping, input_records, heavy, reference_xyz, reference_signs)
        metrics.append({"pose_rank": rank, "score_kcal_mol": float(match.group(1)),
                        "heavy_atom_RMSD_A": symmetry_rmsd(xyz, reference_xyz, automorphisms)})
    if any(metrics[i]["score_kcal_mol"] < metrics[i-1]["score_kcal_mol"] - 0.01
           for i in range(1, len(metrics))):
        raise RuntimeError("Vina pose energies unexpectedly out of rank order")
    return metrics


def csv_bytes(rows, fields) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue().encode("utf-8")


def main() -> int:
    hashes = verify_hashes()
    ref_mol, heavy, reference_xyz, mapping, input_records, automorphisms = get_reference()
    reference_signs = stereo_signs(heavy, reference_xyz)
    center, dimensions, box_detail = derive_box(reference_xyz)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    manifest = {
        "experiment": "3L4W neutral-MIG no-NAG dry/HID diagnostic redocking",
        "created_UTC": datetime.now(timezone.utc).isoformat(),
        "vina_version": "1.2.7", "rdkit_version": Chem.rdBase.rdkitVersion,
        "meeko_version": __import__("meeko").__version__,
        "python_version": sys.version.split()[0],
        "vina_executable": str(VINA), "vina_executable_SHA256": digest(VINA),
        "audited_input_SHA256": hashes,
        "seeds": list(SEEDS), "cpu_per_run": CPU,
        "exhaustiveness": EXHAUSTIVENESS, "num_modes": NUM_MODES,
        "energy_range_kcal_mol": ENERGY_RANGE, "min_rmsd_A": MIN_RMSD,
        "box_center_A": center.tolist(), "box_dimensions_A": dimensions,
        "box_derivation": box_detail,
        "atom_map": str(MAP.relative_to(ROOT)), "chemical_automorphisms": len(automorphisms),
        "starting_pose_method": "NumPy PCG64(seed) rigid quaternion rotation plus +/-3 A translation; at least 5 A unmapped crystal-frame RMSD; no atom/bond/protonation changes. Vina global_search uses independent randomized starts and --seed controls its Monte Carlo generator.",
        "receptor_model_limit": "No-NAG local sensitivity; Asp443-/Asp542-/His600-HID; dry; not full glycan-retaining validation",
    }
    manifest_path = OUTDIR / "experiment_manifest.json"
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text())
        manifest["created_UTC"] = previous["created_UTC"]
    write_immutable(manifest_path, (json.dumps(manifest, indent=2) + "\n").encode("utf-8"))
    print("Verified audited receptor, neutral MIG, CCD-backed maps and deposited coordinates", flush=True)
    print("Box center", center.tolist(), "dimensions", dimensions, "graph automorphisms", len(automorphisms), flush=True)
    rows = []
    poses = []
    ligand_text = LIGAND.read_text()
    for seed in SEEDS:
        base = f"seed_{seed}"
        start = OUTDIR / f"{base}_randomized_start.pdbqt"
        output = OUTDIR / f"{base}_docked.pdbqt"
        log = OUTDIR / f"{base}_vina.log"
        start_text, start_rmsd, shift, attempt = deterministic_start(
            seed, ligand_text, center, dimensions, reference_xyz, mapping,
            input_records, heavy, automorphisms, reference_signs)
        write_immutable(start, start_text.encode("ascii"))
        command = [str(VINA), "--receptor", str(RECEPTOR), "--ligand", str(start),
                   "--center_x", f"{center[0]:.8f}", "--center_y", f"{center[1]:.8f}",
                   "--center_z", f"{center[2]:.8f}", "--size_x", str(dimensions[0]),
                   "--size_y", str(dimensions[1]), "--size_z", str(dimensions[2]),
                   "--scoring", "vina", "--exhaustiveness", str(EXHAUSTIVENESS),
                   "--num_modes", str(NUM_MODES), "--energy_range", str(ENERGY_RANGE),
                   "--min_rmsd", str(MIN_RMSD), "--seed", str(seed), "--cpu", str(CPU),
                   "--out", str(output), "--verbosity", "1"]
        command_string = subprocess.list2cmdline(command)
        start_time = time.monotonic()
        if output.exists() and log.exists():
            stdout = log.read_text()
            code_match = re.search(r"(?m)^EXIT_CODE: (\d+)\s*$", stdout)
            elapsed_match = re.search(r"(?m)^ELAPSED_SECONDS: ([\d.]+)\s*$", stdout)
            if code_match is None or elapsed_match is None:
                raise RuntimeError(f"Existing log lacks completion metadata for seed {seed}")
            run_code = int(code_match.group(1))
            elapsed = float(elapsed_match.group(1))
            if command_string not in stdout:
                raise RuntimeError(f"Existing log does not match command for seed {seed}")
            print(f"Seed {seed}: reanalyzing existing immutable output", flush=True)
        elif output.exists() or log.exists():
            raise RuntimeError(f"Partial existing seed artifacts require manual review: {seed}")
        else:
            print(f"Seed {seed}: running Vina 1.2.7, start RMSD {start_rmsd:.3f} A", flush=True)
            run = subprocess.run(command, cwd=str(ROOT), capture_output=True, text=True, timeout=3600)
            elapsed = time.monotonic() - start_time
            run_code = run.returncode
            stdout = (f"COMMAND: {command_string}\nSEED: {seed}\n"
                      f"RECEPTOR_SHA256: {hashes[str(RECEPTOR.relative_to(ROOT))]}\n"
                      f"RANDOMIZED_LIGAND_SHA256: {digest(start)}\n"
                      f"EXIT_CODE: {run_code}\nELAPSED_SECONDS: {elapsed:.3f}\n"
                      f"--- STDOUT ---\n{run.stdout}\n--- STDERR ---\n{run.stderr}\n")
            write_immutable(log, stdout.encode("utf-8"))
        row = {"PDB_ID": "3L4W", "seed": seed, "status": "failed", "vina_version": "1.2.7",
               "receptor_model": "dry_noNAG_D443m_D542m_H600HID_pH6p5",
               "receptor_path": str(RECEPTOR.relative_to(ROOT)),
               "receptor_SHA256": digest(RECEPTOR),
               "crystal_reference_path": str(REFERENCE.relative_to(ROOT)),
               "crystal_reference_SHA256": digest(REFERENCE),
               "randomized_input_path": str(start.relative_to(ROOT)),
               "randomized_input_SHA256": digest(start),
               "randomized_start_RMSD_A": f"{start_rmsd:.4f}",
               "randomized_translation_xyz_A": ";".join(f"{x:.4f}" for x in shift),
               "randomization_attempt": attempt,
               "box_center_xyz_A": ";".join(f"{x:.8f}" for x in center),
               "box_size_xyz_A": ";".join(map(str, dimensions)),
               "exhaustiveness": EXHAUSTIVENESS, "num_modes_requested": NUM_MODES,
               "energy_range_kcal_mol": ENERGY_RANGE, "min_rmsd_A": MIN_RMSD, "cpu": CPU,
               "command": command_string, "log_path": str(log.relative_to(ROOT)),
               "log_SHA256": digest(log), "output_path": str(output.relative_to(ROOT)) if output.exists() else "",
               "output_SHA256": digest(output) if output.exists() else "",
               "exit_code": run_code, "elapsed_seconds": f"{elapsed:.3f}" if elapsed is not None else "",
               "n_modes_returned": "", "top_score_kcal_mol": "", "top_RMSD_A": "",
               "best_RMSD_A": "", "best_RMSD_pose_rank": "", "best_RMSD_pose_score_kcal_mol": "",
               "top_pose_within_2A": "", "any_pose_within_2A": "", "error": ""}
        try:
            if run_code or not output.is_file() or not output.stat().st_size:
                raise RuntimeError(f"Vina process failed with exit code {run_code}; see saved log")
            if "Performing docking (random seed: " + str(seed) not in stdout:
                raise RuntimeError("Saved Vina log lacks requested search seed confirmation")
            modes = parse_output(output, mapping, input_records, heavy, reference_xyz,
                                 automorphisms, reference_signs)
            if not modes:
                raise RuntimeError("No valid docked poses returned")
            best = min(modes, key=lambda p: p["heavy_atom_RMSD_A"])
            top = modes[0]
            row.update({"status": "completed", "n_modes_returned": len(modes),
                        "top_score_kcal_mol": f"{top['score_kcal_mol']:.3f}",
                        "top_RMSD_A": f"{top['heavy_atom_RMSD_A']:.4f}",
                        "best_RMSD_A": f"{best['heavy_atom_RMSD_A']:.4f}",
                        "best_RMSD_pose_rank": best["pose_rank"],
                        "best_RMSD_pose_score_kcal_mol": f"{best['score_kcal_mol']:.3f}",
                        "top_pose_within_2A": top["heavy_atom_RMSD_A"] <= BENCHMARK,
                        "any_pose_within_2A": best["heavy_atom_RMSD_A"] <= BENCHMARK})
            poses.extend({"seed": seed, **pose, "within_2A": pose["heavy_atom_RMSD_A"] <= BENCHMARK}
                         for pose in modes)
            print(f"Seed {seed}: {len(modes)} modes; top RMSD={top['heavy_atom_RMSD_A']:.3f} A; "
                  f"best={best['heavy_atom_RMSD_A']:.3f} A rank {best['pose_rank']}", flush=True)
        except Exception as error:
            row["error"] = str(error)
            print(f"Seed {seed}: FAILED: {error}", flush=True)
        rows.append(row)
        write_immutable(OUTDIR / f"seed_{seed}_pose_metrics.csv",
                        csv_bytes([p for p in poses if p["seed"] == seed],
                                  ["seed", "pose_rank", "score_kcal_mol", "heavy_atom_RMSD_A", "within_2A"]))
        if row["status"] != "completed":
            break
    result_fields = list(rows[0])
    write_immutable(RESULTS, csv_bytes(rows, result_fields))
    write_immutable(OUTDIR / "all_pose_metrics.csv", csv_bytes(
        poses, ["seed", "pose_rank", "score_kcal_mol", "heavy_atom_RMSD_A", "within_2A"]))
    return 0 if len(rows) == len(SEEDS) and all(row["status"] == "completed" for row in rows) else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"FATAL PREFLIGHT OR ANALYSIS ERROR: {exc}", file=sys.stderr)
        sys.exit(2)

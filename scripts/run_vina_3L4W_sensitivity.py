"""Audited 3L4W/MIG search and charge-state sensitivity redocking.

Run with .venv/Scripts/python.exe scripts/run_vina_3L4W_sensitivity.py
Existing diagnostic outputs and all source inputs remain immutable.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import importlib.util
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
BASE_SPEC = importlib.util.spec_from_file_location("redock_base", ROOT / "scripts/run_vina_3L4W_redocking.py")
base = importlib.util.module_from_spec(BASE_SPEC)
BASE_SPEC.loader.exec_module(base)
OUT = ROOT / "docking/validation/vina/3L4W_sensitivity"
RESULTS = ROOT / "data/metadata/vina_3L4W_sensitivity_results.csv"
CHARGED_SDF = ROOT / "docking/validation/vina/3L4W_MIG_NHplus_crystal_pose_sensitivity.sdf"
CHARGED_PDBQT = ROOT / "docking/validation/vina/3L4W_MIG_NHplus_crystal_pose_sensitivity.pdbqt"
SEEDS = base.SEEDS
EXHAUSTIVENESS = 64


def audited_variant(name: str) -> dict:
    with base.REVIEW.open(newline="", encoding="utf-8") as f:
        row = next((r for r in csv.DictReader(f) if r["component"] == name), None)
    if row is None:
        raise RuntimeError(f"Missing chemistry-review variant {name}")
    paths = row["artifact_paths"].split(";")
    hashes = row["artifact_SHA256"].split(";")
    if len(paths) != len(hashes) or not paths:
        raise RuntimeError(f"Invalid audit paths/hashes for {name}")
    for p, h in zip(paths, hashes):
        path = ROOT / p
        if not path.is_file() or base.digest(path) != h:
            raise RuntimeError(f"Missing or changed audited variant {path}")
    return row


def smiles_without_nitrogen_tag(mol: Chem.Mol) -> str:
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        if atom.GetSymbol() == "N":
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            if atom.HasProp("_CIPCode"):
                atom.ClearProp("_CIPCode")
    return Chem.MolToSmiles(copy, isomericSmiles=True)


def nitrogen_geometry(records, mapping, mol: Chem.Mol) -> float:
    """Signed tetrahedron for the assigned SDF N and its three C plus one H."""
    nitrogen = mol.GetAtomWithIdx(13)
    neighbors = [atom.GetIdx() for atom in nitrogen.GetNeighbors()]
    if [mol.GetAtomWithIdx(i).GetSymbol() for i in neighbors] != ["C", "C", "C", "H"]:
        raise RuntimeError("Unexpected charged MIG nitrogen neighbors")
    inverse = {idx: serial for serial, idx in mapping.items()}
    n_serial = inverse[13]
    n_xyz = records[n_serial][2]
    h_candidates = [serial for serial, (symbol, atype, xyz) in records.items()
                    if symbol == "H" and atype == "HD" and 0.8 < np.linalg.norm(xyz-n_xyz) < 1.2]
    if len(h_candidates) != 1:
        raise RuntimeError("Protonated MIG nitrogen hydrogen missing/ambiguous")
    vectors = np.asarray([records[inverse[i]][2]-n_xyz for i in neighbors[:3]] +
                         [records[h_candidates[0]][2]-n_xyz])
    determinant = float(np.linalg.det(vectors[:3] - vectors[3]))
    if abs(determinant) < 2.0:
        raise RuntimeError("Protonated MIG nitrogen tetrahedron collapsed")
    return determinant


def ligand_case(path: Path, sdf_path: Path, reference_xyz, heavy, reference_signs):
    mols = [m for m in Chem.SDMolSupplier(str(sdf_path), removeHs=False) if m]
    if len(mols) != 1:
        raise RuntimeError(f"Expected one SDF molecule: {sdf_path}")
    mol = mols[0]
    expected_charge = 0 if path == base.LIGAND else 1
    expected_formula = "C8H17NO5" if expected_charge == 0 else "C8H18NO5+"
    if Chem.GetFormalCharge(mol) != expected_charge or rdMolDescriptors.CalcMolFormula(mol) != expected_formula:
        raise RuntimeError(f"Ligand formula/charge mismatch: {sdf_path}")
    candidate_heavy = Chem.RemoveHs(mol)
    if candidate_heavy.GetNumAtoms() != 14:
        raise RuntimeError("MIG must have 14 heavy atoms")
    if dict(Chem.FindMolChiralCenters(candidate_heavy, includeUnassigned=True)).items() < dict(Chem.FindMolChiralCenters(heavy, includeUnassigned=True)).items():
        raise RuntimeError("Carbon stereocentres lost")
    candidate_xyz = np.asarray(mol.GetConformer().GetPositions()[:14], dtype=float)
    if np.max(np.abs(candidate_xyz - reference_xyz)) > 0.001:
        raise RuntimeError("Sensitivity ligand heavy coordinates differ from crystal")
    if base.stereo_signs(heavy, candidate_xyz) != reference_signs:
        raise RuntimeError("Sensitivity ligand stereochemistry changed")
    records = base.atom_records(path.read_text())
    if sum(r[1] != "HD" for r in records.values()) != 14:
        raise RuntimeError("Ligand PDBQT heavy-atom count differs from 14")
    # Deposited heavy coordinates provide an unambiguous serial-to-SDF map for either charge state.
    mapping = {}
    for serial, (element, atype, xyz) in records.items():
        if atype == "HD":
            continue
        matches = [i for i, ref in enumerate(reference_xyz) if np.linalg.norm(ref - xyz) <= 0.001]
        if len(matches) != 1 or element != heavy.GetAtomWithIdx(matches[0]).GetSymbol():
            raise RuntimeError(f"Ambiguous or wrong MIG atom identity at serial {serial}")
        mapping[serial] = matches[0]
    if len(mapping) != 14 or set(mapping.values()) != set(range(14)):
        raise RuntimeError("Ligand heavy atom map incomplete")
    if path == base.LIGAND:
        _, _, _, original_map, _, _ = base.get_reference()
        if mapping != original_map:
            raise RuntimeError("Neutral MIG map changed from previous experiment")
    rebuilt = [m for m in RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(path.read_text())) if m]
    if len(rebuilt) != 1 or smiles_without_nitrogen_tag(Chem.RemoveHs(rebuilt[0])) != smiles_without_nitrogen_tag(candidate_heavy):
        raise RuntimeError("Ligand Meeko graph/stereochemistry readback failed")
    if expected_charge == 1:
        if dict(Chem.FindMolChiralCenters(candidate_heavy, includeUnassigned=True)).get(13) != "R":
            raise RuntimeError("Charged SDF nitrogen reference stereochemistry changed")
        sdf_xyz = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
        n = sdf_xyz[13]
        neighbor_ids = [a.GetIdx() for a in mol.GetAtomWithIdx(13).GetNeighbors()]
        vectors = sdf_xyz[neighbor_ids] - n
        sdf_det = float(np.linalg.det(vectors[:3]-vectors[3]))
        if sdf_det * nitrogen_geometry(records, mapping, mol) <= 0:
            raise RuntimeError("Charged PDBQT nitrogen geometry disagrees with SDF")
    # Meeko PDBQT partial charges are bookkeeping checks, not Vina binding free energies.
    partial = sum(float(line[70:76]) for line in path.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM")))
    if abs(partial - expected_charge) > 0.1:
        raise RuntimeError(f"Ligand PDBQT charge sum {partial} differs from {expected_charge}")
    return records, mapping, partial, mol


def random_start(seed, source, records, mapping, heavy, reference_xyz, center, dims, automorphisms, signs, target_smiles, ligand_mol):
    rng = np.random.default_rng(seed)
    for attempt in range(100):
        q = rng.normal(size=4)
        q /= np.linalg.norm(q)
        w, x, y, z = q
        rotation = np.array([[1-2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w)],
                             [2*(x*y+z*w), 1-2*(x*x+z*z), 2*(y*z-x*w)],
                             [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x*x+y*y)]])
        shift = rng.uniform(-3, 3, size=3)
        if np.linalg.norm(shift) < 3:
            continue
        transformed = {serial: (rec[2]-center) @ rotation.T + center + shift for serial, rec in records.items()}
        trial = {s: (records[s][0], records[s][1], xyz) for s, xyz in transformed.items()}
        xyz = base.heavy_coordinates(trial, mapping)
        rmsd = base.symmetry_rmsd(xyz, reference_xyz, automorphisms)
        if rmsd < 5 or np.any(np.abs(xyz-center) > np.asarray(dims)/2-0.5):
            continue
        if base.stereo_signs(heavy, xyz) != signs:
            raise RuntimeError("Rigid transform inverted stereochemistry")
        lines = []
        for line in source.splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                p = transformed[int(line[6:11])]
                line = line[:30] + "".join(f"{v:8.3f}" for v in p) + line[54:]
            lines.append(line)
        result = "\n".join(lines) + "\n"
        parsed = base.atom_records(result)
        base.verify_pose_records(parsed, mapping, records, heavy, reference_xyz, signs)
        if Chem.GetFormalCharge(ligand_mol) == 1 and nitrogen_geometry(parsed, mapping, ligand_mol) * nitrogen_geometry(records, mapping, ligand_mol) <= 0:
            raise RuntimeError("Randomized start inverted protonated nitrogen")
        rebuilt = [m for m in RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(result)) if m]
        if len(rebuilt) != 1 or smiles_without_nitrogen_tag(Chem.RemoveHs(rebuilt[0])) != target_smiles:
            raise RuntimeError("Randomized ligand chemistry failed readback")
        return result, rmsd, shift.tolist(), attempt+1
    raise RuntimeError(f"Could not generate independent start: {seed}")


def parse_charged_output(path, mapping, input_records, heavy, reference_xyz, automorphisms, signs, ligand_mol):
    text = path.read_text()
    reconstructed = [m for m in RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(text)) if m]
    segments = re.findall(r"(?ms)^MODEL\s+(\d+)\s*\n(.*?)^ENDMDL\s*$", text)
    if not segments or len(reconstructed) != 1 or reconstructed[0].GetNumConformers() != len(segments):
        raise RuntimeError("Charged MIG multi-pose Meeko readback failed")
    target = smiles_without_nitrogen_tag(Chem.RemoveHs(ligand_mol))
    if smiles_without_nitrogen_tag(Chem.RemoveHs(reconstructed[0])) != target:
        raise RuntimeError("Charged MIG graph/carbon stereochemistry differs after docking")
    source_sign = nitrogen_geometry(input_records, mapping, ligand_mol)
    poses = []
    for rank_text, segment in segments:
        rank = int(rank_text)
        if rank != len(poses)+1:
            raise RuntimeError("Nonsequential Vina pose ranks")
        score = re.search(r"(?m)^REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)\s+", segment)
        if score is None:
            raise RuntimeError(f"Missing pose score at rank {rank}")
        records = base.atom_records(segment)
        xyz = base.verify_pose_records(records, mapping, input_records, heavy, reference_xyz, signs)
        if nitrogen_geometry(records, mapping, ligand_mol) * source_sign <= 0:
            raise RuntimeError(f"Charged MIG nitrogen geometry inverted at rank {rank}")
        poses.append({"pose_rank": rank, "score_kcal_mol": float(score.group(1)),
                      "heavy_atom_RMSD_A": base.symmetry_rmsd(xyz, reference_xyz, automorphisms)})
    if any(poses[i]["score_kcal_mol"] < poses[i-1]["score_kcal_mol"]-0.01 for i in range(1,len(poses))):
        raise RuntimeError("Pose score ranks not monotonic")
    return poses


def main() -> int:
    charged_only = len(sys.argv) == 2 and sys.argv[1] == "--charged-only"
    if len(sys.argv) > 1 and not charged_only:
        raise RuntimeError("Only optional --charged-only is supported")
    primary_hashes = base.verify_hashes()
    audited_variant("MIG amine +1")
    audited_variant("3L4W HID sitewaters receptor")
    version = subprocess.run([str(base.VINA), "--version"], capture_output=True, text=True, timeout=20)
    if version.returncode or "1.2.7" not in version.stdout:
        raise RuntimeError("Vina version mismatch")
    with (ROOT / "data/metadata/vina_3L4W_redocking_results.csv").open(newline="", encoding="utf-8") as f:
        previous = list(csv.DictReader(f))
    if len(previous) != 5 or {int(r["seed"]) for r in previous} != set(SEEDS) or any(r["status"] != "completed" for r in previous):
        raise RuntimeError("Previous five-run diagnostic audit incomplete")
    for r in previous:
        for path_field, hash_field in (("output_path", "output_SHA256"), ("log_path", "log_SHA256"), ("randomized_input_path", "randomized_input_SHA256")):
            p = ROOT / r[path_field]
            if not p.is_file() or base.digest(p) != r[hash_field]:
                raise RuntimeError(f"Earlier experiment changed: {p}")
    ref_mol, heavy, xyz, neutral_map, neutral_records, automorphisms = base.get_reference()
    signs = base.stereo_signs(heavy, xyz)
    center, dims, _ = base.derive_box(xyz)
    cases = [("neutral_MIG_dry_HID_E64", base.RECEPTOR, base.LIGAND, base.REFERENCE),
             ("protonated_MIG_dry_HID_E64", base.RECEPTOR, CHARGED_PDBQT, CHARGED_SDF)]
    if charged_only:
        cases = cases[1:]
    # The wet variant is audited but not docked: fixed waters lack a validated exchange/orientation model.
    manifest = {"created_UTC": datetime.now(timezone.utc).isoformat(), "experiment": "3L4W/MIG E64 convergence and amine-protonation sensitivity", "seeds": list(SEEDS),
                "vina_version": "1.2.7", "vina_executable_SHA256": base.digest(base.VINA),
                "rdkit_version": Chem.rdBase.rdkitVersion, "meeko_version": __import__("meeko").__version__,
                "python_version": sys.version.split()[0], "primary_input_SHA256": primary_hashes,
                "charged_SDF_SHA256": base.digest(CHARGED_SDF), "charged_PDBQT_SHA256": base.digest(CHARGED_PDBQT),
                "box_center_A": center.tolist(), "box_size_A": dims, "exhaustiveness": EXHAUSTIVENESS,
                "num_modes": 20, "energy_range_kcal_mol": 5.0, "min_rmsd_A": 1.0, "cpu": 1,
                "excluded_variants": {"sitewaters": "Static-water orientation and exchange unsupported for standard rigid Vina scoring", "glycan_retained": "No validated Meeko glycan template", "Asp542_protonated": "No prepared, audited receptor variant", "His600_HIE_or_HIP": "HIE trial had local hydrogen clash; HIP automatic state not justified"}}
    OUT.mkdir(parents=True, exist_ok=True)
    manifest_path = OUT / "experiment_manifest.json"
    if manifest_path.exists():
        manifest["created_UTC"] = json.loads(manifest_path.read_text())["created_UTC"]
    base.write_immutable(manifest_path, (json.dumps(manifest, indent=2)+"\n").encode())
    rows, all_poses = [], []
    for case_name, receptor, ligand, sdf in cases:
        records, mapping, partial_charge, mol = ligand_case(ligand, sdf, xyz, heavy, signs)
        target_smiles = smiles_without_nitrogen_tag(Chem.RemoveHs(mol))
        source = ligand.read_text()
        for seed in SEEDS:
            prefix = f"{case_name}_seed_{seed}"
            start, output, log = (OUT / f"{prefix}_{suffix}" for suffix in ("start.pdbqt", "docked.pdbqt", "vina.log"))
            start_text, start_rmsd, shift, attempt = random_start(seed, source, records, mapping, heavy, xyz, center, dims, automorphisms, signs, target_smiles, mol)
            base.write_immutable(start, start_text.encode("ascii"))
            cmd = [str(base.VINA), "--receptor", str(receptor), "--ligand", str(start),
                   "--center_x", f"{center[0]:.8f}", "--center_y", f"{center[1]:.8f}", "--center_z", f"{center[2]:.8f}",
                   "--size_x", str(dims[0]), "--size_y", str(dims[1]), "--size_z", str(dims[2]),
                   "--scoring", "vina", "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", "20",
                   "--energy_range", "5.0", "--min_rmsd", "1.0", "--seed", str(seed), "--cpu", "1", "--out", str(output), "--verbosity", "1"]
            command = subprocess.list2cmdline(cmd)
            if output.exists() and log.exists():
                stdout = log.read_text()
                code = re.search(r"(?m)^EXIT_CODE: (\d+)$", stdout)
                elapsed_match = re.search(r"(?m)^ELAPSED_SECONDS: ([\d.]+)$", stdout)
                if command not in stdout or code is None or elapsed_match is None:
                    raise RuntimeError(f"Existing output/log incomplete: {prefix}")
                exit_code, elapsed = int(code.group(1)), float(elapsed_match.group(1))
                print(f"Reanalysing {prefix}", flush=True)
            elif output.exists() or log.exists():
                raise RuntimeError(f"Partial existing run: {prefix}")
            else:
                print(f"Docking {prefix}, start RMSD {start_rmsd:.3f} A", flush=True)
                t0 = time.monotonic()
                run = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=3600)
                elapsed, exit_code = time.monotonic()-t0, run.returncode
                stdout = (f"COMMAND: {command}\nSEED: {seed}\nRECEPTOR_SHA256: {base.digest(receptor)}\nLIGAND_SHA256: {base.digest(ligand)}\nSTART_SHA256: {base.digest(start)}\nEXIT_CODE: {exit_code}\nELAPSED_SECONDS: {elapsed:.3f}\n--- STDOUT ---\n{run.stdout}\n--- STDERR ---\n{run.stderr}\n")
                base.write_immutable(log, stdout.encode())
            row = {"PDB_ID": "3L4W", "case": case_name, "seed": seed, "status": "failed", "receptor_microstate": "Asp443-/Asp542-/His600-HID, dry, no NAG", "ligand_microstate": "neutral" if ligand == base.LIGAND else "amine +1", "vina_version": "1.2.7", "exhaustiveness": EXHAUSTIVENESS,
                   "receptor_path": str(receptor.relative_to(ROOT)), "receptor_SHA256": base.digest(receptor), "ligand_path": str(ligand.relative_to(ROOT)), "ligand_SHA256": base.digest(ligand), "ligand_partial_charge_sum": f"{partial_charge:.3f}", "crystal_reference_path": str(base.REFERENCE.relative_to(ROOT)), "crystal_reference_SHA256": base.digest(base.REFERENCE),
                   "start_path": str(start.relative_to(ROOT)), "start_SHA256": base.digest(start), "start_RMSD_A": f"{start_rmsd:.4f}", "randomization_attempt": attempt, "randomized_translation_xyz_A": ";".join(f"{n:.4f}" for n in shift),
                   "box_center_xyz_A": ";".join(f"{v:.8f}" for v in center), "box_size_xyz_A": ";".join(str(v) for v in dims), "num_modes_requested": 20, "energy_range_kcal_mol": 5.0, "min_rmsd_A": 1.0, "cpu": 1,
                   "command": command, "log_path": str(log.relative_to(ROOT)), "log_SHA256": base.digest(log), "output_path": str(output.relative_to(ROOT)) if output.exists() else "", "output_SHA256": base.digest(output) if output.exists() else "", "exit_code": exit_code, "elapsed_seconds": f"{elapsed:.3f}",
                   "n_modes_returned": "", "top_score_kcal_mol": "", "top_RMSD_A": "", "best_RMSD_A": "", "best_RMSD_rank": "", "best_RMSD_score_kcal_mol": "", "top_within_2A": "", "any_within_2A": "", "error": ""}
            try:
                if exit_code != 0 or not output.is_file() or output.stat().st_size == 0:
                    raise RuntimeError(f"Vina exit {exit_code}; inspect log")
                if f"Performing docking (random seed: {seed}" not in stdout:
                    raise RuntimeError("Vina log does not confirm seed")
                poses = (base.parse_output(output, mapping, records, heavy, xyz, automorphisms, signs)
                         if ligand == base.LIGAND else parse_charged_output(output, mapping, records, heavy, xyz, automorphisms, signs, mol))
                top, best = poses[0], min(poses, key=lambda p: p["heavy_atom_RMSD_A"])
                row.update(status="completed", n_modes_returned=len(poses), top_score_kcal_mol=f"{top['score_kcal_mol']:.3f}", top_RMSD_A=f"{top['heavy_atom_RMSD_A']:.4f}", best_RMSD_A=f"{best['heavy_atom_RMSD_A']:.4f}", best_RMSD_rank=best["pose_rank"], best_RMSD_score_kcal_mol=f"{best['score_kcal_mol']:.3f}", top_within_2A=top["heavy_atom_RMSD_A"] <= 2.0, any_within_2A=best["heavy_atom_RMSD_A"] <= 2.0)
                case_poses = [{"case": case_name, "seed": seed, **p, "within_2A": p["heavy_atom_RMSD_A"] <= 2.0} for p in poses]
                all_poses.extend(case_poses)
                base.write_immutable(OUT / f"{prefix}_pose_metrics.csv", base.csv_bytes(case_poses, ["case", "seed", "pose_rank", "score_kcal_mol", "heavy_atom_RMSD_A", "within_2A"]))
                print(f"{prefix}: top {top['heavy_atom_RMSD_A']:.3f} A, best {best['heavy_atom_RMSD_A']:.3f} A rank {best['pose_rank']}", flush=True)
            except Exception as exc:
                row["error"] = str(exc)
                print(f"FAILED {prefix}: {exc}", flush=True)
            rows.append(row)
            if row["status"] != "completed":
                break
        if rows[-1]["status"] != "completed":
            break
    if rows and not charged_only:
        base.write_immutable(RESULTS, base.csv_bytes(rows, list(rows[0])))
    if all_poses and not charged_only:
        base.write_immutable(OUT / "all_pose_metrics.csv", base.csv_bytes(all_poses, ["case", "seed", "pose_rank", "score_kcal_mol", "heavy_atom_RMSD_A", "within_2A"]))
    return 0 if len(rows) == len(cases)*len(SEEDS) and all(r["status"] == "completed" for r in rows) else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"FATAL PREFLIGHT/ANALYSIS ERROR: {exc}", file=sys.stderr)
        sys.exit(2)

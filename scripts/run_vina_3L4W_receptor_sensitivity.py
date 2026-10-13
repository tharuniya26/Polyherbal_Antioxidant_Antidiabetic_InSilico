"""Prepare an audited 3L4W Asp542H/His600-HID local receptor and redock MIG.

Run .venv/Scripts/python.exe scripts/run_vina_3L4W_receptor_sensitivity.py --prepare-only
then inspect the quality-control output. Without that flag, run or reanalyse the
five prespecified Vina seeds for neutral and protonated MIG. All earlier inputs,
poses and logs remain unchanged; new outputs are immutable on rerun.
"""
from __future__ import annotations

from collections import Counter
import csv
from datetime import datetime, timezone
import importlib.util
import io
import json
import math
from pathlib import Path
import random
import re
import subprocess
import sys
import tempfile
import time

import numpy as np
import openmm
from openmm.app import Modeller, PDBFile
from openmm import unit
from meeko import PDBQTReceptor
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("mig_sens", ROOT / "scripts/run_vina_3L4W_sensitivity.py")
sens = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(sens)
base = sens.base
SOURCE_PDB = base.RECEPTOR_PDB
SOURCE_PDBQT = base.RECEPTOR
STEM = "3L4W_noNAG_D443m_D542H_OD2_H600HID_pH6p5_dry_sensitivity"
NEW_PDB = ROOT / "structures/proteins/prepared" / f"{STEM}.pdb"
NEW_PDBQT = ROOT / "structures/proteins/vina_prepared" / f"{STEM}_Meeko0.8.0.pdbqt"
OUT = ROOT / "docking/validation/vina/3L4W_receptor_sensitivity"
RESULTS = ROOT / "data/metadata/vina_3L4W_receptor_sensitivity_results.csv"
CLI = ROOT / ".venv/Scripts/mk_prepare_receptor.exe"
TEMPLATE = ROOT / ".venv/Lib/site-packages/meeko/data/residue_chem_templates.json"
SEEDS = base.SEEDS


def pos(line: str) -> np.ndarray:
    return np.array([float(line[i:i+8]) for i in (30,38,46)], dtype=float)


def atom_lines(text: str):
    return [line for line in text.splitlines() if line.startswith(("ATOM  ", "HETATM"))]


def key(line: str):
    return (line[21], line[22:26], line[26], line[12:16].strip())


def heavy_map(text: str):
    return {key(line): tuple(pos(line)) for line in atom_lines(text) if line[76:78].strip() not in ("H", "D")}


def _variant_proton(source_text: str):
    # OpenMM proposes orientation without the crystal ligand. Its provisional
    # long O-H vector is normalized to an ordinary 0.98 A bond; other H remain.
    random.seed(20261013)
    np.random.seed(20261013)
    pdb = PDBFile(io.StringIO(source_text))
    modeller = Modeller(pdb.topology, pdb.positions)
    residues = list(modeller.topology.residues())
    aspartate = [i for i, r in enumerate(residues) if r.chain.id == "A" and r.id == "542"]
    histidine = [i for i, r in enumerate(residues) if r.chain.id == "A" and r.id == "600"]
    if len(aspartate) != 1 or len(histidine) != 1:
        raise RuntimeError("Asp542/His600 residue selection ambiguous")
    variants = [None] * len(residues)
    variants[aspartate[0]] = "ASH"
    variants[histidine[0]] = "HID"
    chosen = modeller.addHydrogens(pH=6.5, variants=variants,
                                   platform=openmm.Platform.getPlatformByName("Reference"))
    if chosen[aspartate[0]] != "ASH" or chosen[histidine[0]] != "HID":
        raise RuntimeError("OpenMM did not apply specified catalytic variants")
    coordinates = {atom.name: np.array(position.value_in_unit(unit.angstrom), dtype=float)
                   for atom, position in zip(modeller.topology.atoms(), modeller.positions)
                   if atom.residue.chain.id == "A" and atom.residue.id == "542"}
    if not {"CG", "OD1", "OD2", "HD2"} <= set(coordinates):
        raise RuntimeError("OpenMM ASH atoms incomplete")
    oxygen = coordinates["OD2"]
    direction = coordinates["HD2"] - oxygen
    hydrogen = oxygen + 0.98 * direction / np.linalg.norm(direction)
    axis = coordinates["CG"] - oxygen
    angle = float(np.degrees(np.arccos(np.clip(np.dot(axis, hydrogen-oxygen) / np.linalg.norm(axis) / 0.98, -1, 1))))
    if not 95.0 <= angle <= 125.0:
        raise RuntimeError(f"Asp542 CG-OD2-HD2 angle implausible: {angle}")
    return hydrogen, angle


def _pdb_variant(source_text: str, hydrogen: np.ndarray) -> str:
    lines = source_text.splitlines()
    if any(line.startswith("ATOM  ") and line[21] == "A" and line[22:26].strip() == "542"
           and line[12:16].strip() in ("HD1", "HD2") for line in lines):
        raise RuntimeError("Source Asp542 already has acid hydrogen")
    serial = max(int(line[6:11]) for line in atom_lines(source_text)) + 1
    if serial > 99999:
        raise RuntimeError("PDB serial overflow")
    out = ["REMARK 900 3L4W NO-NAG DRY SENSITIVITY D443-/D542H(OD2)/H600-HID PH 6.5",
           "REMARK 901 OPENMM 8.6.1 LIGAND-FREE ASH DIRECTION; OH NORMALIZED 0.98 A; HEAVY FIXED"]
    inserted = 0
    for line in lines:
        if line == "END":
            continue
        out.append(line)
        if (line.startswith("ATOM  ") and line[21] == "A" and line[22:26].strip() == "542"
                and line[12:16].strip() == "HB3"):
            new = (line[:6] + f"{serial:5d}" + " " + " HD2" + line[16:30]
                   + "".join(f"{float(v):8.3f}" for v in hydrogen)
                   + f"{1.00:6.2f}{0.00:6.2f}" + "           H  ")
            out.append(new)
            inserted += 1
    out.append("END")
    result = "\n".join(out) + "\n"
    if inserted != 1 or heavy_map(result) != heavy_map(source_text):
        raise RuntimeError("Asp542 variant did not preserve all deposited heavy atoms")
    old_lines = {int(line[6:11]): line for line in atom_lines(source_text)}
    new_lines = {int(line[6:11]): line for line in atom_lines(result)}
    if len(new_lines) != len(old_lines)+1 or any(new_lines.get(s) != line for s, line in old_lines.items()):
        raise RuntimeError("Existing receptor atoms changed in sensitivity PDB")
    return result


def _contacts(source_text: str, hydrogen: np.ndarray):
    protein = []
    for line in atom_lines(source_text):
        if line[21] == "A" and line[22:26].strip() == "542":
            continue
        protein.append((float(np.linalg.norm(pos(line)-hydrogen)), line[76:78].strip()))
    if not protein:
        raise RuntimeError("Empty protein contact screen")
    min_protein = min(d for d, _ in protein)
    if min_protein < 1.4:
        raise RuntimeError(f"New Asp542 hydrogen clashes with protein: {min_protein:.3f} A")
    min_ligand = {}
    for name, path in (("neutral", base.REFERENCE), ("charged", sens.CHARGED_SDF)):
        mol = Chem.SDMolSupplier(str(path), removeHs=False)[0]
        if mol is None:
            raise RuntimeError(f"Could not read {path}")
        xyz = np.asarray(mol.GetConformer().GetPositions())
        distances = np.linalg.norm(xyz-hydrogen, axis=1)
        min_ligand[name] = float(np.min(distances))
        if min_ligand[name] < 1.2:
            raise RuntimeError(f"New Asp542 H clashes with {name} MIG modeled hydrogens")
    return min_protein, min_ligand


def prepare_receptor():
    base.verify_hashes()
    if openmm.__version__ != "8.6.1":
        raise RuntimeError(f"OpenMM version changed: {openmm.__version__}")
    if not CLI.is_file() or not TEMPLATE.is_file():
        raise RuntimeError("Meeko receptor CLI/template unavailable")
    template = json.loads(TEMPLATE.read_text())["residue_templates"]["ASH"]
    if "HD2" not in template["atom_name"] or "OD2" not in template["atom_name"]:
        raise RuntimeError("ASH template lacks OD2 hydroxyl")
    source = SOURCE_PDB.read_text()
    hydrogen, angle = _variant_proton(source)
    model = _pdb_variant(source, hydrogen)
    min_protein, min_ligand = _contacts(source, hydrogen)
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        temp_pdb = tmp / "variant.pdb"
        temp_pdbqt = tmp / "variant.pdbqt"
        temp_pdb.write_text(model)
        cmd = [str(CLI), "--read_pdb", str(temp_pdb), "--write_pdbqt", str(temp_pdbqt),
               "--set_template", "A:542=ASH"]
        run = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if run.returncode or not temp_pdbqt.is_file():
            raise RuntimeError(f"Meeko ASH conversion failed: {run.stdout[-800:]} {run.stderr[-800:]}")
        receptor = PDBQTReceptor.from_pdbqt_filename(str(temp_pdbqt))
        atoms = receptor.atoms()
        if not len(atoms) or any(not np.isfinite(float(a["partial_charges"])) for a in atoms):
            raise RuntimeError("Receptor PDBQT empty or nonfinite charges")
        protein_heavy = Counter(tuple(round(float(v), 3) for v in pos(line)) for line in atom_lines(model)
                                if line[76:78].strip() not in ("H", "D"))
        qt_heavy = Counter(tuple(round(float(v), 3) for v in a["xyz"]) for a in atoms if a["atom_type"] != "HD")
        if qt_heavy != protein_heavy:
            raise RuntimeError(f"Receptor PDBQT heavy mismatch: expected {sum(protein_heavy.values())}, "
                               f"found {sum(qt_heavy.values())}; missing {(protein_heavy-qt_heavy).most_common(3)}, "
                               f"extra {(qt_heavy-protein_heavy).most_common(3)}")
        d542 = [a for a in atoms if a["chain"] == "A" and a["resid"] == 542]
        typed = {str(a["name"]): str(a["atom_type"]) for a in d542}
        if typed.get("HD2") != "HD" or typed.get("OD2") != "OA" or typed.get("OD1") != "OA":
            raise RuntimeError(f"Asp542 atom typing incorrect: {typed}")
        charge = float(sum(float(a["partial_charges"]) for a in d542))
        if abs(charge) > 0.05:
            raise RuntimeError(f"Protonated Asp542 charge differs from neutral: {charge:.3f}")
        def residue_atoms(number):
            return [a for a in atoms if a["chain"] == "A" and a["resid"] == number]
        d443 = residue_atoms(443)
        h600 = residue_atoms(600)
        q443 = sum(float(a["partial_charges"]) for a in d443)
        q600 = sum(float(a["partial_charges"]) for a in h600)
        h600types = {str(a["name"]): str(a["atom_type"]) for a in h600}
        if abs(q443+1.0) > 0.05 or abs(q600) > 0.05 or h600types.get("HD1") != "HD" or "HE2" in h600types:
            raise RuntimeError("Asp443-/His600-HID baseline states changed")
        original_atoms = PDBQTReceptor.from_pdbqt_filename(str(SOURCE_PDBQT)).atoms()
        atom_key = lambda a: (str(a["chain"]), int(a["resid"]), str(a["resname"]), str(a["name"]))
        old = {atom_key(a): a for a in original_atoms}
        new = {atom_key(a): a for a in atoms}
        if len(old) != len(original_atoms) or len(new) != len(atoms) or set(old)-set(new) or set(new)-set(old) != {("A",542,"ASP","HD2")}:
            raise RuntimeError("Receptor PDBQT atom identities changed outside added Asp542 HD2")
        for atom_id, prior in old.items():
            if atom_id[1] == 542:
                continue
            current = new[atom_id]
            if (str(prior["atom_type"]) != str(current["atom_type"])
                    or abs(float(prior["partial_charges"])-float(current["partial_charges"])) > 0.001
                    or np.max(np.abs(prior["xyz"]-current["xyz"])) > 0.001):
                raise RuntimeError(f"Receptor PDBQT atom changed outside Asp542: {atom_id}")
        pqt_text = temp_pdbqt.read_bytes()
    base.write_immutable(NEW_PDB, model.encode("ascii"))
    base.write_immutable(NEW_PDBQT, pqt_text)
    for ligand in (base.LIGAND, sens.CHARGED_PDBQT):
        probe = subprocess.run([str(base.VINA), "--receptor", str(NEW_PDBQT), "--ligand", str(ligand),
                                "--center_x", "45.38664286", "--center_y", "92.35157143", "--center_z", "34.822",
                                "--size_x", "0", "--size_y", "1", "--size_z", "1", "--score_only"],
                               cwd=ROOT, capture_output=True, text=True, timeout=30)
        if probe.returncode != 1 or "Grid box dimensions must be greater than 0 Angstrom" not in probe.stdout + probe.stderr:
            raise RuntimeError(f"Vina 1.2.7 receptor/ligand parse probe failed for {ligand.name}")
    return {"H_coordinate_A": hydrogen.tolist(), "CG_OD2_HD2_angle_deg": angle,
            "min_new_H_to_other_protein_A": min_protein, "min_new_H_to_MIG_A": min_ligand,
            "D542_PDBQT_charge": charge, "receptor_atom_count": len(atoms),
            "protein_heavy_atom_count": sum(protein_heavy.values()),
            "ASH_template_SHA256": base.digest(TEMPLATE),
            "D443_PDBQT_charge": q443, "H600_PDBQT_charge": q600,
            "Vina_parse_probe": "neutral and charged MIG reached deliberate zero-grid guard; no scoring/docking",
            "Meeko_command": "mk_prepare_receptor.exe --read_pdb <temporary variant.pdb> --write_pdbqt <temporary variant.pdbqt> --set_template A:542=ASH",
            "new_PDB_SHA256": base.digest(NEW_PDB), "new_PDBQT_SHA256": base.digest(NEW_PDBQT)}


def main():
    prepare_only = len(sys.argv) == 2 and sys.argv[1] == "--prepare-only"
    if len(sys.argv) > 1 and not prepare_only:
        raise RuntimeError("Only optional --prepare-only is supported")
    qc = prepare_receptor()
    print("Prepared and checked Asp542H(OD2) receptor", json.dumps(qc, indent=2), flush=True)
    if prepare_only:
        return 0
    with (ROOT / "data/metadata/vina_3L4W_sensitivity_results.csv").open(newline="", encoding="utf-8") as f:
        baseline = {(r["ligand_microstate"], int(r["seed"])): r for r in csv.DictReader(f)}
    if len(baseline) != 10 or any(r["status"] != "completed" for r in baseline.values()):
        raise RuntimeError("Ten baseline E64 runs not complete")
    ref_mol, heavy, xyz, _, _, automorphisms = base.get_reference()
    signs = base.stereo_signs(heavy, xyz)
    center, dims, _ = base.derive_box(xyz)
    cases = (("neutral", base.LIGAND, base.REFERENCE),
             ("amine +1", sens.CHARGED_PDBQT, sens.CHARGED_SDF))
    manifest = {"created_UTC": datetime.now(timezone.utc).isoformat(),
                "experiment": "3L4W D542H(OD2) dry no-NAG receptor microstate sensitivity",
                "receptor_model": str(NEW_PDBQT.relative_to(ROOT)), "receptor_QC": qc,
                "baseline_receptor_SHA256": base.digest(SOURCE_PDBQT),
                "native_reference_SHA256": base.digest(base.REFERENCE),
                "Vina_executable_SHA256": base.digest(base.VINA), "Vina_version": "1.2.7",
                "OpenMM_version": openmm.__version__, "Meeko_version": __import__("meeko").__version__,
                "RDKit_version": Chem.rdBase.rdkitVersion, "seeds": list(SEEDS),
                "box_center_A": center.tolist(), "box_size_A": list(dims),
                "exhaustiveness": 64, "num_modes": 20, "energy_range_kcal_mol": 5.0,
                "min_rmsd_A": 1.0, "cpu": 1,
                "model_limit": "No NAG, no waters, Asp542 OH orientation ligand-free provisional; catalytic chemistry not experimentally assigned"}
    OUT.mkdir(parents=True, exist_ok=True)
    manifest_path = OUT / "experiment_manifest.json"
    if manifest_path.exists():
        manifest["created_UTC"] = json.loads(manifest_path.read_text())["created_UTC"]
    base.write_immutable(manifest_path, (json.dumps(manifest, indent=2)+"\n").encode())
    rows, all_poses = [], []
    for microstate, ligand, sdf in cases:
        records, mapping, partial, ligand_mol = sens.ligand_case(ligand, sdf, xyz, heavy, signs)
        target_smiles = sens.smiles_without_nitrogen_tag(Chem.RemoveHs(ligand_mol))
        for seed in SEEDS:
            tag = "neutral" if microstate == "neutral" else "NHplus"
            prefix = f"D542H_OD2_{tag}_seed_{seed}"
            start = OUT / f"{prefix}_start.pdbqt"
            output = OUT / f"{prefix}_docked.pdbqt"
            log = OUT / f"{prefix}_vina.log"
            start_text, start_rmsd, shift, attempt = sens.random_start(seed, ligand.read_text(), records, mapping, heavy, xyz, center, dims, automorphisms, signs, target_smiles, ligand_mol)
            base.write_immutable(start, start_text.encode("ascii"))
            cmd = [str(base.VINA), "--receptor", str(NEW_PDBQT), "--ligand", str(start),
                   "--center_x", f"{center[0]:.8f}", "--center_y", f"{center[1]:.8f}", "--center_z", f"{center[2]:.8f}",
                   "--size_x", str(dims[0]), "--size_y", str(dims[1]), "--size_z", str(dims[2]),
                   "--scoring", "vina", "--exhaustiveness", "64", "--num_modes", "20", "--energy_range", "5.0",
                   "--min_rmsd", "1.0", "--seed", str(seed), "--cpu", "1", "--out", str(output), "--verbosity", "1"]
            command = subprocess.list2cmdline(cmd)
            if output.exists() and log.exists():
                stdout = log.read_text()
                code = re.search(r"(?m)^EXIT_CODE: (\d+)$", stdout)
                elapsed_match = re.search(r"(?m)^ELAPSED_SECONDS: ([\d.]+)$", stdout)
                if command not in stdout or code is None or elapsed_match is None:
                    raise RuntimeError(f"Existing run metadata incomplete: {prefix}")
                exit_code, elapsed = int(code.group(1)), float(elapsed_match.group(1))
                print("Reanalysing", prefix, flush=True)
            elif output.exists() or log.exists():
                raise RuntimeError(f"Partial existing run: {prefix}")
            else:
                print("Docking", prefix, "start RMSD", round(start_rmsd,3), flush=True)
                t0 = time.monotonic()
                run = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=3600)
                elapsed, exit_code = time.monotonic()-t0, run.returncode
                stdout = (f"COMMAND: {command}\nSEED: {seed}\nRECEPTOR_SHA256: {base.digest(NEW_PDBQT)}\nLIGAND_SHA256: {base.digest(ligand)}\nSTART_SHA256: {base.digest(start)}\nEXIT_CODE: {exit_code}\nELAPSED_SECONDS: {elapsed:.3f}\n--- STDOUT ---\n{run.stdout}\n--- STDERR ---\n{run.stderr}\n")
                base.write_immutable(log, stdout.encode())
            previous = baseline[(microstate, seed)]
            row = {"PDB_ID": "3L4W", "receptor_variant": "D443-/D542H(OD2)/H600-HID, dry, no NAG", "ligand_microstate": microstate,
                   "seed": seed, "status": "failed", "vina_version": "1.2.7", "exhaustiveness": 64,
                   "receptor_path": str(NEW_PDBQT.relative_to(ROOT)), "receptor_SHA256": base.digest(NEW_PDBQT),
                   "ligand_path": str(ligand.relative_to(ROOT)), "ligand_SHA256": base.digest(ligand),
                   "reference_path": str(base.REFERENCE.relative_to(ROOT)), "reference_SHA256": base.digest(base.REFERENCE),
                   "start_path": str(start.relative_to(ROOT)), "start_SHA256": base.digest(start), "start_RMSD_A": f"{start_rmsd:.4f}",
                   "box_center_xyz_A": ";".join(f"{v:.8f}" for v in center), "box_size_xyz_A": ";".join(map(str,dims)),
                   "num_modes_requested": 20, "energy_range_kcal_mol": 5.0, "min_rmsd_A": 1.0, "cpu": 1,
                   "command": command, "log_path": str(log.relative_to(ROOT)), "log_SHA256": base.digest(log),
                   "output_path": str(output.relative_to(ROOT)) if output.exists() else "",
                   "output_SHA256": base.digest(output) if output.exists() else "", "exit_code": exit_code,
                   "elapsed_seconds": f"{elapsed:.3f}", "n_modes_returned": "", "top_score_kcal_mol": "",
                   "top_RMSD_A": "", "best_RMSD_A": "", "best_RMSD_rank": "", "best_RMSD_score_kcal_mol": "",
                   "baseline_E64_top_RMSD_A": previous["top_RMSD_A"], "baseline_E64_best_RMSD_A": previous["best_RMSD_A"],
                   "baseline_E64_top_score_kcal_mol": previous["top_score_kcal_mol"], "top_within_2A": "", "any_within_2A": "", "error": ""}
            try:
                if exit_code != 0 or not output.is_file() or output.stat().st_size == 0:
                    raise RuntimeError(f"Vina exit {exit_code}; inspect log")
                if f"Performing docking (random seed: {seed}" not in stdout:
                    raise RuntimeError("Vina log did not confirm random seed")
                poses = (base.parse_output(output, mapping, records, heavy, xyz, automorphisms, signs)
                         if microstate == "neutral" else sens.parse_charged_output(output, mapping, records, heavy, xyz, automorphisms, signs, ligand_mol))
                if not poses:
                    raise RuntimeError("No chemically validated modes")
                top, best = poses[0], min(poses, key=lambda p: p["heavy_atom_RMSD_A"])
                row.update(status="completed", n_modes_returned=len(poses), top_score_kcal_mol=f"{top['score_kcal_mol']:.3f}",
                           top_RMSD_A=f"{top['heavy_atom_RMSD_A']:.4f}", best_RMSD_A=f"{best['heavy_atom_RMSD_A']:.4f}",
                           best_RMSD_rank=best["pose_rank"], best_RMSD_score_kcal_mol=f"{best['score_kcal_mol']:.3f}",
                           top_within_2A=top["heavy_atom_RMSD_A"] <= 2.0, any_within_2A=best["heavy_atom_RMSD_A"] <= 2.0)
                saved = [{"ligand_microstate": microstate, "seed": seed, **p, "within_2A": p["heavy_atom_RMSD_A"] <= 2.0} for p in poses]
                all_poses.extend(saved)
                base.write_immutable(OUT / f"{prefix}_pose_metrics.csv", base.csv_bytes(saved,
                                     ["ligand_microstate", "seed", "pose_rank", "score_kcal_mol", "heavy_atom_RMSD_A", "within_2A"]))
                print(prefix, "top", row["top_RMSD_A"], "best", row["best_RMSD_A"], "rank", row["best_RMSD_rank"], flush=True)
            except Exception as exc:
                row["error"] = str(exc)
                print("FAILED", prefix, exc, flush=True)
            rows.append(row)
            if row["status"] != "completed":
                break
        if rows[-1]["status"] != "completed":
            break
    if rows:
        base.write_immutable(RESULTS, base.csv_bytes(rows, list(rows[0])))
    if all_poses:
        base.write_immutable(OUT / "all_pose_metrics.csv", base.csv_bytes(all_poses,
                             ["ligand_microstate", "seed", "pose_rank", "score_kcal_mol", "heavy_atom_RMSD_A", "within_2A"]))
    return 0 if len(rows) == 10 and all(row["status"] == "completed" for row in rows) else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"FATAL PREFLIGHT/PREPARATION ERROR: {exc}", file=sys.stderr)
        sys.exit(2)

"""Freeze and validate neutral CID 442416 agnuside inputs; no Vina calls.

Run from the project .venv. Source structures are read-only. Derived SDF,
PDBQT, randomized starts, and the CSV are write-once and hash verified.
"""

from __future__ import annotations

import csv
import io
import math
import subprocess
import sys
from pathlib import Path

import meeko
import numpy as np
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolAlign, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import validate_vina_3L4W_pilot_inputs as common  # noqa: E402
import validate_vina_3L4W_gallic_acid_inputs as gallic  # noqa: E402

SOURCE = ROOT / "structures/ligands/prepared/LIG09.sdf"
ORIGINAL_PDBQT = ROOT / "structures/ligands/vina_prepared/LIG09_CID442416_neutral_reference.pdbqt"
STEP65 = ROOT / "data/metadata/vina_remaining_plants_input_manifest.csv"
OUT = ROOT / "data/metadata/vina_3L4W_agnuside_pilot_manifest.csv"
DIR = ROOT / "structures/ligands/vina_prepared/agnuside_3L4W"
START_DIR = ROOT / "docking/pilot/vina/3L4W_agnuside_pose_only/frozen_starts"
CID = 442416
SEEDS = (1101, 2202, 3303, 4404, 5505)
CORE = tuple(range(9))  # zero-based SDF heavy atoms of the fused iridoid bicycle
FLEX = tuple(range(9, 33))
BOX_CENTER = np.array([45.38664286, 92.35157143, 34.82200000])
BOX_SIZE = np.array([30.0, 24.0, 24.0])
ENERGY_WINDOW = 12.0
MIN_HEAVY_RMSD = 1.5
MIN_FLEX_RMSD = 2.0
EXPECTED_RECORDS = (1, 5, 8, 11, 17)
FIELDS = [
    "structure_id", "compound_name", "association_id", "plant_species", "PubChem_CID",
    "botanical_evidence_tier", "pilot_decision", "microstate", "source_sdf_path",
    "source_sdf_SHA256", "source_record_count", "source_record_1based", "source_conformer_number",
    "source_MMFF94s_energy_kcal_mol", "delta_from_min_MMFF94s_kcal_mol",
    "selection_rank", "min_fitted_full_heavy_RMSD_to_prior_A",
    "min_flexible_RMSD_after_core_alignment_to_prior_A",
    "energy_window_kcal_mol", "full_heavy_RMSD_threshold_A", "flex_RMSD_threshold_A",
    "ligand_sdf_path", "ligand_sdf_SHA256", "ligand_pdbqt_path", "ligand_pdbqt_SHA256",
    "original_pdbqt_path", "original_pdbqt_SHA256", "isomeric_SMILES", "formula", "formal_charge",
    "heavy_atoms", "defined_stereocentres", "pdbqt_torsions", "pdbqt_charge_sum",
    "pdbqt_readback", "max_heavy_coordinate_shift_A", "heavy_atom_map_sdf_to_pdbqt",
    "rigid_core_sdf_indices_1based", "rigid_core_pdbqt_serials_1based", "core_automorphisms",
    "heavy_centroid_radius_A", "centered_min_box_clearance_A", "seed", "start_path",
    "start_SHA256", "start_attempt", "start_translation_A", "start_quaternion_wxyz",
    "start_min_box_margin_A", "start_graph_stereo_check", "receptor_state", "receptor_pdbqt_path",
    "receptor_pdbqt_SHA256", "receptor_source_pdb_path", "receptor_source_pdb_SHA256",
    "native_MIG_sdf_path", "native_MIG_sdf_SHA256", "native_MIG_control_path",
    "native_MIG_control_SHA256", "native_MIG_control_status", "Vina_executable_SHA256",
    "Vina_version", "box_center_A", "box_size_A", "scoring", "exhaustiveness", "num_modes",
    "energy_range_kcal_mol", "min_rmsd_A", "cpu", "pocket_rule", "clash_rule",
    "rank1_cluster_rule", "box_face_diagnostic", "deposited_water_diagnostic",
    "failure_policy", "software_versions", "limitations",
]


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def assert_reference(mol: Chem.Mol) -> None:
    if mol is None or mol.GetProp("structure_id") != "LIG09" or mol.GetProp("PubChem_CID") != str(CID):
        raise ValueError("LIG09 record identity mismatch")
    if common.iso(mol) != REFERENCE_SMILES:
        raise ValueError("Agnuside reference graph or stereochemistry differs")
    inferred = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(inferred)
    if common.iso(inferred) != REFERENCE_SMILES:
        raise ValueError("Agnuside 3D-implied stereochemistry differs")
    if rdMolDescriptors.CalcMolFormula(mol) != "C22H26O11" or Chem.GetFormalCharge(mol) != 0:
        raise ValueError("Agnuside formula or charge mismatch")
    if sum(a.GetAtomicNum() > 1 for a in mol.GetAtoms()) != 33 or len(Chem.FindMolChiralCenters(Chem.RemoveHs(mol), includeUnassigned=False)) != 9:
        raise ValueError("Agnuside heavy atoms or nine stereocentres mismatch")
    common.check_3d(mol)
    xyz = mol.GetConformer().GetPositions()
    graph = Chem.GetDistanceMatrix(mol)
    for i in range(33):
        for j in range(i + 1, 33):
            if graph[i, j] >= 4 and np.linalg.norm(xyz[i] - xyz[j]) < 1.6:
                raise ValueError("Agnuside source internal heavy-atom clash")
    if not mol.HasProp("mmff_energy_kcal_mol") or not math.isfinite(float(mol.GetProp("mmff_energy_kcal_mol"))):
        raise ValueError("Agnuside source MMFF energy absent/nonfinite")


def aligned_flex_rmsd(mobile: Chem.Mol, reference: Chem.Mol, automorphisms: tuple) -> float:
    x = mobile.GetConformer().GetPositions()[:33]
    y = reference.GetConformer().GetPositions()[:33]
    best = math.inf
    for mapping in automorphisms:
        if {mapping[i] for i in CORE} != set(CORE):
            continue
        mapped = x[list(mapping)]
        cx = mapped[list(CORE)].mean(axis=0)
        cy = y[list(CORE)].mean(axis=0)
        u, _, vt = np.linalg.svd((mapped[list(CORE)] - cx).T @ (y[list(CORE)] - cy))
        correction = np.eye(3)
        correction[2, 2] = np.linalg.det(u @ vt)
        fitted = (mapped - cx) @ (u @ correction @ vt) + cy
        rmsd = float(np.sqrt(np.mean(np.sum((fitted[list(FLEX)] - y[list(FLEX)]) ** 2, axis=1))))
        best = min(best, rmsd)
    if not math.isfinite(best):
        raise ValueError("No graph automorphism preserves iridoid rigid core")
    return best


def select_records(mols: list[Chem.Mol]) -> tuple[list[int], dict[int, tuple[float, float]], tuple]:
    energies = [float(m.GetProp("mmff_energy_kcal_mol")) for m in mols]
    minimum = min(energies)
    eligible = [i for i, e in enumerate(energies) if e <= minimum + ENERGY_WINDOW + 1e-8]
    first = min(range(len(mols)), key=lambda i: (energies[i], i))
    selected = [first]
    metrics = {first: (math.nan, math.nan)}
    heavy = Chem.RemoveHs(mols[0])
    automorphisms = heavy.GetSubstructMatches(heavy, uniquify=False)
    if len(automorphisms) != 2 or len(CORE) != 9:
        raise ValueError("Agnuside graph automorphism/core definition changed")
    rings = [set(r) for r in mols[0].GetRingInfo().AtomRings()]
    if not any(r == set(CORE[:6]) for r in rings) or not any(r == {4, 5, 6, 7, 8} for r in rings):
        raise ValueError("Source atoms 1-9 no longer define the fused iridoid bicycle")
    while len(selected) < 5:
        candidates = []
        for i in eligible:
            if i in selected:
                continue
            min_full = min(rdMolAlign.GetBestRMS(Chem.RemoveHs(mols[i]), Chem.RemoveHs(mols[j])) for j in selected)
            min_flex = min(aligned_flex_rmsd(mols[i], mols[j], automorphisms) for j in selected)
            if min_full >= MIN_HEAVY_RMSD and min_flex >= MIN_FLEX_RMSD:
                candidates.append((min_flex, min_full, energies[i], i))
        if not candidates:
            raise ValueError("Insufficient qualified diverse agnuside conformers")
        choice = max(candidates, key=lambda v: (round(v[0], 6), round(v[1], 6), -v[2], -v[3]))
        selected.append(choice[3])
        metrics[choice[3]] = choice[1], choice[0]
    if tuple(i + 1 for i in selected) != EXPECTED_RECORDS:
        raise ValueError(f"Deterministic selection drift: {[i+1 for i in selected]}")
    return selected, metrics, automorphisms


def one_record_sdf(mol: Chem.Mol, record: int) -> str:
    copy = Chem.Mol(mol)
    copy.SetProp("selected_parent_record_1based", str(record))
    copy.SetProp("parent_source_sdf_sha256", common.sha(SOURCE))
    copy.SetProp("preparation_scope", "neutral_PubChem_reference_only_not_original_extract_confirmation")
    return (Chem.MolToMolBlock(copy) + "\n" +
            "".join(f">  <{k}>\n{copy.GetProp(k)}\n\n" for k in sorted(copy.GetPropNames())) + "$$$$\n")


def prepare_record(mol: Chem.Mol, record: int) -> tuple[Path, Path, dict]:
    sdf = DIR / f"LIG09_CID442416_neutral_record_{record:04d}.sdf"
    pdbqt = sdf.with_suffix(".pdbqt")
    common.write_once(sdf, one_record_sdf(mol, record))
    readback = list(Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=True))
    if len(readback) != 1 or readback[0] is None:
        raise ValueError(f"Record {record}: SDF roundtrip failed")
    assert_reference(readback[0])
    setup = MoleculePreparation(charge_model="gasteiger", merge_these_atom_types=("H",),
                                add_index_map=True, hydrate=False, rigid_macrocycles=False,
                                flexible_amides=False).prepare(mol)
    if len(setup) != 1:
        raise ValueError(f"Record {record}: Meeko did not return one setup")
    text, ok, error = PDBQTWriterLegacy.write_string(setup[0], add_index_map=True)
    if not ok:
        raise ValueError(f"Record {record}: Meeko write failed: {error}")
    common.write_once(pdbqt, text)
    qc = common.check_pdbqt(pdbqt, mol)
    if qc["heavy"] != 33 or qc["torsions"] != 13 or abs(qc["charge"]) > 0.02:
        raise ValueError(f"Record {record}: PDBQT atom/torsion/charge mismatch")
    return sdf, pdbqt, qc


def frozen_random_start(mol: Chem.Mol, pdbqt: Path, qc: dict, seed: int,
                        record: int) -> tuple[Path, dict]:
    text = pdbqt.read_text(encoding="utf-8")
    atoms = {}
    for line in text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            serial = int(line[6:11])
            atoms[serial] = np.array([float(line[i:i+8]) for i in (30, 38, 46)])
    serials = sorted(atoms)
    original = np.array([atoms[s] for s in serials])
    heavy = mol.GetConformer().GetPositions()[:33]
    center = heavy.mean(axis=0)
    rng = np.random.default_rng(np.random.SeedSequence([seed, CID, record, 0]))
    for attempt in range(1, 10001):
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
        transformed = (original - center) @ rot.T + BOX_CENTER + shift
        point = dict(zip(serials, transformed))
        mapped = np.array([point[qc["map"][i]] for i in range(1, 34)])
        margin = float(np.min(BOX_SIZE / 2 - np.abs(mapped - BOX_CENTER)))
        if margin < 0.5:
            continue
        lines = []
        for line in text.splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                p = point[int(line[6:11])]
                line = line[:30] + "".join(f"{v:8.3f}" for v in p) + line[54:]
            lines.append(line)
        data = "\n".join(lines) + "\n"
        start = START_DIR / f"LIG09_neutral_record_{record:04d}_seed_{seed}_start.pdbqt"
        common.write_once(start, data)
        parsed = PDBQTMolecule(start.read_text(encoding="utf-8"))
        rebuilt = RDKitMolCreate.from_pdbqt_mol(parsed)
        if len(rebuilt) != 1 or rebuilt[0] is None or common.iso(rebuilt[0]) != REFERENCE_SMILES:
            raise ValueError(f"Seed {seed}: start graph/stereochemistry readback failed")
        coords = parsed.atoms()["xyz"]
        backmapped = np.array([coords[qc["map"][i]-1] for i in range(1, 34)])
        if np.max(np.abs(np.linalg.norm(backmapped[:, None, :] - backmapped[None, :, :], axis=2) -
                         np.linalg.norm(heavy[:, None, :] - heavy[None, :, :], axis=2))) > 0.002:
            raise ValueError(f"Seed {seed}: rigid randomization distorted intraligand geometry")
        actual_margin = float(np.min(BOX_SIZE / 2 - np.abs(backmapped - BOX_CENTER)))
        if actual_margin < 0.499:
            raise ValueError(f"Seed {seed}: PDBQT rounding crossed box-clearance threshold")
        return start, {"attempt": attempt, "translation": shift, "quaternion": q,
                       "margin": actual_margin}
    raise ValueError(f"Could not place selected record {record} safely inside frozen box")


def main() -> None:
    check = subprocess.run([sys.executable, str(ROOT / "scripts/validate_vina_remaining_plants_inputs.py")],
                           cwd=ROOT, capture_output=True, text=True, timeout=180)
    if check.returncode:
        raise ValueError(f"Step 65 validator failed:\n{check.stdout}\n{check.stderr}")
    gallic.confirm_controls()
    common.check_hash(SOURCE, "5fa53e66be52e08ea9e6561cdd83ce6eb6c1308f31439132e109643395f25c42")
    common.check_hash(ORIGINAL_PDBQT, "2b370df8fb89d528ce0243fce008a2d6e240f3999802ad0735f3aead2fd86875")
    source_rows = [r for r in common.rows(rel(STEP65)) if r["structure_id"] == "LIG09"]
    if len(source_rows) != 2 or any(r["candidate_decision"] != "PASS_FOR_POSE_ONLY_PLANNING" or r["microstate"] != "neutral_agnuside_reference" for r in source_rows):
        raise ValueError("Step 65 agnuside candidate decision changed")
    library = [r for r in common.rows("data/compounds/final_preliminary_ligand_library.csv") if r["structure_id"] == "LIG09"]
    if len(library) != 1 or library[0]["association_id"] != "VN03" or library[0]["plant_species"] != "Vitex negundo":
        raise ValueError("LIG09 botanical association changed")
    global REFERENCE_SMILES
    REFERENCE_SMILES = common.iso(Chem.MolFromSmiles(library[0]["isomeric_SMILES"]))
    mols = list(Chem.SDMolSupplier(str(SOURCE), removeHs=False, sanitize=True))
    if len(mols) != 60:
        raise ValueError("LIG09 source record count changed")
    for m in mols:
        assert_reference(m)
    selected, metrics, automorphisms = select_records(mols)
    if common.check_pdbqt(ORIGINAL_PDBQT, mols[0])["torsions"] != 13:
        raise ValueError("Original neutral agnuside PDBQT changed")
    controls_sha = common.sha(gallic.MIG_RESULTS)
    generated = {}
    for rank, i in enumerate(selected, start=1):
        record = i + 1
        sdf, pdbqt, qc = prepare_record(mols[i], record)
        seed = SEEDS[rank - 1]
        start, randomization = frozen_random_start(mols[i], pdbqt, qc, seed, record)
        generated[record] = (sdf, pdbqt, qc, start, randomization, seed, rank)
    rows = []
    minimum = min(float(m.GetProp("mmff_energy_kcal_mol")) for m in mols)
    for i in selected:
        record = i + 1
        mol = mols[i]
        sdf, pdbqt, qc, start, randomization, seed, rank = generated[record]
        energy = float(mol.GetProp("mmff_energy_kcal_mol"))
        xyz = mol.GetConformer().GetPositions()[:33]
        radius = float(np.max(np.linalg.norm(xyz - xyz.mean(axis=0), axis=1)))
        full, flex = metrics[i]
        mapping = qc["map"]
        for receptor_state, (receptor_path, receptor_sha) in common.RECEPTORS.items():
            rows.append({
                "structure_id": "LIG09", "compound_name": "Agnuside", "association_id": "VN03",
                "plant_species": "Vitex negundo", "PubChem_CID": CID,
                "botanical_evidence_tier": "Tier B leaf methanol-first; direct water-only evidence abstract-level only",
                "pilot_decision": "READY_CONDITIONAL_REFERENCE_POSE_ONLY",
                "microstate": "neutral_agnuside_reference", "source_sdf_path": rel(SOURCE),
                "source_sdf_SHA256": common.sha(SOURCE), "source_record_count": 60,
                "source_record_1based": record, "source_conformer_number": mol.GetProp("conformer_number"),
                "source_MMFF94s_energy_kcal_mol": f"{energy:.8f}",
                "delta_from_min_MMFF94s_kcal_mol": f"{energy-minimum:.8f}",
                "selection_rank": rank,
                "min_fitted_full_heavy_RMSD_to_prior_A": "" if rank == 1 else f"{full:.4f}",
                "min_flexible_RMSD_after_core_alignment_to_prior_A": "" if rank == 1 else f"{flex:.4f}",
                "energy_window_kcal_mol": ENERGY_WINDOW, "full_heavy_RMSD_threshold_A": MIN_HEAVY_RMSD,
                "flex_RMSD_threshold_A": MIN_FLEX_RMSD,
                "ligand_sdf_path": rel(sdf), "ligand_sdf_SHA256": common.sha(sdf),
                "ligand_pdbqt_path": rel(pdbqt), "ligand_pdbqt_SHA256": common.sha(pdbqt),
                "original_pdbqt_path": rel(ORIGINAL_PDBQT), "original_pdbqt_SHA256": common.sha(ORIGINAL_PDBQT),
                "isomeric_SMILES": REFERENCE_SMILES, "formula": "C22H26O11", "formal_charge": 0,
                "heavy_atoms": 33, "defined_stereocentres": 9, "pdbqt_torsions": qc["torsions"],
                "pdbqt_charge_sum": f"{qc['charge']:.3f}", "pdbqt_readback": "PASS Meeko/RDKit exact graph and stereo",
                "max_heavy_coordinate_shift_A": f"{qc['shift']:.4f}",
                "heavy_atom_map_sdf_to_pdbqt": ";".join(f"{j}:{mapping[j]}" for j in range(1, 34)),
                "rigid_core_sdf_indices_1based": "1;2;3;4;5;6;7;8;9",
                "rigid_core_pdbqt_serials_1based": ";".join(str(mapping[j]) for j in range(1, 10)),
                "core_automorphisms": len(automorphisms), "heavy_centroid_radius_A": f"{radius:.4f}",
                "centered_min_box_clearance_A": f"{12-radius:.4f}", "seed": seed,
                "start_path": rel(start), "start_SHA256": common.sha(start),
                "start_attempt": randomization["attempt"],
                "start_translation_A": ";".join(f"{v:.6f}" for v in randomization["translation"]),
                "start_quaternion_wxyz": ";".join(f"{v:.9f}" for v in randomization["quaternion"]),
                "start_min_box_margin_A": f"{randomization['margin']:.4f}",
                "start_graph_stereo_check": "PASS Meeko/RDKit exact graph and stereo",
                "receptor_state": receptor_state, "receptor_pdbqt_path": receptor_path,
                "receptor_pdbqt_SHA256": receptor_sha,
                "receptor_source_pdb_path": common.RECEPTOR_PDB[receptor_state][0],
                "receptor_source_pdb_SHA256": common.RECEPTOR_PDB[receptor_state][1],
                "native_MIG_sdf_path": common.NATIVE[0], "native_MIG_sdf_SHA256": common.NATIVE[1],
                "native_MIG_control_path": rel(gallic.MIG_RESULTS), "native_MIG_control_SHA256": controls_sha,
                "native_MIG_control_status": "5/5 rank1 <=2 A per receptor in enlarged box; 10/10 total",
                "Vina_executable_SHA256": common.VINA[1], "Vina_version": "1.2.7",
                "box_center_A": "45.38664286;92.35157143;34.82200000", "box_size_A": "30;24;24",
                "scoring": "vina", "exhaustiveness": 64, "num_modes": 20,
                "energy_range_kcal_mol": 5.0, "min_rmsd_A": 1.0, "cpu": 1,
                "pocket_rule": "all heavy atoms inside box; core centroid <=6 A from MIG; >=50% core atoms <=4 A from MIG heavy",
                "clash_rule": "protein-heavy distance <1.8 A OR vdW overlap >0.8 A; nonbonded internal heavy distance <1.6 A",
                "rank1_cluster_rule": "unfitted symmetry-aware 9-heavy-core RMSD <=2 A to within-receptor rank1 medoid; >=4/5 seeds full geometric passes",
                "box_face_diagnostic": "record every face margin; flag any ligand heavy atom <1.0 A from a face",
                "deposited_water_diagnostic": "read-only compare all returned poses with original 3L4W water oxygens; flag any <1.8 A; no change to PASS rule",
                "failure_policy": "retain every run/attempt and exit code; retry technical memory failure only with identical frozen input/settings and new attempt log; unresolved exit is unconfirmed; invalid mode is failure, not omitted",
                "software_versions": f"RDKit {rdBase.rdkitVersion}; Meeko {meeko.__version__}; Vina 1.2.7; NumPy {np.__version__}",
                "limitations": "Reference PubChem stereochemistry not source-resolved; water-only leaf chromatogram unavailable; original extracts untested; 13 torsions; dry/no-NAG receptor; omitted waters and provisional protonation",
            })
    if len(rows) != 10 or len({(r["seed"], r["receptor_state"]) for r in rows}) != 10:
        raise ValueError("Expected five seeds for both receptor states")
    for seed in SEEDS:
        pair = [r for r in rows if r["seed"] == seed]
        if len(pair) != 2 or pair[0]["start_SHA256"] != pair[1]["start_SHA256"]:
            raise ValueError("Receptor states do not use the same frozen seed start")
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    common.write_once(OUT, stream.getvalue())
    print("READY_CONDITIONAL_REFERENCE_POSE_ONLY; selected records", ",".join(map(str, EXPECTED_RECORDS)))
    print("Five selected conformers; ten receptor-seed input pairs; no Vina scoring/docking")
    for record in EXPECTED_RECORDS:
        sdf, pdbqt, qc, start, info, seed, rank = generated[record]
        print(f"seed={seed} record={record} torsions={qc['torsions']} start_margin={info['margin']:.3f} A")
    print(rel(OUT), common.sha(OUT))


if __name__ == "__main__":
    main()

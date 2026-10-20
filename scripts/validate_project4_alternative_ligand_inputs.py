"""Step 72: audit reference ligands and freeze new, score-blind 3L4W starts.

Run with .venv/Scripts/python.exe. Existing structures/results are read-only.
Derived rutin conformers and starts are write-once; Vina is never called.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import sys
import urllib.request
from pathlib import Path

import meeko
import numpy as np
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolAlign, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import validate_vina_3L4W_pilot_inputs as common  # noqa: E402

OUT = ROOT / "data/metadata/project4_alternative_ligand_input_manifest.csv"
DERIVED = ROOT / "structures/ligands/vina_prepared/alternative_3L4W"
STARTS = ROOT / "docking/pilot/vina/3L4W_alternative_input_starts"
SEEDS = (1101, 2202, 3303, 4404, 5505)
CENTER = np.array([45.38664286, 92.35157143, 34.82200000])
SIZE = np.array([30.0, 24.0, 24.0])
CONTROL = ROOT / "data/metadata/vina_3L4W_large_box_mig_results.csv"
SPECS = {
    "LIG02": ("Aegle marmelos", "Rutin", 5280805, "C27H30O16", 60, 43, 10, 16),
    "LIG08": ("Vitex negundo", "Luteolin", 5280445, "C15H10O6", 24, 21, 0, 5),
    "LIG06": ("Boerhavia diffusa", "Quercetin", 5280343, "C15H10O7", 24, 22, 0, 6),
}
FIELDS = [
    "structure_id", "compound_name", "plant_species", "PubChem_CID", "association_scope",
    "botanical_tier", "microstate", "formal_charge", "decision", "decision_scope", "source_sdf",
    "source_sdf_sha256", "source_pdbqt", "source_pdbqt_sha256", "source_record_count",
    "selected_record_1based", "selection_rule", "MMFF94s_energy_kcal_mol", "energy_above_min_kcal_mol",
    "source_formula", "isomeric_SMILES", "heavy_atoms", "stereocentres", "torsions",
    "qualified_sdf", "qualified_sdf_sha256", "qualified_pdbqt", "qualified_pdbqt_sha256",
    "pdbqt_charge_sum", "heavy_atom_map", "rigid_core_sdf_indices_1based", "rigid_core_pdbqt_serials",
    "core_rule", "start_seed", "start_pdbqt", "start_sha256", "start_box_margin_A",
    "centered_conformer_radius_A", "receptor_state", "receptor_pdbqt", "receptor_sha256",
    "MIG_reference_sha256", "MIG_control_results_sha256", "MIG_control_pass_count", "box_center_A",
    "box_size_A", "vina_version", "scoring", "exhaustiveness", "num_modes", "energy_range",
    "min_rmsd_A", "cpu", "source_evidence_DOI", "limitations",
]


def rel(p: Path) -> str:
    return p.relative_to(ROOT).as_posix()


def digest(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write_csv_once(path: Path, rows: list[dict]) -> None:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    common.write_once(path, stream.getvalue())


def official_pubchem() -> dict[int, dict]:
    cids = ",".join(str(SPECS[s][2]) for s in SPECS)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids}/property/MolecularFormula,IsomericSMILES/JSON"
    with urllib.request.urlopen(url, timeout=25) as handle:
        data = json.load(handle)
    return {int(x["CID"]): x for x in data["PropertyTable"]["Properties"]}


def reference_core(mol: Chem.Mol, sid: str) -> list[int]:
    if sid == "LIG06":
        return common.core_indices(mol, sid)
    rings = [r for r in mol.GetRingInfo().AtomRings()
             if len(r) == 6 and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    core = set().union(*map(set, rings))
    carbonyl = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 8 and any(
        b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(a).GetIdx() in core
        for b in a.GetBonds())]
    core.update(carbonyl)
    if len(rings) != 3 or len(core) < 15 or len(core) > 17:
        raise ValueError(f"{sid}: rigid flavonoid core is not the expected three-ring scaffold")
    return sorted(i + 1 for i in core)


def check_mol(mol: Chem.Mol, ref_smiles: str, formula: str, heavy: int, stereo: int) -> None:
    if mol is None or common.iso(mol) != ref_smiles:
        raise ValueError("Source molecular graph or stereochemistry differs from PubChem")
    inferred = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(inferred)
    if common.iso(inferred) != ref_smiles:
        raise ValueError("3D-implied stereochemistry differs from reference")
    if rdMolDescriptors.CalcMolFormula(mol) != formula or Chem.GetFormalCharge(mol) != 0:
        raise ValueError("Formula or formal charge mismatch")
    if sum(a.GetAtomicNum() > 1 for a in mol.GetAtoms()) != heavy:
        raise ValueError("Heavy-atom count mismatch")
    if len(Chem.FindMolChiralCenters(Chem.RemoveHs(mol), includeUnassigned=True)) != stereo:
        raise ValueError("Reference stereocentre count mismatch")
    common.check_3d(mol)
    xyz = mol.GetConformer().GetPositions()
    graph = Chem.GetDistanceMatrix(mol)
    for i in range(heavy):
        for j in range(i + 1, heavy):
            if graph[i, j] >= 4 and np.linalg.norm(xyz[i] - xyz[j]) < 1.6:
                raise ValueError("Severe nonbonded internal heavy-atom clash")
    if not mol.HasProp("mmff_energy_kcal_mol") or not math.isfinite(float(mol.GetProp("mmff_energy_kcal_mol"))):
        raise ValueError("Missing/nonfinite stored MMFF energy")


def select_rutin(pool: list[Chem.Mol]) -> list[int]:
    """Score-blind max-min heavy RMSD among conformers within 12 kcal/mol."""
    energy = [float(m.GetProp("mmff_energy_kcal_mol")) for m in pool]
    eligible = [i for i, e in enumerate(energy) if e <= min(energy) + 12.0]
    selected = [min(eligible, key=lambda i: (energy[i], i))]
    heavy = [Chem.RemoveHs(m) for m in pool]
    while len(selected) < 5:
        candidates = []
        for i in eligible:
            if i not in selected:
                d = min(rdMolAlign.GetBestRMS(heavy[i], heavy[j]) for j in selected)
                candidates.append((d, -energy[i], -i, i))
        if not candidates:
            raise ValueError("Rutin energy window has fewer than five conformers")
        best = max(candidates)
        if best[0] < 1.5:
            raise ValueError("Insufficient diverse rutin conformers at 1.5 A fitted RMSD")
        selected.append(best[3])
    return selected


def make_rutin_input(mol: Chem.Mol, record: int) -> tuple[Path, Path]:
    sdf = DERIVED / f"LIG02_CID5280805_neutral_record_{record:04d}.sdf"
    pdbqt = sdf.with_suffix(".pdbqt")
    copy = Chem.Mol(mol)
    copy.SetProp("selected_parent_record_1based", str(record))
    copy.SetProp("parent_source_sdf_sha256", digest(ROOT / "structures/ligands/prepared/LIG02.sdf"))
    body = Chem.MolToMolBlock(copy) + "\n" + "".join(
        f">  <{k}>\n{copy.GetProp(k)}\n\n" for k in sorted(copy.GetPropNames())) + "$$$$\n"
    common.write_once(sdf, body)
    reread = list(Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=True))
    if len(reread) != 1 or reread[0] is None or common.iso(reread[0]) != common.iso(mol):
        raise ValueError("Derived rutin SDF graph/stereochemistry readback failed")
    if not np.allclose(reread[0].GetConformer().GetPositions(), mol.GetConformer().GetPositions(), atol=0.0001):
        raise ValueError("Derived rutin SDF shifted parent conformer coordinates")
    setup = MoleculePreparation(charge_model="gasteiger", merge_these_atom_types=("H",),
                                add_index_map=True, hydrate=False, rigid_macrocycles=False,
                                flexible_amides=False).prepare(mol)
    if len(setup) != 1:
        raise ValueError("Rutin Meeko setup count differs from one")
    text, ok, error = PDBQTWriterLegacy.write_string(setup[0], add_index_map=True)
    if not ok:
        raise ValueError(f"Rutin Meeko conversion failed: {error}")
    common.write_once(pdbqt, text)
    return sdf, pdbqt


def make_start(pdbqt: Path, mol: Chem.Mol, mapping: dict[int, int], sid: str,
               seed: int, record: int, expected_smiles: str) -> tuple[Path, float]:
    text = pdbqt.read_text(encoding="utf-8")
    lines = text.splitlines()
    atom_lines = [x for x in lines if x.startswith(("ATOM  ", "HETATM"))]
    serials = [int(x[6:11]) for x in atom_lines]
    original = np.array([[float(x[i:i + 8]) for i in (30, 38, 46)] for x in atom_lines])
    point_by_serial = dict(zip(serials, original))
    nheavy = sum(a.GetAtomicNum() > 1 for a in mol.GetAtoms())
    heavy = np.array([point_by_serial[mapping[i]] for i in range(1, nheavy + 1)])
    centroid = heavy.mean(axis=0)
    rng = np.random.default_rng(np.random.SeedSequence([seed, int(sid[-2:]), record, 72]))
    for _ in range(10000):
        q = rng.normal(size=4)
        q /= np.linalg.norm(q)
        w, x, y, z = q
        rot = np.array([
            [1 - 2 * (y*y + z*z), 2 * (x*y - z*w), 2 * (x*z + y*w)],
            [2 * (x*y + z*w), 1 - 2 * (x*x + z*z), 2 * (y*z - x*w)],
            [2 * (x*z - y*w), 2 * (y*z + x*w), 1 - 2 * (x*x + y*y)],
        ])
        shift = rng.uniform(-3.0, 3.0, 3)
        if np.linalg.norm(shift) < 2.0:
            continue
        coords = (original - centroid) @ rot.T + CENTER + shift
        transformed = dict(zip(serials, coords))
        hxyz = np.array([transformed[mapping[i]] for i in range(1, nheavy + 1)])
        margin = float(np.min(SIZE / 2 - np.abs(hxyz - CENTER)))
        if margin < 0.5:
            continue
        out = []
        for line in lines:
            if line.startswith(("ATOM  ", "HETATM")):
                p = transformed[int(line[6:11])]
                line = line[:30] + "".join(f"{v:8.3f}" for v in p) + line[54:]
            out.append(line)
        path = STARTS / f"{sid}_neutral_record_{record:04d}_seed_{seed}_start.pdbqt"
        common.write_once(path, "\n".join(out) + "\n")
        back = PDBQTMolecule(path.read_text(encoding="utf-8"))
        rebuilt = RDKitMolCreate.from_pdbqt_mol(back)
        if len(rebuilt) != 1 or rebuilt[0] is None or common.iso(rebuilt[0]) != expected_smiles:
            raise ValueError("Frozen start failed graph/stereo readback")
        atom_xyz = back.atoms()["xyz"]
        mapped = np.array([atom_xyz[mapping[i] - 1] for i in range(1, nheavy + 1)])
        original_dist = np.linalg.norm(heavy[:, None, :] - heavy[None, :, :], axis=2)
        start_dist = np.linalg.norm(mapped[:, None, :] - mapped[None, :, :], axis=2)
        if np.max(np.abs(original_dist - start_dist)) > 0.002:
            raise ValueError("Frozen start changed intraligand heavy-atom geometry")
        actual = float(np.min(SIZE / 2 - np.abs(mapped - CENTER)))
        if actual < 0.499:
            raise ValueError("Frozen start crossed box face after rounding")
        return path, actual
    raise ValueError("No non-native start with 0.5 A clearance in frozen box")


def verify_controls() -> tuple[int, str]:
    for pair in [*common.RECEPTORS.values(), common.NATIVE]:
        common.check_hash(ROOT / pair[0], pair[1])
    rows = common.rows(rel(CONTROL))
    if len(rows) != 10 or {(r["receptor_state"], int(r["seed"])) for r in rows} != {
        (state, seed) for state in common.RECEPTORS for seed in SEEDS}:
        raise ValueError("Native MIG enlarged-box control row set incomplete")
    for r in rows:
        if r["status"] != "completed" or r["exit_code"] != "0" or float(r["rank1_RMSD_A"]) > 2:
            raise ValueError("Native MIG control status or RMSD failed")
        if r["box_size_xyz_A"] != "30;24;24":
            raise ValueError("Native MIG control box differs")
        for path, sha in (("receptor_path", "receptor_SHA256"), ("ligand_path", "ligand_SHA256"),
                          ("reference_path", "reference_SHA256"), ("output_path", "output_SHA256"),
                          ("log_path", "log_SHA256")):
            common.check_hash(ROOT / r[path], r[sha])
    return len(rows), digest(CONTROL)


def verify_prior_quercetin_outputs() -> dict[str, int]:
    counts = {}
    for label, filename in (
        ("original_box", "data/metadata/vina_3L4W_phytochemical_pilot_results.csv"),
        ("enlarged_box", "data/metadata/vina_3L4W_box_adequacy_results.csv"),
    ):
        rows = [r for r in common.rows(filename) if r["structure_id"] == "LIG06"]
        if len(rows) != 20 or len({(r["microstate"], r["receptor_state"], r["seed"]) for r in rows}) != 20:
            raise ValueError(f"Incomplete prior quercetin runs: {filename}")
        for r in rows:
            for path, sha in (("receptor_path", "receptor_sha256"),
                              ("ligand_path", "ligand_sha256"),
                              ("start_path", "start_sha256"),
                              ("output_path", "output_sha256"),
                              ("log_path", "log_sha256")):
                common.check_hash(ROOT / r[path], r[sha])
        counts[label] = len(rows)
    return counts


def main() -> None:
    control_count, control_hash = verify_controls()
    prior_counts = verify_prior_quercetin_outputs()
    library = {r["structure_id"]: r for r in common.rows("data/compounds/final_preliminary_ligand_library.csv")
               if r["structure_id"] in SPECS and r["structure_id"] != ""}
    prep = {r["structure_id"]: r for r in common.rows("data/metadata/ligand_preparation_results.csv")}
    vina = {r["structure_id"]: r for r in common.rows("data/metadata/vina_ligand_preparation_results.csv")}
    pubchem = official_pubchem()
    prior = [r for r in common.rows("data/metadata/vina_3L4W_pilot_input_manifest.csv")
             if r["structure_id"] == "LIG06"]
    result = []
    for sid, (plant, name, cid, formula, nrec, nheavy, nstereo, ntors) in SPECS.items():
        source = ROOT / prep[sid]["sdf_path"]
        original = ROOT / vina[sid]["pdbqt_path"]
        common.check_hash(source, prep[sid]["sdf_sha256"])
        common.check_hash(original, vina[sid]["pdbqt_sha256"])
        if str(cid) != library[sid]["PubChem_CID"] or str(cid) != prep[sid]["PubChem_CID"]:
            raise ValueError(f"{sid}: CID metadata conflict")
        ref = common.iso(Chem.MolFromSmiles(library[sid]["isomeric_SMILES"]))
        official = pubchem[cid]
        if official["MolecularFormula"] != formula or common.iso(Chem.MolFromSmiles(official["SMILES"])) != ref:
            raise ValueError(f"{sid}: live official PubChem identity mismatch")
        pool = list(Chem.SDMolSupplier(str(source), removeHs=False, sanitize=True))
        if len(pool) != nrec:
            raise ValueError(f"{sid}: source record count mismatch")
        for mol in pool:
            check_mol(mol, ref, formula, nheavy, nstereo)
        original_qc = common.check_pdbqt(original, pool[0])
        if original_qc["torsions"] != ntors or original_qc["heavy"] != nheavy:
            raise ValueError(f"{sid}: original PDBQT count mismatch")
        if sid == "LIG02":
            chosen = select_rutin(pool)
        elif sid == "LIG08":
            chosen = [0] * 5  # planar flavone; five independent rotations/translations
        else:
            chosen = [0]  # previous neutral/anion inputs and results reused, no new starts
        core = reference_core(pool[0], sid)
        energy_min = min(float(x.GetProp("mmff_energy_kcal_mol")) for x in pool)
        for slot, idx in enumerate(chosen):
            mol = pool[idx]
            record = idx + 1
            if sid == "LIG02":
                input_sdf, input_pdbqt = make_rutin_input(mol, record)
            else:
                input_sdf, input_pdbqt = source, original
            qc = common.check_pdbqt(input_pdbqt, mol if sid == "LIG02" else pool[0])
            if qc["heavy"] != nheavy or qc["torsions"] != ntors:
                raise ValueError(f"{sid} record {record}: prepared input QC failed")
            mapped_core = [qc["map"][i] for i in core]
            heavy_xyz = np.array([mol.GetConformer().GetAtomPosition(i) for i in range(nheavy)])
            radius = float(np.linalg.norm(heavy_xyz - heavy_xyz.mean(axis=0), axis=1).max())
            if radius + 0.5 > min(SIZE / 2):
                raise ValueError(f"{sid} record {record}: ligand cannot fit frozen box with clearance")
            if sid != "LIG06":
                start, margin = make_start(input_pdbqt, mol, qc["map"], sid, SEEDS[slot], record, ref)
            else:
                start, margin = None, None
            for receptor, (recpath, recsha) in common.RECEPTORS.items():
                row = dict.fromkeys(FIELDS, "")
                row.update({
                    "structure_id": sid, "compound_name": name, "plant_species": plant,
                    "PubChem_CID": cid, "association_scope": "prospective_new_association" if sid != "LIG02" else "existing_AM02",
                    "botanical_tier": "A_qualified_HPLC_DAD" if sid == "LIG06" else "B_organic_leaf",
                    "microstate": "neutral_reference", "formal_charge": 0, "decision": "PASS",
                    "decision_scope": "chemical_input_only_not_docking_validation",
                    "source_sdf": rel(source), "source_sdf_sha256": digest(source),
                    "source_pdbqt": rel(original), "source_pdbqt_sha256": digest(original),
                    "source_record_count": nrec, "selected_record_1based": record,
                    "selection_rule": "12_kcal_window_maxmin_fitted_heavy_RMSD_1p5A" if sid == "LIG02" else
                        "preexisting_frozen_record_1_reused" if sid == "LIG06" else "lowest_MMFF_record_1_five_independent_starts",
                    "MMFF94s_energy_kcal_mol": mol.GetProp("mmff_energy_kcal_mol"),
                    "energy_above_min_kcal_mol": f"{float(mol.GetProp('mmff_energy_kcal_mol'))-energy_min:.6f}",
                    "source_formula": formula, "isomeric_SMILES": ref, "heavy_atoms": nheavy,
                    "stereocentres": nstereo, "torsions": ntors, "qualified_sdf": rel(input_sdf),
                    "qualified_sdf_sha256": digest(input_sdf), "qualified_pdbqt": rel(input_pdbqt),
                    "qualified_pdbqt_sha256": digest(input_pdbqt), "pdbqt_charge_sum": f"{qc['charge']:.3f}",
                    "heavy_atom_map": ";".join(f"{i}:{qc['map'][i]}" for i in range(1, nheavy+1)),
                    "rigid_core_sdf_indices_1based": ";".join(map(str, core)),
                    "rigid_core_pdbqt_serials": ";".join(map(str, mapped_core)),
                    "core_rule": "three_flavonoid_rings_plus_carbonyl_O" if sid != "LIG06" else "prior_frozen_chromone_core",
                    "start_seed": SEEDS[slot] if sid != "LIG06" else "prior_frozen_seeds_1101_2202_3303_4404_5505",
                    "start_pdbqt": rel(start) if start else "prior_3L4W_pilot_inputs",
                    "start_sha256": digest(start) if start else "see_prior_pilot_manifest",
                    "start_box_margin_A": f"{margin:.3f}" if margin is not None else "see_prior_pilot_manifest",
                    "centered_conformer_radius_A": f"{radius:.3f}", "receptor_state": receptor,
                    "receptor_pdbqt": recpath, "receptor_sha256": recsha,
                    "MIG_reference_sha256": common.NATIVE[1],
                    "MIG_control_results_sha256": control_hash, "MIG_control_pass_count": control_count,
                    "box_center_A": ";".join(f"{v:.8f}" for v in CENTER),
                    "box_size_A": "30;24;24", "vina_version": "1.2.7", "scoring": "vina",
                    "exhaustiveness": 64, "num_modes": 20, "energy_range": 5.0,
                    "min_rmsd_A": 1.0, "cpu": 1,
                    "source_evidence_DOI": "10.3390/plants10122573" if sid == "LIG02" else
                        "10.1016/S1995-7645(14)60248-0" if sid == "LIG08" else "10.1016/j.fct.2009.05.033",
                    "limitations": "flexible_glycoside_sampling_and_source_sugar_stereo" if sid == "LIG02" else
                        "water_only_leaf_occurrence_unverified" if sid == "LIG08" else
                        "shared_LIG06_prior_poses_not_independent_Boerhavia_experiment",
                })
                result.append(row)
    # Reuse and recheck prior quercetin neutral/4'-O anion input hashes, without cloning files.
    if len(prior) != 4:
        raise ValueError("Prior frozen quercetin manifest is incomplete")
    for p in prior:
        common.check_hash(ROOT / p["prepared_variant_pdbqt_path"], p["prepared_variant_pdbqt_sha256"])
        if p["microstate"] == "neutral_keto_reference":
            continue
        common.check_hash(ROOT / p["prepared_variant_sdf_path"], p["prepared_variant_sdf_sha256"])
        mol = next(iter(Chem.SDMolSupplier(str(ROOT / p["prepared_variant_sdf_path"]), removeHs=False)))
        qc = common.check_pdbqt(ROOT / p["prepared_variant_pdbqt_path"], mol)
        if Chem.GetFormalCharge(mol) != -1 or qc["heavy"] != 22 or abs(qc["charge"] + 1) > 0.02:
            raise ValueError("Prior quercetin 4'-O phenolate failed readback")
        base = next(x for x in result if x["structure_id"] == "LIG06" and x["receptor_state"] == p["receptor_state"])
        row = base.copy()
        row.update({"microstate": "4prime_O_monoanion_hypothesis", "formal_charge": -1,
                    "source_pdbqt": p["source_pdbqt_path"], "source_pdbqt_sha256": p["source_pdbqt_sha256"],
                    "qualified_sdf": p["prepared_variant_sdf_path"],
                    "qualified_sdf_sha256": p["prepared_variant_sdf_sha256"],
                    "qualified_pdbqt": p["prepared_variant_pdbqt_path"],
                    "qualified_pdbqt_sha256": p["prepared_variant_pdbqt_sha256"],
                    "pdbqt_charge_sum": f"{qc['charge']:.3f}",
                    "heavy_atom_map": ";".join(f"{i}:{qc['map'][i]}" for i in range(1, 23)),
                    "rigid_core_pdbqt_serials": ";".join(str(qc["map"][int(i)]) for i in base["rigid_core_sdf_indices_1based"].split(";")),
                    "selection_rule": "reused_previously_frozen_4prime_O_phenolate_no_new_preparation"})
        result.append(row)
    # Explicitly retain unprepared phenolate hypotheses in the manifest. These
    # rows cannot be mistaken for qualified inputs because no derived file or
    # frozen starting pose is assigned.
    for sid, hypotheses in (
        ("LIG02", ("phenolate_site_unassigned",)),
        ("LIG08", ("7_O_monoanion_candidate", "4prime_O_monoanion_candidate")),
    ):
        for state in hypotheses:
            for receptor in common.RECEPTORS:
                base = next(x for x in result if x["structure_id"] == sid and x["receptor_state"] == receptor)
                hold = base.copy()
                hold.update({"microstate": state, "formal_charge": -1, "decision": "HOLD",
                             "decision_scope": "not_prepared_not_approved_for_search",
                             "qualified_sdf": "", "qualified_sdf_sha256": "",
                             "qualified_pdbqt": "", "qualified_pdbqt_sha256": "",
                             "pdbqt_charge_sum": "", "heavy_atom_map": "",
                             "rigid_core_pdbqt_serials": "", "start_seed": "",
                             "start_pdbqt": "", "start_sha256": "", "start_box_margin_A": "",
                             "limitations": "Phenolate site/population at pH 6.5 not sufficiently qualified; no derived input made"})
                result.append(hold)
    write_csv_once(OUT, result)
    print(json.dumps({"manifest": rel(OUT), "rows": len(result),
                      "rutin_records": sorted({int(r["selected_record_1based"]) for r in result if r["structure_id"] == "LIG02"}),
                      "new_start_files": len({r["start_pdbqt"] for r in result if r["structure_id"] != "LIG06" and r["decision"] == "PASS"}),
                      "MIG_controls": control_count, "prior_quercetin_runs_checked": prior_counts,
                      "status": "chemical_input_audit_pass"}, indent=2))


if __name__ == "__main__":
    main()

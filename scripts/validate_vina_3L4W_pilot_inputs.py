"""Audit and freeze 3L4W pilot inputs; generate only justified ionization variants.

Run from the repository root with .venv/Scripts/python.exe. This script never
invokes Vina and refuses to replace an existing file with different contents.
Atom indices in the manifest are one-based indices in the original SDF record.
"""

from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path

import meeko
import numpy as np
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/metadata/vina_3L4W_pilot_input_manifest.csv"
VARIANT_DIR = ROOT / "structures/ligands/vina_prepared/pilot_3L4W"
RECEPTORS = {
    "Asp542-": ("structures/proteins/vina_prepared/3L4W_noNAG_D443m_D542m_H600HID_pH6p5_dry_sensitivity_Meeko0.8.0.pdbqt", "7a76e0fa3b7bf2197714a57a88c60d2f4eefedccb18636951f3f643473afadf3"),
    "Asp542H(OD2)": ("structures/proteins/vina_prepared/3L4W_noNAG_D443m_D542H_OD2_H600HID_pH6p5_dry_sensitivity_Meeko0.8.0.pdbqt", "220e3585558e4bc261e6e76ddb036738dea827543e1b135399506e95ea9f93fd"),
}
RECEPTOR_PDB = {
    "Asp542-": ("structures/proteins/prepared/3L4W_noNAG_D443m_D542m_H600HID_pH6p5_dry_sensitivity.pdb", "3afebaeab24dfea3407446b653711cbc40ed1875eedc83e7e877d0512fe91f8e"),
    "Asp542H(OD2)": ("structures/proteins/prepared/3L4W_noNAG_D443m_D542H_OD2_H600HID_pH6p5_dry_sensitivity.pdb", "8e0872b236cb6975dc81ef131ffd3aacddbd8edd3bbc178e45eb1a18710192fd"),
}
NATIVE = ("docking/validation/vina/3L4W_MIG_native_pose.sdf", "ea673a59043bea05abbe7d66d8bd29c3786021f564c394027a8b0201385503f7")
NATIVE_PDBQT = ("docking/validation/vina/3L4W_MIG_crystal_pose_reference.pdbqt", "c7417f6bcb661e0f2217c6ba3619b6fd6c03ffb029fe5840d4b0d13c1bf16169")
CONTROL_SUMMARY = ("data/metadata/vina_target_validation_summary.csv", "39a7a5a36366e7ac2f05bc28c28ae2afe2d5539ba47748947310f6ceb69a867c")
VINA = (Path(r"E:\docking tools\vina_1.2.7_win.exe"), "e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5")
BOX_CENTER = np.array([45.38664286, 92.35157143, 34.82200000])
BOX_SIZE = np.array([24.0, 18.0, 18.0])
LIGANDS = {
    "LIG04": {"name": "Mangiferin", "cid": "5281647", "association": "SR01", "plant": "Salacia reticulata", "tier": "A", "deprot_o": 7, "deprot_bonded_carbon": 4, "deprot_label": "xanthone_6_OH", "expected_records": 60, "expected_torsions": 10},
    "LIG06": {"name": "Quercetin", "cid": "5280343", "association": "SA02", "plant": "Senna auriculata", "tier": "B", "deprot_o": 21, "deprot_bonded_carbon": 2, "deprot_label": "flavonol_4prime_OH", "expected_records": 24, "expected_torsions": 6},
}
FIELDS = [
    "record_type", "structure_id", "compound_name", "PubChem_CID", "plant_species", "association_id", "botanical_evidence_tier", "microstate", "formal_charge", "formula", "decision", "decision_reason",
    "source_sdf_path", "source_sdf_sha256", "source_record_count", "records_graph_stereo_3d_pass", "selected_record_1based", "selected_source_conformer_number", "selected_energy_kcal_mol", "conformer_energy_range_kcal_mol", "unique_3d_conformers_0p5A", "source_pdbqt_path", "source_pdbqt_sha256", "prepared_variant_sdf_path", "prepared_variant_sdf_sha256", "prepared_variant_pdbqt_path", "prepared_variant_pdbqt_sha256",
    "heavy_atoms", "stereocentres", "pdbqt_heavy_atoms", "pdbqt_charge_sum", "torsions", "pdbqt_graph_stereo_readback", "max_heavy_coordinate_shift_A", "heavy_atom_map_sdf_to_pdbqt", "rigid_core_sdf_indices_1based", "rigid_core_pdbqt_serials_1based", "ionization_site_sdf_index_1based", "ionization_evidence", "tautomer_status", "conformer_limitation",
    "receptor_state", "receptor_source_pdb_path", "receptor_source_pdb_sha256", "receptor_pdbqt_path", "receptor_pdbqt_sha256", "receptor_heavy_atoms", "receptor_state_heavy_coordinate_match", "native_mig_sdf_path", "native_mig_sdf_sha256", "native_mig_pdbqt_path", "native_mig_pdbqt_sha256", "native_control_summary_path", "native_control_summary_sha256", "native_control_status", "vina_executable_sha256", "vina_version", "box_center_A", "box_size_A", "exhaustiveness", "num_modes", "energy_range_kcal_mol", "min_rmsd_A", "cpu", "seeds", "software_versions", "limitations",
]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_hash(path: Path, expected: str) -> None:
    if not path.is_file() or sha(path) != expected:
        raise ValueError(f"Missing file or SHA256 mismatch: {path}")


def rows(path: str) -> list[dict[str, str]]:
    with (ROOT / path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def iso(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)


def check_3d(mol: Chem.Mol) -> None:
    if mol.GetNumConformers() != 1 or not mol.GetConformer().Is3D():
        raise ValueError("Expected exactly one 3D conformer per SDF record")
    xyz = np.asarray(mol.GetConformer().GetPositions())
    if not np.isfinite(xyz).all():
        raise ValueError("Nonfinite coordinates")
    if not any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()) or any(a.GetNumImplicitHs() for a in mol.GetAtoms()):
        raise ValueError("Missing explicit hydrogens")
    for bond in mol.GetBonds():
        d = np.linalg.norm(xyz[bond.GetBeginAtomIdx()] - xyz[bond.GetEndAtomIdx()])
        if not 0.65 < d < 2.2:
            raise ValueError(f"Invalid bond length {d:.3f} A")


def index_map(text: str) -> dict[int, int]:
    result: dict[int, int] = {}
    for line in text.splitlines():
        if line.startswith("REMARK INDEX MAP "):
            nums = list(map(int, line[len("REMARK INDEX MAP "):].split()))
            if len(nums) % 2:
                raise ValueError("Incomplete Meeko index map")
            result.update(zip(nums[::2], nums[1::2]))
    if not result or len(set(result.values())) != len(result):
        raise ValueError("Missing or duplicate Meeko atom map")
    return result


def check_pdbqt(path: Path, mol: Chem.Mol) -> dict[str, object]:
    text = path.read_text(encoding="utf-8")
    p = PDBQTMolecule(text)
    atoms = p.atoms()
    if not np.isfinite(atoms["xyz"]).all() or not np.isfinite(atoms["partial_charges"]).all():
        raise ValueError("Nonfinite PDBQT coordinates or charges")
    rebuilt = RDKitMolCreate.from_pdbqt_mol(p)
    if len(rebuilt) != 1 or rebuilt[0] is None or iso(rebuilt[0]) != iso(mol):
        raise ValueError(f"Meeko readback graph/stereochemistry mismatch: {path}")
    mapping = index_map(text)
    heavy = [a for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    if any(a.GetIdx() + 1 not in mapping for a in heavy):
        raise ValueError("Heavy-atom map incomplete")
    xyz = mol.GetConformer().GetPositions()
    shifts = [np.linalg.norm(xyz[a.GetIdx()] - atoms["xyz"][mapping[a.GetIdx() + 1] - 1]) for a in heavy]
    charge = float(atoms["partial_charges"].sum())
    formal = sum(a.GetFormalCharge() for a in mol.GetAtoms())
    if max(shifts) > 0.0011 or abs(charge - formal) > 0.02:
        raise ValueError("PDBQT coordinates or net charge disagree with SDF")
    torsion = [int(line.split()[1]) for line in text.splitlines() if line.startswith("TORSDOF ")]
    if len(torsion) != 1:
        raise ValueError("No unique torsion count")
    return {"map": mapping, "heavy": len(heavy), "charge": charge, "shift": max(shifts), "torsions": torsion[0]}


def core_indices(mol: Chem.Mol, sid: str) -> list[int]:
    # Xanthone's three fused aromatic rings (LIG04), or quercetin's two
    # fused chromone rings (LIG06); append their directly attached carbonyl O.
    rings = mol.GetRingInfo().AtomRings()
    chosen = rings[:3] if sid == "LIG04" else rings[1:]
    core = set().union(*map(set, chosen))
    carbonyl = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 8 and any(b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(a).GetIdx() in core for b in a.GetBonds())]
    core.update(carbonyl)
    if len(core) != (15 if sid == "LIG04" else 11):
        raise ValueError("Unexpected fused-ring rigid-core definition")
    return sorted(i + 1 for i in core)


def deprotonate(mol: Chem.Mol, oxygen_idx: int) -> Chem.Mol:
    o = mol.GetAtomWithIdx(oxygen_idx)
    hs = [a.GetIdx() for a in o.GetNeighbors() if a.GetAtomicNum() == 1]
    if o.GetSymbol() != "O" or o.GetFormalCharge() != 0 or len(hs) != 1:
        raise ValueError("Specified phenolic oxygen is not singly protonated")
    edit = Chem.RWMol(mol)
    edit.GetAtomWithIdx(oxygen_idx).SetFormalCharge(-1)
    edit.RemoveAtom(hs[0])
    out = edit.GetMol()
    Chem.SanitizeMol(out)
    Chem.AssignStereochemistryFrom3D(out)
    check_3d(out)
    old_xyz = mol.GetConformer().GetPositions()
    new_xyz = out.GetConformer().GetPositions()
    for i, atom in enumerate(mol.GetAtoms()):
        if atom.GetAtomicNum() > 1 and np.linalg.norm(old_xyz[i] - new_xyz[i]) > 1e-9:
            raise ValueError("Deprotonation changed a heavy-atom coordinate")
    return out


def write_once(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_text(encoding="utf-8") != content:
            raise FileExistsError(f"Refusing to overwrite changed output: {path}")
        return
    path.write_text(content, encoding="utf-8", newline="\n")


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def receptor_heavy_atoms(path: Path) -> list[tuple[str, np.ndarray]]:
    found = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atom_type = line.split()[-1]
            if atom_type not in ("H", "HD"):
                identity = line[12:26].strip()
                xyz = np.array([float(line[i:i+8]) for i in (30, 38, 46)])
                found.append((identity, xyz))
    return found


def main() -> None:
    for rel, expected in [*RECEPTORS.values(), *RECEPTOR_PDB.values(), NATIVE, NATIVE_PDBQT, CONTROL_SUMMARY]:
        check_hash(ROOT / rel, expected)
    check_hash(*VINA)
    mig = next(iter(Chem.SDMolSupplier(str(ROOT / NATIVE[0]), removeHs=False, sanitize=True)))
    if mig is None:
        raise ValueError("Deposited MIG SDF could not be parsed")
    mig_xyz = np.array([mig.GetConformer().GetAtomPosition(a.GetIdx()) for a in mig.GetAtoms() if a.GetAtomicNum() > 1], dtype=float)
    if len(mig_xyz) != 14 or np.max(np.abs(mig_xyz.mean(axis=0) - BOX_CENTER)) > 1e-7:
        raise ValueError("Frozen box center does not equal deposited MIG heavy-atom centroid")
    if np.any(np.abs(mig_xyz - BOX_CENTER) > BOX_SIZE / 2 - 4.0):
        raise ValueError("Frozen box does not retain a four-Angstrom MIG coordinate margin")
    controls = [r for r in rows(CONTROL_SUMMARY[0]) if r["PDB_ID"] == "3L4W"]
    if len(controls) != 5 or sum(int(r["completed_runs"]) for r in controls) != 25 or any(int(r["failed_runs"]) or int(r["rank1_passes_at_2A"]) != 5 for r in controls):
        raise ValueError("Archived 3L4W native-MIG positive-control summary is inconsistent")
    receptor_atoms = {state: receptor_heavy_atoms(ROOT / rel) for state, (rel, _) in RECEPTORS.items()}
    a, b = receptor_atoms.values()
    if len(a) != 6912 or len(b) != len(a) or any(x[0] != y[0] or np.linalg.norm(x[1] - y[1]) > 0.0011 for x, y in zip(a, b)):
        raise ValueError("3L4W receptor heavy atoms/coordinates differ between frozen microstates")
    prep = {r["structure_id"]: r for r in rows("data/metadata/ligand_preparation_results.csv")}
    vina_prep = {r["structure_id"]: r for r in rows("data/metadata/vina_ligand_preparation_results.csv")}
    library = rows("data/compounds/final_preliminary_ligand_library.csv")
    output: list[dict[str, object]] = []
    for sid, info in LIGANDS.items():
        old = prep[sid]
        old_vina = vina_prep[sid]
        matching = [r for r in library if r["structure_id"] == sid]
        if len(matching) != 1 or matching[0]["association_id"] != info["association"] or matching[0]["PubChem_CID"] != info["cid"]:
            raise ValueError(f"Association/CID mismatch for {sid}")
        sdf = ROOT / old["sdf_path"]
        pdbqt = ROOT / old_vina["pdbqt_path"]
        check_hash(sdf, old["sdf_sha256"])
        check_hash(pdbqt, old_vina["pdbqt_sha256"])
        ms = list(Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=True))
        if len(ms) != info["expected_records"] or len(ms) != int(old["conformers_exported"]):
            raise ValueError(f"Incorrect SDF record count for {sid}")
        canonical = Chem.MolToSmiles(Chem.MolFromSmiles(old["source_isomeric_SMILES"]), isomericSmiles=True)
        library_iso = Chem.MolToSmiles(Chem.MolFromSmiles(matching[0]["isomeric_SMILES"]), isomericSmiles=True)
        library_canonical = Chem.MolToSmiles(Chem.MolFromSmiles(matching[0]["canonical_SMILES"]), isomericSmiles=False)
        if library_iso != canonical or library_canonical != Chem.MolToSmiles(Chem.MolFromSmiles(canonical), isomericSmiles=False):
            raise ValueError(f"Library canonical/isomeric structure differs from prepared parent: {sid}")
        energies: list[float] = []
        for m in ms:
            if m is None or m.GetProp("structure_id") != sid or m.GetProp("PubChem_CID") != info["cid"] or iso(m) != canonical:
                raise ValueError(f"Identity/graph/stereo mismatch in {sid} SDF pool")
            if Chem.MolToSmiles(Chem.MolFromSmiles(m.GetProp("source_isomeric_SMILES")), isomericSmiles=True) != canonical:
                raise ValueError(f"Stored source isomeric SMILES mismatch in {sid} SDF pool")
            if rdMolDescriptors.CalcMolFormula(m) != old["neutral_parent_formula"] or sum(a.GetFormalCharge() for a in m.GetAtoms()) != 0:
                raise ValueError(f"Stored formula/charge mismatch in {sid} SDF pool")
            t = Chem.Mol(m)
            Chem.AssignStereochemistryFrom3D(t)
            if iso(t) != canonical:
                raise ValueError(f"3D stereo mismatch in {sid} SDF pool")
            check_3d(m)
            energies.append(float(m.GetProp("mmff_energy_kcal_mol")))
        if not all(math.isfinite(e) for e in energies):
            raise ValueError(f"Nonfinite conformer energy: {sid}")
        chosen_idx = min(range(len(ms)), key=lambda i: (energies[i], i))
        if chosen_idx + 1 != int(old_vina["selected_sdf_record_index_1based"]):
            raise ValueError(f"Existing PDBQT does not use lowest recorded energy: {sid}")
        selected = ms[chosen_idx]
        assert selected is not None
        ionization_oxygen = selected.GetAtomWithIdx(info["deprot_o"])
        if sorted(a.GetIdx() for a in ionization_oxygen.GetNeighbors() if a.GetAtomicNum() > 1) != [info["deprot_bonded_carbon"]]:
            raise ValueError(f"Ionization site does not match the frozen phenolic graph: {sid}")
        oldcheck = check_pdbqt(pdbqt, selected)
        if oldcheck["torsions"] != info["expected_torsions"]:
            raise ValueError(f"Unexpected torsion count: {sid}")
        # An approximate full-heavy conformer-diversity count in the source
        # coordinate frame; distinct means >=0.5 A after pairwise fitting.
        representatives: list[Chem.Mol] = []
        for m in ms:
            if not representatives or all(Chem.rdMolAlign.GetBestRMS(Chem.RemoveHs(m), Chem.RemoveHs(ref)) >= 0.5 for ref in representatives):
                representatives.append(m)
        core = core_indices(selected, sid)
        parent_chirality = Chem.FindMolChiralCenters(Chem.RemoveHs(selected), includeUnassigned=True)
        stereo = len(parent_chirality)
        states: list[tuple[str, Chem.Mol, Path, Path | None, str]] = [
            ("neutral_keto_reference", selected, pdbqt, None, "preexisting_reference")
        ]
        anion = deprotonate(selected, info["deprot_o"])
        if Chem.FindMolChiralCenters(Chem.RemoveHs(anion), includeUnassigned=True) != parent_chirality:
            raise ValueError(f"Deprotonation changed stereocentre identities or configurations: {sid}")
        variant_sdf = VARIANT_DIR / f"{sid}_{info['deprot_label']}_monoanion_pH6p5_hypothesis.sdf"
        variant_pdbqt = variant_sdf.with_suffix(".pdbqt")
        anion.SetProp("microstate", f"{info['deprot_label']}_monoanion_pH6p5_hypothesis")
        anion.SetProp("parent_source_sdf_sha256", sha(sdf))
        anion.SetProp("selected_parent_record_1based", str(chosen_idx + 1))
        anion.SetProp("ionization_atom_index_1based", str(info["deprot_o"] + 1))
        block = Chem.MolToMolBlock(anion, forceV3000=False)
        props = "".join(f">  <{k}>\n{anion.GetProp(k)}\n\n" for k in sorted(anion.GetPropNames()))
        write_once(variant_sdf, block + "\n" + props + "$$$$\n")
        reread = next(iter(Chem.SDMolSupplier(str(variant_sdf), removeHs=False, sanitize=True)))
        if reread is None or iso(reread) != iso(anion):
            raise ValueError(f"Monoanion SDF readback failed: {sid}")
        setup = MoleculePreparation(charge_model="gasteiger", merge_these_atom_types=("H",), add_index_map=True, hydrate=False, rigid_macrocycles=False, flexible_amides=False).prepare(anion)
        if len(setup) != 1:
            raise ValueError(f"Meeko did not produce one setup: {sid}")
        qtext, ok, err = PDBQTWriterLegacy.write_string(setup[0], add_index_map=True)
        if not ok:
            raise ValueError(f"Meeko anion conversion failed: {sid}: {err}")
        write_once(variant_pdbqt, qtext)
        states.append((f"{info['deprot_label']}_monoanion_hypothesis", anion, variant_pdbqt, variant_sdf, "new_pH6p5_sensitivity"))
        for state_name, mol, state_pdbqt, state_sdf, provenance in states:
            qc = check_pdbqt(state_pdbqt, mol)
            expected_state_torsions = info["expected_torsions"] - (1 if state_sdf else 0)
            if qc["torsions"] != expected_state_torsions:
                raise ValueError(f"Unexpected state torsions: {sid}/{state_name}")
            # Only heavy-atom indices are stable across removal of a phenolic
            # hydrogen. Hydrogen indices after that point are not original-SDF
            # indices and must never be presented as a frozen reference map.
            mapstr = ";".join(f"{a.GetIdx() + 1}:{qc['map'][a.GetIdx() + 1]}" for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
            core_serials = ";".join(str(qc["map"][i]) for i in core)
            for receptor_name, (receptor_rel, receptor_hash) in RECEPTORS.items():
                record: dict[str, object] = {f: "" for f in FIELDS}
                record.update({
                    "record_type": "pilot_input_pair", "structure_id": sid, "compound_name": info["name"], "PubChem_CID": info["cid"], "plant_species": info["plant"], "association_id": info["association"], "botanical_evidence_tier": info["tier"],
                    "microstate": state_name, "formal_charge": sum(a.GetFormalCharge() for a in mol.GetAtoms()), "formula": rdMolDescriptors.CalcMolFormula(mol), "decision": "PASS", "decision_reason": "Exact graph/stereo/coordinates and Meeko readback pass; explicitly labelled pH6.5 microstate hypothesis",
                    "source_sdf_path": old["sdf_path"], "source_sdf_sha256": sha(sdf), "source_record_count": len(ms), "records_graph_stereo_3d_pass": len(ms), "selected_record_1based": chosen_idx + 1, "selected_source_conformer_number": selected.GetProp("conformer_number"), "selected_energy_kcal_mol": f"{energies[chosen_idx]:.8f}", "conformer_energy_range_kcal_mol": f"{min(energies):.8f}..{max(energies):.8f}", "unique_3d_conformers_0p5A": len(representatives),
                    "source_pdbqt_path": old_vina["pdbqt_path"], "source_pdbqt_sha256": sha(pdbqt), "prepared_variant_sdf_path": relative(state_sdf) if state_sdf else "", "prepared_variant_sdf_sha256": sha(state_sdf) if state_sdf else "", "prepared_variant_pdbqt_path": relative(state_pdbqt), "prepared_variant_pdbqt_sha256": sha(state_pdbqt),
                    "heavy_atoms": sum(a.GetAtomicNum() > 1 for a in mol.GetAtoms()), "stereocentres": stereo, "pdbqt_heavy_atoms": qc["heavy"], "pdbqt_charge_sum": f"{qc['charge']:.3f}", "torsions": qc["torsions"], "pdbqt_graph_stereo_readback": "pass", "max_heavy_coordinate_shift_A": f"{qc['shift']:.4f}", "heavy_atom_map_sdf_to_pdbqt": mapstr, "rigid_core_sdf_indices_1based": ";".join(map(str, core)), "rigid_core_pdbqt_serials_1based": core_serials, "ionization_site_sdf_index_1based": info["deprot_o"] + 1 if state_sdf else "", "ionization_evidence": "aqueous_mangiferin_pKa1_6p52_at_6OH_DOI_10.1016/j.saa.2005.09.009" if sid == "LIG04" else "aqueous_quercetin_pKa_near_7p3_DOI_10.1002/chem.202001956_4prime_site_hypothesis", "tautomer_status": "keto_parent_and_phenolate_hypothesis_only; other tautomers not enumerated", "conformer_limitation": "MMFF94s isolated-molecule pool; selected minimum is not asserted bioactive; receptor-bound coverage unproven",
                    "receptor_state": receptor_name, "receptor_source_pdb_path": RECEPTOR_PDB[receptor_name][0], "receptor_source_pdb_sha256": RECEPTOR_PDB[receptor_name][1], "receptor_pdbqt_path": receptor_rel, "receptor_pdbqt_sha256": receptor_hash, "receptor_heavy_atoms": len(receptor_atoms[receptor_name]), "receptor_state_heavy_coordinate_match": "pass_0p0011A", "native_mig_sdf_path": NATIVE[0], "native_mig_sdf_sha256": NATIVE[1], "native_mig_pdbqt_path": NATIVE_PDBQT[0], "native_mig_pdbqt_sha256": NATIVE_PDBQT[1], "native_control_summary_path": CONTROL_SUMMARY[0], "native_control_summary_sha256": CONTROL_SUMMARY[1], "native_control_status": "25_of_25_completed;25_of_25_rank1_within_2A;no_rerun", "vina_executable_sha256": VINA[1], "vina_version": "1.2.7", "box_center_A": "45.38664286;92.35157143;34.82200000", "box_size_A": "24;18;18", "exhaustiveness": 64, "num_modes": 20, "energy_range_kcal_mol": 5.0, "min_rmsd_A": 1.0, "cpu": 1, "seeds": "1101;2202;3303;4404;5505", "software_versions": f"RDKit {rdBase.rdkitVersion}; Meeko {meeko.__version__}; Vina 1.2.7", "limitations": f"{provenance}; dry/no-NAG model; missing active-site waters; His600 HID only; literature occurrence does not establish detection in original aqueous extract; phenolic ionization/tautomer population unresolved",
                })
                output.append(record)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"Validated {len(output)} ligand-state/receptor-state pairs; wrote {relative(OUT)}")


if __name__ == "__main__":
    main()

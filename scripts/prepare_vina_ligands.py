"""Prepare and audit one neutral-reference PDBQT per LIG01-LIG10 with Meeko 0.8.0.

Run with .venv/Scripts/python.exe from the repository root. Every SDF record
is checked, while only the lowest stored MMFF94s-energy record is converted.
This script does not determine receptor-specific protonation or run Vina.
"""

from __future__ import annotations

import csv
import hashlib
import math
import os
import sys
from collections import defaultdict
from pathlib import Path

import meeko
import numpy as np
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors


ROOT = Path(__file__).resolve().parents[1]
SOURCE_RESULTS = ROOT / "data/metadata/ligand_preparation_results.csv"
SOURCE_VERIFY = ROOT / "data/metadata/ligand_structure_verification.csv"
SOURCE_LIBRARY = ROOT / "data/compounds/final_preliminary_ligand_library.csv"
OUTPUT_DIR = ROOT / "structures/ligands/vina_prepared"
OUTPUT_CSV = ROOT / "data/metadata/vina_ligand_preparation_results.csv"
FIELDS = [
    "structure_id", "compound_name", "PubChem_CID", "association_ids",
    "source_sdf_path", "source_sdf_sha256", "sdf_record_count", "valid_record_count",
    "selected_sdf_record_index_1based", "source_conformer_number", "selection_method",
    "selected_mmff94s_energy_kcal_mol", "source_embedding_seed", "source_force_field",
    "reference_formula", "reference_formal_charge", "selected_formula", "selected_formal_charge",
    "source_atom_count", "source_heavy_atom_count", "source_explicit_H_count",
    "source_3D_check", "source_connectivity_check", "source_stereochemistry_check",
    "PubChem_reference_match", "meeko_version", "rdkit_version", "meeko_settings",
    "pdbqt_path", "pdbqt_sha256", "pdbqt_atom_count", "pdbqt_heavy_atom_count",
    "pdbqt_polar_H_count", "nonpolar_H_merged_count", "pdbqt_charge_sum",
    "pdbqt_atom_types", "torsdof", "pdbqt_parser_check", "pdbqt_graph_stereo_roundtrip",
    "pdbqt_heavy_coordinate_check", "status", "protonation_tautomer_concern",
    "source_identity_caveat", "error_message",
]
SETTINGS = (
    "MoleculePreparation(charge_model='gasteiger', "
    "merge_these_atom_types=('H',), add_index_map=True, hydrate=False, "
    "rigid_macrocycles=False, flexible_amides=False); "
    "PDBQTWriterLegacy.write_string(add_index_map=True)"
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_iso_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)


def check_coordinates(mol: Chem.Mol) -> None:
    if mol.GetNumConformers() != 1 or not mol.GetConformer().Is3D():
        raise ValueError("SDF record does not contain exactly one labelled 3D conformer")
    xyz = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    if xyz.shape != (mol.GetNumAtoms(), 3) or not np.isfinite(xyz).all():
        raise ValueError("missing or nonfinite 3D atom coordinates")
    heavy = np.asarray([xyz[a.GetIdx()] for a in mol.GetAtoms() if a.GetAtomicNum() > 1])
    if heavy.shape[0] < 2 or np.linalg.norm(heavy.max(axis=0) - heavy.min(axis=0)) < 1.0:
        raise ValueError("degenerate 3D geometry")
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        distance = float(np.linalg.norm(xyz[a] - xyz[b]))
        if distance < 0.65 or distance > 2.2:
            raise ValueError(f"implausible covalent distance {distance:.3f} A for atoms {a + 1}-{b + 1}")


def check_record(
    mol: Chem.Mol | None,
    record_index: int,
    identity: dict[str, str],
    reference_smiles: str,
) -> float:
    if mol is None:
        raise ValueError(f"SDF record {record_index} did not parse or sanitize")
    sid = identity["structure_id"]
    for prop in ("structure_id", "PubChem_CID", "source_isomeric_SMILES", "conformer_number", "mmff_energy_kcal_mol"):
        if not mol.HasProp(prop):
            raise ValueError(f"SDF record {record_index} lacks property {prop}")
    if mol.GetProp("structure_id") != sid or mol.GetProp("PubChem_CID") != identity["PubChem_CID"]:
        raise ValueError(f"SDF record {record_index} structure ID or PubChem CID mismatch")
    if Chem.MolToSmiles(Chem.MolFromSmiles(mol.GetProp("source_isomeric_SMILES")), isomericSmiles=True) != reference_smiles:
        raise ValueError(f"SDF record {record_index} stored source SMILES differs from verified reference")
    if canonical_iso_smiles(mol) != reference_smiles:
        raise ValueError(f"SDF record {record_index} molecular connectivity or stereo differs from verified reference")
    from_3d = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(from_3d)
    if canonical_iso_smiles(from_3d) != reference_smiles:
        raise ValueError(f"SDF record {record_index} 3D-implied stereochemistry differs from reference")
    formula = rdMolDescriptors.CalcMolFormula(mol)
    if formula != identity["formula"]:
        raise ValueError(f"SDF record {record_index} formula {formula} differs from {identity['formula']}")
    formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if formal_charge != int(identity["formal_charge"]):
        raise ValueError(f"SDF record {record_index} formal charge mismatch")
    explicit_h = sum(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms())
    if explicit_h == 0 or sum(atom.GetNumImplicitHs() for atom in mol.GetAtoms()):
        raise ValueError(f"SDF record {record_index} lacks complete explicit hydrogens")
    check_coordinates(mol)
    energy = float(mol.GetProp("mmff_energy_kcal_mol"))
    if not math.isfinite(energy):
        raise ValueError(f"SDF record {record_index} has nonfinite MMFF energy")
    return energy


def parse_index_map(pdbqt: str) -> dict[int, int]:
    mapping = {}
    for line in pdbqt.splitlines():
        if line.startswith("REMARK INDEX MAP "):
            values = [int(x) for x in line[len("REMARK INDEX MAP "):].split()]
            if len(values) % 2:
                raise ValueError("Meeko index map has an odd number of entries")
            for source, output in zip(values[::2], values[1::2]):
                if source in mapping:
                    raise ValueError("Meeko index map has duplicate source atom")
                mapping[source] = output
    if not mapping or len(set(mapping.values())) != len(mapping):
        raise ValueError("Meeko index map is absent or duplicates output atoms")
    return mapping


def validate_pdbqt(pdbqt: str, mol: Chem.Mol) -> dict[str, object]:
    if not pdbqt.strip() or "ROOT" not in pdbqt or "TORSDOF " not in pdbqt:
        raise ValueError("PDBQT is empty or lacks ligand torsion-tree records")
    parsed = PDBQTMolecule(pdbqt)
    atoms = parsed.atoms()
    if len(atoms) == 0 or not np.isfinite(atoms["xyz"]).all() or not np.isfinite(atoms["partial_charges"]).all():
        raise ValueError("Meeko could not parse finite PDBQT atoms and charges")
    recovered = RDKitMolCreate.from_pdbqt_mol(parsed)
    if len(recovered) != 1 or recovered[0] is None:
        raise ValueError("Meeko PDBQT-to-RDKit reconstruction did not yield one molecule")
    if canonical_iso_smiles(recovered[0]) != canonical_iso_smiles(mol):
        raise ValueError("PDBQT read-back changes molecular graph or stereochemistry")
    mapping = parse_index_map(pdbqt)
    original = mol.GetConformer().GetPositions()
    max_heavy_shift = 0.0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        idx = atom.GetIdx() + 1
        if idx not in mapping:
            raise ValueError(f"PDBQT omits input heavy atom {idx}")
        target = atoms["xyz"][mapping[idx] - 1]
        max_heavy_shift = max(max_heavy_shift, float(np.linalg.norm(original[idx - 1] - target)))
    if max_heavy_shift > 0.0011:
        raise ValueError(f"PDBQT heavy coordinates differ from selected SDF by {max_heavy_shift:.4f} A")
    input_heavy = sum(a.GetAtomicNum() > 1 for a in mol.GetAtoms())
    polar_h = sum(x["atom_type"] == "HD" for x in atoms)
    if len(atoms) != input_heavy + polar_h:
        raise ValueError("PDBQT output heavy atom or polar H count is inconsistent")
    torsion_lines = [line for line in pdbqt.splitlines() if line.startswith("TORSDOF ")]
    if len(torsion_lines) != 1:
        raise ValueError("PDBQT must contain exactly one TORSDOF")
    torsdof = int(torsion_lines[0].split()[1])
    branches = sum(line.startswith("BRANCH ") for line in pdbqt.splitlines())
    ends = sum(line.startswith("ENDBRANCH ") for line in pdbqt.splitlines())
    if branches != ends or torsdof < 0:
        raise ValueError("PDBQT branch tree is unbalanced")
    charge = float(atoms["partial_charges"].sum())
    if abs(charge - sum(a.GetFormalCharge() for a in mol.GetAtoms())) > 0.02:
        raise ValueError(f"PDBQT net partial charge {charge:.3f} disagrees with formal charge")
    return {
        "pdbqt_atom_count": len(atoms),
        "pdbqt_heavy_atom_count": input_heavy,
        "pdbqt_polar_H_count": polar_h,
        "nonpolar_H_merged_count": sum(a.GetAtomicNum() == 1 for a in mol.GetAtoms()) - polar_h,
        "pdbqt_charge_sum": f"{charge:.3f}",
        "pdbqt_atom_types": ";".join(sorted(set(atoms["atom_type"]))),
        "torsdof": torsdof,
        "pdbqt_parser_check": "pass_Meeko_PDBQTMolecule",
        "pdbqt_graph_stereo_roundtrip": "pass_Meeko_RDKitMolCreate_exact_isomeric_graph",
        "pdbqt_heavy_coordinate_check": f"pass_max_shift_{max_heavy_shift:.4f}_A",
    }


def concern(structure_id: str) -> str:
    notes = {
        "LIG01": "Neutral gallic acid selected; carboxylate microstate at docking pH unresolved.",
        "LIG02": "Neutral rutin reference; phenolic microstates and flexible-glycoside sampling remain receptor-specific.",
        "LIG03": "Neutral E/trans caffeic acid reference; carboxylate microstate and botanical E/Z assignment unresolved.",
        "LIG04": "Neutral mangiferin reference; phenolic microstates and source-level sugar stereo need review.",
        "LIG05": "Neutral PubChem (2R,3R) (-)-epicatechin reference; published non-chiral HPLC did not establish optical isomer.",
        "LIG06": "Neutral quercetin keto/phenolic reference; alternate tautomer/protonation states need receptor-specific review.",
        "LIG07": "Neutral kaempferol keto/phenolic reference; alternate tautomer/protonation states need receptor-specific review.",
        "LIG08": "Neutral luteolin keto/phenolic reference; alternate tautomer/protonation states need receptor-specific review.",
        "LIG09": "Neutral agnuside PubChem reference; full water-only botanical evidence and source-level stereo remain conditional.",
        "LIG10": "Neutral nicotiflorin PubChem reference; source-level glycoside identity and flexible conformers need review.",
    }
    return notes[structure_id]


def process(row: dict[str, str], identity: dict[str, str], associations: list[dict[str, str]]) -> dict[str, object]:
    sid = row["structure_id"]
    result: dict[str, object] = {field: "" for field in FIELDS}
    result.update({
        "structure_id": sid,
        "compound_name": row["compound_name"],
        "PubChem_CID": row["PubChem_CID"],
        "association_ids": ";".join(a["association_id"] for a in associations),
        "source_sdf_path": row["sdf_path"],
        "source_sdf_sha256": row["sdf_sha256"],
        "reference_formula": identity["formula"],
        "reference_formal_charge": identity["formal_charge"],
        "meeko_version": meeko.__version__,
        "rdkit_version": rdBase.rdkitVersion,
        "meeko_settings": SETTINGS,
        "protonation_tautomer_concern": concern(sid),
        "source_identity_caveat": identity["notes"],
        "status": "failed",
    })
    try:
        if identity["PubChem_CID"] != row["PubChem_CID"]:
            raise ValueError("structure-verification and preparation CID mismatch")
        if not identity["source_URL"].rstrip("/").endswith("/" + row["PubChem_CID"]):
            raise ValueError("structure-verification PubChem URL/CID mismatch")
        if {a["PubChem_CID"] for a in associations} != {row["PubChem_CID"]}:
            raise ValueError("plant-association PubChem CID mismatch")
        reference = Chem.MolFromSmiles(row["source_isomeric_SMILES"])
        if reference is None:
            raise ValueError("verified source isomeric SMILES cannot be parsed")
        reference_smiles = Chem.MolToSmiles(reference, isomericSmiles=True)
        for association in associations:
            ref2 = Chem.MolFromSmiles(association["isomeric_SMILES"])
            if ref2 is None or Chem.MolToSmiles(ref2, isomericSmiles=True) != reference_smiles:
                raise ValueError("plant-association reference structure differs")
        sdf = ROOT / row["sdf_path"]
        if not sdf.is_file() or sha256(sdf) != row["sdf_sha256"]:
            raise ValueError("source SDF missing or SHA-256 differs from preparation record")
        records = list(Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=True))
        result["sdf_record_count"] = len(records)
        if len(records) != int(row["conformers_exported"]):
            raise ValueError("SDF record count differs from ligand-preparation metadata")
        reference_check = {**row, **identity}
        energies = [check_record(mol, i, reference_check, reference_smiles) for i, mol in enumerate(records, 1)]
        result.update({
            "valid_record_count": len(energies),
            "source_3D_check": "pass_all_records",
            "source_connectivity_check": "pass_all_records_exact_reference_graph",
            "source_stereochemistry_check": "pass_all_records_including_3D_implied_stereo",
            "PubChem_reference_match": "pass_CID_formula_charge_and_preverified_isomeric_SMILES",
        })
        chosen = min(range(len(records)), key=lambda i: (energies[i], i))
        mol = records[chosen]
        assert mol is not None
        result.update({
            "selected_sdf_record_index_1based": chosen + 1,
            "source_conformer_number": mol.GetProp("conformer_number"),
            "selection_method": "minimum_stored_MMFF94s_energy_within_same_neutral_reference; tie_lowest_SDF_record_index",
            "selected_mmff94s_energy_kcal_mol": f"{energies[chosen]:.8f}",
            "source_embedding_seed": mol.GetProp("random_seed"),
            "source_force_field": mol.GetProp("minimization_method"),
            "selected_formula": rdMolDescriptors.CalcMolFormula(mol),
            "selected_formal_charge": sum(a.GetFormalCharge() for a in mol.GetAtoms()),
            "source_atom_count": mol.GetNumAtoms(),
            "source_heavy_atom_count": sum(a.GetAtomicNum() > 1 for a in mol.GetAtoms()),
            "source_explicit_H_count": sum(a.GetAtomicNum() == 1 for a in mol.GetAtoms()),
        })
        setups = MoleculePreparation(
            charge_model="gasteiger", merge_these_atom_types=("H",),
            add_index_map=True, hydrate=False, rigid_macrocycles=False,
            flexible_amides=False,
        ).prepare(mol)
        if len(setups) != 1:
            raise ValueError(f"Meeko returned {len(setups)} molecule setups, expected one")
        pdbqt, ok, message = PDBQTWriterLegacy.write_string(setups[0], add_index_map=True)
        if not ok:
            raise ValueError(f"Meeko writer failed: {message}")
        result.update(validate_pdbqt(pdbqt, mol))
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        destination = OUTPUT_DIR / f"{sid}_CID{row['PubChem_CID']}_neutral_reference.pdbqt"
        temp = destination.with_suffix(".pdbqt.tmp")
        try:
            temp.write_text(pdbqt, encoding="ascii", newline="\n")
            reread = PDBQTMolecule.from_file(str(temp))
            if len(reread.atoms()) != int(result["pdbqt_atom_count"]):
                raise ValueError("PDBQT file read-back atom count mismatch")
            if not temp.stat().st_size:
                raise ValueError("PDBQT output file is empty")
            os.replace(temp, destination)
        finally:
            if temp.exists():
                temp.unlink()
        if sha256(sdf) != row["sdf_sha256"]:
            raise ValueError("source SDF changed during conversion")
        result.update({
            "pdbqt_path": destination.relative_to(ROOT).as_posix(),
            "pdbqt_sha256": sha256(destination),
            "status": "preliminary_PDBQT_QC_pass; docking_microstate_review_pending",
        })
    except Exception as exc:
        result["status"] = "failed"
        result["error_message"] = f"{type(exc).__name__}: {exc}"
    return result


def main() -> None:
    if Path(sys.prefix).resolve() != (ROOT / ".venv").resolve():
        raise SystemExit("Use the project's .venv Python")
    if meeko.__version__ != "0.8.0":
        raise SystemExit(f"Expected Meeko 0.8.0; found {meeko.__version__}")
    source = read_csv(SOURCE_RESULTS)
    verified = {r["structure_id"]: r for r in read_csv(SOURCE_VERIFY)}
    library = defaultdict(list)
    for row in read_csv(SOURCE_LIBRARY):
        if row["structure_id"]:
            library[row["structure_id"]].append(row)
    if len(source) != 10 or {r["structure_id"] for r in source} != set(verified) or set(verified) != set(library):
        raise SystemExit("Expected exactly ten matched assigned structures")
    results = []
    for row in source:
        result = process(row, verified[row["structure_id"]], library[row["structure_id"]])
        results.append(result)
        print(row["structure_id"], result["status"], result["error_message"])
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(results)


if __name__ == "__main__":
    main()

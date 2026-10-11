r"""Prepare neutral PubChem reference conformers; botanical occurrence is separate.

Run from the repository root with::

    .venv\Scripts\python.exe scripts\prepare_ligands.py

One multi-record SDF is written per unique structure_id. Each record is a distinct
minimized conformer of the same reference chemical, not a new plant association.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LIBRARY = ROOT / "data/compounds/final_preliminary_ligand_library.csv"
DEFAULT_VERIFICATION = ROOT / "data/metadata/ligand_structure_verification.csv"
DEFAULT_OUTPUT = ROOT / "structures/ligands/prepared"
DEFAULT_RESULTS = ROOT / "data/metadata/ligand_preparation_results.csv"
FLEXIBLE = {"LIG02", "LIG04", "LIG09", "LIG10"}
RESULT_FIELDS = [
    "structure_id", "compound_name", "PubChem_CID", "association_ids",
    "plant_species", "botanical_evidence_tiers", "botanical_evidence_status",
    "reference_stereochemistry_caveat", "source_isomeric_SMILES", "neutral_parent_formula",
    "formal_charge", "rdkit_version", "embedding_method", "random_seed",
    "conformers_requested", "conformers_embedded", "embedding_failures",
    "force_field", "mmff_parameters_supported", "conformers_minimized",
    "minimization_converged", "minimization_nonconverged", "minimization_failed",
    "conformers_exported", "roundtrip_qc_passed", "connectivity_check",
    "stereochemistry_check", "coordinate_check", "file_integrity_check",
    "preparation_status", "docking_readiness", "sdf_path", "sdf_sha256",
    "unresolved_issues",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def chemical_key(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True, isomericSmiles=True)


def valid_coordinates(mol: Chem.Mol) -> bool:
    conf = mol.GetConformer()
    if not conf.Is3D():
        return False
    coords = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    if not all(math.isfinite(value) for xyz in coords for value in (xyz.x, xyz.y, xyz.z)):
        return False
    heavy = [coords[a.GetIdx()] for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    return max((heavy[0] - xyz).Length() for xyz in heavy[1:]) > 0.5


def verify_3d_stereo(mol: Chem.Mol, reference_key: str) -> bool:
    check = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(check)
    return chemical_key(check) == reference_key


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def prepare_one(
    sid: str,
    associations: list[dict[str, str]],
    verification: dict[str, str],
    output_dir: Path,
    small_count: int,
    flexible_count: int,
    base_seed: int,
) -> dict[str, str | int]:
    source = associations[0]
    row: dict[str, str | int] = {field: "" for field in RESULT_FIELDS}
    association_ids = [r["association_id"] for r in associations]
    tiers = [f'{r["association_id"]}:{r["evidence_tier"]}' for r in associations]
    seed = base_seed + int(sid.removeprefix("LIG"))
    requested = flexible_count if sid in FLEXIBLE else small_count
    row.update({
        "structure_id": sid,
        "compound_name": source["compound_name"],
        "PubChem_CID": source["PubChem_CID"],
        "association_ids": ";".join(association_ids),
        "plant_species": ";".join(r["plant_species"] for r in associations),
        "botanical_evidence_tiers": ";".join(tiers),
        "botanical_evidence_status": ";".join(
            f'{r["association_id"]}:{r["source_evidence_classification"]}' for r in associations
        ),
        "reference_stereochemistry_caveat": source["stereochemical_identity"],
        "source_isomeric_SMILES": source["isomeric_SMILES"],
        "neutral_parent_formula": verification["formula"],
        "formal_charge": verification["formal_charge"],
        "rdkit_version": rdkit.__version__,
        "embedding_method": "ETKDGv3; useRandomCoords=False; numThreads=1; pruneRmsThresh=-1",
        "random_seed": seed,
        "conformers_requested": requested,
        "force_field": "MMFF94s",
        "preparation_status": "failed",
        "docking_readiness": "not_assessed",
        "unresolved_issues": source["limitations"],
    })
    try:
        if any(r["PubChem_CID"] != source["PubChem_CID"] or
               r["isomeric_SMILES"] != source["isomeric_SMILES"] for r in associations):
            raise ValueError("Shared structure_id has inconsistent CID or isomeric SMILES")
        if source["PubChem_CID"] != verification["PubChem_CID"]:
            raise ValueError("Library and verification CID differ")
        neutral = Chem.MolFromSmiles(source["isomeric_SMILES"])
        if neutral is None:
            raise ValueError("RDKit could not parse isomeric SMILES")
        reference_key = chemical_key(neutral)
        if rdMolDescriptors.CalcMolFormula(neutral) != verification["formula"]:
            raise ValueError("RDKit formula differs from PubChem verification")
        if Chem.GetFormalCharge(neutral) != int(verification["formal_charge"]):
            raise ValueError("RDKit charge differs from PubChem verification")
        row["connectivity_check"] = "input_formula_charge_and_graph_pass"

        molecule = Chem.AddHs(neutral)
        if not AllChem.MMFFHasAllMoleculeParams(molecule):
            row["mmff_parameters_supported"] = "no"
            raise ValueError("MMFF parameters unavailable")
        row["mmff_parameters_supported"] = "yes"
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        params.numThreads = 1
        params.pruneRmsThresh = -1.0
        params.useRandomCoords = False
        conf_ids = list(AllChem.EmbedMultipleConfs(molecule, numConfs=requested, params=params))
        row["conformers_embedded"] = len(conf_ids)
        row["embedding_failures"] = requested - len(conf_ids)
        if not conf_ids:
            raise ValueError("No conformer embedded")

        optimization = AllChem.MMFFOptimizeMoleculeConfs(
            molecule, numThreads=1, maxIters=1000, mmffVariant="MMFF94s"
        )
        if len(optimization) != len(conf_ids):
            raise ValueError("Optimization result count differs from embedded count")
        row["conformers_minimized"] = len(optimization)
        row["minimization_converged"] = sum(status == 0 for status, _ in optimization)
        row["minimization_nonconverged"] = sum(status == 1 for status, _ in optimization)
        row["minimization_failed"] = sum(status not in (0, 1) for status, _ in optimization)

        output_dir.mkdir(parents=True, exist_ok=True)
        final_path = output_dir / f"{sid}.sdf"
        staging_path = output_dir / f".{sid}.sdf.tmp"
        accepted: list[tuple[int, float, Chem.Mol]] = []
        rejected_stereo = rejected_geometry = 0
        for index, (conf_id, (status, energy)) in enumerate(zip(conf_ids, optimization), 1):
            if status != 0 or not math.isfinite(energy):
                continue
            conformer = Chem.Mol(molecule)
            conformer.RemoveAllConformers()
            conformer.AddConformer(molecule.GetConformer(conf_id), assignId=True)
            if not valid_coordinates(conformer):
                rejected_geometry += 1
                continue
            if not verify_3d_stereo(conformer, reference_key):
                rejected_stereo += 1
                continue
            conformer.SetProp("_Name", f"{sid}_neutral_reference_c{index:03d}")
            for key, value in {
                "structure_id": sid,
                "PubChem_CID": source["PubChem_CID"],
                "source_isomeric_SMILES": source["isomeric_SMILES"],
                "reference_type": "neutral_PubChem_reference_not_extract_confirmation",
                "association_ids": ";".join(association_ids),
                "botanical_evidence_tiers": ";".join(tiers),
                "conformer_number": str(index),
                "random_seed": str(seed),
                "embedding_method": "ETKDGv3",
                "minimization_method": "MMFF94s",
                "mmff_energy_kcal_mol": f"{energy:.8f}",
                "rdkit_version": rdkit.__version__,
                "formal_charge": verification["formal_charge"],
            }.items():
                conformer.SetProp(key, value)
            accepted.append((index, energy, conformer))
        row["stereochemistry_check"] = f"pre_export_pass={len(accepted)};rejected={rejected_stereo}"
        row["coordinate_check"] = f"pre_export_pass={len(accepted)};rejected={rejected_geometry}"
        if not accepted:
            raise ValueError("No converged conformer passed 3D and stereo QC")
        accepted.sort(key=lambda item: (item[1], item[0]))
        try:
            writer = Chem.SDWriter(str(staging_path))
            for _, _, conformer in accepted:
                writer.write(conformer)
            writer.close()
            loaded = [m for m in Chem.SDMolSupplier(str(staging_path), removeHs=False, sanitize=True)]
            if len(loaded) != len(accepted) or any(m is None for m in loaded):
                raise ValueError("SDF round-trip record count or sanitization failed")
            for written, loaded_mol in zip((item[2] for item in accepted), loaded):
                if chemical_key(loaded_mol) != reference_key:
                    raise ValueError("SDF round-trip connectivity or stereochemistry mismatch")
                if not verify_3d_stereo(loaded_mol, reference_key):
                    raise ValueError("SDF round-trip geometry-derived stereochemistry mismatch")
                if not valid_coordinates(loaded_mol):
                    raise ValueError("SDF round-trip coordinate check failed")
                if loaded_mol.GetNumAtoms() != written.GetNumAtoms():
                    raise ValueError("SDF round-trip atom count mismatch")
                if rdMolDescriptors.CalcMolFormula(loaded_mol) != verification["formula"]:
                    raise ValueError("SDF round-trip formula mismatch")
                if Chem.GetFormalCharge(loaded_mol) != int(verification["formal_charge"]):
                    raise ValueError("SDF round-trip formal charge mismatch")
                if loaded_mol.GetProp("PubChem_CID") != source["PubChem_CID"]:
                    raise ValueError("SDF round-trip CID property mismatch")
            staging_path.replace(final_path)
        finally:
            if staging_path.exists():
                staging_path.unlink()
        row["conformers_exported"] = len(accepted)
        row["roundtrip_qc_passed"] = len(loaded)
        row["connectivity_check"] = "pass_input_and_all_exported_records"
        row["stereochemistry_check"] = f"pass_all_exported_records;pre_export_rejected={rejected_stereo}"
        row["coordinate_check"] = f"pass_all_exported_records;pre_export_rejected={rejected_geometry}"
        row["file_integrity_check"] = "pass_SDF_roundtrip_and_SHA256"
        row["preparation_status"] = "prepared_neutral_reference_3D_QC_pass"
        row["sdf_path"] = final_path.relative_to(ROOT).as_posix()
        row["sdf_sha256"] = file_hash(final_path)
    except Exception as exc:  # keep a result row for every assigned structure
        row["preparation_status"] = "failed"
        row["file_integrity_check"] = f"failed:{type(exc).__name__}"
        row["unresolved_issues"] = f'{row["unresolved_issues"]} PREPARATION_FAILURE: {exc}'
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", type=Path, default=DEFAULT_LIBRARY)
    parser.add_argument("--verification", type=Path, default=DEFAULT_VERIFICATION)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--small-conformers", type=int, default=24)
    parser.add_argument("--flexible-conformers", type=int, default=60)
    parser.add_argument("--base-seed", type=int, default=20261011)
    args = parser.parse_args()
    if min(args.small_conformers, args.flexible_conformers) < 2:
        parser.error("conformer counts must be at least two")
    library = read_csv(args.library)
    checks = {r["structure_id"]: r for r in read_csv(args.verification)}
    groups: dict[str, list[dict[str, str]]] = {}
    for association in library:
        if association["structure_id"]:
            groups.setdefault(association["structure_id"], []).append(association)
    if len(library) != 12 or len(groups) != 10 or set(groups) != set(checks):
        parser.error("Expected 12 original associations and exactly 10 verified structures")
    results = [
        prepare_one(sid, groups[sid], checks[sid], args.output_dir,
                    args.small_conformers, args.flexible_conformers, args.base_seed)
        for sid in sorted(groups, key=lambda sid: int(sid.removeprefix("LIG")))
    ]
    args.results.parent.mkdir(parents=True, exist_ok=True)
    with args.results.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        writer.writerows(results)
    for row in results:
        print(f'{row["structure_id"]}: {row["preparation_status"]}; '
              f'{row["conformers_exported"]}/{row["conformers_embedded"]} exported/embedded')
    return 1 if any(row["preparation_status"] == "failed" for row in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())

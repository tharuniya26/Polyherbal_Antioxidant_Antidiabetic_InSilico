"""Build deposited native reference ligands from official RCSB definitions.

The script does not protonate or minimize receptors and never runs docking.
It fails closed: an SDF is not written if the reconstructed graph, formula,
stereochemistry, atom list, or crystallographic heavy-atom coordinates differ
from the RCSB reference definition and original 2QMJ/1XD0 coordinates.
"""

from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path

import gemmi
import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


ROOT = Path(__file__).resolve().parents[1]
ORIGINAL = ROOT / "structures/proteins/original"
DEFINITIONS = ROOT / "docking/validation/definitions"
VALIDATION = ROOT / "docking/validation"
PREPARED = ROOT / "structures/proteins/prepared"
RESULTS = ROOT / "data/metadata/receptor_preparation_results.csv"
AUDIT = ROOT / "data/metadata/receptor_structure_audit.csv"
CONFIG = {
    "2QMJ": {"reference_id": "PRD_900007", "reference_name": "alpha-acarbose", "chemistry": "PRDCC_900007", "components": ("GLC", "GLC", "AC1"), "author_chain": "B"},
    "1XD0": {"reference_id": "ARE", "reference_name": "acarbose-derived pentasaccharide", "chemistry": "ARE", "components": ("ARE",), "author_chain": "A"},
}
FIELDS = [
    "PDB_ID", "reference_ligand_id", "reference_ligand_name", "official_coordinate_URL",
    "official_definition_URL", "original_PDB_SHA256", "original_mmCIF_SHA256",
    "reference_definition_SHA256", "PDB_heavy_atoms", "expected_heavy_atoms",
    "full_connectivity_recovered", "formula", "formula_verified", "stereochemistry_verified",
    "atom_names_complete", "crystallographic_pose_preserved", "ligand_graph_InChIKey",
    "site_center_x", "site_center_y", "site_center_z", "working_receptor_path",
    "working_receptor_SHA256", "native_ligand_SDF", "native_ligand_SHA256",
    "atom_mapping_CSV", "atom_mapping_SHA256", "receptor_hydrogens_added",
    "ions_retained", "waters_retained", "other_hetero_retained", "GOLD_ligand_format_compatible",
    "GOLD_receptor_ready", "native_redocking_ready", "remaining_issues", "status",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def block_for(name: str) -> gemmi.cif.Block:
    return gemmi.cif.read_file(str(DEFINITIONS / f"{name}.cif")).sole_block()


def official_chemistry(name: str):
    block = block_for(name)
    atom = block.get_mmcif_category("_chem_comp_atom.")
    bond = block.get_mmcif_category("_chem_comp_bond.")
    atoms = {
        atom["atom_id"][i]: atom["type_symbol"][i]
        for i in range(len(atom["atom_id"])) if atom["type_symbol"][i] != "H"
    }
    bonds = [
        (bond["atom_id_1"][i], bond["atom_id_2"][i], bond["value_order"][i])
        for i in range(len(bond["atom_id_1"]))
    ]
    descriptor = block.get_mmcif_category("_pdbx_chem_comp_descriptor.")
    smiles = [
        "".join(descriptor["descriptor"][i].split())
        for i in range(len(descriptor["type"]))
        if descriptor["type"][i] == "SMILES_CANONICAL"
        and descriptor["program"][i] == "OpenEye OEToolkits"
    ]
    if len(smiles) != 1:
        raise ValueError(f"{name}: exactly one official stereodefined descriptor required")
    reference = Chem.MolFromSmiles(smiles[0])
    if reference is None:
        raise ValueError(f"{name}: official stereodefined descriptor did not parse")
    formula = "".join(block.find_value("_chem_comp.formula").strip('"').split())
    return atoms, bonds, reference, formula


def ligand_atom_rows(pdb_id: str):
    rows = []
    for line in (ORIGINAL / f"{pdb_id}.pdb").read_text().splitlines():
        if not line.startswith("HETATM"):
            continue
        comp, chain, seq = line[17:20].strip(), line[21].strip(), int(line[22:26])
        is_native = (pdb_id == "2QMJ" and chain == "B" and comp in ("GLC", "AC1")) or (
            pdb_id == "1XD0" and chain == "A" and comp == "ARE" and seq == 501
        )
        if is_native:
            rows.append({
                "serial": int(line[6:11]), "component": comp, "chain": chain,
                "residue": seq, "atom_name": line[12:16].strip(),
                "element": line[76:78].strip().title(),
                "coords": tuple(float(line[pos:pos + 8]) for pos in (30, 38, 46)),
                "occupancy": float(line[54:60]),
            })
    return rows


def verify_against_original_cif(pdb_id: str, pdb_atoms: list[dict]):
    block = gemmi.cif.read_file(str(ORIGINAL / f"{pdb_id}.cif")).sole_block()
    table = block.get_mmcif_category("_atom_site.")
    cif_atoms = {}
    for i in range(len(table["id"])):
        key = (table["auth_comp_id"][i], table["auth_asym_id"][i],
               table["auth_seq_id"][i], table["auth_atom_id"][i])
        cif_atoms[key] = tuple(float(table[field][i]) for field in ("Cartn_x", "Cartn_y", "Cartn_z"))
    for atom in pdb_atoms:
        key = (atom["component"], atom["chain"], str(atom["residue"]), atom["atom_name"])
        if key not in cif_atoms or any(abs(a - b) > 0.0011 for a, b in zip(atom["coords"], cif_atoms[key])):
            raise ValueError(f"{pdb_id}: PDB/mmCIF coordinate mismatch at {key}")


def build_reference(pdb_id: str, pdb_atoms: list[dict]):
    config = CONFIG[pdb_id]
    mutable = Chem.RWMol()
    positions = Chem.Conformer(len(pdb_atoms))
    index = {}
    for i, atom in enumerate(pdb_atoms):
        key = (atom["component"], atom["chain"], atom["residue"], atom["atom_name"])
        if key in index:
            raise ValueError(f"Duplicate deposited atom name: {key}")
        index[key] = i
        rd_atom = Chem.Atom(atom["element"])
        rd_atom.SetProp("PDB_atom_identifier", f'{atom["component"]}:{atom["chain"]}:{atom["residue"]}:{atom["atom_name"]}')
        mutable.AddAtom(rd_atom)
        positions.SetAtomPosition(i, atom["coords"])
    positions.Set3D(True)
    mutable.AddConformer(positions)
    residues = sorted({(atom["component"], atom["chain"], atom["residue"]) for atom in pdb_atoms}, key=lambda x: x[2])
    if tuple(res[0] for res in residues) != config["components"]:
        raise ValueError(f"{pdb_id}: unexpected component order or count")

    for component, chain, seq in residues:
        definition_atoms, definition_bonds, _, _ = official_chemistry(component)
        deposited = {atom["atom_name"] for atom in pdb_atoms if (atom["component"], atom["chain"], atom["residue"]) == (component, chain, seq)}
        if not deposited <= set(definition_atoms):
            raise ValueError(f"{pdb_id}: deposited atom outside CCD component {component} {seq}")
        expected_missing = set() if pdb_id == "1XD0" or seq == 1 else {"O1"}
        if set(definition_atoms) - deposited != expected_missing:
            raise ValueError(f"{pdb_id}: unexpected missing heavy atoms in {component} {seq}: {set(definition_atoms)-deposited}")
        for name in deposited:
            observed = pdb_atoms[index[(component, chain, seq, name)]]["element"].upper()
            if observed != definition_atoms[name].upper():
                raise ValueError(f"{pdb_id}: element mismatch at {component} {seq} {name}")
        for atom1, atom2, order in definition_bonds:
            key1 = (component, chain, seq, atom1)
            key2 = (component, chain, seq, atom2)
            if key1 in index and key2 in index:
                mutable.AddBond(index[key1], index[key2], {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE, "TRIP": Chem.BondType.TRIPLE}[order])

    if pdb_id == "2QMJ":
        bird = block_for("PRD_900007")
        entity = bird.get_mmcif_category("_pdbx_reference_molecule.")
        if entity["name"][0] != "alpha-acarbose":
            raise ValueError("Unexpected BIRD identity")
        link = bird.get_mmcif_category("_pdbx_reference_entity_poly_link.")
        pdb_lines = (ORIGINAL / "2QMJ.pdb").read_text().splitlines()
        pdb_links = [line for line in pdb_lines if line.startswith("LINK") and line[21] == "B" and line[51] == "B"]
        if len(pdb_links) != 2 or len(link["link_id"]) != 2:
            raise ValueError("2QMJ: expected two matching deposited and BIRD inter-component links")
        for i in range(len(link["link_id"])):
            if link["value_order"][i] != "SING":
                raise ValueError("2QMJ: unexpected BIRD glycosidic bond order")
            key1 = (link["comp_id_1"][i], "B", int(link["entity_seq_num_1"][i]), link["atom_id_1"][i])
            key2 = (link["comp_id_2"][i], "B", int(link["entity_seq_num_2"][i]), link["atom_id_2"][i])
            if key1 not in index or key2 not in index:
                raise ValueError(f"2QMJ: BIRD link atom missing {key1} {key2}")
            if not any(line[12:16].strip() == key1[3] and line[17:20].strip() == key1[0]
                       and int(line[22:26]) == key1[2] and line[42:46].strip() == key2[3]
                       and line[47:50].strip() == key2[0] and int(line[52:56]) == key2[2]
                       for line in pdb_links):
                raise ValueError(f"2QMJ: original PDB LINK not found for {key1} to {key2}")
            mutable.AddBond(index[key1], index[key2], Chem.BondType.SINGLE)

    molecule = mutable.GetMol()
    Chem.SanitizeMol(molecule)
    Chem.AssignStereochemistryFrom3D(molecule)
    _, _, official, formula = official_chemistry(config["chemistry"])
    if len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError(f"{pdb_id}: reconstructed ligand is disconnected")
    if molecule.GetNumHeavyAtoms() != official.GetNumHeavyAtoms():
        raise ValueError(f"{pdb_id}: heavy-atom count mismatch with official whole-ligand definition")
    if rdMolDescriptors.CalcMolFormula(molecule) != formula:
        raise ValueError(f"{pdb_id}: formula mismatch with official whole-ligand definition")
    if Chem.MolToSmiles(molecule, isomericSmiles=True) != Chem.MolToSmiles(official, isomericSmiles=True):
        raise ValueError(f"{pdb_id}: whole-ligand connectivity or crystallographic stereochemistry differs from official definition")
    return molecule, official, formula


def write_working_receptor(pdb_id: str, removed_serials: set[int]) -> Path:
    source = ORIGINAL / f"{pdb_id}.pdb"
    destination = PREPARED / f"{pdb_id}_working_no_native_ligand.pdb"
    PREPARED.mkdir(parents=True, exist_ok=True)
    output = [
        f"REMARK 900 DERIVED FROM RCSB {pdb_id}; ORIGINAL SHA256 {sha256(source)}",
        "REMARK 900 NATIVE INHIBITOR OMITTED; NO PROTONATION OR MINIMIZATION APPLIED",
    ]
    keep_prefixes = ("HEADER", "TITLE ", "COMPND", "SOURCE", "EXPDTA", "CRYST1", "DBREF ",
                     "SEQRES", "MODRES", "SSBOND", "LINK  ", "ATOM  ", "HETATM", "TER   ", "CONECT", "END   ")
    for line in source.read_text().splitlines():
        if not line.startswith(keep_prefixes):
            continue
        if line.startswith(("ATOM  ", "HETATM")) and int(line[6:11]) in removed_serials:
            continue
        if line.startswith("LINK  "):
            # The only removed LINK records are the two internal acarbose links.
            if pdb_id == "2QMJ" and line[21] == "B" and line[51] == "B":
                continue
        if line.startswith("CONECT"):
            numbers = [int(line[i:i + 5]) for i in range(6, len(line), 5) if line[i:i + 5].strip().isdigit()]
            if not numbers or numbers[0] in removed_serials:
                continue
            numbers = [n for n in numbers if n not in removed_serials]
            if len(numbers) < 2:
                continue
            line = "CONECT" + "".join(f"{n:5d}" for n in numbers)
        output.append(line)
    destination.write_text("\n".join(output) + "\n", encoding="ascii", newline="\n")
    remaining = destination.read_text().splitlines()
    if any(line.startswith(("ATOM  ", "HETATM")) and int(line[6:11]) in removed_serials for line in remaining):
        raise ValueError(f"{pdb_id}: native ligand atoms remain in receptor copy")
    if pdb_id == "1XD0" and not all(any(line.startswith("HETATM") and line[17:20].strip() == ion for line in remaining) for ion in ("CA", "CL", "PCA")):
        raise ValueError("1XD0: required CA/CL/PCA missing in receptor copy")
    return destination


def write_mapping(pdb_id: str, reference_id: str, atoms: list[dict]) -> Path:
    path = VALIDATION / f"{pdb_id}_{reference_id}_atom_map.csv"
    fields = ["sdf_atom_index_1based", "PDB_serial", "auth_chain", "auth_residue_name", "auth_residue_number",
              "auth_atom_name", "element", "deposited_x", "deposited_y", "deposited_z"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for i, atom in enumerate(atoms, 1):
            writer.writerow({
                "sdf_atom_index_1based": i,
                "PDB_serial": atom["serial"],
                "auth_chain": atom["chain"],
                "auth_residue_name": atom["component"],
                "auth_residue_number": atom["residue"],
                "auth_atom_name": atom["atom_name"],
                "element": atom["element"],
                "deposited_x": f'{atom["coords"][0]:.3f}',
                "deposited_y": f'{atom["coords"][1]:.3f}',
                "deposited_z": f'{atom["coords"][2]:.3f}',
            })
    return path


def write_native_sdf(pdb_id: str, reference_id: str, molecule: Chem.Mol, atoms: list[dict], mapping: Path) -> Path:
    ligand = Chem.AddHs(molecule, addCoords=True)
    if ligand.GetNumHeavyAtoms() != len(atoms):
        raise ValueError(f"{pdb_id}: hydrogen addition changed heavy-atom count")
    ligand.SetProp("_Name", f"{pdb_id}_{reference_id}_deposited_reference")
    ligand.SetProp("reference_id", reference_id)
    ligand.SetProp("source_PDB", pdb_id)
    ligand.SetProp("original_PDB_SHA256", sha256(ORIGINAL / f"{pdb_id}.pdb"))
    ligand.SetProp("atom_mapping_file", mapping.relative_to(ROOT).as_posix())
    ligand.SetProp("coordinate_note", "Deposited heavy-atom coordinates unchanged; hydrogen coordinates generated by RDKit")
    ligand.SetProp("graph_source", "RCSB CCD/BIRD definitions and deposited PDB LINK records")
    path = VALIDATION / f"{pdb_id}_{reference_id}_native_pose.sdf"
    writer = Chem.SDWriter(str(path))
    writer.write(ligand)
    writer.close()
    loaded = list(Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True))
    if len(loaded) != 1 or loaded[0] is None:
        raise ValueError(f"{pdb_id}: SDF read-back failed")
    rebuilt = loaded[0]
    if Chem.MolToSmiles(Chem.RemoveHs(rebuilt)) != Chem.MolToSmiles(molecule):
        raise ValueError(f"{pdb_id}: SDF graph/stereo mismatch")
    for i, atom in enumerate(atoms):
        pos = rebuilt.GetConformer().GetAtomPosition(i)
        if any(abs(v - target) > 0.0011 for v, target in zip((pos.x, pos.y, pos.z), atom["coords"])):
            raise ValueError(f"{pdb_id}: SDF heavy-atom reference pose changed at {atom['atom_name']}")
    if not rebuilt.GetConformer().Is3D() or not all(
        math.isfinite(v) for atom in rebuilt.GetAtoms() for v in
        (rebuilt.GetConformer().GetAtomPosition(atom.GetIdx()).x,
         rebuilt.GetConformer().GetAtomPosition(atom.GetIdx()).y,
         rebuilt.GetConformer().GetAtomPosition(atom.GetIdx()).z)
    ):
        raise ValueError(f"{pdb_id}: nonfinite or non-3D coordinates")
    return path


def run_one(pdb_id: str, audit: dict[str, str]) -> dict[str, str]:
    cfg = CONFIG[pdb_id]
    pdb_path = ORIGINAL / f"{pdb_id}.pdb"
    cif_path = ORIGINAL / f"{pdb_id}.cif"
    row = {field: "" for field in FIELDS}
    row.update({
        "PDB_ID": pdb_id, "reference_ligand_id": cfg["reference_id"],
        "reference_ligand_name": cfg["reference_name"],
        "official_coordinate_URL": f"https://files.rcsb.org/download/{pdb_id}.cif",
        "official_definition_URL": "https://files.rcsb.org/birds/download/PRD_900007.cif;https://files.rcsb.org/birds/download/PRDCC_900007.cif;https://files.rcsb.org/ligands/download/GLC.cif;https://files.rcsb.org/ligands/download/AC1.cif" if pdb_id == "2QMJ" else "https://files.rcsb.org/ligands/download/ARE.cif",
        "original_PDB_SHA256": sha256(pdb_path), "original_mmCIF_SHA256": sha256(cif_path),
        "reference_definition_SHA256": sha256(DEFINITIONS / f'{cfg["chemistry"]}.cif'),
        "receptor_hydrogens_added": "no", "GOLD_receptor_ready": "no",
        "native_redocking_ready": "no", "status": "failed",
    })
    try:
        if row["original_PDB_SHA256"] != audit["SHA256"]:
            raise ValueError("Original PDB hash differs from prior audit")
        atoms = ligand_atom_rows(pdb_id)
        if not atoms:
            raise ValueError("Native ligand atoms missing")
        verify_against_original_cif(pdb_id, atoms)
        molecule, official, formula = build_reference(pdb_id, atoms)
        row.update({
            "PDB_heavy_atoms": len(atoms), "expected_heavy_atoms": official.GetNumHeavyAtoms(),
            "full_connectivity_recovered": "yes", "formula": formula,
            "formula_verified": "yes", "stereochemistry_verified": "yes; deposited 3D and whole-ligand RCSB reference agree",
            "atom_names_complete": "yes; component/whole-reference heavy atom list matches",
            "ligand_graph_InChIKey": Chem.inchi.MolToInchiKey(molecule),
        })
        center = tuple(sum(atom["coords"][axis] for atom in atoms) / len(atoms) for axis in range(3))
        row["site_center_x"], row["site_center_y"], row["site_center_z"] = (f"{v:.3f}" for v in center)
        VALIDATION.mkdir(parents=True, exist_ok=True)
        map_path = write_mapping(pdb_id, cfg["reference_id"], atoms)
        sdf_path = write_native_sdf(pdb_id, cfg["reference_id"], molecule, atoms, map_path)
        receptor_path = write_working_receptor(pdb_id, {atom["serial"] for atom in atoms})
        row.update({
            "crystallographic_pose_preserved": "yes; every heavy atom within 0.0011 A of original PDB and mmCIF",
            "working_receptor_path": receptor_path.relative_to(ROOT).as_posix(),
            "working_receptor_SHA256": sha256(receptor_path),
            "native_ligand_SDF": sdf_path.relative_to(ROOT).as_posix(),
            "native_ligand_SHA256": sha256(sdf_path),
            "atom_mapping_CSV": map_path.relative_to(ROOT).as_posix(),
            "atom_mapping_SHA256": sha256(map_path),
            "ions_retained": "SO4 retained" if pdb_id == "2QMJ" else "CA A498 and CL A500 retained",
            "waters_retained": "all 626 deposited waters retained pending site review" if pdb_id == "2QMJ" else "all 142 deposited waters retained including calcium-coordinating HOH536/546",
            "other_hetero_retained": "NAG glycosylation and nine GOL retained" if pdb_id == "2QMJ" else "PCA A1 and NAG A499 retained",
            "GOLD_ligand_format_compatible": "yes; one connected explicit-H 3D SDF record; GOLD import not tested",
            "remaining_issues": "Receptor protonation/atom typing; site NAG A2005 and waters; alternate locations; inhibitor-density review; GOLD import and redocking untested" if pdb_id == "2QMJ" else "Receptor protonation/atom typing; CA/CL and coordinating-water treatment; site waters; ARE crystal-mate contact; GOLD import and redocking untested",
            "status": "reference_ligand_validated; working_receptor_coordinate_copy_created",
        })
    except Exception as exc:
        row["remaining_issues"] = f"PREPARATION_BLOCKER: {type(exc).__name__}: {exc}"
    return row


def main() -> int:
    with AUDIT.open(newline="", encoding="utf-8-sig") as handle:
        audits = {row["PDB_ID"]: row for row in csv.DictReader(handle)}
    rows = [run_one(pdb_id, audits[pdb_id]) for pdb_id in ("2QMJ", "1XD0")]
    RESULTS.parent.mkdir(parents=True, exist_ok=True)
    with RESULTS.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(row["PDB_ID"], row["status"], row["remaining_issues"])
    return 1 if any(row["status"] == "failed" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())

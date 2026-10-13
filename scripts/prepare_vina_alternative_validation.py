"""Prepare audited 4GQR/MYC and 3L4W/MIG crystal references and local receptor probes.

Run: .venv/Scripts/python.exe scripts/prepare_vina_alternative_validation.py
Only downloads official RCSB originals if absent. Refuses to replace changed outputs.
No docking, scoring, minimization, or pKa prediction is performed.
"""
from __future__ import annotations

import csv
import hashlib
import io
import math
import subprocess
import tempfile
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import gemmi
import meeko
import numpy as np
import openmm
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTReceptor, PDBQTWriterLegacy, RDKitMolCreate
from openmm.app import Modeller, PDBFile
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from resolve_vina_receptor_templates import check_deposited_pca_stereo

ROOT = Path(__file__).resolve().parents[1]
ORIGINAL = ROOT / "structures/proteins/original"
DEFINITIONS = ROOT / "docking/validation/definitions"
VALIDATION = ROOT / "docking/validation/vina"
PREPARED = ROOT / "structures/proteins/prepared"
VINA_PREPARED = ROOT / "structures/proteins/vina_prepared"
RESULTS = ROOT / "data/metadata/vina_alternative_receptor_preparation_results.csv"
CLI = ROOT / ".venv/Scripts/mk_prepare_receptor.exe"
PCA_TEMPLATE = ROOT / "docking/configs/vina_templates/PCA_C_terminal_Meeko_0.8.0.json"

CONFIG = {
    "4GQR": {"ligand": "MYC", "seq": 504, "pH": 7.5, "expected_heavy": 23,
             "pdb_sha": "e3d36d41f8d59498cb24073997d5ddacf407269b00eef9b6b50c6530ca36b584",
             "cif_sha": "bae3af37a5f6e69298fcf941c2b737d08d30d6a14fe6f370e4ff06b914aaf333",
             "ccd_sha": "fe429a0e62755f5b3296fb359e31df4dad3397d5909135b851074341c7a97324"},
    "3L4W": {"ligand": "MIG", "seq": 1001, "pH": 6.5, "expected_heavy": 14,
             "pdb_sha": "57638e45666e2f848e37691a914a978d85c377c13eaed9e5f94af24d2c52c95e",
             "cif_sha": "81f5582ab3356b1e6318c3474367268712c5b8c2886193598564fe40e0218c46",
             "ccd_sha": "4420e2e083b312422cc2e08833a297c68ae25fdebcde01d180daeef56df28b28"},
}
FIELDS = ["PDB_ID", "protein", "native_ligand", "original_PDB_path", "original_PDB_SHA256",
          "original_mmCIF_path", "original_mmCIF_SHA256", "CCD_path", "CCD_SHA256",
          "native_SDF_path", "native_SDF_SHA256", "native_PDBQT_path", "native_PDBQT_SHA256",
          "native_atom_map_path", "native_atom_map_SHA256", "native_bridge_map_path", "native_bridge_map_SHA256",
          "native_heavy_atoms", "native_torsions",
          "native_graph_stereo_pose_check", "native_microstate_status", "site_center_xyz_A",
          "working_receptor_path", "working_receptor_SHA256", "local_receptor_PDB_path",
          "local_receptor_PDB_SHA256", "local_receptor_PDBQT_path", "local_receptor_PDBQT_SHA256",
          "receptor_conversion_status", "receptor_heavy_coordinate_check", "receptor_heavy_atoms",
          "receptor_PDBQT_atoms", "ion_PDBQT_fields", "water_OH_range_A", "catalytic_microstates", "prep_pH",
          "glycan_decision", "glycan_closest_A", "glycerol_closest_A", "ions_retained", "waters_retained",
          "PCA_link_stereo_check", "missing_residue_status",
          "altloc_rule", "template_path", "template_SHA256", "missing_atom_QC",
          "software_versions", "redocking_ready", "blockers"]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def safe_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != data:
            raise RuntimeError(f"refusing to overwrite differing output: {path}")
    else:
        path.write_bytes(data)


def get_official(url: str, path: Path, expected_sha: str) -> None:
    if not path.exists():
        with urllib.request.urlopen(url, timeout=60) as response:
            safe_write(path, response.read())
    if sha(path) != expected_sha:
        raise RuntimeError(f"official source checksum mismatch: {path}")


def source_line_atoms(lines: list[str], resname: str, resnum: int) -> list[dict]:
    out = []
    for line in lines:
        if line.startswith("HETATM") and line[17:20].strip() == resname and line[21] == "A" and int(line[22:26]) == resnum:
            out.append({"serial": int(line[6:11]), "name": line[12:16].strip(),
                        "element": line[76:78].strip().upper(), "occ": float(line[54:60]),
                        "xyz": tuple(float(line[i:i+8]) for i in (30, 38, 46))})
    return out


def ccd_reference(path: Path):
    block = gemmi.cif.read_file(str(path)).sole_block()
    atom = block.get_mmcif_category("_chem_comp_atom.")
    bond = block.get_mmcif_category("_chem_comp_bond.")
    desc = block.get_mmcif_category("_pdbx_chem_comp_descriptor.")
    heavy = {atom["atom_id"][i]: atom["type_symbol"][i].upper()
             for i in range(len(atom["atom_id"])) if atom["type_symbol"][i].upper() != "H"}
    bonds = [(bond["atom_id_1"][i], bond["atom_id_2"][i], bond["value_order"][i])
             for i in range(len(bond["atom_id_1"]))]
    smiles = [desc["descriptor"][i] for i in range(len(desc["type"]))
              if desc["type"][i] == "SMILES_CANONICAL" and desc["program"][i] == "OpenEye OEToolkits"]
    if len(smiles) != 1:
        raise RuntimeError(f"CCD stereochemical descriptor count {len(smiles)} for {path.name}")
    ref = Chem.MolFromSmiles(smiles[0])
    if ref is None:
        raise RuntimeError("CCD stereochemical SMILES failed RDKit parsing")
    formula = "".join(block.find_value("_chem_comp.formula").strip('"').split())
    return heavy, bonds, ref, formula


def verify_cif_coords(path: Path, atoms: list[dict], ligand: str, seq: int) -> None:
    table = gemmi.cif.read_file(str(path)).sole_block().get_mmcif_category("_atom_site.")
    found = {}
    for i in range(len(table["id"])):
        if (table["auth_comp_id"][i], table["auth_asym_id"][i], table["auth_seq_id"][i]) == (ligand, "A", str(seq)):
            found[table["auth_atom_id"][i]] = tuple(float(table[c][i]) for c in ("Cartn_x", "Cartn_y", "Cartn_z"))
    if len(found) != len(atoms):
        raise RuntimeError("PDB/mmCIF native atom count disagrees")
    for a in atoms:
        if a["name"] not in found or math.dist(a["xyz"], found[a["name"]]) > 0.0011:
            raise RuntimeError(f"PDB/mmCIF native coordinate mismatch: {a['name']}")


def reference_molecule(atoms: list[dict], ccd_path: Path):
    expected, bonds, ref, formula = ccd_reference(ccd_path)
    if len(atoms) != len(expected) or {a["name"]: a["element"] for a in atoms} != expected:
        raise RuntimeError("deposited native atoms are incomplete or disagree with CCD")
    mol = Chem.RWMol()
    conf = Chem.Conformer(len(atoms))
    index = {}
    for i, a in enumerate(atoms):
        index[a["name"]] = i
        mol.AddAtom(Chem.Atom(a["element"].title()))
        conf.SetAtomPosition(i, a["xyz"])
    conf.Set3D(True)
    mol.AddConformer(conf)
    for a, b, order in bonds:
        if a in index and b in index:
            mol.AddBond(index[a], index[b], {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
                                            "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}[order])
    native = mol.GetMol()
    Chem.SanitizeMol(native)
    Chem.AssignStereochemistryFrom3D(native)
    canon = lambda m: Chem.MolToSmiles(m, isomericSmiles=True)
    if canon(native) != canon(ref) or rdMolDescriptors.CalcMolFormula(native) != formula:
        raise RuntimeError(f"deposited native graph/stereo/formula differs from CCD: {canon(native)} versus {canon(ref)}")
    return native, formula


def csv_bytes(fields: list[str], rows: list[dict]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def parse_index_map(text: str) -> dict[int, int]:
    pairs = {}
    for line in text.splitlines():
        if line.startswith("REMARK INDEX MAP "):
            nums = [int(x) for x in line[17:].split()]
            pairs.update(zip(nums[::2], nums[1::2]))
    if not pairs or len(pairs) != len(set(pairs.values())):
        raise RuntimeError("Meeko index map missing or ambiguous")
    return pairs


def export_native(code: str, cfg: dict, atoms: list[dict], mol: Chem.Mol):
    ligand = cfg["ligand"]
    name = f"{code}_{ligand}"
    map_path = VALIDATION / f"{name}_atom_map.csv"
    map_fields = ["sdf_atom_index_1based", "PDB_serial", "auth_chain", "auth_residue_name",
                  "auth_residue_number", "auth_atom_name", "element", "deposited_x", "deposited_y", "deposited_z"]
    maprows = [dict(zip(map_fields, [i, a["serial"], "A", ligand, cfg["seq"], a["name"], a["element"],
                                      *[f"{v:.3f}" for v in a["xyz"]]])) for i, a in enumerate(atoms, 1)]
    safe_write(map_path, csv_bytes(map_fields, maprows))
    with_h = Chem.AddHs(mol, addCoords=True)
    with_h.SetProp("_Name", f"{name}_deposited_crystal_pose")
    with_h.SetProp("source_PDB", code)
    with_h.SetProp("source_CCD", ligand)
    with_h.SetProp("microstate", "neutral CCD reference; not experimentally assigned")
    with_h.SetProp("coordinate_note", "Deposited heavy coordinates unchanged; H positions generated by RDKit")
    buffer = io.StringIO()
    writer = Chem.SDWriter(buffer)
    writer.write(with_h)
    writer.close()
    sdf = VALIDATION / f"{name}_native_pose.sdf"
    safe_write(sdf, buffer.getvalue().encode("utf-8"))
    loaded = list(Chem.SDMolSupplier(str(sdf), removeHs=False))
    if len(loaded) != 1 or loaded[0] is None:
        raise RuntimeError("SDF readback failed")
    for i, a in enumerate(atoms):
        if math.dist(tuple(loaded[0].GetConformer().GetAtomPosition(i)), a["xyz"]) > 0.0011:
            raise RuntimeError("crystallographic native heavy coordinate moved")
    prep = MoleculePreparation(add_index_map=True, charge_model="gasteiger")
    setups = prep.prepare(loaded[0])
    if len(setups) != 1:
        raise RuntimeError("Meeko ligand setup count differs from one")
    text, ok, err = PDBQTWriterLegacy.write_string(setups[0], add_index_map=True)
    if not ok:
        raise RuntimeError(f"Meeko ligand conversion failed: {err}")
    parsed = PDBQTMolecule(text)
    made = RDKitMolCreate.from_pdbqt_mol(parsed)
    canon = lambda m: Chem.MolToSmiles(Chem.RemoveHs(m), isomericSmiles=True)
    if len(made) != 1 or made[0] is None or canon(made[0]) != canon(mol):
        raise RuntimeError("PDBQT graph/stereo roundtrip differs from CCD")
    index = parse_index_map(text)
    readback = parsed.atoms()
    bridge = []
    for i, a in enumerate(atoms, 1):
        if i not in index:
            raise RuntimeError("PDBQT heavy atom absent in Meeko bridge map")
        serial = index[i]
        if math.dist(readback["xyz"][serial-1], a["xyz"]) > 0.0011:
            raise RuntimeError("PDBQT crystal heavy coordinate differs")
        bridge.append({"sdf_atom_index_1based": i, "PDB_serial": a["serial"],
                       "auth_atom_name": a["name"], "pdbqt_atom_serial": serial})
    bridge_path = VALIDATION / f"{name}_sdf_to_pdbqt_atom_map.csv"
    safe_write(bridge_path, csv_bytes(list(bridge[0]), bridge))
    pdbqt = VALIDATION / f"{name}_crystal_pose_reference.pdbqt"
    safe_write(pdbqt, text.encode("utf-8"))
    tors = int(next(line.split()[1] for line in text.splitlines() if line.startswith("TORSDOF ")))
    return sdf, map_path, pdbqt, bridge_path, tors


def atom_pos(line: str) -> tuple[float, float, float]:
    return tuple(float(line[i:i+8]) for i in (30, 38, 46))


def residue_key(line: str):
    return line[21], line[22:26].strip(), line[26], line[17:20].strip()


def choose_altlocs(lines: list[str]):
    scores = defaultdict(lambda: defaultdict(list))
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and line[16].strip():
            scores[residue_key(line)][line[16]].append(float(line[54:60]))
    return {res: sorted(choices, key=lambda a: (-sum(choices[a])/len(choices[a]), a))[0]
            for res, choices in scores.items()}


def working_copy(code: str, lines: list[str], native_serials: set[int]) -> Path:
    path = PREPARED / f"{code}_vina_validation_working_no_native.pdb"
    kept = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and int(line[6:11]) in native_serials:
            continue
        if line.startswith("CONECT"):
            nums = [int(line[i:i+5]) for i in range(6, len(line), 5) if line[i:i+5].strip().isdigit()]
            nums = [n for n in nums if n not in native_serials]
            if len(nums) < 2:
                continue
            line = "CONECT" + "".join(f"{n:5d}" for n in nums)
        kept.append(line)
    safe_write(path, ("\n".join(kept) + "\n").encode("ascii"))
    return path


def make_local_coordinate_model(code: str, lines: list[str], native: list[dict]):
    choices = choose_altlocs(lines)
    native_xyz = [a["xyz"] for a in native]
    chloride = next((atom_pos(s) for s in lines if s.startswith("HETATM") and s[17:20].strip() == "CL"), None)
    keep_water = set()
    for line in lines:
        if not line.startswith("HETATM") or line[17:20] != "HOH":
            continue
        xyz = atom_pos(line)
        if min(math.dist(xyz, p) for p in native_xyz) <= 4.0 or (
                chloride is not None and math.dist(xyz, chloride) <= 3.6):
            keep_water.add(residue_key(line))
    kept = []
    prior_protein_num = None
    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        res = line[17:20].strip()
        if res in ("NAG", "GOL", "SO4") or (res == "HOH" and residue_key(line) not in keep_water):
            continue
        if res == CONFIG[code]["ligand"]:
            continue
        if line[16].strip() and line[16] != choices[residue_key(line)]:
            continue
        if line[16].strip():
            line = line[:16] + " " + line[17:]
        if code == "3L4W" and line.startswith("ATOM  "):
            num = int(line[22:26])
            if prior_protein_num == 836 and num == 838:
                kept.append("TER")
            prior_protein_num = num
        kept.append(line)
    header = [line for line in lines if line.startswith("CRYST1")]
    text = "\n".join([*header, *kept, "END"]) + "\n"
    return text, choices, keep_water


def add_provisional_hydrogens(pdb_text: str, pH: float) -> str:
    import random
    random.seed(20261012)
    np.random.seed(20261012)
    pdb = PDBFile(io.StringIO(pdb_text))
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(pH=pH, platform=openmm.Platform.getPlatformByName("Reference"))
    hydrogens = defaultdict(list)
    for atom, position in zip(modeller.topology.atoms(), modeller.positions):
        if atom.element.symbol != "H":
            continue
        r = atom.residue
        key = (r.chain.id, r.id, getattr(r, "insertionCode", " ") or " ", r.name)
        xyz = position.value_in_unit(openmm.unit.angstrom)
        hydrogens[key].append((atom.name, xyz))
    lines = pdb_text.splitlines()
    by_res = defaultdict(list)
    last_index = {}
    for i, line in enumerate(lines):
        if line.startswith(("ATOM  ", "HETATM")):
            key = residue_key(line)
            by_res[key].append(line)
            last_index[key] = i
    if set(hydrogens) - set(by_res):
        raise RuntimeError("added hydrogen without deposited heavy-atom residue")
    serial = max(int(s[6:11]) for s in lines if s.startswith(("ATOM  ", "HETATM"))) + 1
    out = [f"REMARK 901 OPENMM {openmm.__version__} PROVISIONAL HYDROGENS PH {pH:.1f}; HEAVY COORDINATES UNCHANGED"]
    for i, line in enumerate(lines):
        if line == "END":
            continue
        out.append(line)
        if line.startswith(("ATOM  ", "HETATM")) and last_index[residue_key(line)] == i:
            key = residue_key(line)
            for name, xyz in hydrogens.get(key, []):
                if serial > 99999:
                    raise RuntimeError("PDB atom serial overflow")
                template = by_res[key][0]
                hline = template[:6] + f"{serial:5d} {name:4s} " + template[17:30]
                hline += "".join(f"{float(v):8.3f}" for v in xyz)
                hline += f"{1.00:6.2f}{0.00:6.2f}          {'H':>2s}  "
                out.append(hline)
                serial += 1
    out.append("END")
    return "\n".join(out) + "\n"


def receptor_probe(code: str, source: Path) -> tuple[str, str, str]:
    command = [str(CLI), "--read_pdb", str(source)]
    if code == "4GQR":
        command += ["--add_templates", str(PCA_TEMPLATE)]
    with tempfile.TemporaryDirectory() as tmp:
        target = Path(tmp) / "probe.pdbqt"
        result = subprocess.run(command + ["--write_pdbqt", str(target)], capture_output=True,
                                text=True, timeout=240)
        if result.returncode or not target.is_file():
            return "blocked_meeko", "", " ".join((result.stderr + " " + result.stdout).split())[:1500]
        parsed = PDBQTReceptor.from_pdbqt_filename(str(target))
        arr = parsed.atoms()
        if not len(arr) or not np.isfinite(arr["xyz"]).all() or not np.isfinite(arr["partial_charges"]).all():
            return "failed_QC", "", "PDBQT nonfinite or empty"
        source_xyz = Counter(tuple(round(float(s[i:i+8]), 3) for i in (30, 38, 46))
                             for s in source.read_text().splitlines()
                             if s.startswith(("ATOM  ", "HETATM")) and s[76:78].strip().upper() not in ("H", "D"))
        output_xyz = Counter(tuple(round(float(v), 3) for v in atom["xyz"])
                             for atom in arr if atom["atom_type"] != "HD")
        if source_xyz != output_xyz:
            return "failed_QC", "", f"heavy-coordinate multiset differs: source={sum(source_xyz.values())}, PDBQT={sum(output_xyz.values())}"
        path = VINA_PREPARED / f"{code}_validation_local_noNAG_Meeko0.8.0.pdbqt"
        safe_write(path, target.read_bytes())
        return "preliminary_local_model", rel(path), "all heavy coordinates preserved to 0.001 A; finite PDBQT atoms/charges"


def receptor_details(code: str, pdb: Path, pdbqt: Path):
    lines = pdb.read_text().splitlines()
    qlines = pdbqt.read_text().splitlines()
    heavy = [s for s in lines if s.startswith(("ATOM  ", "HETATM")) and s[76:78].strip().upper() not in ("H", "D")]
    qatoms = [s for s in qlines if s.startswith(("ATOM  ", "HETATM"))]
    ions = [f"{s[17:20].strip()} A{s[22:26].strip()}:{s.split()[-1]}:{s.split()[-2]}"
            for s in qatoms if s[17:20].strip() in ("CA", "CL")]
    water_bonds = []
    for oxygen in (s for s in lines if s.startswith("HETATM") and s[17:20].strip() == "HOH" and s[76:78].strip() == "O"):
        for hydrogen in (s for s in lines if s.startswith("HETATM") and s[17:20].strip() == "HOH"
                         and s[22:26] == oxygen[22:26] and s[76:78].strip() == "H"):
            water_bonds.append(math.dist(atom_pos(oxygen), atom_pos(hydrogen)))
    if not water_bonds:
        raise RuntimeError("selected receptor waters lost hydrogens")
    cat = [("ASP", "197"), ("GLU", "233"), ("ASP", "300"), ("HIS", "201")] if code == "4GQR" else [
        ("ASP", "443"), ("ASP", "542"), ("HIS", "600")]
    micro = []
    for res, seq in cat:
        h = {s[12:16].strip() for s in lines if s.startswith(("ATOM  ", "HETATM"))
             and s[17:20].strip() == res and s[22:26].strip() == seq and s[76:78].strip() == "H"}
        if res in ("ASP", "GLU"):
            status = "no sidechain acidic H" if not (h & {"HD1", "HD2", "HE1", "HE2"}) else "sidechain H present"
        else:
            status = "HD1=" + str("HD1" in h) + ";HE2=" + str("HE2" in h)
        micro.append(f"{res}{seq}:{status}")
    fixer = PDBFixer(filename=str(pdb))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    missing = f"missing standard sidechain residues={len(fixer.missingAtoms)}; terminal={[(r.name,r.id,names) for r,names in fixer.missingTerminals.items()]}"
    return {"receptor_heavy_atoms": len(heavy), "receptor_PDBQT_atoms": len(qatoms),
            "ion_PDBQT_fields": ";".join(ions),
            "water_OH_range_A": f"{min(water_bonds):.3f}-{max(water_bonds):.3f}",
            "catalytic_microstates": ";".join(micro), "missing_atom_QC": missing}


def run_one(code: str, cfg: dict) -> dict:
    ligand = cfg["ligand"]
    pdb = ORIGINAL / f"{code}.pdb"
    cif = ORIGINAL / f"{code}.cif"
    ccd = DEFINITIONS / f"{ligand}.cif"
    get_official(f"https://files.rcsb.org/download/{code}.pdb", pdb, cfg["pdb_sha"])
    get_official(f"https://files.rcsb.org/download/{code}.cif", cif, cfg["cif_sha"])
    get_official(f"https://files.rcsb.org/ligands/download/{ligand}.cif", ccd, cfg["ccd_sha"])
    lines = pdb.read_text().splitlines()
    native = source_line_atoms(lines, ligand, cfg["seq"])
    if len(native) != cfg["expected_heavy"] or any(a["occ"] <= 0 for a in native):
        raise RuntimeError(f"{code} incomplete/zero-occupancy native inhibitor")
    verify_cif_coords(cif, native, ligand, cfg["seq"])
    def closest_component(component: str) -> float | None:
        positions = [atom_pos(s) for s in lines if s.startswith("HETATM") and s[17:20].strip() == component]
        return min((math.dist(p, a["xyz"]) for p in positions for a in native), default=None)
    nag_distance = closest_component("NAG")
    gol_distance = closest_component("GOL")
    if nag_distance is None or nag_distance < (30.0 if code == "4GQR" else 10.0):
        raise RuntimeError("NAG too close or absent; local glycan exclusion needs re-assessment")
    if gol_distance is not None and gol_distance < 10.0:
        raise RuntimeError("glycerol too close to native site for this local exclusion rule")
    mol, formula = reference_molecule(native, ccd)
    sdf, atommap, native_pdbqt, bridge, tors = export_native(code, cfg, native, mol)
    working = working_copy(code, lines, {a["serial"] for a in native})
    coordinate, choices, water_set = make_local_coordinate_model(code, lines, native)
    local = PREPARED / f"{code}_validation_local_noNAG_pH{str(cfg['pH']).replace('.', 'p')}.pdb"
    status, path, qc = "blocked", "", ""
    issue = ""
    try:
        protonated = add_provisional_hydrogens(coordinate, cfg["pH"])
        safe_write(local, protonated.encode("ascii"))
        pca_check = "not deposited"
        if code == "4GQR":
            pca_check = check_deposited_pca_stereo(PCA_TEMPLATE.read_bytes(), local)
            pca_c = next(atom_pos(s) for s in lines if s.startswith("HETATM") and s[17:20] == "PCA"
                         and s[22:26].strip() == "1" and s[12:16].strip() == "C")
            tyr_n = next(atom_pos(s) for s in lines if s.startswith("ATOM  ") and s[17:20] == "TYR"
                         and s[22:26].strip() == "2" and s[12:16].strip() == "N")
            link_distance = math.dist(pca_c, tyr_n)
            if not 1.2 < link_distance < 1.5:
                raise RuntimeError("PCA-Tyr2 peptide geometry invalid")
            pca_check += f"; PCA C-Tyr2 N {link_distance:.3f} A"
        status, path, qc = receptor_probe(code, local)
        if status != "preliminary_local_model":
            issue = qc
    except Exception as exc:
        status = "blocked_preparation"
        issue = f"{type(exc).__name__}: {exc}"
        pca_check = "not completed" if code == "4GQR" else "not deposited"
    center = tuple(sum(a["xyz"][i] for a in native) / len(native) for i in range(3))
    glycan = ("linked NAG A501 omitted only from local model; nearest MYC ~50.9 A" if code == "4GQR" else
              "six protein-linked NAG units omitted only from local model; nearest MIG ~16.4 A")
    blocker = ("Standard Vina Ca/Cl scoring is not physically validated; catalytic/ligand protonation and selected waters provisional. " if code == "4GQR" else
               "NAG-free receptor is a local-pocket sensitivity approximation; Gln837 is missing; miglitol amine and catalytic protonation provisional. ") + issue
    row = {key: "" for key in FIELDS}
    details = receptor_details(code, local, ROOT / path) if path else {}
    row.update({"PDB_ID": code, "protein": "human pancreatic alpha-amylase" if code == "4GQR" else "human ntMGAM",
                "native_ligand": ligand, "original_PDB_path": rel(pdb), "original_PDB_SHA256": sha(pdb),
                "original_mmCIF_path": rel(cif), "original_mmCIF_SHA256": sha(cif),
                "CCD_path": rel(ccd), "CCD_SHA256": sha(ccd),
                "native_SDF_path": rel(sdf), "native_SDF_SHA256": sha(sdf),
                "native_PDBQT_path": rel(native_pdbqt), "native_PDBQT_SHA256": sha(native_pdbqt),
                "native_atom_map_path": rel(atommap), "native_atom_map_SHA256": sha(atommap),
                "native_bridge_map_path": rel(bridge), "native_bridge_map_SHA256": sha(bridge),
                "native_heavy_atoms": len(native), "native_torsions": tors,
                "native_graph_stereo_pose_check": f"pass CCD graph/stereo/formula {formula}; PDB/mmCIF and PDBQT heavy coordinates preserved",
                "native_microstate_status": "neutral CCD reference; not a measured protonation assignment",
                "site_center_xyz_A": ";".join(f"{v:.3f}" for v in center),
                "working_receptor_path": rel(working), "working_receptor_SHA256": sha(working),
                "local_receptor_PDB_path": rel(local) if local.exists() else "",
                "local_receptor_PDB_SHA256": sha(local) if local.exists() else "",
                "local_receptor_PDBQT_path": path,
                "local_receptor_PDBQT_SHA256": sha(ROOT / path) if path else "",
                "receptor_conversion_status": status, "receptor_heavy_coordinate_check": qc if path else "",
                "prep_pH": cfg["pH"], "glycan_decision": glycan,
                "glycan_closest_A": f"{nag_distance:.3f}",
                "glycerol_closest_A": f"{gol_distance:.3f}" if gol_distance is not None else "not deposited",
                "ions_retained": "Ca A502 and Cl A503" if code == "4GQR" else "none deposited",
                "waters_retained": ";".join(sorted((k[1] for k in water_set), key=int)),
                "PCA_link_stereo_check": pca_check,
                "missing_residue_status": "496/496 modeled; PCA A1 retained" if code == "4GQR" else
                                          "863/875 modeled; residues 1-6 and 837 and 871-875 missing; TER at 836/838 gap",
                "altloc_rule": f"highest mean occupancy per residue then alphabetic tie; {len(choices)} residues affected",
                "template_path": rel(PCA_TEMPLATE) if code == "4GQR" else "none",
                "template_SHA256": sha(PCA_TEMPLATE) if code == "4GQR" else "",
                "software_versions": f"RDKit {Chem.rdBase.rdkitVersion}; Gemmi {gemmi.__version__}; OpenMM {openmm.__version__}; Meeko {meeko.__version__}",
                "redocking_ready": "no", "blockers": blocker.strip()})
    row.update(details)
    return row


def main() -> None:
    if meeko.__version__ != "0.8.0" or not CLI.is_file():
        raise RuntimeError("expected project Meeko 0.8.0 CLI")
    rows = []
    for code, cfg in CONFIG.items():
        try:
            row = run_one(code, cfg)
            rows.append(row)
            print(code, row["native_ligand"], row["receptor_conversion_status"], row["blockers"])
        except Exception as exc:
            row = {key: "" for key in FIELDS}
            row.update({"PDB_ID": code, "native_ligand": cfg["ligand"], "receptor_conversion_status": "blocked",
                        "redocking_ready": "no", "blockers": f"{type(exc).__name__}: {exc}"})
            rows.append(row)
            print(code, "BLOCKED", row["blockers"])
    safe_write(RESULTS, csv_bytes(FIELDS, rows))


if __name__ == "__main__":
    main()

"""Reproduce labelled 3L4W local-pocket chemistry sensitivity inputs.

Run with .venv/Scripts/python.exe scripts/review_3l4w_chemistry.py.
This script does not dock, score, alter the deposited coordinates, or replace
the original 3L4W and MIG files. It fails if audited source hashes change.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from pathlib import Path
import math
import subprocess
import tempfile

from meeko import MoleculePreparation, PDBQTMolecule, PDBQTReceptor, PDBQTWriterLegacy, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PDB = ROOT / "structures/proteins/prepared/3L4W_validation_local_noNAG_pH6p5.pdb"
SOURCE_SDF = ROOT / "docking/validation/vina/3L4W_MIG_native_pose.sdf"
RECEPTOR_DIR = ROOT / "structures/proteins/vina_prepared"
PREPARED_DIR = ROOT / "structures/proteins/prepared"
LIGAND_DIR = ROOT / "docking/validation/vina"
MEEKO_CLI = ROOT / ".venv/Scripts/mk_prepare_receptor.exe"
VINA_EXE = Path(r"E:\docking tools\vina_1.2.7_win.exe")
NEUTRAL_MIG_PDBQT = ROOT / "docking/validation/vina/3L4W_MIG_crystal_pose_reference.pdbqt"
EXPECTED = {
    ROOT / "structures/proteins/original/3L4W.pdb": "57638e45666e2f848e37691a914a978d85c377c13eaed9e5f94af24d2c52c95e",
    ROOT / "structures/proteins/original/3L4W.cif": "81f5582ab3356b1e6318c3474367268712c5b8c2886193598564fe40e0218c46",
    ROOT / "docking/validation/definitions/MIG.cif": "4420e2e083b312422cc2e08833a297c68ae25fdebcde01d180daeef56df28b28",
    ROOT / "docking/validation/vina/3L4W_MIG_atom_map.csv": "ff1c77265befde60237c5a1fc129172b40620c9f45beeb264679f25180729df9",
    ROOT / "docking/validation/vina/3L4W_MIG_sdf_to_pdbqt_atom_map.csv": "68c9254203234c64557de148d64e161e2a4cf972c3f1e156e3b9658a0c8fcd50",
    NEUTRAL_MIG_PDBQT: "c7417f6bcb661e0f2217c6ba3619b6fd6c03ffb029fe5840d4b0d13c1bf16169",
    SOURCE_PDB: "f4b176042ed743d957dccb2c6b1206069b4d9d50f7f56c393ecbac8982ee437d",
    SOURCE_SDF: "ea673a59043bea05abbe7d66d8bd29c3786021f564c394027a8b0201385503f7",
}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def write_new(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != data:
            raise RuntimeError(f"Refusing to overwrite an existing variant: {path}")
        return
    path.write_bytes(data)


def atom_xyz(line: str) -> tuple[float, float, float]:
    return tuple(float(line[i : i + 8]) for i in (30, 38, 46))


def prepare_receptor(wet: bool) -> tuple[Path, Path]:
    tag = "sitewaters" if wet else "dry"
    stem = f"3L4W_noNAG_D443m_D542m_H600HID_pH6p5_{tag}_sensitivity"
    pdb_path = PREPARED_DIR / f"{stem}.pdb"
    pdbqt_path = RECEPTOR_DIR / f"{stem}_Meeko0.8.0.pdbqt"
    source = SOURCE_PDB.read_text().splitlines()
    out = [
        "REMARK 900 3L4W LOCAL NO-NAG SENSITIVITY; D443-/D542-/H600-HID; PH 6.5",
        f"REMARK 901 {'HOH4002/4004 RETAINED' if wet else 'ALL HOH EXCLUDED'}; DEPOSITED HEAVY COORDINATES FIXED",
    ]
    removed = Counter()
    for line in source:
        if line.startswith("REMARK 901"):
            continue
        if line.startswith(("ATOM  ", "HETATM")):
            if line[17:20] == "HIS" and line[22:26].strip() == "600" and line[12:16].strip() == "HE2":
                removed["HE2"] += 1
                continue
            if not wet and line[17:20] == "HOH":
                removed["water atoms"] += 1
                continue
        out.append(line)
    if removed["HE2"] != 1 or (not wet and removed["water atoms"] != 6):
        raise RuntimeError(f"Unexpected H600/water source composition: {removed}")
    text = "\n".join(out) + "\n"
    original_heavy = Counter(atom_xyz(s) for s in source if s.startswith("ATOM  ") and s[76:78].strip() != "H")
    variant_heavy = Counter(atom_xyz(s) for s in out if s.startswith("ATOM  ") and s[76:78].strip() != "H")
    if original_heavy != variant_heavy:
        raise RuntimeError("Protein heavy coordinates changed")
    with tempfile.TemporaryDirectory() as td:
        tmp_pdb = Path(td) / "input.pdb"
        tmp_qt = Path(td) / "input.pdbqt"
        tmp_pdb.write_text(text)
        cmd = [str(MEEKO_CLI), "--read_pdb", str(tmp_pdb), "--write_pdbqt", str(tmp_qt)]
        run = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if run.returncode or not tmp_qt.exists():
            raise RuntimeError(f"Meeko receptor conversion failed: {run.stdout[-1000:]} {run.stderr[-1000:]}")
        receptor = PDBQTReceptor.from_pdbqt_filename(str(tmp_qt))
        atoms = receptor.atoms()
        if not len(atoms):
            raise RuntimeError("Empty receptor PDBQT")
        qt_heavy = Counter(tuple(round(float(v), 3) for v in a["xyz"])
                           for a in atoms if a["atom_type"] != "HD")
        expected_heavy = Counter(atom_xyz(s) for s in out if s.startswith(("ATOM  ", "HETATM"))
                                 and s[76:78].strip() not in ("H", "D"))
        if qt_heavy != expected_heavy:
            raise RuntimeError(f"Receptor heavy coordinate mismatch: {sum(qt_heavy.values())} vs {sum(expected_heavy.values())}")
        if not all(float(a["partial_charges"]) == float(a["partial_charges"]) for a in atoms):
            raise RuntimeError("Nonfinite receptor charge")
        write_new(pdb_path, text.encode("ascii"))
        write_new(pdbqt_path, tmp_qt.read_bytes())
    return pdb_path, pdbqt_path


def prepare_protonated_mig() -> tuple[Path, Path]:
    mol = Chem.SDMolSupplier(str(SOURCE_SDF), removeHs=False)[0]
    if mol is None:
        raise RuntimeError("Native MIG SDF cannot be read")
    neutral = Chem.RemoveHs(mol)
    original_cips = Chem.FindMolChiralCenters(neutral, includeUnassigned=True)
    original_xyz = [tuple(mol.GetConformer().GetAtomPosition(i)) for i in range(14)]
    amines = [a for a in neutral.GetAtoms() if a.GetSymbol() == "N"]
    if len(amines) != 1 or amines[0].GetIdx() != 13 or amines[0].GetFormalCharge() != 0:
        raise RuntimeError("Unexpected MIG amine identity")
    edit = Chem.RWMol(neutral)
    edit.GetAtomWithIdx(13).SetFormalCharge(1)
    charged = edit.GetMol()
    Chem.SanitizeMol(charged)
    charged = Chem.AddHs(charged, addCoords=True)
    if rdMolDescriptors.CalcMolFormula(charged) != "C8H18NO5+":
        raise RuntimeError("Incorrect protonated MIG formula")
    charged_cips = Chem.FindMolChiralCenters(Chem.RemoveHs(charged), includeUnassigned=True)
    if [x for x in charged_cips if x[0] != 13] != original_cips:
        raise RuntimeError("MIG carbon stereocentres changed")
    # Protonating tertiary N creates an unassigned pyramidal N stereocentre in
    # RDKit. Do not invent a fixed N configuration from the X-ray heavy atoms.
    if (13, "?") not in charged_cips:
        raise RuntimeError("Unexpected protonated amine stereochemical handling")
    coords = [tuple(charged.GetConformer().GetAtomPosition(i)) for i in range(14)]
    if any(sum((a - b) ** 2 for a, b in zip(p, q)) > 1e-8 for p, q in zip(original_xyz, coords)):
        raise RuntimeError("MIG deposited heavy coordinates changed")
    charged.SetProp("_Name", "3L4W_MIG_NHplus_crystal_heavy_pose_pH6p5_sensitivity")
    charged.SetProp("MICROSTATE", "tertiary amine +1; hypothesis, not CCD-deposited state")
    prep = MoleculePreparation(add_index_map=True, charge_model="gasteiger")
    setups = prep.prepare(charged)
    if len(setups) != 1:
        raise RuntimeError("MIG Meeko preparation did not yield exactly one setup")
    pdbqt, ok, error = PDBQTWriterLegacy.write_string(setups[0], add_index_map=True)
    if not ok:
        raise RuntimeError(f"Meeko MIG conversion failed: {error}")
    parsed = PDBQTMolecule(pdbqt)
    roundtrips = [m for m in RDKitMolCreate.from_pdbqt_mol(parsed) if m is not None]
    if len(roundtrips) != 1 or Chem.GetFormalCharge(roundtrips[0]) != 1:
        raise RuntimeError("MIG PDBQT charge readback failed")
    if sorted(label for _, label in Chem.FindMolChiralCenters(Chem.RemoveHs(roundtrips[0]), includeUnassigned=True)
              if label != "?") != sorted(label for _, label in original_cips):
        raise RuntimeError("MIG PDBQT carbon stereochemistry readback failed")
    heavy_lines = [line for line in pdbqt.splitlines()
                   if line.startswith(("ATOM  ", "HETATM")) and line.split()[-1] != "HD"]
    if len(heavy_lines) != 14:
        raise RuntimeError(f"MIG PDBQT heavy-atom count mismatch: {len(heavy_lines)}")
    qt_coords = Counter(tuple(float(line[i : i + 8]) for i in (30, 38, 46)) for line in heavy_lines)
    if qt_coords != Counter(tuple(round(v, 3) for v in xyz) for xyz in original_xyz):
        raise RuntimeError("MIG PDBQT deposited heavy coordinates changed")
    sdf_path = LIGAND_DIR / "3L4W_MIG_NHplus_crystal_pose_sensitivity.sdf"
    pdbqt_path = LIGAND_DIR / "3L4W_MIG_NHplus_crystal_pose_sensitivity.pdbqt"
    write_new(sdf_path, (Chem.MolToMolBlock(charged) + "\n$$$$\n").encode("ascii"))
    write_new(pdbqt_path, pdbqt.encode("ascii"))
    return sdf_path, pdbqt_path


def check_intercomponent_clashes(receptor: Path, ligand_sdf: Path) -> float:
    mol = Chem.SDMolSupplier(str(ligand_sdf), removeHs=False)[0]
    if mol is None:
        raise RuntimeError(f"Cannot read ligand for clash check: {ligand_sdf}")
    atoms = [line for line in receptor.read_text().splitlines()
             if line.startswith(("ATOM  ", "HETATM"))]
    heavy_min = float("inf")
    h_min = float("inf")
    conformer = mol.GetConformer()
    for atom in mol.GetAtoms():
        xyz = tuple(conformer.GetAtomPosition(atom.GetIdx()))
        ligand_h = atom.GetSymbol() == "H"
        for line in atoms:
            distance = math.dist(xyz, atom_xyz(line))
            if ligand_h or line[76:78].strip() == "H":
                h_min = min(h_min, distance)
            else:
                heavy_min = min(heavy_min, distance)
    if heavy_min < 2.0 or h_min < 1.2:
        raise RuntimeError(f"Intercomponent clash in {receptor.name} / {ligand_sdf.name}: "
                           f"heavy={heavy_min:.3f}, H-containing={h_min:.3f} A")
    print("Intercomponent distance screen:", receptor.name, ligand_sdf.name,
          f"heavy={heavy_min:.3f} A", f"H-containing={h_min:.3f} A")
    return h_min


def main() -> None:
    for source, expected in EXPECTED.items():
        if digest(source) != expected:
            raise RuntimeError(f"Audited source hash changed: {source}")
    outputs = [*prepare_receptor(False), *prepare_receptor(True), *prepare_protonated_mig()]
    for output in outputs:
        print(output.relative_to(ROOT), digest(output))
    for receptor in (outputs[0], outputs[2]):
        for ligand_sdf in (SOURCE_SDF, outputs[4]):
            check_intercomponent_clashes(receptor, ligand_sdf)
    if VINA_EXE.is_file():
        # In Vina 1.2.7 main.cpp, set_receptor and set_ligand_from_file run
        # before compute_vina_maps checks grid dimensions. An intentionally
        # invalid zero-size box therefore tests *parsing only* and aborts
        # before map generation, score calculation, or global search.
        for receptor in (outputs[1], outputs[3]):
            for ligand in (NEUTRAL_MIG_PDBQT, outputs[5]):
                with tempfile.TemporaryDirectory() as td:
                    forbidden_output = Path(td) / "must_not_exist.pdbqt"
                    cmd = [str(VINA_EXE), "--receptor", str(receptor), "--ligand", str(ligand),
                           "--center_x", "45.387", "--center_y", "92.352", "--center_z", "34.822",
                           "--size_x", "0", "--size_y", "0", "--size_z", "0",
                           "--out", str(forbidden_output), "--cpu", "1", "--verbosity", "0"]
                    run = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
                    message = run.stdout + run.stderr
                    if run.returncode != 1 or "Grid box dimensions must be greater than 0 Angstrom" not in message or forbidden_output.exists():
                        raise RuntimeError(f"Unexpected Vina parse-only probe result: {message[:1200]}")
                print("Vina 1.2.7 parse reached grid guard:", receptor.name, ligand.name)
    else:
        print("Vina executable absent; parser probe not performed")


if __name__ == "__main__":
    main()

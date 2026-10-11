"""Stage receptor protonation/water variants; never claim GOLD import validation.

Requires OpenMM 8.6.1, PDBFixer 1.12.0, RDKit, and the unchanged files
recorded in receptor_preparation_results.csv. Coordinates of retained heavy
atoms are copied exactly from the previous ligand-free working PDB files.
"""

from __future__ import annotations

import csv
import hashlib
import math
import random
from collections import defaultdict
from pathlib import Path

import numpy as np
import openmm
from openmm.app import Modeller, PDBFile
from pdbfixer import PDBFixer
from rdkit import Chem


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/metadata/receptor_preparation_results.csv"
OUTPUT = ROOT / "data/metadata/gold_receptor_preparation_results.csv"
PDB_OUT = ROOT / "structures/proteins/prepared"
SEED = 20261011
WATER_CUTOFF_A = 4.0
VARIANTS = {
    "2QMJ": [
        ("allwaters_keepNAG", "all", True),
        ("sitewaters_keepNAG", "site", True),
        ("sitewaters_noNAG", "site", False),
    ],
    "1XD0": [
        ("allwaters_CA_CL", "all", False),
        ("sitewaters_CA_CL", "site", False),
    ],
}
FIELDS = [
    "PDB_ID", "variant", "input_receptor_path", "input_SHA256", "native_reference_path",
    "native_reference_SHA256", "atom_map_path", "atom_map_SHA256", "output_path",
    "output_SHA256", "prep_tool", "prep_pH", "altloc_rule", "altloc_residues",
    "standard_sidechain_missing_heavy_atoms", "missing_terminal_atoms", "sequence_gaps",
    "water_rule", "waters_retained", "waters_removed", "site_waters_retained",
    "calcium_waters_retained", "ions_retained", "glycan_decision", "glycerol_decision",
    "sulfate_decision", "PCA_decision", "protein_hydrogens_added", "water_hydrogens_added",
    "hetero_hydrogen_limitation", "catalytic_acid_state", "site_histidine_state",
    "heavy_coordinates_preserved", "reference_mapping_valid", "GOLD_import_test",
    "chemical_preparation_status", "redocking_ready", "issues",
]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atom_key(line: str) -> tuple[str, str, str, str, str]:
    return line[21], line[22:26].strip(), line[26], line[17:20].strip(), line[12:16].strip()


def residue_key(line: str) -> tuple[str, str, str, str]:
    return atom_key(line)[:4]


def pos(line: str) -> tuple[float, float, float]:
    return tuple(float(line[i:i + 8]) for i in (30, 38, 46))


def reference_atoms(row: dict[str, str]) -> list[tuple[float, float, float]]:
    ref = ROOT / row["native_ligand_SDF"]
    mapping = ROOT / row["atom_mapping_CSV"]
    if digest(ref) != row["native_ligand_SHA256"] or digest(mapping) != row["atom_mapping_SHA256"]:
        raise ValueError(f'{row["PDB_ID"]}: native ligand or map changed')
    mol = Chem.SDMolSupplier(str(ref), removeHs=False)[0]
    if mol is None:
        raise ValueError(f'{row["PDB_ID"]}: native SDF cannot be read')
    atoms = [mol.GetConformer().GetAtomPosition(i) for i, a in enumerate(mol.GetAtoms()) if a.GetAtomicNum() > 1]
    maprows = list(csv.DictReader(mapping.open(newline="", encoding="utf-8")))
    if len(atoms) != len(maprows) or len(atoms) != int(row["expected_heavy_atoms"]):
        raise ValueError(f'{row["PDB_ID"]}: native map length changed')
    for i, (p, m) in enumerate(zip(atoms, maprows), 1):
        if int(m["sdf_atom_index_1based"]) != i:
            raise ValueError("Native mapping index is not sequential")
        if math.dist(p, (float(m["deposited_x"]), float(m["deposited_y"]), float(m["deposited_z"]))) > 0.0011:
            raise ValueError("Native mapping coordinates disagree with reference SDF")
    return [tuple(p) for p in atoms]


def choose_altlocs(lines: list[str]) -> dict[tuple[str, str, str, str], str]:
    scores = defaultdict(lambda: defaultdict(list))
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and line[16].strip():
            scores[residue_key(line)][line[16]].append(float(line[54:60]))
    return {k: sorted(v, key=lambda label: (-sum(v[label]) / len(v[label]), label))[0] for k, v in scores.items()}


def filter_receptor(code: str, source: list[str], ligand: list[tuple[float, float, float]], water_rule: str, keep_nag: bool):
    choices = choose_altlocs(source)
    keep_waters = set()
    for line in source:
        if line.startswith("HETATM") and line[17:20] == "HOH":
            d = min(math.dist(pos(line), p) for p in ligand)
            if water_rule == "all" or d <= WATER_CUTOFF_A or (code == "1XD0" and line[22:26].strip() in ("536", "546")):
                keep_waters.add((line[21], line[22:26].strip(), line[26]))
    filtered = []
    original_serials = set()
    expected_heavy = {}
    for line in source:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        rk = residue_key(line)
        alt = line[16]
        if alt.strip() and alt != choices[rk]:
            continue
        name = line[17:20].strip()
        if name == "HOH" and (line[21], line[22:26].strip(), line[26]) not in keep_waters:
            continue
        if name in ("GOL", "SO4") or (name == "NAG" and not keep_nag):
            continue
        if alt.strip():
            line = line[:16] + " " + line[17:]
        filtered.append(line)
        original_serials.add(int(line[6:11]))
        expected_heavy[atom_key(line)] = pos(line)
    return filtered, original_serials, expected_heavy, choices, keep_waters


def build_coordinate_pdb(source: list[str], filtered: list[str], old_serials: set[int], choices: dict) -> str:
    kept = {int(line[6:11]): line for line in filtered}
    atom_ids = {atom_key(line) for line in filtered}
    result = []
    previous_protein_residue = None
    for line in source:
        if line.startswith(("ATOM  ", "HETATM")):
            serial = int(line[6:11])
            if serial in kept:
                if line.startswith("ATOM  ") and line[21] == "A":
                    this_residue = line[22:26].strip()
                    if previous_protein_residue == "836" and this_residue == "838":
                        # Construct residue 837 is absent.  Without a chain
                        # break OpenMM invents a 4.13 A peptide bond.
                        result.append("TER")
                    previous_protein_residue = this_residue
                result.append(kept[serial])
            continue
        if line.startswith("LINK  "):
            endpoints = [(line[21], line[22:26].strip(), line[26], line[17:20].strip(), line[12:16].strip(), line[16]),
                         (line[51], line[52:56].strip(), line[56], line[47:50].strip(), line[42:46].strip(), line[46])]
            if any(e[:5] not in atom_ids or (e[5].strip() and choices.get(e[:4]) != e[5]) for e in endpoints):
                continue
            line = line[:16] + " " + line[17:46] + " " + line[47:]
        if line.startswith("CONECT"):
            ids = [int(line[i:i + 5]) for i in range(6, len(line), 5) if line[i:i + 5].strip().isdigit()]
            ids = [n for n in ids if n in old_serials]
            if len(ids) < 2:
                continue
            line = "CONECT" + "".join(f"{n:5d}" for n in ids)
        if line.startswith("END"):
            continue
        result.append(line)
    result.append("END")
    return "\n".join(result) + "\n"


def protonated_pdb(pdb_text: str, expected_heavy: dict, code: str):
    import io

    p = PDBFile(io.StringIO(pdb_text))
    modeller = Modeller(p.topology, p.positions)
    random.seed(SEED)
    np.random.seed(SEED)
    modeller.addHydrogens(pH=7.0, platform=openmm.Platform.getPlatformByName("Reference"))
    hydrogen_lines = defaultdict(list)
    h_counts = defaultdict(int)
    protonation = {}
    for atom, xyz in zip(modeller.topology.atoms(), modeller.positions):
        if atom.element.symbol != "H":
            continue
        r = atom.residue
        rk = (r.chain.id, r.id, getattr(r, "insertionCode", " ") or " ", r.name)
        h_counts[r.name] += 1
        hydrogen_lines[rk].append((atom.name, xyz.value_in_unit(__import__("openmm").unit.angstrom)))
        if r.name in ("ASP", "GLU", "HIS"):
            protonation.setdefault((r.name, r.id), []).append(atom.name)
    lines = pdb_text.splitlines()
    atom_lines = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
    actual_heavy = {atom_key(line): pos(line) for line in atom_lines}
    if actual_heavy != expected_heavy or len(atom_lines) != len(expected_heavy):
        raise ValueError(f"{code}: heavy atom identity or coordinates changed")
    grouped = defaultdict(list)
    for line in atom_lines:
        grouped[residue_key(line)].append(line)
    if set(hydrogen_lines) - set(grouped):
        raise ValueError(f"{code}: OpenMM added H to unknown residue")
    serial = max(int(line[6:11]) for line in atom_lines) + 1
    h_output = []
    for rk, atoms in hydrogen_lines.items():
        template = grouped[rk][0]
        for name, xyz in atoms:
            if serial > 99999:
                raise ValueError("PDB five-column atom serial overflow")
            record = template[:6] + f"{serial:5d} {name:4s} " + template[17:30]
            record += "".join(f"{float(v):8.3f}" for v in xyz)
            record += f"{1.00:6.2f}{0.00:6.2f}          {'H':>2s}  "
            h_output.append(record)
            serial += 1
    # Insert H immediately after its parent residue: appending all H after
    # CONECT would cause PDB readers to interpret duplicate residue records.
    h_by_residue = defaultdict(list)
    for line in h_output:
        h_by_residue[residue_key(line)].append(line)
    last_index = {}
    for i, line in enumerate(lines):
        if line.startswith(("ATOM  ", "HETATM")):
            last_index[residue_key(line)] = i
    out = [f"REMARK 901 OPENMM 8.6.1 HYDROGENS PH 7.0; SEED {SEED}; GOLD IMPORT UNTESTED"]
    for i, line in enumerate(lines):
        if line == "END":
            continue
        out.append(line)
        if line.startswith(("ATOM  ", "HETATM")) and last_index[residue_key(line)] == i:
            out.extend(h_by_residue.pop(residue_key(line), []))
    if h_by_residue:
        raise ValueError(f"{code}: some hydrogen residues had no heavy-atom parent")
    out.append("END")
    return "\n".join(out) + "\n", h_counts, protonation


def missing_atom_report(pdb_text: str):
    import io

    fixer = PDBFixer(pdbfile=io.StringIO(pdb_text))
    fixer.findMissingResidues()
    gaps = [(k, v) for k, v in fixer.missingResidues.items()]
    fixer.findMissingAtoms()
    return len(fixer.missingAtoms), len(fixer.missingTerminals), str(gaps)


def main() -> None:
    rows = list(csv.DictReader(SOURCE.open(newline="", encoding="utf-8")))
    if {r["PDB_ID"] for r in rows} != set(VARIANTS):
        raise ValueError("Expected 2QMJ and 1XD0 source records")
    results = []
    for row in rows:
        code = row["PDB_ID"]
        source_path = ROOT / row["working_receptor_path"]
        if not source_path.is_file() or digest(source_path) != row["working_receptor_SHA256"]:
            raise ValueError(f"{code}: working receptor missing or changed")
        ligand = reference_atoms(row)
        source = source_path.read_text(encoding="ascii").splitlines()
        for label, water_rule, keep_nag in VARIANTS[code]:
            filtered, old_serials, expected_heavy, choices, keep_waters = filter_receptor(code, source, ligand, water_rule, keep_nag)
            coordinate = build_coordinate_pdb(source, filtered, old_serials, choices)
            missing, terminal, gaps = missing_atom_report(coordinate)
            if keep_nag:
                # OpenMM does not add NAG hydrogens and mis-protonates the
                # covalently linked Asn amide. Preserve a coordinate-only
                # glycan variant rather than writing invalid chemistry.
                prepared = coordinate
                h_counts = {}
                output = PDB_OUT / f"{code}_gold_stage_{label}_unprotonated.pdb"
            else:
                prepared, h_counts, protonation = protonated_pdb(coordinate, expected_heavy, code)
                output = PDB_OUT / f"{code}_gold_stage_{label}_pH7.pdb"
            output.write_text(prepared, encoding="ascii", newline="\n")
            # Validate OpenMM can parse the final PDB and that retained heavy coordinates are exact.
            check = PDBFile(str(output))
            actual = {atom_key(line): pos(line) for line in output.read_text(encoding="ascii").splitlines()
                      if line.startswith(("ATOM  ", "HETATM")) and line[76:78].strip().upper() not in ("H", "D")}
            if actual != expected_heavy or len(actual) != len(expected_heavy):
                raise ValueError(f"{code}/{label}: heavy coordinates changed after export")
            if not keep_nag and not any(a.element.symbol == "H" for a in check.topology.atoms()):
                raise ValueError(f"{code}/{label}: exported hydrogens not parsable")
            total_waters = sum(line.startswith("HETATM") and line[17:20] == "HOH" for line in source)
            retained_waters = len(keep_waters)
            site_water_count = sum(min(math.dist(pos(line), p) for p in ligand) <= WATER_CUTOFF_A
                                   for line in filtered if line.startswith("HETATM") and line[17:20] == "HOH")
            missing_hetero_h = "NAG and linked Asn protonation unresolved; variant intentionally has no hydrogens" if keep_nag else "none among retained residue classes by OpenMM templates; GOLD atom typing untested"
            if code == "2QMJ":
                acid = "Asp443/Asp542 carboxylates: no side-chain H at pH 7; proton-transfer state unvalidated"
                his = "His600 ND1-H (HD1 present; HE2 absent), OpenMM protein-only estimate"
            else:
                acid = "Asp197/Glu233/Asp300 carboxylates: no side-chain H at pH 7; proton-transfer state unvalidated"
                his = "His201 ND1-H (HD1 present; HE2 absent), OpenMM protein-only estimate"
            if keep_nag:
                acid = "unassigned; coordinate-only variant has no added H"
                his = "unassigned; coordinate-only variant has no added H"
            results.append({
                "PDB_ID": code,
                "variant": label,
                "input_receptor_path": row["working_receptor_path"],
                "input_SHA256": row["working_receptor_SHA256"],
                "native_reference_path": row["native_ligand_SDF"],
                "native_reference_SHA256": row["native_ligand_SHA256"],
                "atom_map_path": row["atom_mapping_CSV"],
                "atom_map_SHA256": row["atom_mapping_SHA256"],
                "output_path": output.relative_to(ROOT).as_posix(),
                "output_SHA256": digest(output),
                "prep_tool": "PDBFixer 1.12.0 missing-atom audit; OpenMM 8.6.1 Modeller.addHydrogens only for NAG-free variants",
                "prep_pH": "not assigned; coordinate-only" if keep_nag else "7.0 assumed, not experimental",
                "altloc_rule": "highest mean occupancy per residue; alphabetical label on tie; chosen A for all six 2QMJ residues",
                "altloc_residues": len(choices),
                "standard_sidechain_missing_heavy_atoms": missing,
                "missing_terminal_atoms": terminal,
                "sequence_gaps": gaps,
                "water_rule": f"all deposited" if water_rule == "all" else f"within {WATER_CUTOFF_A:.1f} A of native ligand plus 1XD0 calcium-coordinating HOH536/546",
                "waters_retained": retained_waters,
                "waters_removed": total_waters - retained_waters,
                "site_waters_retained": site_water_count,
                "calcium_waters_retained": "HOH536 and HOH546" if code == "1XD0" else "not applicable",
                "ions_retained": "CA A498 and CL A500" if code == "1XD0" else "none; SO4 A4001 removed as remote crystallization anion",
                "glycan_decision": "all deposited NAG retained; no H added because linked Asn/NAG protonation is unsupported" if keep_nag else ("all NAG removed for sensitivity; NAG A2005 was 3.895 A from native inhibitor" if code == "2QMJ" else "NAG A499 removed; nearest native ligand distance 47.4 A"),
                "glycerol_decision": "all nine GOL removed; nearest native ligand distance 13.86 A" if code == "2QMJ" else "not deposited",
                "sulfate_decision": "SO4 A4001 removed; nearest native ligand distance 49.37 A" if code == "2QMJ" else "not deposited",
                "PCA_decision": "PCA A1 retained and protonated by OpenMM; peptide linkage retained" if code == "1XD0" else "not deposited",
                "protein_hydrogens_added": sum(v for k, v in h_counts.items() if k not in ("HOH", "NAG", "GOL", "SO4", "CA", "CL")),
                "water_hydrogens_added": h_counts.get("HOH", 0),
                "hetero_hydrogen_limitation": missing_hetero_h,
                "catalytic_acid_state": acid,
                "site_histidine_state": his,
                "heavy_coordinates_preserved": "yes; all retained heavy atom names and coordinates unchanged",
                "reference_mapping_valid": "yes; reference SDF and map checksums and coordinates match source record",
                "GOLD_import_test": "not performed; GOLD/Hermes executable and CCDC API unavailable in this workspace",
                "chemical_preparation_status": "coordinate-only; NAG protonation unresolved" if keep_nag else "provisionally protonated; GOLD atom typing and histidine/acid-state review pending",
                "redocking_ready": "no",
                "issues": "No GOLD import; NAG-linked chemistry unresolved where retained; catalytic protonation provisional; missing construct residues 1-6 and 837, explicit chain break after 836, terminal OXT absent; deposited Cys26-Cys44 SG-SG 2.916 A; pH7 NAG-free model HOH4627 H1 to native ligand heavy atom 1.553 A" if code == "2QMJ" else "No GOLD import; PCA/ion typing and metal coordination in GOLD unverified; catalytic protonation provisional; nearest pH7 generated H to native ARE heavy atom about 1.96 A; ARE crystal-mate contact",
            })
            print(code, label, "waters", retained_waters, "H", sum(h_counts.values()), output.name)
    with OUTPUT.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(results)


if __name__ == "__main__":
    main()

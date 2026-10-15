"""Qualify LIG01 chemical inputs for a future 3L4W pose-only pilot.

Run with .venv/Scripts/python.exe from the repository root. No Vina search or
scoring is invoked. Existing inputs are read-only; new outputs are write-once.
All atom indices in the manifest are one-based SDF atom indices.
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
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolAlign, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import validate_vina_3L4W_pilot_inputs as common  # noqa: E402

SOURCE_SDF = ROOT / "structures/ligands/prepared/LIG01.sdf"
NEUTRAL_PDBQT = ROOT / "structures/ligands/vina_prepared/LIG01_CID370_neutral_reference.pdbqt"
VARIANT_SDF = ROOT / "structures/ligands/vina_prepared/pilot_3L4W/LIG01_CID370_carboxylate_pH6p5_hypothesis.sdf"
VARIANT_PDBQT = VARIANT_SDF.with_suffix(".pdbqt")
OUT = ROOT / "data/metadata/vina_3L4W_gallic_acid_input_manifest.csv"
MIG_RESULTS = ROOT / "data/metadata/vina_3L4W_large_box_mig_results.csv"
BOX_CENTER = common.BOX_CENTER
BOX_SIZE = np.array([30.0, 24.0, 24.0])
SOURCE_SMILES = "C1=C(C=C(C(=C1O)O)O)C(=O)O"
SEEDS = "1101;2202;3303;4404;5505"
CORE = tuple(range(1, 11))  # ring C1-C6, phenolic O7-O9, carboxyl C10
FIELDS = [
    "record_type", "structure_id", "compound_name", "PubChem_CID", "association_id", "plant_species", "plant_part", "literature_extraction_solvent", "evidence_tier", "literature_DOI", "microstate", "decision", "decision_scope", "decision_reason", "formula", "formal_charge", "source_isomeric_SMILES", "verified_isomeric_SMILES", "source_sdf_path", "source_sdf_sha256", "source_record_count", "all_source_records_QC", "selected_record_1based", "source_conformer_number", "selected_mmff94s_energy_kcal_mol", "variant_mmff94s_energy_kcal_mol", "variant_mmff94s_converged", "carboxyl_co_distances_A", "source_energy_range_kcal_mol", "unique_heavy_conformers_0p5A", "selected_heavy_span_xyz_A", "maximum_source_heavy_centroid_radius_A", "minimum_centered_box_clearance_A", "source_pdbqt_path", "source_pdbqt_sha256", "prepared_sdf_path", "prepared_sdf_sha256", "prepared_pdbqt_path", "prepared_pdbqt_sha256", "heavy_atoms", "stereocentres", "pdbqt_charge_sum", "torsions", "pdbqt_graph_stereo_readback", "max_heavy_coordinate_shift_A", "heavy_atom_map_sdf_to_pdbqt", "rigid_core_sdf_indices_1based", "rigid_core_pdbqt_serials_1based", "rigid_core_symmetry_automorphisms", "ionization_site_sdf_index_1based", "ionization_evidence", "tautomer_limit", "receptor_state", "receptor_source_pdb_path", "receptor_source_pdb_sha256", "receptor_pdbqt_path", "receptor_pdbqt_sha256", "receptor_heavy_atoms", "native_mig_sdf_path", "native_mig_sdf_sha256", "native_mig_pdbqt_path", "native_mig_pdbqt_sha256", "large_box_mig_results_path", "large_box_mig_results_sha256", "large_box_mig_control", "vina_executable_sha256", "vina_version", "box_center_A", "box_size_A", "scoring", "exhaustiveness", "num_modes", "energy_range_kcal_mol", "min_rmsd_A", "cpu", "seeds", "pocket_rule", "clash_rule", "cluster_rule", "software_versions", "limitations",
]


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def keyed_rows(path: str, key: str) -> dict[str, dict[str, str]]:
    rows = common.rows(path)
    if len(rows) != len({r[key] for r in rows}):
        raise ValueError(f"Nonunique {key} in {path}")
    return {r[key]: r for r in rows}


def assert_graph(mol: Chem.Mol, smiles: str, formula: str, charge: int) -> None:
    if common.iso(mol) != smiles or rdMolDescriptors.CalcMolFormula(mol) != formula:
        raise ValueError("Graph or formula differs from reference")
    if sum(a.GetFormalCharge() for a in mol.GetAtoms()) != charge:
        raise ValueError("Formal charge differs from reference")
    if Chem.FindMolChiralCenters(Chem.RemoveHs(mol), includeUnassigned=True):
        raise ValueError("Unexpected gallic-acid stereocentre")
    if sum(a.GetAtomicNum() > 1 for a in mol.GetAtoms()) != 12:
        raise ValueError("Expected twelve heavy atoms")
    common.check_3d(mol)
    positions = mol.GetConformer().GetPositions()
    graph_distances = Chem.GetDistanceMatrix(mol)
    heavy_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    for n, i in enumerate(heavy_indices):
        for j in heavy_indices[n + 1:]:
            if graph_distances[i, j] >= 4 and np.linalg.norm(positions[i] - positions[j]) < 1.6:
                raise ValueError("Severe internal nonbonded heavy-atom clash")


def confirm_controls() -> tuple[str, dict[str, list[tuple[str, np.ndarray]]]]:
    for item in [*common.RECEPTORS.values(), *common.RECEPTOR_PDB.values(), common.NATIVE, common.NATIVE_PDBQT]:
        common.check_hash(ROOT / item[0], item[1])
    common.check_hash(*common.VINA)
    version = subprocess.run([str(common.VINA[0]), "--version"], capture_output=True, text=True, timeout=15)
    if version.returncode != 0 or "1.2.7" not in version.stdout + version.stderr:
        raise ValueError("Frozen Vina executable does not identify as version 1.2.7")
    frozen = common.rows("data/metadata/vina_3L4W_pilot_input_manifest.csv")
    for state, (receptor, digest) in common.RECEPTORS.items():
        matching = [r for r in frozen if r["receptor_state"] == state]
        if not matching or any(r["receptor_pdbqt_path"] != receptor or r["receptor_pdbqt_sha256"] != digest for r in matching):
            raise ValueError(f"Step 57 receptor freeze mismatch: {state}")
    mig = next(iter(Chem.SDMolSupplier(str(ROOT / common.NATIVE[0]), removeHs=False)))
    if mig is None:
        raise ValueError("Native MIG reference SDF cannot be parsed")
    mig_heavy = np.asarray([mig.GetConformer().GetAtomPosition(a.GetIdx()) for a in mig.GetAtoms() if a.GetAtomicNum() > 1])
    if len(mig_heavy) != 14 or not np.allclose(mig_heavy.mean(axis=0), BOX_CENTER, atol=1e-7):
        raise ValueError("Frozen center is not the deposited MIG heavy-atom centroid")
    native_rows = common.rows(rel(MIG_RESULTS))
    if len(native_rows) != 10 or {(r["receptor_state"], r["seed"]) for r in native_rows} != {(s, v) for s in common.RECEPTORS for v in SEEDS.split(";")}:
        raise ValueError("Incomplete enlarged-box native-MIG controls")
    for row in native_rows:
        if row["status"] != "completed" or row["exit_code"] != "0" or float(row["rank1_RMSD_A"]) > 2.0 or row["box_size_xyz_A"] != "30;24;24":
            raise ValueError("Enlarged-box MIG positive control failed")
        for pcol, hcol in [("receptor_path", "receptor_SHA256"), ("ligand_path", "ligand_SHA256"), ("reference_path", "reference_SHA256"), ("output_path", "output_SHA256"), ("log_path", "log_SHA256")]:
            common.check_hash(ROOT / row[pcol], row[hcol])
    atoms = {s: common.receptor_heavy_atoms(ROOT / p) for s, (p, _) in common.RECEPTORS.items()}
    a, b = atoms.values()
    if len(a) != 6912 or len(b) != len(a) or any(x[0] != y[0] or np.linalg.norm(x[1] - y[1]) > 0.0011 for x, y in zip(a, b)):
        raise ValueError("Receptor heavy atoms differ across frozen microstates")
    return common.sha(MIG_RESULTS), atoms


def main() -> None:
    control_hash, receptor_atoms = confirm_controls()
    prep = keyed_rows("data/metadata/ligand_preparation_results.csv", "structure_id")["LIG01"]
    vina_prep = keyed_rows("data/metadata/vina_ligand_preparation_results.csv", "structure_id")["LIG01"]
    library = [r for r in common.rows("data/compounds/final_preliminary_ligand_library.csv") if r["structure_id"] == "LIG01"]
    if len(library) != 1:
        raise ValueError("LIG01 association count changed")
    assoc = library[0]
    if (assoc["association_id"], assoc["PubChem_CID"], assoc["plant_species"], assoc["plant_part_reported"], assoc["extraction_solvent_reported"], assoc["evidence_tier"]) != ("AM01", "370", "Aegle marmelos", "leaves", "ethanol", "B"):
        raise ValueError("Botanical evidence fields changed")
    if assoc["primary_literature_DOI"] != "https://doi.org/10.3390/plants10122573":
        raise ValueError("Botanical source DOI changed")
    canonical = Chem.MolToSmiles(Chem.MolFromSmiles(SOURCE_SMILES), isomericSmiles=True)
    if any(Chem.MolToSmiles(Chem.MolFromSmiles(s), isomericSmiles=True) != canonical for s in [prep["source_isomeric_SMILES"], assoc["isomeric_SMILES"], assoc["canonical_SMILES"]]):
        raise ValueError("Source and library SMILES disagree with PubChem CID 370")
    common.check_hash(SOURCE_SDF, prep["sdf_sha256"])
    common.check_hash(NEUTRAL_PDBQT, vina_prep["pdbqt_sha256"])
    if prep["sdf_path"] != rel(SOURCE_SDF) or vina_prep["pdbqt_path"] != rel(NEUTRAL_PDBQT):
        raise ValueError("Metadata source paths disagree")
    mols = list(Chem.SDMolSupplier(str(SOURCE_SDF), removeHs=False, sanitize=True))
    if len(mols) != 24 or len(mols) != int(prep["conformers_exported"]):
        raise ValueError("Incorrect LIG01 conformer count")
    energies: list[float] = []
    radii: list[float] = []
    for m in mols:
        if m is None or m.GetProp("structure_id") != "LIG01" or m.GetProp("PubChem_CID") != "370":
            raise ValueError("Conformer identity mismatch")
        if Chem.MolToSmiles(Chem.MolFromSmiles(m.GetProp("source_isomeric_SMILES")), isomericSmiles=True) != canonical:
            raise ValueError("Conformer source SMILES mismatch")
        assert_graph(m, canonical, "C7H6O5", 0)
        t = Chem.Mol(m)
        Chem.AssignStereochemistryFrom3D(t)
        if common.iso(t) != canonical:
            raise ValueError("Unexpected 3D stereochemistry")
        energy = float(m.GetProp("mmff_energy_kcal_mol"))
        if not math.isfinite(energy):
            raise ValueError("Nonfinite MMFF energy")
        energies.append(energy)
        xyz = np.asarray([m.GetConformer().GetAtomPosition(a.GetIdx()) for a in m.GetAtoms() if a.GetAtomicNum() > 1])
        radii.append(float(np.max(np.linalg.norm(xyz - xyz.mean(axis=0), axis=1))))
    idx = min(range(len(mols)), key=lambda i: (energies[i], i))
    if idx + 1 != int(vina_prep["selected_sdf_record_index_1based"]) or not math.isclose(energies[idx], float(vina_prep["selected_mmff94s_energy_kcal_mol"]), abs_tol=1e-7):
        raise ValueError("Existing PDBQT conformer selection mismatch")
    selected = mols[idx]
    assert selected is not None
    if selected.GetAtomWithIdx(11).GetSymbol() != "O" or selected.GetBondBetweenAtoms(9, 11) is None or selected.GetBondBetweenAtoms(9, 10) is None:
        raise ValueError("Carboxyl-oxygen atom assignment is not the expected source graph")
    neutral_qc = common.check_pdbqt(NEUTRAL_PDBQT, selected)
    if neutral_qc["torsions"] != 5:
        raise ValueError("Original neutral PDBQT torsion count changed")
    # Original atom 12 is the carboxylic OH oxygen; the phenolic oxygens are 7-9.
    anion = common.deprotonate(selected, 11)
    anion_smiles = common.iso(anion)
    assert_graph(anion, anion_smiles, "C7H5O5-", -1)
    for i in range(12):
        if selected.GetAtomWithIdx(i).GetSymbol() != anion.GetAtomWithIdx(i).GetSymbol():
            raise ValueError("Deprotonation changed heavy-atom index")
    if any(a.GetFormalCharge() for a in anion.GetAtoms() if a.GetIdx() != 11):
        raise ValueError("Unexpected additional ionization")
    if not AllChem.MMFFHasAllMoleculeParams(anion):
        raise ValueError("MMFF94s lacks parameters for gallic-acid carboxylate")
    mmff_props = AllChem.MMFFGetMoleculeProperties(anion, mmffVariant="MMFF94s")
    force_field = AllChem.MMFFGetMoleculeForceField(anion, mmff_props)
    if force_field.Minimize(maxIts=1000) != 0:
        raise ValueError("MMFF94s carboxylate minimization did not converge")
    anion_energy = force_field.CalcEnergy()
    if not math.isfinite(anion_energy):
        raise ValueError("Nonfinite carboxylate MMFF94s energy")
    assert_graph(anion, anion_smiles, "C7H5O5-", -1)
    carboxyl_xyz = anion.GetConformer().GetPositions()
    co_distances = [float(np.linalg.norm(carboxyl_xyz[9] - carboxyl_xyz[i])) for i in (10, 11)]
    if any(not 1.15 < d < 1.40 for d in co_distances) or abs(co_distances[0] - co_distances[1]) > 0.05:
        raise ValueError("Carboxylate bond geometry is inconsistent after minimization")
    anion.SetProp("microstate", "carboxylate_pH6p5_hypothesis")
    anion.SetProp("variant_mmff94s_energy_kcal_mol", f"{anion_energy:.8f}")
    anion.SetProp("variant_mmff94s_converged", "yes")
    anion.SetProp("parent_source_sdf_sha256", common.sha(SOURCE_SDF))
    anion.SetProp("selected_parent_record_1based", str(idx + 1))
    anion.SetProp("deprotonated_parent_atom_index_1based", "12")
    block = Chem.MolToMolBlock(anion)
    props = "".join(f">  <{k}>\n{anion.GetProp(k)}\n\n" for k in sorted(anion.GetPropNames()))
    common.write_once(VARIANT_SDF, block + "\n" + props + "$$$$\n")
    check_variant = list(Chem.SDMolSupplier(str(VARIANT_SDF), removeHs=False, sanitize=True))
    if len(check_variant) != 1 or check_variant[0] is None:
        raise ValueError("Carboxylate SDF roundtrip failed")
    assert_graph(check_variant[0], anion_smiles, "C7H5O5-", -1)
    setup = MoleculePreparation(charge_model="gasteiger", merge_these_atom_types=("H",), add_index_map=True, hydrate=False, rigid_macrocycles=False, flexible_amides=False).prepare(anion)
    if len(setup) != 1:
        raise ValueError("Meeko did not produce one anion setup")
    qtext, ok, err = PDBQTWriterLegacy.write_string(setup[0], add_index_map=True)
    if not ok:
        raise ValueError(f"Meeko carboxylate conversion failed: {err}")
    common.write_once(VARIANT_PDBQT, qtext)
    anion_qc = common.check_pdbqt(VARIANT_PDBQT, anion)
    if anion_qc["torsions"] not in (4, 5):
        raise ValueError("Unexpected carboxylate torsion count")
    representatives: list[Chem.Mol] = []
    for m in mols:
        if not representatives or all(rdMolAlign.GetBestRMS(Chem.RemoveHs(m), Chem.RemoveHs(ref)) >= 0.5 for ref in representatives):
            representatives.append(m)
    xyz = np.asarray([selected.GetConformer().GetAtomPosition(i) for i in range(12)])
    span = np.ptp(xyz, axis=0)
    anion_xyz = np.asarray([anion.GetConformer().GetAtomPosition(i) for i in range(12)])
    anion_radius = float(np.max(np.linalg.norm(anion_xyz - anion_xyz.mean(axis=0), axis=1)))
    min_clearance = float(np.min(BOX_SIZE / 2 - max(max(radii), anion_radius)))
    if min_clearance < 4.0:
        raise ValueError("Insufficient centered-box ligand clearance")
    heavy = Chem.RemoveHs(selected)
    automorphisms = len(heavy.GetSubstructMatches(heavy, uniquify=False))
    output: list[dict[str, object]] = []
    for state, molecule, pfile, sfile, qc, scope in [
        ("neutral_acid_reference", selected, NEUTRAL_PDBQT, SOURCE_SDF, neutral_qc, "PASS: minority aqueous neutral-acid sensitivity hypothesis"),
        ("carboxylate_pH6p5_hypothesis", anion, VARIANT_PDBQT, VARIANT_SDF, anion_qc, "PASS: primary bulk-water carboxylate hypothesis"),
    ]:
        atom_map = qc["map"]
        for receptor_state, (receptor_path, receptor_sha) in common.RECEPTORS.items():
            output.append({
                "record_type": "gallic_acid_input_pair", "structure_id": "LIG01", "compound_name": "Gallic acid", "PubChem_CID": "370", "association_id": "AM01", "plant_species": "Aegle marmelos", "plant_part": "leaves", "literature_extraction_solvent": "ethanol", "evidence_tier": "B", "literature_DOI": assoc["primary_literature_DOI"], "microstate": state, "decision": "PASS", "decision_scope": scope, "decision_reason": "CID-backed graph; complete 3D/SDF and Meeko PDBQT readback; exact source hashes; frozen receptor and native-MIG control verified", "formula": rdMolDescriptors.CalcMolFormula(molecule), "formal_charge": sum(a.GetFormalCharge() for a in molecule.GetAtoms()), "source_isomeric_SMILES": SOURCE_SMILES, "verified_isomeric_SMILES": common.iso(molecule), "source_sdf_path": rel(SOURCE_SDF), "source_sdf_sha256": common.sha(SOURCE_SDF), "source_record_count": len(mols), "all_source_records_QC": "24/24", "selected_record_1based": idx + 1, "source_conformer_number": selected.GetProp("conformer_number"), "selected_mmff94s_energy_kcal_mol": f"{energies[idx]:.8f}", "variant_mmff94s_energy_kcal_mol": f"{anion_energy:.8f}" if state.startswith("carboxylate") else "", "variant_mmff94s_converged": "yes" if state.startswith("carboxylate") else "not_applicable", "carboxyl_co_distances_A": ";".join(f"{d:.4f}" for d in co_distances) if state.startswith("carboxylate") else "", "source_energy_range_kcal_mol": f"{min(energies):.8f}..{max(energies):.8f}", "unique_heavy_conformers_0p5A": len(representatives), "selected_heavy_span_xyz_A": ";".join(f"{x:.4f}" for x in span), "maximum_source_heavy_centroid_radius_A": f"{max(radii):.4f}", "minimum_centered_box_clearance_A": f"{min_clearance:.4f}", "source_pdbqt_path": rel(NEUTRAL_PDBQT), "source_pdbqt_sha256": common.sha(NEUTRAL_PDBQT), "prepared_sdf_path": rel(sfile), "prepared_sdf_sha256": common.sha(sfile), "prepared_pdbqt_path": rel(pfile), "prepared_pdbqt_sha256": common.sha(pfile), "heavy_atoms": 12, "stereocentres": 0, "pdbqt_charge_sum": f"{qc['charge']:.3f}", "torsions": qc["torsions"], "pdbqt_graph_stereo_readback": "pass", "max_heavy_coordinate_shift_A": f"{qc['shift']:.4f}", "heavy_atom_map_sdf_to_pdbqt": ";".join(f"{i}:{atom_map[i]}" for i in range(1, 13)), "rigid_core_sdf_indices_1based": ";".join(map(str, CORE)), "rigid_core_pdbqt_serials_1based": ";".join(str(atom_map[i]) for i in CORE), "rigid_core_symmetry_automorphisms": automorphisms, "ionization_site_sdf_index_1based": "12" if state.startswith("carboxylate") else "", "ionization_evidence": "aqueous carboxyl pKa about 4.4 (DOI 10.1080/00032718508068758); pH 6.5 favors carboxylate in bulk water; bound state unassigned", "tautomer_limit": "Phenolic deprotonation and bound-state proton shifts not enumerated; no tautomer selected by score", "receptor_state": receptor_state, "receptor_source_pdb_path": common.RECEPTOR_PDB[receptor_state][0], "receptor_source_pdb_sha256": common.RECEPTOR_PDB[receptor_state][1], "receptor_pdbqt_path": receptor_path, "receptor_pdbqt_sha256": receptor_sha, "receptor_heavy_atoms": len(receptor_atoms[receptor_state]), "native_mig_sdf_path": common.NATIVE[0], "native_mig_sdf_sha256": common.NATIVE[1], "native_mig_pdbqt_path": common.NATIVE_PDBQT[0], "native_mig_pdbqt_sha256": common.NATIVE_PDBQT[1], "large_box_mig_results_path": rel(MIG_RESULTS), "large_box_mig_results_sha256": control_hash, "large_box_mig_control": "5/5 rank1 RMSD <=2A in this receptor; 10/10 combined; no rerun", "vina_executable_sha256": common.VINA[1], "vina_version": "1.2.7", "box_center_A": ";".join(f"{x:.8f}" for x in BOX_CENTER), "box_size_A": "30;24;24", "scoring": "vina", "exhaustiveness": 64, "num_modes": 20, "energy_range_kcal_mol": 5.0, "min_rmsd_A": 1.0, "cpu": 1, "seeds": SEEDS, "pocket_rule": "all heavy atoms inside box; core centroid <=6A from MIG centroid; >=50% core atoms <=4A from MIG heavy atom", "clash_rule": "receptor-ligand heavy pair <1.8A OR VdW overlap >0.8A; nonbonded internal heavy pair >=4 bonds apart <1.6A", "cluster_rule": "unfitted symmetry-aware core RMSD <=2.0A to rank1 medoid; >=4/5 rank1 poses pass all rules per receptor state", "software_versions": f"RDKit {rdBase.rdkitVersion}; Meeko {meeko.__version__}; Vina 1.2.7", "limitations": "Tier B ethanol leaf evidence; not measured in project aqueous extract; dry/no-NAG receptor; His600-HID; missing binding-site waters; bound ligand microstate unknown; isolated MMFF minimum not proven bioactive; no docking performed",
            })
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(output)
    common.write_once(OUT, stream.getvalue())
    print(f"LIG01 source records: {len(mols)}; states: neutral PASS, carboxylate PASS; pairs: {len(output)}; no docking")
    for p in [SOURCE_SDF, NEUTRAL_PDBQT, VARIANT_SDF, VARIANT_PDBQT, MIG_RESULTS, OUT]:
        print(f"{rel(p)} {common.sha(p)}")


if __name__ == "__main__":
    main()



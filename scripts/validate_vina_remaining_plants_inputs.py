"""Qualify LIG03 and LIG09 reference inputs for future 3L4W pose-only planning.

Run from this repository with .venv/Scripts/python.exe. No Vina scoring or
docking is called. Source files are read-only; derived files are write-once.
Use --online-pubchem to repeat the live PubChem CID property check.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import sys
import urllib.request
from pathlib import Path

import meeko
import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolAlign, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import validate_vina_3L4W_pilot_inputs as common  # noqa: E402
import validate_vina_3L4W_gallic_acid_inputs as gallic  # noqa: E402

OUT = ROOT / "data/metadata/vina_remaining_plants_input_manifest.csv"
VARIANT_SDF = ROOT / "structures/ligands/vina_prepared/remaining_plants_3L4W/LIG03_CID689043_E_carboxylate_pH6p5_hypothesis.sdf"
VARIANT_PDBQT = VARIANT_SDF.with_suffix(".pdbqt")
PUBCHEM = {
    "LIG03": {"cid": 689043, "formula": "C9H8O4", "charge": 0,
              "smiles": "C1=CC(=C(C=C1/C=C/C(=O)O)O)O", "heavy": 13, "records": 24,
              "association": "BD01", "species": "Boerhavia diffusa", "core": tuple(range(1, 10))},
    "LIG09": {"cid": 442416, "formula": "C22H26O11", "charge": 0,
              "smiles": "C1=CO[C@H]([C@H]2[C@@H]1[C@@H](C=C2COC(=O)C3=CC=C(C=C3)O)O)O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O",
              "heavy": 33, "records": 60, "association": "VN03", "species": "Vitex negundo",
              "core": tuple(range(1, 10))},
}
FIELDS = [
    "structure_id", "compound_name", "association_id", "plant_species", "PubChem_CID", "PubChem_URL",
    "botanical_evidence_tier", "botanical_evidence_scope", "primary_DOI", "corroborating_DOI",
    "candidate_decision", "candidate_decision_reason", "microstate", "microstate_input_decision",
    "docking_execution_status", "formula", "formal_charge", "reference_isomeric_SMILES",
    "prepared_isomeric_SMILES", "source_sdf_path", "source_sdf_SHA256", "source_records",
    "all_source_records_QC", "selected_source_record_1based", "selected_source_conformer_number",
    "selected_source_MMFF94s_energy_kcal_mol", "source_energy_range_kcal_mol",
    "distinct_heavy_conformers_at_0p5A", "source_pdbqt_path", "source_pdbqt_SHA256",
    "prepared_sdf_path", "prepared_sdf_SHA256", "prepared_pdbqt_path", "prepared_pdbqt_SHA256",
    "variant_MMFF94s_energy_kcal_mol", "variant_MMFF94s_converged", "ionized_source_atom_1based",
    "heavy_atoms", "defined_stereocentres", "double_bond_stereo", "pdbqt_charge_sum",
    "pdbqt_torsions", "pdbqt_readback", "max_heavy_coordinate_shift_A",
    "heavy_atom_map_sdf_to_pdbqt", "rigid_core_sdf_indices_1based",
    "rigid_core_pdbqt_serials_1based", "rigid_core_automorphisms",
    "maximum_source_heavy_centroid_radius_A", "centered_box_clearance_A",
    "receptor_state", "receptor_pdbqt_path", "receptor_pdbqt_SHA256",
    "receptor_source_pdb_path", "receptor_source_pdb_SHA256", "MIG_reference_sdf_path",
    "MIG_reference_sdf_SHA256", "MIG_control_csv_path", "MIG_control_csv_SHA256",
    "MIG_control_status", "Vina_executable_SHA256", "Vina_version", "box_center_A",
    "box_size_A", "future_settings", "future_seeds", "pocket_rule", "clash_rule",
    "proposed_cluster_rule", "software_versions", "limitations",
]


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def online_pubchem_check() -> None:
    for sid, reference in PUBCHEM.items():
        cid = reference["cid"]
        url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/"
               "MolecularFormula,CanonicalSMILES,IsomericSMILES,Charge/JSON")
        with urllib.request.urlopen(url, timeout=30) as response:
            data = json.load(response)["PropertyTable"]["Properties"]
        if len(data) != 1:
            raise ValueError(f"PubChem returned nonunique record for {sid}")
        item = data[0]
        smiles = item.get("SMILES", item.get("IsomericSMILES"))
        if (item["CID"], item["MolecularFormula"], item["Charge"], smiles) != (
            cid, reference["formula"], reference["charge"], reference["smiles"]
        ):
            raise ValueError(f"Official PubChem property mismatch for {sid}: {item}")


def verify_molecule(mol: Chem.Mol, sid: str, smiles: str, formula: str, charge: int,
                    heavy: int) -> None:
    if common.iso(mol) != smiles or rdMolDescriptors.CalcMolFormula(mol) != formula:
        raise ValueError(f"{sid}: molecular graph, stereo or formula mismatch")
    if Chem.GetFormalCharge(mol) != charge or sum(a.GetAtomicNum() > 1 for a in mol.GetAtoms()) != heavy:
        raise ValueError(f"{sid}: charge/heavy-atom mismatch")
    common.check_3d(mol)
    inferred = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(inferred)
    if common.iso(inferred) != smiles:
        raise ValueError(f"{sid}: coordinate-implied stereochemistry mismatch")
    graph = Chem.GetDistanceMatrix(mol)
    xyz = mol.GetConformer().GetPositions()
    indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    for n, i in enumerate(indices):
        for j in indices[n + 1:]:
            if graph[i, j] >= 4 and np.linalg.norm(xyz[i] - xyz[j]) < 1.6:
                raise ValueError(f"{sid}: severe intramolecular heavy-atom contact")


def conformer_audit(sid: str, prep: dict, vp: dict, library: dict, reference: dict) -> dict:
    sdf = ROOT / prep["sdf_path"]
    pdbqt = ROOT / vp["pdbqt_path"]
    common.check_hash(sdf, prep["sdf_sha256"])
    common.check_hash(pdbqt, vp["pdbqt_sha256"])
    if vp["source_sdf_path"] != prep["sdf_path"] or vp["source_sdf_sha256"] != prep["sdf_sha256"]:
        raise ValueError(f"{sid}: source metadata disagreement")
    if (int(library["PubChem_CID"]), library["association_id"], library["plant_species"]) != (
        reference["cid"], reference["association"], reference["species"]
    ):
        raise ValueError(f"{sid}: library botanical association/CID mismatch")
    expected = common.iso(Chem.MolFromSmiles(reference["smiles"]))
    if common.iso(Chem.MolFromSmiles(library["isomeric_SMILES"])) != expected or common.iso(Chem.MolFromSmiles(prep["source_isomeric_SMILES"])) != expected:
        raise ValueError(f"{sid}: reference/source isomeric SMILES mismatch")
    mols = list(Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=True))
    if len(mols) != reference["records"] or len(mols) != int(prep["conformers_exported"]) or len(mols) != int(vp["sdf_record_count"]):
        raise ValueError(f"{sid}: SDF conformer count mismatch")
    energies, radii, distinct = [], [], []
    for i, m in enumerate(mols):
        if m is None or m.GetProp("structure_id") != sid or int(m.GetProp("PubChem_CID")) != reference["cid"]:
            raise ValueError(f"{sid}: record {i+1} parse/identity failure")
        if common.iso(Chem.MolFromSmiles(m.GetProp("source_isomeric_SMILES"))) != expected:
            raise ValueError(f"{sid}: record {i+1} source SMILES mismatch")
        verify_molecule(m, sid, expected, reference["formula"], reference["charge"], reference["heavy"])
        energy = float(m.GetProp("mmff_energy_kcal_mol"))
        if not math.isfinite(energy):
            raise ValueError(f"{sid}: nonfinite stored MMFF energy")
        energies.append(energy)
        xyz = np.array([m.GetConformer().GetAtomPosition(a.GetIdx()) for a in m.GetAtoms() if a.GetAtomicNum() > 1])
        radii.append(float(np.max(np.linalg.norm(xyz - xyz.mean(axis=0), axis=1))))
        if all(rdMolAlign.GetBestRMS(Chem.RemoveHs(m), Chem.RemoveHs(x)) >= 0.5 for x in distinct):
            distinct.append(m)
    selected_idx = min(range(len(mols)), key=lambda i: (energies[i], i))
    selected = mols[selected_idx]
    if selected_idx + 1 != int(vp["selected_sdf_record_index_1based"]) or selected.GetProp("conformer_number") != vp["source_conformer_number"] or not math.isclose(energies[selected_idx], float(vp["selected_mmff94s_energy_kcal_mol"]), abs_tol=1e-7):
        raise ValueError(f"{sid}: selected conformer/energy provenance mismatch")
    qc = common.check_pdbqt(pdbqt, selected)
    if qc["heavy"] != reference["heavy"] or qc["torsions"] != int(vp["torsdof"]):
        raise ValueError(f"{sid}: PDBQT heavy atoms/torsions mismatch")
    if max(radii) > 10.0:
        raise ValueError(f"{sid}: at least one source conformer has less than 2 A centered box clearance")
    return {"selected": selected, "sdf": sdf, "pdbqt": pdbqt, "qc": qc,
            "energies": energies, "radii": radii, "distinct": len(distinct),
            "selected_index": selected_idx + 1, "source_smiles": expected}


def prepare_caffeate(source: Chem.Mol) -> tuple[Chem.Mol, dict, float]:
    # Source atom 11 (zero-based 10) is the carboxylic OH oxygen; 12/13 are phenolic.
    o = source.GetAtomWithIdx(10)
    if o.GetSymbol() != "O" or source.GetBondBetweenAtoms(8, 10) is None:
        raise ValueError("Caffeic acid source carboxyl assignment changed")
    anion = common.deprotonate(source, 10)
    expected = common.iso(anion)
    verify_molecule(anion, "LIG03-carboxylate", expected, "C9H7O4-", -1, 13)
    if not AllChem.MMFFHasAllMoleculeParams(anion):
        raise ValueError("MMFF94s lacks trans-caffeate parameters")
    props = AllChem.MMFFGetMoleculeProperties(anion, mmffVariant="MMFF94s")
    force = AllChem.MMFFGetMoleculeForceField(anion, props)
    if force is None or force.Minimize(maxIts=1000) != 0:
        raise ValueError("MMFF94s trans-caffeate minimization failed")
    energy = float(force.CalcEnergy())
    if not math.isfinite(energy):
        raise ValueError("MMFF94s trans-caffeate energy invalid")
    verify_molecule(anion, "LIG03-carboxylate", expected, "C9H7O4-", -1, 13)
    xyz = anion.GetConformer().GetPositions()
    co = [float(np.linalg.norm(xyz[8] - xyz[i])) for i in (9, 10)]
    if any(not 1.15 < d < 1.40 for d in co) or abs(co[0] - co[1]) > 0.05:
        raise ValueError("Trans-caffeate carboxylate C-O geometry invalid")
    anion.SetProp("microstate", "E_carboxylate_pH6p5_hypothesis")
    anion.SetProp("parent_source_sdf_sha256", common.sha(ROOT / "structures/ligands/prepared/LIG03.sdf"))
    anion.SetProp("selected_parent_record_1based", "1")
    anion.SetProp("deprotonated_parent_atom_index_1based", "11")
    anion.SetProp("variant_mmff94s_energy_kcal_mol", f"{energy:.8f}")
    block = Chem.MolToMolBlock(anion)
    text = block + "\n" + "".join(f">  <{k}>\n{anion.GetProp(k)}\n\n" for k in sorted(anion.GetPropNames())) + "$$$$\n"
    common.write_once(VARIANT_SDF, text)
    readback = list(Chem.SDMolSupplier(str(VARIANT_SDF), removeHs=False, sanitize=True))
    if len(readback) != 1 or readback[0] is None:
        raise ValueError("Carboxylate SDF roundtrip failed")
    verify_molecule(readback[0], "LIG03-carboxylate", expected, "C9H7O4-", -1, 13)
    setup = MoleculePreparation(charge_model="gasteiger", merge_these_atom_types=("H",),
                                add_index_map=True, hydrate=False, rigid_macrocycles=False,
                                flexible_amides=False).prepare(anion)
    if len(setup) != 1:
        raise ValueError("Meeko did not produce exactly one trans-caffeate setup")
    qtext, ok, error = PDBQTWriterLegacy.write_string(setup[0], add_index_map=True)
    if not ok:
        raise ValueError(f"Meeko trans-caffeate conversion failed: {error}")
    common.write_once(VARIANT_PDBQT, qtext)
    qc = common.check_pdbqt(VARIANT_PDBQT, anion)
    return anion, qc, energy


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--online-pubchem", action="store_true")
    args = parser.parse_args()
    if args.online_pubchem:
        online_pubchem_check()
    control_sha, receptors = gallic.confirm_controls()
    prep = {r["structure_id"]: r for r in common.rows("data/metadata/ligand_preparation_results.csv")}
    vp = {r["structure_id"]: r for r in common.rows("data/metadata/vina_ligand_preparation_results.csv")}
    library = {r["structure_id"]: r for r in common.rows("data/compounds/final_preliminary_ligand_library.csv") if r["structure_id"] in PUBCHEM}
    audited = {sid: conformer_audit(sid, prep[sid], vp[sid], library[sid], ref) for sid, ref in PUBCHEM.items()}
    caffeate, caffeate_qc, caffeate_energy = prepare_caffeate(audited["LIG03"]["selected"])
    states = {
        "LIG03": [("E_neutral_reference", audited["LIG03"]["selected"], audited["LIG03"]["sdf"], audited["LIG03"]["pdbqt"], audited["LIG03"]["qc"], ""),
                  ("E_carboxylate_pH6p5_hypothesis", caffeate, VARIANT_SDF, VARIANT_PDBQT, caffeate_qc, f"{caffeate_energy:.8f}")],
        "LIG09": [("neutral_agnuside_reference", audited["LIG09"]["selected"], audited["LIG09"]["sdf"], audited["LIG09"]["pdbqt"], audited["LIG09"]["qc"], "")],
    }
    rows = []
    for sid, variants in states.items():
        ref, source, meta = PUBCHEM[sid], audited[sid], library[sid]
        for state, mol, sdf_path, q_path, qc, variant_energy in variants:
            verify_molecule(mol, sid, common.iso(mol), rdMolDescriptors.CalcMolFormula(mol), Chem.GetFormalCharge(mol), ref["heavy"])
            xyz = np.array([mol.GetConformer().GetAtomPosition(a.GetIdx()) for a in mol.GetAtoms() if a.GetAtomicNum() > 1])
            radius = float(np.max(np.linalg.norm(xyz - xyz.mean(axis=0), axis=1)))
            centered_clearance = float(min(15.0, 12.0, 12.0) - max(max(source["radii"]), radius))
            if centered_clearance < 2.0:
                raise ValueError(f"{sid}: insufficient centered-box clearance")
            core = ref["core"]
            if len(core) != 9 or any(mol.GetAtomWithIdx(i-1).GetAtomicNum() == 1 for i in core):
                raise ValueError(f"{sid}: rigid-core atom mapping invalid")
            graph = Chem.RemoveHs(mol)
            auto = len(graph.GetSubstructMatches(graph, uniquify=False))
            maptext = ";".join(f"{a.GetIdx()+1}:{qc['map'][a.GetIdx()+1]}" for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
            botanical = ("C: species water infusion of mixed aerial parts; methanolic leaf tri-TMS derivative only"
                         if sid == "LIG03" else
                         "B/abstract-only aqueous: leaf HPLC-PDA after methanol-first partition; water-only study abstract accessible")
            decision = "HOLD" if sid == "LIG03" else "PASS_FOR_POSE_ONLY_PLANNING"
            decision_reason = ("Leaf-only parent identity and source E/Z geometry unverified" if sid == "LIG03" else
                               "Leaf agnuside association corroborated by primary HPLC; water-only chromatogram and source stereochemistry unresolved")
            limitations = ("PubChem E/trans reference does not assign source E/Z; no isolated aqueous-leaf parent confirmation; bound microstate unknown; dry/no-NAG receptor"
                           if sid == "LIG03" else
                           "Water-only leaf full paper inaccessible; PubChem 9-center stereochemistry not source-resolved; 13 torsions and glycoside sampling; dry/no-NAG receptor")
            for receptor_state, (receptor_path, receptor_sha) in common.RECEPTORS.items():
                rows.append({
                    "structure_id": sid, "compound_name": meta["compound_name"], "association_id": ref["association"],
                    "plant_species": ref["species"], "PubChem_CID": ref["cid"],
                    "PubChem_URL": f"https://pubchem.ncbi.nlm.nih.gov/compound/{ref['cid']}",
                    "botanical_evidence_tier": botanical, "botanical_evidence_scope": meta["source_evidence_classification"],
                    "primary_DOI": meta["primary_literature_DOI"], "corroborating_DOI": meta["corroborating_DOI"],
                    "candidate_decision": decision, "candidate_decision_reason": decision_reason,
                    "microstate": state, "microstate_input_decision": "PASS_REFERENCE_CHEMISTRY",
                    "docking_execution_status": "NOT_AUTHORIZED_BY_THIS_INPUT_AUDIT",
                    "formula": rdMolDescriptors.CalcMolFormula(mol), "formal_charge": Chem.GetFormalCharge(mol),
                    "reference_isomeric_SMILES": ref["smiles"], "prepared_isomeric_SMILES": common.iso(mol),
                    "source_sdf_path": rel(source["sdf"]), "source_sdf_SHA256": common.sha(source["sdf"]),
                    "source_records": ref["records"], "all_source_records_QC": f"{ref['records']}/{ref['records']}",
                    "selected_source_record_1based": source["selected_index"],
                    "selected_source_conformer_number": source["selected"].GetProp("conformer_number"),
                    "selected_source_MMFF94s_energy_kcal_mol": f"{source['energies'][source['selected_index']-1]:.8f}",
                    "source_energy_range_kcal_mol": f"{min(source['energies']):.8f}..{max(source['energies']):.8f}",
                    "distinct_heavy_conformers_at_0p5A": source["distinct"],
                    "source_pdbqt_path": rel(source["pdbqt"]), "source_pdbqt_SHA256": common.sha(source["pdbqt"]),
                    "prepared_sdf_path": rel(sdf_path), "prepared_sdf_SHA256": common.sha(sdf_path),
                    "prepared_pdbqt_path": rel(q_path), "prepared_pdbqt_SHA256": common.sha(q_path),
                    "variant_MMFF94s_energy_kcal_mol": variant_energy,
                    "variant_MMFF94s_converged": "yes" if variant_energy else "not_applicable",
                    "ionized_source_atom_1based": "11" if variant_energy else "",
                    "heavy_atoms": ref["heavy"],
                    "defined_stereocentres": len(Chem.FindMolChiralCenters(Chem.RemoveHs(mol), includeUnassigned=False)),
                    "double_bond_stereo": "E/trans retained" if sid == "LIG03" else "no specified alkene E/Z in reference",
                    "pdbqt_charge_sum": f"{qc['charge']:.3f}", "pdbqt_torsions": qc["torsions"],
                    "pdbqt_readback": "Meeko/RDKit graph and stereo PASS",
                    "max_heavy_coordinate_shift_A": f"{qc['shift']:.4f}",
                    "heavy_atom_map_sdf_to_pdbqt": maptext,
                    "rigid_core_sdf_indices_1based": ";".join(map(str, core)),
                    "rigid_core_pdbqt_serials_1based": ";".join(str(qc["map"][i]) for i in core),
                    "rigid_core_automorphisms": auto,
                    "maximum_source_heavy_centroid_radius_A": f"{max(source['radii']):.4f}",
                    "centered_box_clearance_A": f"{centered_clearance:.4f}",
                    "receptor_state": receptor_state, "receptor_pdbqt_path": receptor_path,
                    "receptor_pdbqt_SHA256": receptor_sha,
                    "receptor_source_pdb_path": common.RECEPTOR_PDB[receptor_state][0],
                    "receptor_source_pdb_SHA256": common.RECEPTOR_PDB[receptor_state][1],
                    "MIG_reference_sdf_path": common.NATIVE[0], "MIG_reference_sdf_SHA256": common.NATIVE[1],
                    "MIG_control_csv_path": rel(gallic.MIG_RESULTS), "MIG_control_csv_SHA256": control_sha,
                    "MIG_control_status": "5/5 rank1 RMSD <=2 A for this receptor; no rerun",
                    "Vina_executable_SHA256": common.VINA[1], "Vina_version": "1.2.7",
                    "box_center_A": ";".join(f"{x:.8f}" for x in common.BOX_CENTER),
                    "box_size_A": "30;24;24",
                    "future_settings": "vina;exhaustiveness=64;num_modes=20;energy_range=5.0;min_rmsd=1.0;cpu=1",
                    "future_seeds": "1101;2202;3303;4404;5505",
                    "pocket_rule": "all heavy inside box; rigid-core centroid <=6 A from MIG centroid; >=50% core atoms <=4 A from MIG heavy",
                    "clash_rule": "receptor-ligand heavy <1.8 A or VdW overlap >0.8 A; nonbonded internal heavy <1.6 A",
                    "proposed_cluster_rule": "unfitted symmetry-aware rigid-core RMSD <=2 A to within-pair rank1 medoid; >=4/5 rank1 seeds pass full geometry rule",
                    "software_versions": f"RDKit {rdBase.rdkitVersion}; Meeko {meeko.__version__}; Vina 1.2.7",
                    "limitations": limitations,
                })
    if len(rows) != 6 or len({(r["structure_id"], r["microstate"], r["receptor_state"]) for r in rows}) != 6:
        raise ValueError("Expected six unique frozen chemical/receptor input pairs")
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    common.write_once(OUT, stream.getvalue())
    print("PASS: 24/24 LIG03 and 60/60 LIG09 conformers; six input pairs; no docking")
    for path in [audited["LIG03"]["sdf"], audited["LIG03"]["pdbqt"], VARIANT_SDF,
                 VARIANT_PDBQT, audited["LIG09"]["sdf"], audited["LIG09"]["pdbqt"], OUT]:
        print(rel(path), common.sha(path))


if __name__ == "__main__":
    main()

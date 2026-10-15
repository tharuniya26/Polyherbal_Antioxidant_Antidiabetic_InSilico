"""Audited 4GQR/MYC dry-pocket diagnostic redocking with Vina 1.2.7.

Run from the repository root with .venv/Scripts/python.exe. Existing sources and
completed seed artifacts are immutable; --preflight-only stops before docking.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from hashlib import sha256
import io
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import time

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from meeko import PDBQTMolecule, RDKitMolCreate


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "data/metadata/vina_alternative_receptor_preparation_results.csv"
ORIGINAL = ROOT / "structures/proteins/original/4GQR.pdb"
CCD = ROOT / "docking/validation/definitions/MYC.cif"
SOURCE_RECEPTOR = ROOT / "structures/proteins/vina_prepared/4GQR_validation_local_noNAG_Meeko0.8.0.pdbqt"
RECEPTOR = ROOT / "structures/proteins/vina_prepared/4GQR_dry_noNAG_CA_CL_PCA_pH7p5_diagnostic_Meeko0.8.0.pdbqt"
LIGAND = ROOT / "docking/validation/vina/4GQR_MYC_crystal_pose_reference.pdbqt"
REFERENCE = ROOT / "docking/validation/vina/4GQR_MYC_native_pose.sdf"
BRIDGE = ROOT / "docking/validation/vina/4GQR_MYC_sdf_to_pdbqt_atom_map.csv"
DEPOSITED = ROOT / "docking/validation/vina/4GQR_MYC_atom_map.csv"
VINA = Path(r"E:\docking tools\vina_1.2.7_win.exe")
OUT = ROOT / "docking/validation/vina/4GQR_diagnostic_redocking"
RESULTS = ROOT / "data/metadata/vina_4GQR_redocking_results.csv"
SEEDS = (1101, 2202, 3303, 4404, 5505)
EXHAUSTIVENESS = 32
NUM_MODES = 20
ENERGY_RANGE = 5.0
MIN_RMSD = 1.0
CPU = 1


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def immutable(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != data:
            raise RuntimeError(f"Refusing to overwrite changed artifact: {path}")
    else:
        path.write_bytes(data)


def csv_data(rows: list[dict], fields: list[str]) -> bytes:
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
    return out.getvalue().encode("utf-8")


def audit_inputs() -> dict:
    with AUDIT.open(newline="", encoding="utf-8") as handle:
        row = next(r for r in csv.DictReader(handle) if r["PDB_ID"] == "4GQR")
    pairs = {
        ORIGINAL: "original_PDB_SHA256",
        ROOT / "structures/proteins/original/4GQR.cif": "original_mmCIF_SHA256",
        CCD: "CCD_SHA256",
        SOURCE_RECEPTOR: "local_receptor_PDBQT_SHA256",
        ROOT / row["local_receptor_PDB_path"]: "local_receptor_PDB_SHA256",
        REFERENCE: "native_SDF_SHA256",
        LIGAND: "native_PDBQT_SHA256",
        BRIDGE: "native_bridge_map_SHA256",
        DEPOSITED: "native_atom_map_SHA256",
        ROOT / row["template_path"]: "template_SHA256",
    }
    for path, key in pairs.items():
        if not path.is_file() or digest(path).lower() != row[key].lower():
            raise RuntimeError(f"Missing or changed audited input: {path} ({key})")
    if not VINA.is_file():
        raise RuntimeError(f"Vina executable missing: {VINA}")
    version = subprocess.run([str(VINA), "--version"], capture_output=True, text=True, timeout=20)
    if version.returncode or "1.2.7" not in version.stdout:
        raise RuntimeError(f"Vina version mismatch: {version.stdout} {version.stderr}")
    if row["native_ligand"] != "MYC" or row["native_heavy_atoms"] != "23":
        raise RuntimeError("4GQR/MYC preparation audit identity mismatch")
    return {str(p.relative_to(ROOT)): digest(p) for p in pairs}


def atom_records(text: str) -> dict[int, tuple[str, str, np.ndarray]]:
    records = {}
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        serial = int(line[6:11])
        if serial in records:
            raise RuntimeError(f"Duplicate PDBQT atom serial {serial}")
        records[serial] = (line[12:16].strip()[0], line.split()[-1],
                           np.array([float(line[i:i+8]) for i in (30, 38, 46)]))
    return records


def reference_data():
    molecules = [m for m in Chem.SDMolSupplier(str(REFERENCE), removeHs=False) if m is not None]
    if len(molecules) != 1:
        raise RuntimeError("MYC reference SDF does not contain exactly one valid molecule")
    mol = molecules[0]
    heavy = Chem.RemoveHs(mol)
    if (rdMolDescriptors.CalcMolFormula(mol) != "C15H10O8" or
            Chem.GetFormalCharge(mol) != 0 or heavy.GetNumAtoms() != 23 or
            Chem.FindMolChiralCenters(heavy, includeUnassigned=True)):
        raise RuntimeError("MYC formula, formal charge, heavy atoms or stereochemistry mismatch")
    xyz = np.asarray(mol.GetConformer().GetPositions()[:23], dtype=float)
    if not np.isfinite(xyz).all() or not mol.GetConformer().Is3D():
        raise RuntimeError("Invalid crystal MYC 3D coordinates")
    with BRIDGE.open(newline="", encoding="utf-8") as handle:
        mapping = {int(r["pdbqt_atom_serial"]): int(r["sdf_atom_index_1based"])-1
                   for r in csv.DictReader(handle)}
    if len(mapping) != 23 or set(mapping.values()) != set(range(23)):
        raise RuntimeError("Incomplete or nonbijective MYC heavy-atom bridge")
    with DEPOSITED.open(newline="", encoding="utf-8") as handle:
        deposit_rows = list(csv.DictReader(handle))
    if len(deposit_rows) != 23:
        raise RuntimeError("Incomplete deposited MYC atom map")
    pdb_serials = {int(line[6:11]): line for line in ORIGINAL.read_text().splitlines()
                   if line.startswith("HETATM")}
    for r in deposit_rows:
        idx, serial = int(r["sdf_atom_index_1based"])-1, int(r["PDB_serial"])
        line = pdb_serials.get(serial)
        if line is None or line[17:20] != "MYC" or line[12:16].strip() != r["auth_atom_name"]:
            raise RuntimeError(f"MYC deposited atom identity mismatch: {serial}")
        if max(abs(float(line[i:i+8])-xyz[idx,j]) for j,i in enumerate((30,38,46))) > 0.0011:
            raise RuntimeError(f"MYC crystal coordinate mismatch: {serial}")
    source = LIGAND.read_text()
    input_records = atom_records(source)
    if len(input_records) != 29 or sum(t != "HD" for _,t,_ in input_records.values()) != 23:
        raise RuntimeError("MYC PDBQT atom counts differ from audited neutral reference")
    for serial, idx in mapping.items():
        symbol, atom_type, pos = input_records[serial]
        if symbol != heavy.GetAtomWithIdx(idx).GetSymbol() or atom_type == "HD":
            raise RuntimeError(f"MYC PDBQT element/map mismatch: {serial}")
        if max(abs(pos-xyz[idx])) > 0.0011:
            raise RuntimeError(f"MYC PDBQT native-pose mismatch: {serial}")
    built = [m for m in RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(source)) if m]
    if len(built) != 1 or Chem.MolToSmiles(Chem.RemoveHs(built[0]), True) != Chem.MolToSmiles(heavy, True):
        raise RuntimeError("Meeko MYC graph/stereochemistry readback failed")
    automorphisms = heavy.GetSubstructMatches(heavy, uniquify=False, useChirality=True, maxMatches=1000)
    if not automorphisms:
        raise RuntimeError("No valid MYC graph automorphism")
    return heavy, xyz, mapping, input_records, automorphisms


def prepare_dry_receptor() -> dict:
    lines = SOURCE_RECEPTOR.read_text().splitlines()
    removed = [line for line in lines if line.startswith(("ATOM  ","HETATM")) and line[17:20] == "HOH"]
    keep = [line for line in lines if line not in removed]
    if len(removed) != 24 or {int(line[22:26]) for line in removed} != {905,1083,1101,1142,1144,1177,1204,1244}:
        raise RuntimeError("Unexpected local receptor water composition")
    if any(line[17:20] in ("CA ","CL ","PCA") for line in removed):
        raise RuntimeError("Ion or PCA unexpectedly removed")
    contents = ("\n".join(keep)+"\n").encode("ascii")
    immutable(RECEPTOR, contents)
    records = [line for line in keep if line.startswith(("ATOM  ","HETATM"))]
    if len(records) != 4864-24 or any(line[17:20] == "HOH" for line in records):
        raise RuntimeError("Dry receptor atom-count/water QC failed")
    ions = {(line[17:20].strip(), int(line[22:26])): (line.split()[-1], float(line.split()[-2]))
            for line in records if line[17:20].strip() in ("CA","CL")}
    if ions != {("CA",502): ("Ca",2.0), ("CL",503): ("Cl",-1.0)}:
        raise RuntimeError(f"Ca/Cl retention or PDBQT fields changed: {ions}")
    if not any(line[17:20] == "PCA" and int(line[22:26]) == 1 for line in records):
        raise RuntimeError("PCA A1 missing from dry receptor")
    if any(not np.isfinite(np.array([float(line[i:i+8]) for i in (30,38,46)])).all() for line in records):
        raise RuntimeError("Nonfinite receptor coordinate")
    return {"removed_waters": sorted({int(line[22:26]) for line in removed}),
            "receptor_atoms": len(records), "receptor_SHA256": digest(RECEPTOR),
            "ion_PDBQT_fields": {f"{name}_A{number}": value for (name, number), value in ions.items()}}


def box(xyz: np.ndarray):
    center = xyz.mean(axis=0)
    side = []
    for line in ORIGINAL.read_text().splitlines():
        if line.startswith("ATOM  ") and line[21] == "A" and line[22:26].strip() in ("197","233","300","201"):
            if line[12:16].strip() in ("CG","OD1","OD2","CD","OE1","OE2","ND1","NE2","CE1","CD2"):
                side.append(np.array([float(line[i:i+8]) for i in (30,38,46)]))
    if len(side) != 15:
        raise RuntimeError(f"Unexpected catalytic side-chain atom coverage: {len(side)}")
    extent = np.max(abs(np.vstack([xyz]+side)-center), axis=0)
    dims = tuple(int(2*math.ceil(e+5.0)) for e in extent)
    if dims != (26,24,28):
        raise RuntimeError(f"Unexpected 4GQR box after source check: {dims}")
    return center, dims, extent


def mapped_xyz(records, mapping):
    pos = np.empty((23,3))
    for serial, idx in mapping.items():
        pos[idx] = records[serial][2]
    return pos


def symmetry_rmsd(pos, reference, automorphisms):
    return min(float(np.sqrt(np.mean(np.sum((pos[list(match)]-reference)**2,axis=1))))
               for match in automorphisms)


def validate_pose(records, mapping, input_records, heavy, crystal):
    if set(records) != set(input_records):
        raise RuntimeError("Docked MYC atom serial set changed")
    for serial, rec in records.items():
        if rec[:2] != input_records[serial][:2] or not np.isfinite(rec[2]).all():
            raise RuntimeError(f"Docked MYC atom identity/type/coordinate failure: {serial}")
    xyz = mapped_xyz(records,mapping)
    for bond in heavy.GetBonds():
        a,b=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        if abs(np.linalg.norm(xyz[a]-xyz[b])-np.linalg.norm(crystal[a]-crystal[b]))>0.08:
            raise RuntimeError(f"Docked MYC bond geometry changed: {a}-{b}")
    return xyz


def randomized_start(seed, source, center, dims, crystal, mapping, input_records, heavy, autos):
    rng=np.random.default_rng(seed)
    for attempt in range(100):
        q=rng.normal(size=4); q/=np.linalg.norm(q)
        w,x,y,z=q
        rot=np.array([[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],
                      [2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],
                      [2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]])
        shift=rng.uniform(-3,3,size=3)
        if np.linalg.norm(shift)<3: continue
        moved={serial:(rec[2]-center)@rot.T+center+shift for serial,rec in input_records.items()}
        trial={serial:(input_records[serial][0],input_records[serial][1],pos) for serial,pos in moved.items()}
        pos=validate_pose(trial,mapping,input_records,heavy,crystal)
        rmsd=symmetry_rmsd(pos,crystal,autos)
        if rmsd<5 or np.any(abs(pos-center)>np.asarray(dims)/2-0.5): continue
        lines=[]
        for line in source.splitlines():
            if line.startswith(("ATOM  ","HETATM")):
                xyz=moved[int(line[6:11])]
                line=line[:30]+"".join(f"{v:8.3f}" for v in xyz)+line[54:]
            lines.append(line)
        text="\n".join(lines)+"\n"
        parsed=atom_records(text)
        validate_pose(parsed,mapping,input_records,heavy,crystal)
        rebuilt=[m for m in RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(text)) if m]
        if len(rebuilt)!=1 or Chem.MolToSmiles(Chem.RemoveHs(rebuilt[0]),True)!=Chem.MolToSmiles(heavy,True):
            raise RuntimeError("Randomized MYC Meeko graph readback failed")
        return text,symmetry_rmsd(mapped_xyz(parsed,mapping),crystal,autos),shift.tolist(),attempt+1
    raise RuntimeError(f"Unable to randomize MYC away from crystal for seed {seed}")


def parse_output(path,mapping,input_records,heavy,crystal,autos):
    text=path.read_text()
    segments=re.findall(r"(?ms)^MODEL\s+(\d+)\s*\n(.*?)^ENDMDL\s*$",text)
    rebuilt=[m for m in RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(text)) if m]
    if not segments or len(rebuilt)!=1 or rebuilt[0].GetNumConformers()!=len(segments):
        raise RuntimeError("Meeko could not reconstruct all docked MYC modes")
    if Chem.MolToSmiles(Chem.RemoveHs(rebuilt[0]),True)!=Chem.MolToSmiles(heavy,True):
        raise RuntimeError("Docked MYC graph/stereochemistry differs from native CCD reference")
    result=[]
    for rank_text,segment in segments:
        rank=int(rank_text)
        if rank!=len(result)+1: raise RuntimeError("Nonsequential pose ranks")
        match=re.search(r"(?m)^REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)\s+",segment)
        if not match: raise RuntimeError("Missing Vina score")
        pos=validate_pose(atom_records(segment),mapping,input_records,heavy,crystal)
        result.append({"pose_rank":rank,"score_kcal_mol":float(match.group(1)),
                       "heavy_atom_RMSD_A":symmetry_rmsd(pos,crystal,autos)})
    if any(result[i]["score_kcal_mol"]<result[i-1]["score_kcal_mol"]-0.01 for i in range(1,len(result))):
        raise RuntimeError("Unexpected Vina score order")
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight-only",action="store_true")
    args=parser.parse_args()
    hashes=audit_inputs()
    heavy,crystal,mapping,input_records,autos=reference_data()
    receptor_qc=prepare_dry_receptor()
    center,dims,extent=box(crystal)
    manifest={"experiment":"4GQR/MYC dry no-NAG local-pocket diagnostic redocking",
              "created_UTC":datetime.now(timezone.utc).isoformat(),
              "vina_version":"1.2.7","vina_executable":str(VINA),"vina_SHA256":digest(VINA),
              "rdkit_version":Chem.rdBase.rdkitVersion,"meeko_version":__import__('meeko').__version__,
              "python_version":sys.version.split()[0],"audited_input_SHA256":hashes,
              "receptor_QC":receptor_qc,"reference_heavy_atoms":23,"reference_charge":0,
              "reference_CCD":"MYC","automorphisms":len(autos),
              "box_center_A":center.tolist(),"box_dimensions_A":dims,
              "box_max_site_deviation_A":extent.tolist(),
              "box_rule":"deposited MYC heavy-atom centroid; >=5 A beyond MYC and Asp197/Glu233/Asp300/His201 catalytic side-chain heavy atoms; rounded half-extents",
              "seeds":list(SEEDS),"exhaustiveness":EXHAUSTIVENESS,"num_modes":NUM_MODES,
              "energy_range_kcal_mol":ENERGY_RANGE,"min_rmsd_A":MIN_RMSD,"cpu":CPU,
              "starting_pose":"NumPy PCG64(seed) rigid quaternion plus +/-3 A translation, >=5 A crystal-frame RMSD; Vina uses its independent seeded global search",
              "model_limit":"Dry no-NAG local model; retained Ca/Cl/PCA. Standard Vina ion scoring, pH7.5 catalytic/phenol states, omitted waters, MYC 0.4 occupancy and geometry outliers unresolved."}
    OUT.mkdir(parents=True,exist_ok=True)
    manifest_file=OUT/'experiment_manifest.json'
    if manifest_file.exists(): manifest['created_UTC']=json.loads(manifest_file.read_text())['created_UTC']
    immutable(manifest_file,(json.dumps(manifest,indent=2,default=str)+'\n').encode())
    print('Preflight passed:',len(input_records),'ligand atoms,',receptor_qc['receptor_atoms'],'receptor atoms; box',center,dims,flush=True)
    if args.preflight_only: return 0
    rows=[]; all_poses=[]; source=LIGAND.read_text()
    for seed in SEEDS:
        prefix=f'seed_{seed}'
        start=OUT/f'{prefix}_randomized_start.pdbqt'
        output=OUT/f'{prefix}_docked.pdbqt'
        log=OUT/f'{prefix}_vina.log'
        start_text,start_rmsd,shift,attempt=randomized_start(seed,source,center,dims,crystal,mapping,input_records,heavy,autos)
        immutable(start,start_text.encode('ascii'))
        cmd=[str(VINA),'--receptor',str(RECEPTOR),'--ligand',str(start),
             '--center_x',f'{center[0]:.8f}','--center_y',f'{center[1]:.8f}','--center_z',f'{center[2]:.8f}',
             '--size_x',str(dims[0]),'--size_y',str(dims[1]),'--size_z',str(dims[2]),
             '--scoring','vina','--exhaustiveness','32','--num_modes','20','--energy_range','5.0',
             '--min_rmsd','1.0','--seed',str(seed),'--cpu','1','--out',str(output),'--verbosity','1']
        command=subprocess.list2cmdline(cmd)
        if output.exists() and log.exists():
            stdout=log.read_text(); m=re.search(r'(?m)^EXIT_CODE: (\d+)$',stdout); elapsed_m=re.search(r'(?m)^ELAPSED_SECONDS: ([\d.]+)$',stdout)
            if command not in stdout or not m or not elapsed_m: raise RuntimeError(f'Existing run metadata incomplete: {seed}')
            code=int(m.group(1)); elapsed=float(elapsed_m.group(1)); print('Reanalysing',seed,flush=True)
        elif output.exists() or log.exists():
            raise RuntimeError(f'Partial existing run requires review: {seed}')
        else:
            print('Docking seed',seed,'start RMSD',round(start_rmsd,3),flush=True)
            t0=time.monotonic(); run=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,timeout=3600)
            elapsed=time.monotonic()-t0;code=run.returncode
            stdout=(f'COMMAND: {command}\nSEED: {seed}\nRECEPTOR_SHA256: {digest(RECEPTOR)}\n'
                    f'REFERENCE_LIGAND_SHA256: {digest(LIGAND)}\nSTART_SHA256: {digest(start)}\n'
                    f'EXIT_CODE: {code}\nELAPSED_SECONDS: {elapsed:.3f}\n--- STDOUT ---\n{run.stdout}\n--- STDERR ---\n{run.stderr}\n')
            immutable(log,stdout.encode())
        row={'PDB_ID':'4GQR','ligand':'MYC neutral CCD reference','seed':seed,'status':'failed',
             'receptor_model':'dry_noNAG_CA_CL_PCA_pH7p5_local_diagnostic','vina_version':'1.2.7',
             'receptor_path':str(RECEPTOR.relative_to(ROOT)),'receptor_SHA256':digest(RECEPTOR),
             'ligand_path':str(LIGAND.relative_to(ROOT)),'ligand_SHA256':digest(LIGAND),
             'crystal_reference_path':str(REFERENCE.relative_to(ROOT)),'crystal_reference_SHA256':digest(REFERENCE),
             'start_path':str(start.relative_to(ROOT)),'start_SHA256':digest(start),'start_RMSD_A':f'{start_rmsd:.4f}',
             'start_translation_xyz_A':';'.join(f'{x:.4f}' for x in shift),'randomization_attempt':attempt,
             'box_center_xyz_A':';'.join(f'{x:.8f}' for x in center),'box_size_xyz_A':';'.join(map(str,dims)),
             'exhaustiveness':32,'num_modes_requested':20,'energy_range_kcal_mol':5.0,'min_rmsd_A':1.0,'cpu':1,
             'command':command,'log_path':str(log.relative_to(ROOT)),'log_SHA256':digest(log),
             'output_path':str(output.relative_to(ROOT)) if output.exists() else '',
             'output_SHA256':digest(output) if output.exists() else '',
             'exit_code':code,'elapsed_seconds':f'{elapsed:.3f}','n_modes_returned':'',
             'top_score_kcal_mol':'','top_RMSD_A':'','best_RMSD_A':'','best_RMSD_rank':'',
             'best_RMSD_score_kcal_mol':'','top_within_2A':'','any_within_2A':'','error':''}
        try:
            if code or not output.is_file() or not output.stat().st_size: raise RuntimeError(f'Vina exit {code}; see log')
            if f'Performing docking (random seed: {seed}' not in stdout: raise RuntimeError('Vina seed confirmation absent')
            poses=parse_output(output,mapping,input_records,heavy,crystal,autos)
            if not poses: raise RuntimeError('No validated returned poses')
            top=poses[0];best=min(poses,key=lambda p:p['heavy_atom_RMSD_A'])
            row.update(status='completed',n_modes_returned=len(poses),top_score_kcal_mol=f"{top['score_kcal_mol']:.3f}",
                       top_RMSD_A=f"{top['heavy_atom_RMSD_A']:.4f}",best_RMSD_A=f"{best['heavy_atom_RMSD_A']:.4f}",
                       best_RMSD_rank=best['pose_rank'],best_RMSD_score_kcal_mol=f"{best['score_kcal_mol']:.3f}",
                       top_within_2A=top['heavy_atom_RMSD_A']<=2.0,any_within_2A=best['heavy_atom_RMSD_A']<=2.0)
            pose_rows=[{'seed':seed,**p,'within_2A':p['heavy_atom_RMSD_A']<=2.0} for p in poses]
            all_poses.extend(pose_rows)
            immutable(OUT/f'{prefix}_pose_metrics.csv',csv_data(pose_rows,['seed','pose_rank','score_kcal_mol','heavy_atom_RMSD_A','within_2A']))
            print(seed,'completed; top RMSD',row['top_RMSD_A'],'best',row['best_RMSD_A'],'rank',row['best_RMSD_rank'],flush=True)
        except Exception as exc:
            row['error']=str(exc);print(seed,'FAILED',exc,flush=True)
        rows.append(row)
        if row['status']!='completed':break
    if rows: immutable(RESULTS,csv_data(rows,list(rows[0])))
    if all_poses: immutable(OUT/'all_pose_metrics.csv',csv_data(all_poses,['seed','pose_rank','score_kcal_mol','heavy_atom_RMSD_A','within_2A']))
    return 0 if len(rows)==5 and all(r['status']=='completed' for r in rows) else 1


if __name__=='__main__':
    try: sys.exit(main())
    except Exception as exc:
        print('FATAL PREFLIGHT/ANALYSIS ERROR:',exc,file=sys.stderr);sys.exit(2)

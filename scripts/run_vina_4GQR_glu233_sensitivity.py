"""One-variable 4GQR Glu233H/MYC diagnostic redocking with Vina 1.2.7.

Run from the repository root with .venv/Scripts/python.exe. Use --prepare-only
to construct and audit the receptor without docking. Previously completed E32
and E64 artifacts are read-only. Existing new artifacts are never overwritten.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time

import numpy as np
from meeko import PDBQTReceptor

import run_vina_4GQR_redocking as base
import run_vina_4GQR_convergence as conv


ROOT = base.ROOT
SOURCE_PDB = ROOT / "structures/proteins/prepared/4GQR_validation_local_noNAG_pH7p5.pdb"
OUT = ROOT / "docking/validation/vina/4GQR_Glu233H_sensitivity"
RECEPTOR_PDB = OUT / "4GQR_dry_noNAG_Glu233H_OE2_pH7p5.pdb"
RECEPTOR_PDBQT = OUT / "4GQR_dry_noNAG_Glu233H_OE2_pH7p5_Meeko0.8.0.pdbqt"
RESULTS = ROOT / "data/metadata/vina_4GQR_glu233_sensitivity_results.csv"
CLI = ROOT / ".venv/Scripts/mk_prepare_receptor.exe"
TEMPLATE = ROOT / ".venv/Lib/site-packages/meeko/data/residue_chem_templates.json"
PCA_TEMPLATE = ROOT / "docking/configs/vina_templates/PCA_C_terminal_Meeko_0.8.0.json"
EXPECTED_BASELINE = "71bee218c2cfede049ad8254bc39b33cb7c916055c0355e17236a2227a3e52b4"
WATERS = {905, 1083, 1101, 1142, 1144, 1177, 1204, 1244}
SEEDS = base.SEEDS


def atoms(text):
    return [l for l in text.splitlines() if l.startswith(("ATOM  ", "HETATM"))]


def xyz(line):
    return np.array([float(line[i:i+8]) for i in (30, 38, 46)], dtype=float)


def key(line):
    return line[21], int(line[22:26]), line[17:20].strip(), line[12:16].strip()


def verify_e64():
    if base.digest(base.RECEPTOR) != EXPECTED_BASELINE:
        raise RuntimeError("Existing E64 receptor hash changed")
    with conv.RESULTS.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 5 or [int(r["seed"]) for r in rows] != list(SEEDS):
        raise RuntimeError("Expected five E64 records")
    for row in rows:
        if (row["status"] != "completed" or row["receptor_SHA256"] != EXPECTED_BASELINE
                or row["exhaustiveness"] != "64" or row["box_center_xyz_A"] !=
                "13.44317391;15.13530435;39.43613043" or row["box_size_xyz_A"] != "26;24;28"):
            raise RuntimeError(f"Unverified E64 baseline seed {row['seed']}")
        for p, h in (("ligand_path", "ligand_SHA256"), ("crystal_reference_path", "crystal_reference_SHA256"),
                     ("output_path", "output_SHA256"), ("log_path", "log_SHA256"),
                     ("start_path", "start_SHA256")):
            path = ROOT / row[p]
            if not path.is_file() or base.digest(path) != row[h]:
                raise RuntimeError(f"Altered E64 artifact: {path}")
    return {int(r["seed"]): r for r in rows}


def acid_hydrogen(source):
    glu = {l[12:16].strip(): xyz(l) for l in atoms(source)
           if l[21] == "A" and l[22:26].strip() == "233" and l[17:20] == "GLU"}
    if not {"CD", "OE1", "OE2"} <= set(glu) or "HE2" in glu:
        raise RuntimeError("Glu233 heavy atoms/initial microstate unexpected")
    cd, oe1, oe2 = (glu[n] for n in ("CD", "OE1", "OE2"))
    d1, d2 = np.linalg.norm(cd-oe1), np.linalg.norm(cd-oe2)
    if not 1.1 < d1 < 1.4 or not 1.1 < d2 < 1.4:
        raise RuntimeError(f"Glu233 deposited carboxyl geometry unexpected: {d1}, {d2}")
    u = (cd-oe2)/d2
    v = (oe1-cd)-np.dot(oe1-cd, u)*u
    v /= np.linalg.norm(v)
    theta = math.radians(108.0)
    # Syn planar O-H arrangement toward the carbonyl oxygen, selected before
    # docking for carboxylic-acid chemistry, not for MYC pose recovery.
    h = oe2 + 0.98*(math.cos(theta)*u + math.sin(theta)*v)
    source_atoms = [l for l in atoms(source) if l[17:20] != "HOH"]
    others = [(np.linalg.norm(h-xyz(l)), l) for l in source_atoms
              if not (l[21] == "A" and l[22:26].strip() == "233")]
    nearest = min(others, key=lambda t: t[0])
    if nearest[0] < 1.4:
        raise RuntimeError(f"Added HE2 clashes with receptor: {nearest[0]:.3f} A {key(nearest[1])}")
    mol = base.Chem.SDMolSupplier(str(base.REFERENCE), removeHs=False)[0]
    if mol is None: raise RuntimeError("Missing native MYC crystal reference")
    min_ligand = float(np.min(np.linalg.norm(np.asarray(mol.GetConformer().GetPositions())-h, axis=1)))
    if min_ligand < 1.2: raise RuntimeError(f"Added HE2 clashes with crystal MYC: {min_ligand:.3f} A")
    return h, {"Glu233_CD_OE1_A": float(d1), "Glu233_CD_OE2_A": float(d2),
               "OE2_HE2_A": 0.98, "CD_OE2_HE2_deg": 108.0,
               "HE2_OE1_A": float(np.linalg.norm(h-oe1)),
               "nearest_other_receptor_atom_A": float(nearest[0]),
               "nearest_other_receptor_atom": str(key(nearest[1])),
               "nearest_MYC_atom_A": min_ligand,
               "HE2_xyz_A": h.tolist()}


def build_pdb(source, h):
    lines = source.splitlines()
    existing = atoms(source)
    water_ids = {int(l[22:26]) for l in existing if l[17:20] == "HOH"}
    if water_ids != WATERS or sum(l[17:20] == "HOH" for l in existing) != 24:
        raise RuntimeError(f"Local PDB water set differs from E64 dry-receptor audit: {water_ids}")
    if any(l[21] == "A" and l[22:26].strip() == "233" and l[12:16].strip() == "HE2" for l in existing):
        raise RuntimeError("Source already has Glu233 HE2")
    serial = max(int(l[6:11]) for l in existing)+1
    if serial > 99999: raise RuntimeError("PDB atom serial overflow")
    out = ["REMARK 900 4GQR DRY NO-NAG GLU233H(OE2) PH7.5 SENSITIVITY MODEL",
           "REMARK 901 GLH TEMPLATE; SYN PLANAR HE2; OTHER HEAVY/HYDROGEN ATOMS FIXED"]
    inserted = 0
    for line in lines:
        if line == "END" or (line.startswith(("ATOM  ", "HETATM")) and line[17:20] == "HOH"):
            continue
        out.append(line)
        if (line.startswith("ATOM  ") and line[21] == "A" and line[22:26].strip() == "233"
                and line[12:16].strip() == "HG3"):
            new = line[:6]+f"{serial:5d}"+" " + " HE2" + line[16:30]
            new += "".join(f"{v:8.3f}" for v in h)+f"{1.00:6.2f}{0.00:6.2f}"+"           H  "
            out.append(new); inserted += 1
    out.append("END")
    result = "\n".join(out)+"\n"
    if inserted != 1 or len(atoms(result)) != len(existing)-24+1:
        raise RuntimeError("Glu233H PDB atom count mismatch")
    old = {int(l[6:11]): l for l in existing if l[17:20] != "HOH"}
    new = {int(l[6:11]): l for l in atoms(result)}
    if any(new.get(s) != l for s,l in old.items()):
        raise RuntimeError("An existing receptor atom changed in the new PDB")
    return result


def prep_receptor():
    if not CLI.is_file() or not PCA_TEMPLATE.is_file() or not TEMPLATE.is_file():
        raise RuntimeError("Meeko CLI or residue template missing")
    glh = json.loads(TEMPLATE.read_text())["residue_templates"]["GLH"]
    if not {"OE1", "OE2", "HE2"} <= set(glh["atom_name"]):
        raise RuntimeError("Installed GLH template does not protonate OE2")
    source = SOURCE_PDB.read_text()
    h, geometry = acid_hydrogen(source)
    model = build_pdb(source, h)
    with tempfile.TemporaryDirectory() as td:
        temp = Path(td)
        inp, out = temp/"Glu233H.pdb", temp/"Glu233H.pdbqt"
        inp.write_text(model)
        cmd = [str(CLI), "--read_pdb", str(inp), "--write_pdbqt", str(out),
               "--add_templates", str(PCA_TEMPLATE), "--set_template", "A:233=GLH"]
        run = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=300)
        if run.returncode or not out.is_file():
            raise RuntimeError(f"Meeko GLH conversion failed: {run.stdout[-1200:]} {run.stderr[-1200:]}")
        new_bytes = out.read_bytes()
    parsed = PDBQTReceptor.from_pdbqt_filename(str(base.RECEPTOR)).atoms()
    # The temporary Meeko output is saved only after comprehensive comparison.
    with tempfile.TemporaryDirectory() as td:
        check_file = Path(td)/"check.pdbqt"; check_file.write_bytes(new_bytes)
        new_atoms = PDBQTReceptor.from_pdbqt_filename(str(check_file)).atoms()
    atomkey = lambda a: (str(a["chain"]),int(a["resid"]),str(a["resname"]),str(a["name"]))
    old, new = ({atomkey(a):a for a in group} for group in (parsed,new_atoms))
    if len(old)!=len(parsed) or len(new)!=len(new_atoms):
        raise RuntimeError("Duplicate receptor PDBQT atom identity")
    added = set(new)-set(old)
    if set(old)-set(new) or added != {("A",233,"GLU","HE2")}:
        raise RuntimeError(f"Unexpected receptor PDBQT atom change: removed={set(old)-set(new)}, added={added}")
    if new[("A",233,"GLU","HE2")]["atom_type"] != "HD":
        raise RuntimeError("Glu233 HE2 not typed as donor hydrogen")
    for k,a in old.items():
        b = new[k]
        if np.max(np.abs(np.asarray(a["xyz"])-np.asarray(b["xyz"])))>0.0011:
            raise RuntimeError(f"Existing heavy/hydrogen coordinate changed: {k}")
        if k[1] != 233 and (a["atom_type"] != b["atom_type"] or
                            abs(float(a["partial_charges"])-float(b["partial_charges"]))>0.0011):
            raise RuntimeError(f"Atom type/charge changed outside Glu233: {k}")
        if not np.isfinite(float(b["partial_charges"])):
            raise RuntimeError(f"Nonfinite PDBQT charge: {k}")
    typed = {k[3]:str(a["atom_type"]) for k,a in new.items() if k[1]==233 and k[0]=="A"}
    if typed.get("OE1")!="OA" or typed.get("OE2")!="OA":
        raise RuntimeError(f"Glu233 oxygen atom typing changed unexpectedly: {typed}")
    q_old = sum(float(a["partial_charges"]) for k,a in old.items() if k[0]=="A" and k[1]==233)
    q_new = sum(float(a["partial_charges"]) for k,a in new.items() if k[0]=="A" and k[1]==233)
    if abs((q_new-q_old)-1.0)>0.06 or abs(q_new)>0.06:
        raise RuntimeError(f"GLH charge bookkeeping failed: old={q_old:.3f}, new={q_new:.3f}")
    for name,number,typ,charge in (("CA",502,"Ca",2.0),("CL",503,"Cl",-1.0)):
        match=[a for k,a in new.items() if k[1]==number and k[2]==name]
        if len(match)!=1 or match[0]["atom_type"]!=typ or abs(float(match[0]["partial_charges"])-charge)>0.001:
            raise RuntimeError(f"Retained ion {name} changed")
    if not any(k[2]=="PCA" and k[1]==1 for k in new):
        raise RuntimeError("Terminal PCA missing")
    if any(k[2]=="HOH" for k in new) or any(k[2]=="NAG" for k in new):
        raise RuntimeError("Water/glycan status changed")
    base.immutable(RECEPTOR_PDB,model.encode("ascii"))
    base.immutable(RECEPTOR_PDBQT,new_bytes)
    # Deliberate invalid zero-size grid checks Vina atom parsing without scoring.
    probe = subprocess.run([str(base.VINA),"--receptor",str(RECEPTOR_PDBQT),
                            "--ligand",str(base.LIGAND),"--center_x","13.44317391",
                            "--center_y","15.13530435","--center_z","39.43613043",
                            "--size_x","0","--size_y","24","--size_z","28","--score_only"],
                           capture_output=True,text=True,timeout=30)
    if probe.returncode != 1 or "Grid box dimensions must be greater than 0 Angstrom" not in probe.stdout+probe.stderr:
        raise RuntimeError(f"Vina input parse probe failed: {probe.stdout[-400:]} {probe.stderr[-400:]}")
    return {"geometry":geometry,"baseline_receptor_SHA256":base.digest(base.RECEPTOR),
            "new_PDB_SHA256":base.digest(RECEPTOR_PDB),"new_PDBQT_SHA256":base.digest(RECEPTOR_PDBQT),
            "GLH_template_SHA256":base.digest(TEMPLATE),"PCA_template_SHA256":base.digest(PCA_TEMPLATE),
            "Glu233_charge_baseline":q_old,"Glu233_charge_Glu233H":q_new,
            "PDBQT_atom_count":len(new),"PDBQT_HE2_type":typed["HE2"],
            "Meeko_command":"mk_prepare_receptor --read_pdb Glu233H.pdb --write_pdbqt Glu233H.pdbqt --add_templates PCA_C_terminal_Meeko_0.8.0.json --set_template A:233=GLH",
            "Vina_input_probe":"Reached deliberate zero-size-grid guard without scoring"}


def one_seed(seed, e64, ref):
    heavy, crystal, mapping, input_records, autos = ref
    source = base.LIGAND.read_text()
    center, dims, _ = base.box(crystal)
    start_text, start_rmsd, shift, attempt = base.randomized_start(
        seed,source,center,dims,crystal,mapping,input_records,heavy,autos)
    old_start = ROOT/e64["start_path"]
    if old_start.read_bytes()!=start_text.encode("ascii"):
        raise RuntimeError(f"Seed {seed} randomized start differs from frozen E64 start")
    start=OUT/f"seed_{seed}_randomized_start.pdbqt"
    output=OUT/f"seed_{seed}_docked.pdbqt"
    log=OUT/f"seed_{seed}_vina.log"
    base.immutable(start,start_text.encode("ascii"))
    cmd=[str(base.VINA),"--receptor",str(RECEPTOR_PDBQT),"--ligand",str(start),
         "--center_x","13.44317391","--center_y","15.13530435","--center_z","39.43613043",
         "--size_x","26","--size_y","24","--size_z","28","--scoring","vina",
         "--exhaustiveness","64","--num_modes","20","--energy_range","5.0",
         "--min_rmsd","1.0","--seed",str(seed),"--cpu","1","--out",str(output),"--verbosity","1"]
    command=subprocess.list2cmdline(cmd)
    if output.exists() and log.exists():
        stdout=log.read_text()
        code_match=re.search(r"(?m)^EXIT_CODE: (\d+)\s*$",stdout)
        time_match=re.search(r"(?m)^ELAPSED_SECONDS: ([\d.]+)\s*$",stdout)
        if command not in stdout or not code_match or not time_match:
            raise RuntimeError(f"Existing sensitivity run {seed} metadata incomplete")
        code,elapsed=int(code_match.group(1)),float(time_match.group(1))
        print("Reanalysing sensitivity seed",seed,flush=True)
    elif output.exists() or log.exists():
        raise RuntimeError(f"Partial sensitivity seed {seed} requires manual review")
    else:
        print("Docking sensitivity seed",seed,"start RMSD",round(start_rmsd,4),flush=True)
        t0=time.monotonic()
        try:
            run=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,timeout=3600)
            code,elapsed=run.returncode,time.monotonic()-t0
            captured=f"--- STDOUT ---\n{run.stdout}\n--- STDERR ---\n{run.stderr}\n"
        except subprocess.TimeoutExpired as exc:
            code,elapsed=124,time.monotonic()-t0
            captured=f"ERROR: {exc}\n"
        stdout=(f"COMMAND: {command}\nSEED: {seed}\nRECEPTOR_SHA256: {base.digest(RECEPTOR_PDBQT)}\n"
                f"LIGAND_SHA256: {base.digest(base.LIGAND)}\nSTART_SHA256: {base.digest(start)}\n"
                f"EXIT_CODE: {code}\nELAPSED_SECONDS: {elapsed:.3f}\n{captured}")
        base.immutable(log,stdout.encode())
    row={"PDB_ID":"4GQR","ligand":"MYC neutral CCD reference","seed":seed,"status":"failed",
         "receptor_model":"dry_noNAG_CA_CL_PCA_Glu233H_OE2_pH7p5_sensitivity",
         "vina_version":"1.2.7","meeko_version":"0.8.0",
         "receptor_path":str(RECEPTOR_PDBQT.relative_to(ROOT)),"receptor_SHA256":base.digest(RECEPTOR_PDBQT),
         "ligand_path":str(base.LIGAND.relative_to(ROOT)),"ligand_SHA256":base.digest(base.LIGAND),
         "crystal_reference_path":str(base.REFERENCE.relative_to(ROOT)),"crystal_reference_SHA256":base.digest(base.REFERENCE),
         "start_path":str(start.relative_to(ROOT)),"start_SHA256":base.digest(start),
         "start_RMSD_A":f"{start_rmsd:.4f}","start_translation_xyz_A":";".join(f"{x:.4f}" for x in shift),
         "randomization_attempt":attempt,"box_center_xyz_A":"13.44317391;15.13530435;39.43613043",
         "box_size_xyz_A":"26;24;28","exhaustiveness":64,"num_modes_requested":20,
         "energy_range_kcal_mol":5.0,"min_rmsd_A":1.0,"cpu":1,"command":command,
         "log_path":str(log.relative_to(ROOT)),"log_SHA256":base.digest(log),
         "output_path":str(output.relative_to(ROOT)) if output.exists() else "",
         "output_SHA256":base.digest(output) if output.exists() else "",
         "exit_code":code,"elapsed_seconds":f"{elapsed:.3f}","n_modes_returned":"",
         "top_score_kcal_mol":"","top_RMSD_A":"","best_RMSD_A":"","best_RMSD_rank":"",
         "best_RMSD_score_kcal_mol":"","top_within_2A":"","any_within_2A":"",
         "E64_top_RMSD_A":e64["top_RMSD_A"],"E64_best_RMSD_A":e64["best_RMSD_A"],
         "E64_best_RMSD_rank":e64["best_RMSD_rank"],"E64_top_score_kcal_mol":e64["top_score_kcal_mol"],"error":""}
    poses=[]
    try:
        if code or not output.is_file() or not output.stat().st_size:
            raise RuntimeError(f"Vina exit {code}; inspect log")
        if f"Performing docking (random seed: {seed}" not in stdout:
            raise RuntimeError("Vina log did not confirm seed")
        poses=base.parse_output(output,mapping,input_records,heavy,crystal,autos)
        if not poses: raise RuntimeError("No validated returned MYC modes")
        top,best=poses[0],min(poses,key=lambda p:p["heavy_atom_RMSD_A"])
        row.update(status="completed",n_modes_returned=len(poses),
                   top_score_kcal_mol=f"{top['score_kcal_mol']:.3f}",
                   top_RMSD_A=f"{top['heavy_atom_RMSD_A']:.4f}",
                   best_RMSD_A=f"{best['heavy_atom_RMSD_A']:.4f}",
                   best_RMSD_rank=best["pose_rank"],
                   best_RMSD_score_kcal_mol=f"{best['score_kcal_mol']:.3f}",
                   top_within_2A=top["heavy_atom_RMSD_A"]<=2.0,
                   any_within_2A=best["heavy_atom_RMSD_A"]<=2.0)
        metric=[{"seed":seed,**p,"within_2A":p["heavy_atom_RMSD_A"]<=2.0} for p in poses]
        base.immutable(OUT/f"seed_{seed}_pose_metrics.csv",base.csv_data(
            metric,["seed","pose_rank","score_kcal_mol","heavy_atom_RMSD_A","within_2A"]))
    except Exception as exc:
        row["error"]=str(exc)
    print("Seed",seed,row["status"],"top",row["top_RMSD_A"],"best",row["best_RMSD_A"],row["error"],flush=True)
    return row,[{"seed":seed,**p,"within_2A":p["heavy_atom_RMSD_A"]<=2.0} for p in poses]


def main():
    prepare_only=sys.argv[1:]==["--prepare-only"]
    if sys.argv[1:] and not prepare_only:
        raise RuntimeError("Only optional --prepare-only is supported")
    source_hashes=base.audit_inputs()
    ref=base.reference_data()
    if len(ref[-1])!=2: raise RuntimeError("MYC graph automorphism count differs from audited reference")
    base.prepare_dry_receptor()  # immutable baseline verification
    e64=verify_e64()
    qc=prep_receptor()
    OUT.mkdir(parents=True,exist_ok=True)
    manifest={"experiment":"4GQR/MYC one-variable Glu233H sensitivity", "created_UTC":datetime.now(timezone.utc).isoformat(),
              "source_SHA256":source_hashes,"E64_results_SHA256":base.digest(conv.RESULTS),
              "receptor_QC":qc,"Vina_executable_SHA256":base.digest(base.VINA),
              "Vina_version":"1.2.7","Meeko_version":"0.8.0","RDKit_version":base.Chem.rdBase.rdkitVersion,
              "seeds":list(SEEDS),"exhaustiveness":64,"num_modes":20,"energy_range_kcal_mol":5.0,
              "min_rmsd_A":1.0,"cpu_per_run":1,"box_center_A":[13.44317391,15.13530435,39.43613043],
              "box_size_A":[26,24,28],"MYC_heavy_atoms":23,"MYC_graph_automorphisms":2,
              "starting_pose":"Exact E64 independent randomized start per seed; crystal reference never a search input",
              "limitations":"Provisional Glu233H; dry/no-NAG local model; standard Vina ion and fixed-water limits; MYC low occupancy."}
    manifest_file=OUT/"experiment_manifest.json"
    if manifest_file.exists(): manifest["created_UTC"]=json.loads(manifest_file.read_text())["created_UTC"]
    base.immutable(manifest_file,(json.dumps(manifest,indent=2)+"\n").encode())
    print("PREPARATION QC PASSED",json.dumps(qc,indent=2),flush=True)
    if prepare_only:return 0
    results=[]
    with ThreadPoolExecutor(max_workers=5) as executor:
        jobs={executor.submit(one_seed,seed,e64[seed],ref):seed for seed in SEEDS}
        for job in as_completed(jobs):
            try: results.append(job.result())
            except Exception as exc:
                print("FATAL SEED",jobs[job],exc,file=sys.stderr,flush=True)
                raise
    rows=[r for r,_ in sorted(results,key=lambda pair:int(pair[0]["seed"]))]
    all_poses=[p for _,poses in sorted(results,key=lambda pair:int(pair[0]["seed"])) for p in poses]
    base.immutable(RESULTS,base.csv_data(rows,list(rows[0])))
    if all_poses:
        base.immutable(OUT/"all_pose_metrics.csv",base.csv_data(
            all_poses,["seed","pose_rank","score_kcal_mol","heavy_atom_RMSD_A","within_2A"]))
    return 0 if all(r["status"]=="completed" for r in rows) else 1


if __name__=="__main__":
    try: sys.exit(main())
    except Exception as exc:
        print("FATAL GLU233H PREFLIGHT OR RUN ERROR:",exc,file=sys.stderr,flush=True)
        sys.exit(2)

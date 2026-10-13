"""Read-only chemistry audit plus Vina --score_only input test for 1XD0.

The score-only output tests format/atom-type acceptance; it is not redocking
or a binding-affinity estimate. Existing coordinate files are never edited.
"""
from __future__ import annotations

import csv
import hashlib
import io
import math
import re
import subprocess
from collections import Counter
from pathlib import Path

import numpy as np
from meeko import PDBQTMolecule, PDBQTReceptor, RDKitMolCreate
from rdkit import Chem

ROOT=Path(__file__).resolve().parents[1]
VINA=Path(r'E:\docking tools\vina_1.2.7_win.exe')
AUDIT=ROOT/'data/metadata/vina_receptor_blocker_resolution_results.csv'
LIG_AUDIT=ROOT/'data/metadata/vina_receptor_preparation_results.csv'
SDF=ROOT/'docking/validation/1XD0_ARE_native_pose.sdf'
OUT=ROOT/'data/metadata/vina_1xd0_input_preflight_results.csv'
BOX_MARGIN=6.0
FIELDS=['variant','receptor_path','receptor_sha256','ligand_path','ligand_sha256',
        'vina_version','meeko_version','score_only_exit_code','score_only_parse_status',
        'diagnostic_score_kcal_mol_not_affinity','box_center_xyz_A','box_size_xyz_A',
        'native_heavy_atoms','native_torsdof','native_formal_charge','native_stereo_graph_check',
        'PCA_C_Tyr2_N_distance_A','PCA_charge','calcium_type_charge','chloride_type_charge',
        'catalytic_microstates','water_count','water_atom_types','water_OH_range_A',
        'water_HOH_angle_range_deg','Ca_coordination_waters','min_ligand_protein_heavy_A',
        'min_ligand_water_O_A','min_ligand_chloride_A','min_ligand_calcium_A',
        'heavy_contacts_under_2A','physical_scoring_status','redocking_readiness','notes']

def rows(path):
    with path.open(newline='',encoding='utf-8') as f: return list(csv.DictReader(f))

def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def verify(path, expected):
    if not path.is_file() or digest(path)!=expected: raise RuntimeError(f'missing or SHA256 mismatch: {path}')

def canon(m): return Chem.MolToSmiles(Chem.RemoveHs(m),isomericSmiles=True)

def csv_write_new(path, items):
    f=io.StringIO(newline=''); w=csv.DictWriter(f,fieldnames=FIELDS)
    w.writeheader(); w.writerows(items); data=f.getvalue().encode('utf-8')
    if path.exists():
        if path.read_bytes()!=data: raise RuntimeError(f'refusing to overwrite differing report: {path}')
    else:
        path.parent.mkdir(parents=True,exist_ok=True); path.write_bytes(data)

def ligand():
    mols=list(Chem.SDMolSupplier(str(SDF),removeHs=False))
    if len(mols)!=1 or mols[0] is None: raise RuntimeError('ARE SDF not one valid record')
    mol=mols[0]
    path=ROOT/'docking/validation/vina/1XD0_ARE_crystal_pose_reference.pdbqt'
    pdbqt=PDBQTMolecule(path.read_text())
    made=RDKitMolCreate.from_pdbqt_mol(pdbqt)
    if len(made)!=1 or made[0] is None or canon(made[0])!=canon(mol):
        raise RuntimeError('ARE PDBQT molecular graph/stereo differs from deposited SDF')
    atoms=pdbqt.atoms()
    nheavy=sum(a.GetAtomicNum()>1 for a in mol.GetAtoms())
    if nheavy!=55 or sum(a['atom_type']!='HD' for a in atoms)!=nheavy:
        raise RuntimeError('ARE heavy atom count mismatch')
    tors=[int(l.split()[1]) for l in path.read_text().splitlines() if l.startswith('TORSDOF ')]
    if tors!=[28]: raise RuntimeError('ARE torsion count mismatch')
    xyz=np.array([list(mol.GetConformer().GetAtomPosition(a.GetIdx())) for a in mol.GetAtoms() if a.GetAtomicNum()>1])
    center=xyz.mean(axis=0)
    size=np.ceil(2*np.max(np.abs(xyz-center),axis=0)+2*BOX_MARGIN).astype(int)
    if tuple(size)!=(32,23,26): raise RuntimeError(f'ARE box changed: {size}')
    if np.any((center-size/2)>xyz.min(axis=0)-BOX_MARGIN+1) or np.any((center+size/2)<xyz.max(axis=0)+BOX_MARGIN-1):
        raise RuntimeError('box margin calculation failed')
    charge=sum(a.GetFormalCharge() for a in mol.GetAtoms())
    return path,mol,xyz,center,size,charge

def vina_probe(rec,lig,center,size):
    command=[str(VINA),'--receptor',str(rec),'--ligand',str(lig),
             '--center_x',f'{center[0]:.3f}','--center_y',f'{center[1]:.3f}',
             '--center_z',f'{center[2]:.3f}','--size_x',str(size[0]),
             '--size_y',str(size[1]),'--size_z',str(size[2]),
             '--scoring','vina','--score_only','--cpu','1']
    run=subprocess.run(command,capture_output=True,text=True,timeout=120)
    combined=run.stdout+'\n'+run.stderr
    match=re.search(r'Estimated Free Energy of Binding\s*:\s*([-+]?\d+(?:\.\d+)?)',combined)
    if run.returncode!=0 or match is None: raise RuntimeError(f'Vina score-only input test failed: {combined[-2000:]}')
    return run.returncode,match.group(1)

def coord(a): return np.array(a['xyz'],dtype=float)

def chemistry(rec,lig_xyz):
    atoms=PDBQTReceptor.from_pdbqt_filename(str(rec)).atoms()
    if not np.isfinite(atoms['xyz']).all() or not np.isfinite(atoms['partial_charges']).all():
        raise RuntimeError('receptor coordinates or charges invalid')
    def residue(resname,resid): return [a for a in atoms if a['resname']==resname and a['resid']==resid and a['chain']=='A']
    pca=residue('PCA',1);tyr=residue('TYR',2)
    cp=next(a for a in pca if a['name']=='C'); nt=next(a for a in tyr if a['name']=='N')
    linkage=float(np.linalg.norm(coord(cp)-coord(nt)))
    if not 1.2<linkage<1.5 or len([a for a in pca if a['atom_type']!='HD'])!=8:
        raise RuntimeError('PCA--Tyr2 peptide geometry or PCA heavy atom count invalid')
    pca_charge=float(sum(a['partial_charges'] for a in pca))
    ca=residue('CA',498); cl=residue('CL',500)
    if len(ca)!=1 or len(cl)!=1 or ca[0]['atom_type']!='Ca' or cl[0]['atom_type']!='Cl':
        raise RuntimeError('ion types/identities invalid')
    if abs(float(ca[0]['partial_charges'])-2)>0.001 or abs(float(cl[0]['partial_charges'])+1)>0.001:
        raise RuntimeError('ion charge fields invalid')
    micro=[]
    for name,rid in [('ASP',197),('GLU',233),('ASP',300),('HIS',201)]:
        res=residue(name,rid); names={str(a['name']) for a in res}; q=float(sum(a['partial_charges'] for a in res))
        micro.append(f'{name}{rid}:{q:+.3f}:{"/".join(sorted(names & {"HD1","HE2","HD2","HE1","HE2"})) or "no sidechain polar H"}')
    waters={int(a['resid']) for a in atoms if a['resname']=='HOH'}
    lengths=[];angles=[]; water_oxygen=[]; water_types=Counter()
    for rid in waters:
        group=residue('HOH',rid); water_types.update(str(a['atom_type']) for a in group)
        oxygen=[a for a in group if a['name']=='OW']; hs=[a for a in group if a['atom_type']=='HD']
        if len(oxygen)!=1 or len(hs)!=2: raise RuntimeError(f'water {rid} not OW + two HD')
        o=coord(oxygen[0]); water_oxygen.append(o)
        u=coord(hs[0])-o;v=coord(hs[1])-o
        lengths.extend([float(np.linalg.norm(u)),float(np.linalg.norm(v))])
        angles.append(float(np.degrees(np.arccos(np.clip(np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)),-1,1)))))
    if 536 not in waters or 546 not in waters: raise RuntimeError('calcium waters missing')
    protein=np.array([coord(a) for a in atoms if a['atom_type']!='HD' and a['resname'] not in ('HOH','CA','CL')])
    def minimum(x): return float(np.linalg.norm(lig_xyz[:,None,:]-x[None,:,:],axis=2).min())
    protein_dist=np.linalg.norm(lig_xyz[:,None,:]-protein[None,:,:],axis=2)
    return {'pca_link':linkage,'pca_q':pca_charge,'ca':ca[0],'cl':cl[0],
            'micro':';'.join(micro),'waters':len(waters),'water_types':dict(water_types),
            'oh_min':min(lengths),'oh_max':max(lengths),'angle_min':min(angles),'angle_max':max(angles),
            'min_protein':float(protein_dist.min()),'min_water':minimum(np.array(water_oxygen)),
            'min_cl':minimum(np.array([coord(cl[0])])),'min_ca':minimum(np.array([coord(ca[0])])),
            'under2':int((protein_dist<2).sum())}

def main():
    import meeko
    if meeko.__version__!='0.8.0': raise RuntimeError('Meeko 0.8.0 required')
    version=subprocess.run([str(VINA),'--version'],capture_output=True,text=True,timeout=15)
    if version.returncode or '1.2.7' not in version.stdout: raise RuntimeError('Vina 1.2.7 unavailable')
    audit={r['variant']:r for r in rows(AUDIT) if r['PDB_ID']=='1XD0'}
    old={r['variant']:r for r in rows(LIG_AUDIT) if r['PDB_ID']=='1XD0'}
    lig,mol,lig_xyz,center,size,charge=ligand()
    output=[]
    for variant in ('allwaters_CA_CL','sitewaters_CA_CL'):
        a=audit[variant];b=old[variant]
        rec=ROOT/a['PDBQT_path'];verify(rec,a['PDBQT_SHA256'])
        if lig!=ROOT/b['native_PDBQT_path']: raise RuntimeError('ARE reference path mismatch')
        verify(lig,b['native_PDBQT_SHA256']);verify(SDF,b['native_SDF_SHA256'])
        probe_code,diagnostic=vina_probe(rec,lig,center,size)
        q=chemistry(rec,lig_xyz)
        if q['waters']!=int(a['waters_retained']): raise RuntimeError('audited water count mismatch')
        row={k:'' for k in FIELDS}
        row.update({'variant':variant,'receptor_path':a['PDBQT_path'],'receptor_sha256':a['PDBQT_SHA256'],
            'ligand_path':b['native_PDBQT_path'],'ligand_sha256':b['native_PDBQT_SHA256'],
            'vina_version':'1.2.7','meeko_version':'0.8.0','score_only_exit_code':probe_code,
            'score_only_parse_status':'pass; native crystal pose only; no docking/minimization',
            'diagnostic_score_kcal_mol_not_affinity':diagnostic,
            'box_center_xyz_A':';'.join(f'{x:.3f}' for x in center),
            'box_size_xyz_A':';'.join(str(x) for x in size),
            'native_heavy_atoms':55,'native_torsdof':28,'native_formal_charge':charge,
            'native_stereo_graph_check':'pass; Meeko readback matches reference isomeric graph',
            'PCA_C_Tyr2_N_distance_A':f'{q["pca_link"]:.3f}',
            'PCA_charge':f'{q["pca_q"]:+.3f}',
            'calcium_type_charge':f'{q["ca"]["atom_type"]};{float(q["ca"]["partial_charges"]):+.3f}',
            'chloride_type_charge':f'{q["cl"]["atom_type"]};{float(q["cl"]["partial_charges"]):+.3f}',
            'catalytic_microstates':q['micro'],'water_count':q['waters'],
            'water_atom_types':';'.join(f'{k}:{v}' for k,v in sorted(q['water_types'].items())),
            'water_OH_range_A':f'{q["oh_min"]:.3f}-{q["oh_max"]:.3f}',
            'water_HOH_angle_range_deg':f'{q["angle_min"]:.1f}-{q["angle_max"]:.1f}',
            'Ca_coordination_waters':'HOH A536;HOH A546',
            'min_ligand_protein_heavy_A':f'{q["min_protein"]:.3f}',
            'min_ligand_water_O_A':f'{q["min_water"]:.3f}',
            'min_ligand_chloride_A':f'{q["min_cl"]:.3f}',
            'min_ligand_calcium_A':f'{q["min_ca"]:.3f}',
            'heavy_contacts_under_2A':q['under2'],
            'physical_scoring_status':'blocked; Vina ignores charge, Cl- typed as hydrophobic chlorine, Ca generic metal donor, static water geometry unreviewed',
            'redocking_readiness':'blocked',
            'notes':'Diagnostic score is an input-parsing check only. ARE neutral secondary amine microstate and catalytic states unresolved; ARE has crystal-mate contacts.'})
        output.append(row)
    csv_write_new(OUT,output)
    for row in output: print(row['variant'],row['score_only_parse_status'],row['physical_scoring_status'])

if __name__=='__main__': main()

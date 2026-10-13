"""Audit Meeko receptor conversion and export immutable-pose native ligand PDBQTs.

Run with .venv/Scripts/python.exe. No docking is performed. Existing source files
and previous GOLD outputs are read-only; all receptor CLI probes use temp files.
"""
from __future__ import annotations

import csv
import hashlib
import math
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

import meeko
import numpy as np
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTReceptor, PDBQTWriterLegacy, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
GOLD = ROOT / 'data/metadata/gold_receptor_preparation_results.csv'
PREP = ROOT / 'data/metadata/receptor_preparation_results.csv'
OUT = ROOT / 'data/metadata/vina_receptor_preparation_results.csv'
REC_OUT = ROOT / 'structures/proteins/vina_prepared'
LIG_OUT = ROOT / 'docking/validation/vina'
CLI = ROOT / '.venv/Scripts/mk_prepare_receptor.exe'
FIELDS = ['PDB_ID','variant','input_receptor_path','input_SHA256','Meeko_version',
          'receptor_conversion_status','receptor_pdbqt_path','receptor_pdbqt_SHA256',
          'input_heavy_atoms','pdbqt_heavy_atoms','pdbqt_atoms','coordinate_check',
          'retained_waters','retained_ions','retained_glycans_or_modified_residues',
          'native_ligand','native_SDF_SHA256','native_PDBQT_path','native_PDBQT_SHA256',
          'native_heavy_atoms','native_torsdof','native_roundtrip_status',
          'native_atom_map_path','site_center_xyz_A','proposed_box_xyz_A',
          'redocking_ready','blockers']
BOX = {'2QMJ': (23,19,29), '1XD0': (32,23,26)}

def rows(path):
    with path.open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def verify(path, expected):
    if not path.is_file() or sha(path) != expected:
        raise RuntimeError(f'missing or SHA256 mismatch: {path}')

def rel(path):
    return path.relative_to(ROOT).as_posix()

def safe_write(path, data):
    if path.exists():
        if path.read_bytes() != data:
            raise RuntimeError(f'refusing to overwrite changed output: {path}')
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)

def heavy_pdb_coords(path):
    coords=[]
    for l in path.read_text().splitlines():
        if l.startswith(('ATOM  ','HETATM')) and l[76:78].strip().upper() not in ('H','D'):
            coords.append(tuple(round(float(l[i:i+8]),3) for i in (30,38,46)))
    return coords

def receptor_probe(g):
    source=ROOT/g['output_path']; verify(source,g['output_SHA256'])
    # Meeko 0.8.0 supplies templates for standard protein, water, Ca and Cl;
    # no residue-deletion flags are used here.
    with tempfile.TemporaryDirectory() as td:
        target=Path(td)/'probe.pdbqt'
        run=subprocess.run([str(CLI),'--read_pdb',str(source),'--write_pdbqt',str(target)],
                           capture_output=True,text=True,timeout=180)
        if run.returncode:
            error=' '.join((run.stderr+' '+run.stdout).split())
            return {'status':'blocked_template','error':error[:1800]}
        if not target.is_file() or target.stat().st_size==0:
            return {'status':'failed_empty_output','error':'Meeko returned success without a nonempty PDBQT'}
        parsed=PDBQTReceptor.from_pdbqt_filename(str(target))
        atoms=parsed.atoms()
        if len(atoms)==0 or not np.isfinite(atoms['xyz']).all() or not np.isfinite(atoms['partial_charges']).all():
            return {'status':'failed_QC','error':'invalid atom coordinates or charges'}
        source_heavy=Counter(heavy_pdb_coords(source))
        target_heavy=Counter(tuple(round(float(v),3) for v in a['xyz']) for a in atoms if a['atom_type']!='HD')
        if source_heavy!=target_heavy:
            return {'status':'failed_QC','error':f'heavy-coordinate multiset mismatch: input {sum(source_heavy.values())}, output {sum(target_heavy.values())}'}
        data=target.read_bytes()
        return {'status':'preliminary_sensitivity_only','data':data,'input_heavy':sum(source_heavy.values()),
                'output_heavy':sum(target_heavy.values()),'atoms':len(atoms),
                'coordinate':'all input heavy coordinates retained to 0.001 A; Meeko readback finite'}

def index_map(text):
    mapping={}
    for l in text.splitlines():
        if l.startswith('REMARK INDEX MAP '):
            n=[int(x) for x in l[17:].split()]
            for a,b in zip(n[::2],n[1::2]):
                if a in mapping: raise RuntimeError('duplicate Meeko source index')
                mapping[a]=b
    if not mapping or len(set(mapping.values()))!=len(mapping):
        raise RuntimeError('missing/duplicate index map')
    return mapping

def canon(mol):
    return Chem.MolToSmiles(Chem.RemoveHs(mol),isomericSmiles=True)

def native_prepare(p):
    sdf=ROOT/p['native_ligand_SDF']; verify(sdf,p['native_ligand_SHA256'])
    prior=ROOT/p['atom_mapping_CSV']; verify(prior,p['atom_mapping_SHA256'])
    records=list(Chem.SDMolSupplier(str(sdf),removeHs=False))
    if len(records)!=1 or records[0] is None: raise RuntimeError('native SDF must have one valid record')
    mol=records[0]
    if not mol.GetConformer().Is3D() or not np.isfinite(mol.GetConformer().GetPositions()).all():
        raise RuntimeError('native reference has invalid 3D coordinates')
    prep=MoleculePreparation(add_index_map=True,charge_model='gasteiger')
    setups=prep.prepare(mol)
    if len(setups)!=1: raise RuntimeError('Meeko native ligand setup count differs from one')
    text,ok,error=PDBQTWriterLegacy.write_string(setups[0],add_index_map=True)
    if not ok: raise RuntimeError(str(error))
    parsed=PDBQTMolecule(text)
    atoms=parsed.atoms()
    made=RDKitMolCreate.from_pdbqt_mol(parsed)
    if len(made)!=1 or made[0] is None or canon(made[0])!=canon(mol):
        raise RuntimeError('native ligand graph/stereochemistry roundtrip failed')
    mapping=index_map(text)
    old=rows(prior)
    bridge=[]; maxshift=0.0
    for r in old:
        i=int(r['sdf_atom_index_1based']); atom=mol.GetAtomWithIdx(i-1)
        if atom.GetAtomicNum()==1 or i not in mapping: raise RuntimeError('native heavy atom mapping missing')
        n=mapping[i]; xyz=mol.GetConformer().GetAtomPosition(i-1)
        shift=math.dist((xyz.x,xyz.y,xyz.z),[float(v) for v in atoms['xyz'][n-1]])
        maxshift=max(maxshift,shift)
        bridge.append({'sdf_atom_index_1based':i,'PDB_serial':r['PDB_serial'],
                       'auth_atom_name':r['auth_atom_name'],'pdbqt_atom_serial':n})
    if len(bridge)!=sum(a.GetAtomicNum()>1 for a in mol.GetAtoms()) or maxshift>0.0011:
        raise RuntimeError('native pose atom count or coordinates mismatch')
    stem=p['PDB_ID']+'_'+p['reference_ligand_id']
    out=LIG_OUT/(stem+'_crystal_pose_reference.pdbqt')
    safe_write(out,text.encode())
    mapout=LIG_OUT/(stem+'_sdf_to_pdbqt_atom_map.csv')
    import io
    buf=io.StringIO(newline=''); writer=csv.DictWriter(buf,fieldnames=list(bridge[0]))
    writer.writeheader(); writer.writerows(bridge)
    safe_write(mapout,buf.getvalue().encode())
    tors=int([l.split()[1] for l in text.splitlines() if l.startswith('TORSDOF ')][0])
    return {'path':rel(out),'sha':sha(out),'map':rel(mapout),'heavy':len(bridge),
            'tors':tors,'check':f'pass graph/stereo and {len(bridge)} mapped heavy coordinates; max shift {maxshift:.4f} A; neutral reference, microstate unreviewed'}

def main():
    if meeko.__version__!='0.8.0': raise RuntimeError('expected Meeko 0.8.0')
    if not CLI.is_file(): raise RuntimeError('Meeko receptor CLI missing')
    preps={r['PDB_ID']:r for r in rows(PREP)}
    gold=rows(GOLD)
    for p in preps.values():
        for f,h in [('working_receptor_path','working_receptor_SHA256'),('native_ligand_SDF','native_ligand_SHA256'),('atom_mapping_CSV','atom_mapping_SHA256')]:
            verify(ROOT/p[f],p[h])
    natives={pid:native_prepare(p) for pid,p in preps.items()}
    result=[]
    for g in gold:
        pid=g['PDB_ID']; variant=g['variant']; p=preps[pid]; native=natives[pid]
        probe=receptor_probe(g)
        outpath=''; outhash=''
        if probe['status']=='preliminary_sensitivity_only':
            if pid!='2QMJ' or variant!='sitewaters_noNAG':
                raise RuntimeError('unexpectedly converted unreviewed receptor variant')
            out=REC_OUT/'2QMJ_sitewaters_noNAG_preliminary_sensitivity.pdbqt'
            safe_write(out,probe['data']); outpath=rel(out); outhash=sha(out)
        if pid=='2QMJ':
            blockers='NAG linked glycan lacks validated Meeko template; NAG A2005 is near site; catalytic Asp443/Asp542/His600 protonation provisional; HOH4627 H/native ligand short contact; missing Gln837/terminal issues; no redocking.'
            kept='NAG absent only in sensitivity variant' if variant=='sitewaters_noNAG' else 'NAG retained; template unsupported'
            ions='none'; waters=g['waters_retained']
        else:
            blockers='Peptide-linked PCA A1 lacks validated Meeko template; CA/CL template support alone is insufficient; catalytic Asp197/Glu233/Asp300/His201 protonation provisional; ARE crystal-mate contact; no redocking.'
            kept='PCA A1 retained; template unsupported; remote NAG already removed in GOLD variant'
            ions='CA A498; CL A500'; waters=g['waters_retained']
        row={k:'' for k in FIELDS}
        row.update({'PDB_ID':pid,'variant':variant,'input_receptor_path':g['output_path'],
            'input_SHA256':g['output_SHA256'],'Meeko_version':meeko.__version__,
            'receptor_conversion_status':probe['status'],'receptor_pdbqt_path':outpath,
            'receptor_pdbqt_SHA256':outhash,'input_heavy_atoms':probe.get('input_heavy',''),
            'pdbqt_heavy_atoms':probe.get('output_heavy',''),'pdbqt_atoms':probe.get('atoms',''),
            'coordinate_check':probe.get('coordinate',''), 'retained_waters':waters,
            'retained_ions':ions,'retained_glycans_or_modified_residues':kept,
            'native_ligand':p['reference_ligand_id'],'native_SDF_SHA256':p['native_ligand_SHA256'],
            'native_PDBQT_path':native['path'],'native_PDBQT_SHA256':native['sha'],
            'native_heavy_atoms':native['heavy'],'native_torsdof':native['tors'],
            'native_roundtrip_status':native['check'],'native_atom_map_path':native['map'],
            'site_center_xyz_A':';'.join(p[f'site_center_{axis}'] for axis in 'xyz'),
            'proposed_box_xyz_A':';'.join(map(str,BOX[pid])),'redocking_ready':'no',
            'blockers':blockers+' '+probe.get('error','')})
        result.append(row)
    import io
    buf=io.StringIO(newline=''); writer=csv.DictWriter(buf,fieldnames=FIELDS)
    writer.writeheader(); writer.writerows(result)
    safe_write(OUT,buf.getvalue().encode())
    for r in result: print(r['PDB_ID'],r['variant'],r['receptor_conversion_status'],r['receptor_pdbqt_path'])

if __name__=='__main__': main()

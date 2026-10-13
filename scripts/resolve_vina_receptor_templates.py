"""Validate RCSB PCA chemistry, build its Meeko template and audit 1XD0 PDBQTs.

No receptor/source file is edited. No docking or score calculation is performed.
Meeko 0.8.0 and the project's virtual environment are required.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import subprocess
import tempfile
import urllib.request
from collections import Counter
from pathlib import Path

import gemmi
import meeko
import numpy as np
from meeko import PDBQTReceptor, ResidueChemTemplates
from meeko.chemtempgen import ChemicalComponent, export_chem_templates_to_json
from rdkit import Chem, rdBase

ROOT = Path(__file__).resolve().parents[1]
GOLD = ROOT / 'data/metadata/gold_receptor_preparation_results.csv'
OLD_VINA = ROOT / 'data/metadata/vina_receptor_preparation_results.csv'
ORIG = ROOT / 'structures/proteins/original/1XD0.pdb'
CLI = ROOT / '.venv/Scripts/mk_prepare_receptor.exe'
VINA = Path(r'E:\docking tools\vina_1.2.7_win.exe')
TEMPLATE_DIR = ROOT / 'docking/configs/vina_templates'
OUTPUT_DIR = ROOT / 'structures/proteins/vina_prepared'
OUTPUT_CSV = ROOT / 'data/metadata/vina_receptor_blocker_resolution_results.csv'
CCD_URL = 'https://files.rcsb.org/ligands/download/PCA.cif'
CCD_SHA = '2c2aa7b8ee9b3acc6ad5ec66bb88f4701c7260f3ec60681799a1b8ec68366238'
FIELDS = ['PDB_ID','variant','source_PDB_path','source_PDB_SHA256','CCD_URL','CCD_path','CCD_SHA256',
          'custom_template_path','custom_template_SHA256','Meeko_version','RDKit_version','Gemmi_version','Vina_version',
          'settings','linkage_check','stereochemistry_check','hydrogens_status','waters_retained',
          'glycan_or_PCA_status','ions_status','PDBQT_path','PDBQT_SHA256','source_heavy_atoms',
          'PDBQT_heavy_atoms','PDBQT_total_atoms','identity_coordinate_check','charge_check',
          'water_ion_check','ligand_contact_check','template_status','chemical_preparation_status',
          'controlled_redocking_eligible','unresolved_issues']

def table(path):
    with path.open(newline='',encoding='utf-8') as f: return list(csv.DictReader(f))

def sha(data): return hashlib.sha256(data).hexdigest()

def verify(path, digest):
    if not path.is_file() or sha(path.read_bytes())!=digest:
        raise RuntimeError(f'missing or changed source: {path}')

def save_new(path, data):
    if path.exists():
        if path.read_bytes()!=data: raise RuntimeError(f'will not overwrite differing file: {path}')
    else:
        path.parent.mkdir(parents=True,exist_ok=True); path.write_bytes(data)

def relative(path): return path.relative_to(ROOT).as_posix()

def deposited_pca_links():
    lines=ORIG.read_text().splitlines()
    link=[l for l in lines if l.startswith('LINK') and 'PCA A   1' in l and 'TYR A   2' in l]
    if len(link)!=1 or ' C   PCA A   1' not in link[0] or ' N   TYR A   2' not in link[0]:
        raise RuntimeError('deposited PCA C--Tyr2 N peptide LINK missing or ambiguous')
    pdb_atoms={}
    for l in lines:
        if l.startswith(('ATOM  ','HETATM')) and l[17:20]=='PCA' and l[21]=='A' and l[22:26].strip()=='1':
            pdb_atoms[l[12:16].strip()]=np.array([float(l[i:i+8]) for i in (30,38,46)])
    tyr_n=next(np.array([float(l[i:i+8]) for i in (30,38,46)]) for l in lines if l.startswith('ATOM  ') and l[17:20]=='TYR' and l[21]=='A' and l[22:26].strip()=='2' and l[12:16].strip()=='N')
    dist=float(np.linalg.norm(pdb_atoms['C']-tyr_n))
    if not 1.2<dist<1.5 or set(pdb_atoms)!={'N','CA','CB','CG','CD','OE','C','O'}:
        raise RuntimeError(f'PCA deposited heavy atoms or C--N geometry invalid ({dist:.3f} A)')
    return dist,pdb_atoms

def build_pca_template(ccd, pdb_atoms):
    block=gemmi.cif.read(str(ccd)).find_block('PCA')
    if block is None: raise RuntimeError('PCA CCD block absent')
    ideal=block.find('_chem_comp_atom.',['atom_id','pdbx_model_Cartn_x_ideal','pdbx_model_Cartn_y_ideal','pdbx_model_Cartn_z_ideal'])
    ideal_coords={row[0]:(float(row[1]),float(row[2]),float(row[3])) for row in ideal}
    stereo=dict((row[0],row[1]) for row in block.find('_chem_comp_atom.',['atom_id','pdbx_stereo_config']))
    if stereo.get('CA')!='S': raise RuntimeError('official PCA CA stereo is not S')
    comp=ChemicalComponent.from_cif(str(ccd),'PCA')
    if comp is None: raise RuntimeError('Meeko CCD parser failed PCA')
    mol=comp.rdkit_mol
    conf=Chem.Conformer(mol.GetNumAtoms()); conf.Set3D(True)
    for atom in mol.GetAtoms(): conf.SetAtomPosition(atom.GetIdx(),ideal_coords[atom.GetProp('atom_id')])
    mol.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(mol)
    stereo_ideal=dict((mol.GetAtomWithIdx(i).GetProp('atom_id'), cip) for i,cip in Chem.FindMolChiralCenters(mol,includeUnassigned=True))
    if stereo_ideal.get('CA')!='S': raise RuntimeError('CCD ideal coordinates conflict with declared S chirality')
    # CCD OXT/HXT leave on peptide condensation; PDB LINK confirms C--Tyr2 N.
    comp.make_embedded(allowed_smarts='[#6]',leaving_names={'OXT'}).make_pretty_smiles().make_link_labels_from_names({'C':'C-term'})
    comp.ResidueTemplate_check()
    parser=Chem.SmilesParserParams(); parser.removeHs=False
    tmpl_mol=Chem.MolFromSmiles(comp.smiles_exh,parser)
    if tmpl_mol is None: raise RuntimeError('template SMILES invalid')
    centers=dict(Chem.FindMolChiralCenters(tmpl_mol,includeUnassigned=True))
    ca_index=comp.atom_name.index('CA')
    if centers.get(ca_index)!='S' or comp.link_labels!={str(comp.atom_name.index('C')):'C-term'}:
        raise RuntimeError('PCA template stereo or C-terminal linker invalid')
    if set(n for n in comp.atom_name if not n.startswith('H')) != set(pdb_atoms):
        # H prefix is not sufficient to identify element, so check elements below.
        heavy={comp.atom_name[a.GetIdx()] for a in tmpl_mol.GetAtoms() if a.GetAtomicNum()>1}
        if heavy!=set(pdb_atoms): raise RuntimeError('PCA template heavy atom names differ from deposited residue')
    template_bytes=export_chem_templates_to_json([comp]).encode('utf-8')
    return template_bytes,comp.smiles_exh

def check_deposited_pca_stereo(template_bytes, source):
    item=json.loads(template_bytes)['residue_templates']['PCA']
    parser=Chem.SmilesParserParams(); parser.removeHs=False
    mol=Chem.MolFromSmiles(item['smiles'],parser)
    names=item['atom_name']
    xyz={l[12:16].strip():[float(l[i:i+8]) for i in (30,38,46)]
         for l in source.read_text().splitlines()
         if l.startswith(('ATOM  ','HETATM')) and l[17:20]=='PCA' and l[21]=='A' and l[22:26].strip()=='1'}
    if set(names)!=set(xyz): raise RuntimeError('deposited PCA atom list differs from template')
    conf=Chem.Conformer(mol.GetNumAtoms()); conf.Set3D(True)
    for i,name in enumerate(names): conf.SetAtomPosition(i,xyz[name])
    mol.AddConformer(conf); Chem.AssignStereochemistryFrom3D(mol)
    centers=dict(Chem.FindMolChiralCenters(mol,includeUnassigned=True))
    if centers.get(names.index('CA'))!='S': raise RuntimeError('deposited PCA CA chirality conflicts with official CCD')
    return 'official CCD CA S; ideal 3D S; deposited 3D S; custom template S'

def pdb_heavy_atoms(path):
    result={}
    for l in path.read_text().splitlines():
        if not l.startswith(('ATOM  ','HETATM')) or l[76:78].strip().upper() in ('H','D'): continue
        key=(l[21],int(l[22:26]),l[17:20].strip(),l[12:16].strip())
        if key in result: raise RuntimeError(f'duplicate heavy atom {key}')
        result[key]=np.array([float(l[i:i+8]) for i in (30,38,46)])
    return result

def audit_pdbqt(source, data):
    with tempfile.TemporaryDirectory() as temp:
        probe=Path(temp)/'receptor.pdbqt'; probe.write_bytes(data)
        atoms=PDBQTReceptor.from_pdbqt_filename(str(probe)).atoms()
    if not len(atoms) or not np.isfinite(atoms['xyz']).all() or not np.isfinite(atoms['partial_charges']).all():
        raise RuntimeError('PDBQT readback invalid coordinates or charges')
    source_atoms=pdb_heavy_atoms(source)
    output={}
    for atom in atoms:
        if atom['atom_type']=='HD': continue
        name=atom['name']
        if atom['resname']=='HOH' and name=='OW': name='O'
        key=(atom['chain'],int(atom['resid']),atom['resname'],name)
        if key in output: raise RuntimeError(f'duplicate PDBQT heavy atom {key}')
        output[key]=atom
    if set(source_atoms)!=set(output):
        raise RuntimeError(f'PDBQT heavy identities differ: lost {set(source_atoms)-set(output)}, gained {set(output)-set(source_atoms)}')
    # Meeko may exchange equivalent terminal carboxylate O/OXT atom names while
    # preserving both deposited coordinates. Match these two by position.
    terminal=('A',496,'LEU')
    ko=terminal+('O',); kx=terminal+('OXT',)
    terminal_swap=False
    if ko in source_atoms and kx in source_atoms:
        direct=sum(float(np.linalg.norm(source_atoms[k]-output[k]['xyz'])) for k in (ko,kx))
        swapped=float(np.linalg.norm(source_atoms[ko]-output[kx]['xyz'])+np.linalg.norm(source_atoms[kx]-output[ko]['xyz']))
        if swapped<direct:
            output[ko],output[kx]=output[kx],output[ko]
            terminal_swap=True
    shifts=sorted(((float(np.linalg.norm(source_atoms[k]-output[k]['xyz'])),k) for k in source_atoms),reverse=True)
    shift=shifts[0][0]
    if shift>0.0011: raise RuntimeError(f'heavy coordinate shift {shift:.4f} A; worst {shifts[:8]}')
    pca=[a for a in atoms if a['chain']=='A' and a['resid']==1 and a['resname']=='PCA']
    pca_heavy=[a for a in pca if a['atom_type']!='HD']
    if len(pca_heavy)!=8: raise RuntimeError('PCA heavy atom count mismatch')
    pca_charge=float(sum(a['partial_charges'] for a in pca))
    if abs(pca_charge)>0.08: raise RuntimeError(f'PCA fragment partial charge unexpected {pca_charge:.3f}')
    ions={a['resname']:a for a in atoms if a['resname'] in ('CA','CL')}
    if set(ions)!={'CA','CL'} or ions['CA']['atom_type']!='Ca' or ions['CL']['atom_type']!='Cl':
        raise RuntimeError('calcium/chloride missing or incorrectly atom-typed')
    if abs(float(ions['CA']['partial_charges'])-2)>0.001 or abs(float(ions['CL']['partial_charges'])+1)>0.001:
        raise RuntimeError('calcium/chloride charge mismatch')
    waters={int(a['resid']) for a in atoms if a['resname']=='HOH' and a['name']=='OW'}
    if 536 not in waters or 546 not in waters: raise RuntimeError('Ca coordination waters omitted')
    # The six vs 142 water variants are handled identically; count is checked by caller.
    return {'heavy':len(output),'atoms':len(atoms),'shift':shift,'terminal_swap':terminal_swap,'pca_charge':pca_charge,
            'waters':len(waters),'types':';'.join(sorted(set(atoms['atom_type']))),
            'ion_summary':f"CA A498: Ca +{float(ions['CA']['partial_charges']):.3f}; CL A500: Cl {float(ions['CL']['partial_charges']):.3f}",
            'water_ids':waters}

def min_ligand_distances(source, ligand_sdf):
    mol=list(Chem.SDMolSupplier(str(ligand_sdf),removeHs=False))[0]
    if mol is None: raise RuntimeError('native ligand SDF invalid')
    ligand=np.array([list(mol.GetConformer().GetAtomPosition(a.GetIdx())) for a in mol.GetAtoms() if a.GetAtomicNum()>1])
    mins={'protein_heavy':float('inf'),'water_O':float('inf'),'polar_H':float('inf')}
    for l in source.read_text().splitlines():
        if not l.startswith(('ATOM  ','HETATM')): continue
        pos=np.array([float(l[i:i+8]) for i in (30,38,46)])
        d=float(np.linalg.norm(ligand-pos,axis=1).min())
        element=l[76:78].strip().upper()
        if l[17:20]=='HOH' and element=='O': mins['water_O']=min(mins['water_O'],d)
        elif element=='H' and l[12:16].strip().startswith('H'): mins['polar_H']=min(mins['polar_H'],d)
        elif element not in ('H','D') and l[17:20] not in ('CA ','CL '): mins['protein_heavy']=min(mins['protein_heavy'],d)
    return mins

def main():
    if meeko.__version__!='0.8.0': raise RuntimeError('Meeko 0.8.0 required')
    if not CLI.is_file(): raise RuntimeError('project Meeko CLI missing')
    if not VINA.is_file(): raise RuntimeError('AutoDock Vina 1.2.7 executable missing')
    version=subprocess.run([str(VINA),'--version'],capture_output=True,text=True,timeout=15)
    if version.returncode or '1.2.7' not in version.stdout+version.stderr:
        raise RuntimeError('AutoDock Vina version check failed')
    audit={r['PDB_ID']:r for r in table(ROOT/'data/metadata/receptor_structure_audit.csv')}
    verify(ORIG,audit['1XD0']['SHA256'])
    for r in table(OLD_VINA):
        verify(ROOT/r['input_receptor_path'],r['input_SHA256'])
        if r['receptor_pdbqt_path']: verify(ROOT/r['receptor_pdbqt_path'],r['receptor_pdbqt_SHA256'])
    gold=table(GOLD)
    for g in gold: verify(ROOT/g['output_path'],g['output_SHA256'])
    with urllib.request.urlopen(CCD_URL,timeout=30) as response: ccd_data=response.read()
    if sha(ccd_data)!=CCD_SHA: raise RuntimeError('official RCSB PCA CCD checksum changed; review before use')
    ccd=TEMPLATE_DIR/'PCA_RCSB_CCD_2024-09-27.cif'; save_new(ccd,ccd_data)
    linkdist,pca_atoms=deposited_pca_links()
    tbytes,smiles=build_pca_template(ccd,pca_atoms)
    template=TEMPLATE_DIR/'PCA_C_terminal_Meeko_0.8.0.json'; save_new(template,tbytes)
    residue_templates=ResidueChemTemplates.create_from_defaults(); residue_templates.add_json_file(str(template))
    if 'PCA' not in residue_templates.ambiguous: raise RuntimeError('custom PCA not recognized by Meeko')
    rows=[]
    for g in gold:
        pid=g['PDB_ID']; variant=g['variant']; source=ROOT/g['output_path']
        row={k:'' for k in FIELDS}
        row.update({'PDB_ID':pid,'variant':variant,'source_PDB_path':g['output_path'],
            'source_PDB_SHA256':g['output_SHA256'],'Meeko_version':meeko.__version__,
            'RDKit_version':rdBase.rdkitVersion,'Gemmi_version':gemmi.__version__,
            'Vina_version':'1.2.7 (executable version check only; no docking/scoring)',
            'settings':'mk_prepare_receptor --read_pdb --add_templates PCA_C_terminal_Meeko_0.8.0.json --write_pdbqt; default Gasteiger; no deletion or residue substitution',
            'controlled_redocking_eligible':'no'})
        if pid=='2QMJ':
            row.update({'glycan_or_PCA_status':'all four NAG are in protein-linked glycans; A2005 links Asn209 and is near ligand',
                'waters_retained':g['waters_retained'],'hydrogens_status':'none added in NAG-retaining variants; pH7 provisional only in no-NAG model',
                'template_status':'blocked_no_validated_NAG_linkage_template',
                'chemical_preparation_status':'blocked; earlier no-NAG PDBQT remains sensitivity only',
                'unresolved_issues':'NAG/Asn-linked protonation and custom glycan linker/padder charges/stereo; catalytic Asp443/Asp542/His600; HOH4627 H close contact; Gln837 missing. No new 2QMJ PDBQT.'})
            if variant=='sitewaters_noNAG':
                row['template_status']='no NAG template needed for unchanged NAG-free sensitivity model'
                row['chemical_preparation_status']='unchanged earlier preliminary sensitivity PDBQT; not a glycan-retaining model'
        else:
            deposited_stereo=check_deposited_pca_stereo(tbytes,source)
            row.update({'CCD_URL':CCD_URL,'CCD_path':relative(ccd),'CCD_SHA256':sha(ccd_data),
                'custom_template_path':relative(template),'custom_template_SHA256':sha(tbytes),
                'linkage_check':f'official PCA CCD; deposited PCA A1 C--Tyr2 N {linkdist:.3f} A; C-term linker',
                'stereochemistry_check':deposited_stereo,
                'hydrogens_status':'existing OpenMM pH7 H; Meeko merges nonpolar H and retains polar H',
                'glycan_or_PCA_status':'PCA A1 retained with 8 deposited heavy atoms; remote NAG A499 omitted in prior GOLD variant',
                'waters_retained':g['waters_retained'],'ions_status':'CA A498 and CL A500 retained; HOH536/546 coordinate Ca',
                'template_status':'PCA_template_passes_graph_stereo_linker_and_Meeko_load'})
            with tempfile.TemporaryDirectory() as temp:
                probe=Path(temp)/'output.pdbqt'
                run=subprocess.run([str(CLI),'--read_pdb',str(source),'--add_templates',str(template),'--write_pdbqt',str(probe)],capture_output=True,text=True,timeout=240)
                if run.returncode or not probe.is_file():
                    raise RuntimeError('Meeko PCA-template conversion failed: '+(run.stderr+' '+run.stdout)[-2000:])
                data=probe.read_bytes()
            qc=audit_pdbqt(source,data)
            if qc['waters']!=int(g['waters_retained']): raise RuntimeError('water count differs')
            if any(t not in {'A','C','N','NA','OA','SA','S','HD','Ca','Cl'} for t in qc['types'].split(';')):
                raise RuntimeError('unexpected PDBQT atom type')
            native=ROOT/g['native_reference_path']; verify(native,g['native_reference_SHA256'])
            distances=min_ligand_distances(source,native)
            target=OUTPUT_DIR/f'1XD0_{variant}_PCA_template_Meeko0.8.0.pdbqt'
            save_new(target,data)
            row.update({'PDBQT_path':relative(target),'PDBQT_SHA256':sha(data),
                'source_heavy_atoms':len(pdb_heavy_atoms(source)),'PDBQT_heavy_atoms':qc['heavy'],
                'PDBQT_total_atoms':qc['atoms'],
                'identity_coordinate_check':f'all source heavy atom positions preserved; max shift {qc["shift"]:.4f} A; Leu496 O/OXT names swapped by Meeko={qc["terminal_swap"]}',
                'charge_check':f'finite; PCA fragment sum {qc["pca_charge"]:+.3f}; Gasteiger/template charge model, site suitability untested',
                'water_ion_check':qc['ion_summary']+f'; {qc["waters"]} waters including HOH536/546; atom types {qc["types"]}',
                'ligand_contact_check':'; '.join(f'{k} min {v:.3f} A' for k,v in distances.items()),
                'chemical_preparation_status':'PDBQT structurally/chemically prepared for controlled review; no Vina docking validation',
                'unresolved_issues':'Catalytic Asp197/Glu233/Asp300/His201 microstates provisional; ion/water scoring and hydrogen orientations need review; ARE distal atoms contact crystal mate; no redocking yet.'})
        rows.append(row)
    buf=io.StringIO(newline=''); writer=csv.DictWriter(buf,fieldnames=FIELDS)
    writer.writeheader(); writer.writerows(rows)
    save_new(OUTPUT_CSV,buf.getvalue().encode('utf-8'))
    for r in rows: print(r['PDB_ID'],r['variant'],r['chemical_preparation_status'],r['PDBQT_path'])

if __name__=='__main__': main()

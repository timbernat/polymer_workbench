'''Draw chemical fragments of an mPD-tMC membrane'''

import random
from pathlib import Path
from itertools import combinations

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdMolTransforms

from polymerist.rdutils import rdkdraw
from polymerist.polymers.monomers import specification
from polymerist.rdutils.bonding import portlib

rdkdraw.set_rdkdraw_size(400, 3/2)
rdkdraw.disable_substruct_highlights()


OUTDIR = Path('membrane_fragments')
OUTDIR.mkdir(exist_ok=True)

fragments = {
    'mpd' : [
        'c1ccc(N)cc1N',
        'c1ccc(N[1*])cc1N',
        'c1ccc(N[1*])cc1N[2*]',
    ],
    'tmc' : [
        'c1c(C(=O)Cl)cc(C(=O)Cl)cc1(C(=O)Cl)',
        'c1c(C(=O)[3*])cc(C(=O)Cl)cc1(C(=O)Cl)',
        'c1c(C(=O)[3*])cc(C(=O)[4*])cc1(C(=O)Cl)',
        'c1c(C(=O)[3*])cc(C(=O)[4*])cc1(C(=O)[5*])',
    ],
}
# DIHEDRAL_QUERY = Chem.MolFromSmarts('[#1,#17]-[#7,#8]-[#6]=[#6]')
DIHEDRAL_QUERY = Chem.MolFromSmarts('*~*~*~*')

for mol_class, smis in fragments.items():
    MOLDIR = OUTDIR / mol_class
    MOLDIR.mkdir(exist_ok=True)

    for smi in smis:
        exp_smi = specification.expanded_SMILES(smi, assign_map_nums=False)
        mol = Chem.MolFromSmiles(exp_smi, sanitize=False)
        Chem.SanitizeMol(mol, sanitizeOps=specification.SANITIZE_AS_KEKULE)
        valence = portlib.get_num_ports(mol)
        molname = f'{mol_class}_{valence}-site'

        img_filepath = MOLDIR / f'{molname}.png'
        Draw.MolToImageFile(mol, img_filepath)

        if valence == 0:
            errcode = AllChem.EmbedMolecule(mol, randomSeed=0xf00d, useRandomCoords=True, ETversion=2)
            assert(errcode == 0)
            AllChem.MMFFOptimizeMolecule(mol)
            
            ETKDG_filepath = MOLDIR / f'{molname}_ETKDG.png'
            Draw.MolToImageFile(mol, ETKDG_filepath)

            conf = mol.GetConformer(0)
            for dihed_ids in mol.GetSubstructMatches(DIHEDRAL_QUERY):
                print(dihed_ids)
                try:
                    rdMolTransforms.SetDihedralDeg(conf, *dihed_ids, random.uniform(0, 360))
                except ValueError:
                    print('FAILED')
                    pass
                
            rand_torsion_filepath = MOLDIR / f'{molname}_rand_torsion.png'
            Draw.MolToImageFile(mol, rand_torsion_filepath)

            pdb_filepath = MOLDIR / f'{molname}.pdb'
            Chem.MolToPDBFile(mol, pdb_filepath, flavor=32)
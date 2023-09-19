import re
from pathlib import Path
from rdkit import Chem

pdb_dir = Path('polymer_examples/compatible_pdbs/simple_polymers/') # adjust this as needed
out_dir = Path.cwd() / f'{pdb_dir.name}_cleaned_resnames'
out_dir.mkdir(exist_ok=True)


def adjust_PDB_resname_pos(pdb_block : str, delimiter : str='\n', res_idx : int=17-1) -> str:
    '''For correcting spacing issues with PDB residue names'''
    lines = []
    for line in pdb_block.split(delimiter):
        if re.match('HETATM|ATOM', line) and (line[res_idx] != ' '):
            lines.append(line[:res_idx] + ' ' + line[res_idx:])
        else:
            lines.append(line)

    return delimiter.join(lines)

for path in pdb_dir.iterdir():
    if (path.suffix == '.pdb') and not (re.search(r'\d+', path.stem)): # only read non-squared PDBs
        print(path.stem)
        with path.open('r') as file:
            pdb_block = file.read()

        # initial file load to check for successful RDKit load
        rdmol = Chem.MolFromPDBBlock(pdb_block)
        if rdmol.GetNumAtoms() != 0:
            print('\tWorked the first time! Continuing...')
            with (out_dir / path.name).open('w') as outfile:
                outfile.write(pdb_block)
            continue
        
        print('\tINITIALLY FAILED, attempting PDB Residue Name shift')
        pdb_block = adjust_PDB_resname_pos(pdb_block)

        # attempt shift to see if that fixes things
        rdmol = Chem.MolFromPDBBlock(pdb_block)
        if rdmol.GetNumAtoms() != 0:
            print('\tResidue name shift fixed it! Continuing...')
            # with (out_dir / f'{path.stem}_shifted{path.suffix}').open('w') as outfile:
            with (out_dir / path.name).open('w') as outfile:
                outfile.write(pdb_block)
            continue

        print('\tFAILED AGAIN, issue is something deeper')
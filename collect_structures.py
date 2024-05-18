'''For extracting desired PDB and monomer JSON files from a larger cleaned polymer directory'''

from shutil import copyfile
from pathlib import Path

from polymerist.genutils.fileutils.pathutils import assemble_path
from polymerist.polymers.monomers import MonomerGroup


# CONFIGURE PATHS
pdb_sub = 'simple_polymers'
pdb_src_dir  = Path(f'pdb_test_cleaned/pdbs/{pdb_sub}')
mono_src_dir = Path(f'pdb_test_cleaned/monos/{pdb_sub}')

working_dir = Path('polymer_structures')
working_dir.mkdir(exist_ok=True)

pdb_out_dir = working_dir / 'pdb'
pdb_out_dir.mkdir(exist_ok=True)

mono_out_dir = working_dir / 'monomers'
mono_out_dir.mkdir(exist_ok=True)

# DEFINE TERMINAL GROUPS ORIENTATIONS
term_orients = { 
    'peg_modified' : {
        'peg_TERM2' : 'head',
        'peg_TERM3' : 'tail',
    },
    'paam_modified' : {
        'paam_TERM2' : 'head',
        'paam_TERM3' : 'tail',
    },
    'pnipam_modified' : {
        'pnipam_TERM2' : 'head',
        'pnipam_TERM3' : 'tail',
    },
}

for mol_name, term_orient in term_orients.items():
    # COPY PDB FILE
    pdb_src_path = assemble_path(pdb_src_dir, mol_name, extension='pdb')
    assert(pdb_src_path.exists())
    pdb_path = assemble_path(pdb_out_dir, mol_name, extension='pdb')
    copyfile(pdb_src_path, pdb_path)

    # COPY MONOMER FILE WITH ASSIGNED TERMINAL GROUPS
    mono_src_path = assemble_path(mono_src_dir, mol_name, extension='json')
    assert(mono_src_path.exists())
    mono_path = assemble_path(mono_out_dir, mol_name, extension='json')

    monomers = MonomerGroup.from_file(mono_src_path)
    monomers.term_orient = term_orients[mol_name]
    monomers.to_file(mono_path)
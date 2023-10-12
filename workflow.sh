#!/bin/bash
# ===================================================
# Shell-based wrapper for polymer simulation workflow
# ===================================================


# parameter definitions
# ===================================================
struct_src='polymer_structures'
pdb_src="$struct_src/pdb"
mono_src="$struct_src/monomers"

mol_names=(
    'peg_modified'
    'paam_modified'
    'pnipam_modified'
)

N=150
rct_charge_method='Espaloma-AM1-BCC'
charge_methods=(
    'Espaloma-AM1-BCC'
    'RCT'
)

box_x=5.0
box_y=5.0
box_z=5.0
box_unit='nanometer'
solvent='water_TIP3P'
forcefield='openff-2.0.0'

root_dir='workflow_test'
# TODO : add optional deletion if existing
mkdir -p $root_dir


# actual code
# ===================================================
for mol_name in "${mol_names[@]}"; do
    mol_dir=$root_dir/$mol_name
    python -m components.assign_chem  -wdir $mol_dir -pdb $pdb_src/$mol_name.pdb -mono $mono_src/$mol_name.json -n $mol_name
    python -m components.rct_protocol -wdir $mol_dir -mono $mono_src/$mol_name.json -n $mol_name -cmet $rct_charge_method -N $N
    
    for charge_method in "${charge_methods[@]}"; do
        charge_dir=$mol_dir/$charge_method
        # TODO : also add RCT charge of reduction
        python -m components.charge_mol -wdir $charge_dir -sdf $mol_dir/$mol_name.sdf -cmet $charge_method -lc $mol_dir/${mol_name}_residue_charges.json

        # TODO : add anneal
        # TODO : add solvation
        # TODO : add simulation schedule
        # TODO : add analysis
    done
done


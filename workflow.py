'''Functions for execution automated polymer simulation workflows'''

# LOGGING
from polysaccharide2.genutils.logutils.IOHandlers import LOG_FORMATTER

import logging
logging.basicConfig(
    level=logging.INFO,
    format =LOG_FORMATTER._fmt,
    datefmt=LOG_FORMATTER.datefmt,
    # force=True
)
LOGGER = logging.Logger(__name__)

# CORE IMPORTS
from typing import Optional, Union

from pathlib import Path
from copy import deepcopy

from openff.toolkit import ForceField
from openmm.unit import angstrom
from openmm.app import Topology

# CUSTOM IMPORTS
from polysaccharide2.genutils.unitutils import openmm_to_openff
from polysaccharide2.genutils.fileutils.pathutils import assemble_path
from polysaccharide2.genutils.decorators.functional import allow_string_paths

from polysaccharide2.openmmtools import execution
from polysaccharide2.openmmtools.parameters import SimulationParameters

from polysaccharide2.openfftools import TKREGS, topology
from polysaccharide2.openfftools.pcharge import MolCharger
from polysaccharide2.openfftools.omminter import openff_topology_to_openmm
from polysaccharide2.openfftools.solvation.boxvectors import VectorQuantity, BoxVectorsQuantity

from polysaccharide2.monomers.repr import MonomerGroup
from polysaccharide2.polymers import estimation, building
from polysaccharide2.polymers.exceptions import MorphologyError

from polysaccharide2.residues.partition import partition
from polysaccharide2.residues.rescharge.rctypes import ChargesByResidue
from polysaccharide2.residues.rescharge.calculation import compute_residue_charges


# PARAMETERIZATION FUNCTIONS
@allow_string_paths
def rct_protocol(working_dir : Path, mol_name : str, monogroup : MonomerGroup, term_group_orient : dict[str, str], N : int, charger : MolCharger, delete_pdb : bool=True, save_sdf : bool=False) -> ChargesByResidue:
    '''
    Generates library charges for a monomer group given terminal group head-tail orientations and a maximum chain length
    If delete_pdb=True, will remove the working pdb path after charges are generated
    as currently implemented, only supports monomer groups which constitute a linear homopolymer
    '''
    if not monogroup.is_linear:
        raise MorphologyError('RCT currently only supports linear homopolymers')
    
    DOP = estimation.estimate_DOP_lower(monogroup, max_chain_len=N)
    chain = building.build_linear_polymer(monogroup, DOP, term_orient=term_group_orient)

    pdb_path = assemble_path(working_dir, mol_name, extension='pdb')
    building.mbmol_to_openmm_pdb(pdb_path, chain) # output PDB file to disc 
    offtop = Topology.from_pdb(pdb_path, _custom_substructures=monogroup.monomers, toolkit_registry=TKREGS['The RDKit']) # load custom substructures - raises error if PDB has issues
    if delete_pdb:
        pdb_path.unlink()
        
    was_partitioned = partition(offtop) 
    assert(was_partitioned)
    offmol = topology.get_largest_offmol(offtop)
    offmol.name = mol_name

    cmol = charger.charge_molecule(offmol, in_place=False)
    if save_sdf:
        sdf_path = assemble_path(working_dir, mol_name, extension='sdf')
        topology.topology_to_sdf(sdf_path, cmol.to_topology())
    res_chgs = compute_residue_charges(cmol, monogroup)

    return res_chgs


# SIMULATION FUNCTIONS
@allow_string_paths
def vacuum_anneal(working_dir : Path, offtop : Topology, anneal_params : SimulationParameters, forcefield : Union[ForceField, str, Path],
                   box_vecs : Optional[Union[VectorQuantity, BoxVectorsQuantity]]=None, step_name : str='anneal') -> Topology:
    '''Run a short vacuum simulation of a Topology with a single molecule and reassign its conformer based on the '''
    if offtop.n_molecules != 1:
        raise ValueError(f'Can only run vacuum anneal on Topology with single molecule (Topology contains {offtop.n_molecules} molecules)')
    offmol = topology.get_largest_offmol(offtop)
    mol_name = offmol.name

    # perform conformer anneal simulation in OpenMM
    ommtop, ommsys, ommpos = openff_topology_to_openmm(offtop, forcefield=forcefield, box_vecs=box_vecs)
    anneal_schedule = {step_name : anneal_params}
    anneal_history = execution.run_simulation_schedule(working_dir, anneal_schedule, ommtop, ommsys, ommpos, return_history=True)
    new_conformer = execution.get_simulation_positions(anneal_history[step_name]['simulation'])

    # update conformer with new positions
    new_mol = deepcopy(offmol)
    LOGGER.info(f'Transferring coordinates from anneal to {mol_name} conformer')
    new_mol.conformers[0] = openmm_to_openff(new_conformer.in_units_of(angstrom)) # convert to correct units in the OpenFF format

    return new_mol.to_topology()
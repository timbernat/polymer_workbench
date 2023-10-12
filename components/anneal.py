'''Vacuum anneal a molecule to generate a unique starting configuration'''

import argparse
from pathlib import Path
import logging
LOGGER = logging.Logger(__name__)

from typing import Optional, Union

import numpy as np
from copy import deepcopy

from openff.toolkit import ForceField, Topology
import openmm.unit
from openmm.unit import angstrom, nanometer, Unit

from polysaccharide2.genutils.decorators.functional import allow_string_paths
from polysaccharide2.genutils.fileutils.pathutils import assemble_path
from polysaccharide2.genutils.logutils.IOHandlers import MSFHandlerFlex
from polysaccharide2.genutils.unitutils import openmm_to_openff

from polysaccharide2.openmmtools.parameters import SimulationParameters
from polysaccharide2.openmmtools import execution
from polysaccharide2.openmmtools.parameters import SimulationParameters

from polysaccharide2.openfftools import topology
from polysaccharide2.openfftools.omminter import openff_topology_to_openmm
from polysaccharide2.openfftools.solvation import boxvectors


# accept args from CLI
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-wdir', '--working_directory'  , help='Directory into which files generated should be saved'  , type=Path)
parser.add_argument('-sdf' , '--sdf_path'           , help='Path to SDF file from which to read parameterized molecule structure', type=Path)
parser.add_argument('-bd'  , '--box_dimensions'     , help='XYZ dimensions of the desired periodic box for the system', type=list, nargs=3)
parser.add_argument('-bdu' , '--box_dimension_units', help='Unit to use when assigning box dimensions (default nanometer)', type=str, default='nanometer')
parser.add_argument('-sp'  , '--sim_params_path'    , help='Path to the file storing the simulation parameters for the anneal', type=str)
parser.add_argument('-ff'  , '--forcefield'         , help='Name of the ForceField XML to use for MD', type=str)
parser.add_argument('-a'   , '--affix'              , help='Additional text to attach to the end of the molecule name provided when naming ', type=str)


# post-process args as needed
args = parser.parse_args()
args.working_directory.mkdir(exist_ok=True)
assert(args.sdf_path.suffix == '.sdf')
assert(args.sim_params_paths.suffix == '.json')

box_dim_unit = getattr(openmm.unit, args.box_dimension_unit)
box_dims = np.array(args.box_dimensions) * box_dim_unit
box_vecs = boxvectors.box_vectors_flexible(box_dims)


# utility functions
@allow_string_paths
def vacuum_anneal(working_dir : Path, offtop : Topology, anneal_params : SimulationParameters, forcefield : Union[ForceField, str, Path],
                box_vecs : Optional[Union[boxvectors.VectorQuantity, boxvectors.BoxVectorsQuantity]]=None, step_name : str='anneal') -> Topology:
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


# main code
if __name__ == '__main__':
    # initialize simulation
    with MSFHandlerFlex(args.working_directory, proc_name='vacuum_anneal', loggers='all') as log_handler:
        anneal_params = SimulationParameters.from_file(args.sim_params_path)
        offtop = topology.topology_from_sdf(args.sdf_path)
        offmol = topology.get_largest_offmol(offtop)
        mol_name = offmol.name

        conf_top = vacuum_anneal(args.working_directory, offtop, anneal_params, forcefield=args.forcefield, box_vecs=box_dims)
        conf_top_path = assemble_path(args.working_directory, mol_name, extension='sdf', postfix=args.affix)
        topology.topology_to_sdf(conf_top_path, conf_top)
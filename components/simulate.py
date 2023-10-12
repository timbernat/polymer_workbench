'''Add solvent to a Topology'''

import argparse
from pathlib import Path
import logging
LOGGER = logging.Logger(__name__)

import numpy as np
import openmm.unit

from polysaccharide2.genutils.logutils.IOHandlers import MSFHandlerFlex
from polysaccharide2.openfftools.solvation import boxvectors

from polysaccharide2.genutils.logutils.IOHandlers import MSFHandlerFlex
from polysaccharide2.openmmtools import parameters, execution
from polysaccharide2.openfftools import topology
from polysaccharide2.openfftools.omminter import openff_topology_to_openmm



# accept args from CLI
class ParseKwargs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split('=')
            getattr(namespace, self.dest)[key] = value

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-wdir', '--working_directory'  , help='Directory into which files generated should be saved'  , type=Path)
parser.add_argument('-sdf' , '--sdf_path'           , help='Path to SDF file from which to read parameterized molecule structure', type=Path)
parser.add_argument('-a'   , '--affix'              , help='Additional text to attach to the end of the molecule name provided when naming ', type=str)
parser.add_argument('-sp'  , '--sim_params'         , help='Key-value pairs of step names and simulation parameter paths strings which define the simulation schedule', nargs='*', action=ParseKwargs)

parser.add_argument('-bd'  , '--box_dimensions'     , help='XYZ dimensions of the desired periodic box for the system', type=list, nargs=3)
parser.add_argument('-bdu' , '--box_dimension_units', help='Unit to use when assigning box dimensions (default nanometer)', type=str, default='nanometer')
parser.add_argument('-ff'  , '--forcefield'         , help='Name of the ForceField XML to use for MD', type=str)


# post-process args as needed
args = parser.parse_args()
args.working_directory.mkdir(exist_ok=True)
assert(args.sdf_path.suffix == '.sdf')

box_dim_unit = getattr(openmm.unit, args.box_dimension_unit)
box_dims = np.array(args.box_dimensions) * box_dim_unit
box_vecs = boxvectors.box_vectors_flexible(box_dims)
schedule = {
    step_name : parameters.SimulationParameters.from_file(path_str)
        for step_name, path_str in args.sim_params.items()
}


# main code
if __name__ == '__main__':
    with MSFHandlerFlex(args.working_directory, proc_name=Path(__file__).stem, loggers='all') as log_handler:
        offtop = topology.topology_from_sdf(args.sdf_path)
        ommtop, ommsys, ommpos = openff_topology_to_openmm(offtop, forcefield=args.forcefield, box_vecs=box_dims)
        history = execution.run_simulation_schedule(args.working_directory, schedule, ommtop, ommsys, ommpos, return_history=True)
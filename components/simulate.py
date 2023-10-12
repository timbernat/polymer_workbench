'''Add solvent to a Topology'''

import warnings
warnings.catch_warnings(record=True)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

import argparse
import logging
from polysaccharide2.genutils.logutils.IOHandlers import LOG_FORMATTER
logging.basicConfig(
    level=logging.INFO,
    format =LOG_FORMATTER._fmt,
    datefmt=LOG_FORMATTER.datefmt,
    # force=True
)
LOGGER = logging.Logger(__name__)

import numpy as np
import openmm.unit
from pathlib import Path

from polysaccharide2.genutils.logutils.IOHandlers import MSFHandlerFlex
from polysaccharide2.openfftools.solvation import boxvectors

from polysaccharide2.genutils.logutils.IOHandlers import MSFHandlerFlex
from polysaccharide2.openmmtools import parameters, execution
from polysaccharide2.openfftools import topology
from polysaccharide2.openfftools.omminter import openff_topology_to_openmm


class ParseKwargs(argparse.Action):
    '''Parser action for dict-like key-value pairs'''
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split('=')
            getattr(namespace, self.dest)[key] = value

def parse_args() -> argparse.Namespace:
    '''Accept and pre-process arguments from command line'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-wdir', '--working_directory'  , help='Directory into which files generated should be saved'  , type=Path)
    parser.add_argument('-sdf' , '--sdf_path'           , help='Path to SDF file from which to read parameterized molecule structure', type=Path)
    parser.add_argument('-a'   , '--affix'              , help='Additional text to attach to the end of the molecule name provided when naming ', type=str)

    parser.add_argument('-bd'  , '--box_dimensions'     , help='XYZ dimensions of the desired periodic box for the system', nargs=3, type=float)
    parser.add_argument('-bdu' , '--box_dimension_unit' , help='Unit to use when assigning box dimensions (default nanometer)', type=str, default='nanometer')
    
    parser.add_argument('-sp'  , '--sim_param_paths'    , help='Path to the file storing the simulation parameters for the anneal', type=str, nargs='*', action=ParseKwargs)
    parser.add_argument('-ff'  , '--forcefield'         , help='Name of the ForceField XML to use for MD', type=str)

    # post-process args as needed
    args = parser.parse_args()
    args.working_directory.mkdir(exist_ok=True)
    assert(args.sdf_path.suffix == '.sdf')

    box_dim_unit = getattr(openmm.unit, args.box_dimension_unit)
    box_dims = np.array(args.box_dimensions) * box_dim_unit
    args.box_vecs = boxvectors.box_vectors_flexible(box_dims)

    args.schedule = {
        step_name : parameters.SimulationParameters.from_file(path_str)
            for step_name, path_str in args.sim_param_paths.items()
    }

    return args

def main() -> None:
    '''Define main code body'''
    args = parse_args()
    with MSFHandlerFlex(args.working_directory, proc_name=Path(__file__).stem, loggers='all') as log_handler:
        offtop = topology.topology_from_sdf(args.sdf_path)
        ommtop, ommsys, ommpos = openff_topology_to_openmm(offtop, forcefield=args.forcefield, box_vecs=args.box_vecs)
        history = execution.run_simulation_schedule(args.working_directory, args.schedule, ommtop, ommsys, ommpos, return_history=True)

if __name__ == '__main__':
    main()
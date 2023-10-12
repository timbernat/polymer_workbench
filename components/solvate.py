'''Add solvent to a Topology'''

import argparse
from pathlib import Path
import logging
LOGGER = logging.Logger(__name__)

import numpy as np

import openmm.unit
from openmm.unit import gram, centimeter

from polysaccharide2.genutils.fileutils.pathutils import assemble_path
from polysaccharide2.genutils.logutils.IOHandlers import MSFHandlerFlex
from polysaccharide2.openfftools.solvation import solvents, packing, boxvectors
from polysaccharide2.openfftools import topology


# accept args from CLI
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-wdir', '--working_directory'  , help='Directory into which files generated should be saved'  , type=Path)
parser.add_argument('-sdf' , '--sdf_path'           , help='Path to SDF file from which to read parameterized molecule structure', type=Path)
parser.add_argument('-a'   , '--affix'              , help='Additional text to attach to the end of the molecule name provided when naming ', type=str)

parser.add_argument('-solv', '--solvent'            , help='The name of the solvent to use to pack the box', type=str, default='water_TIP3P')
parser.add_argument('-bd'  , '--box_dimensions'     , help='XYZ dimensions of the desired periodic box for the system', type=list, nargs=3)
parser.add_argument('-bdu' , '--box_dimension_units', help='Unit to use when assigning box dimensions (default nanometer)', type=str, default='nanometer')
parser.add_argument('-rho' , '--density'            , help='Density (in grams/cm**3) up to which the desired box should be packed', type=float, default=1.0)
parser.add_argument('-exc' , '--exclusion'          , help='Name of the ForceField XML to use for MD', type=float, default=1.0)
parser.add_argument('-excu', '--exclusion_units'    , help='Unit to use when assigning box exlusion (default nanometer)', type=str, default='nanometer')


# post-process args as needed
args = parser.parse_args()
args.working_directory.mkdir(exist_ok=True)
assert(args.sdf_path.suffix == '.sdf')
solvent = getattr(solvents, args.solvent)

box_dim_unit = getattr(openmm.unit, args.box_dimension_unit)
box_dims = np.array(args.box_dimensions) * box_dim_unit
box_vecs = boxvectors.box_vectors_flexible(box_dims)

density = args.density * (gram/centimeter**3)
exclusion_unit = getattr(openmm.unit, args.exclusion_unit)
exclusion = args.exclusion * exclusion_unit

# main code
if __name__ == '__main__':
    with MSFHandlerFlex(args.working_dir, proc_name=__name__, loggers='all') as log_handler:
        offtop = topology.topology_from_sdf(args.sdf_path)
        offmol = topology.get_largest_offmol(offtop)
        mol_name = offmol.name

        solv_top = packing.pack_topology_with_solvent(offtop, solvent, box_vecs=box_dims, density=density, exclusion=exclusion)
        solv_top_path = assemble_path(args.working_dir, mol_name, extension='sdf', postfix=f'solv_{solvent.name}')
        topology.topology_to_sdf(solv_top_path, solv_top)
'''For generating chemically-specified SDF files from PDB and monomer JSON structure files'''

import argparse
from pathlib import Path
import logging
LOGGER = logging.Logger(__name__)

from openff.toolkit import Topology

from polysaccharide2.genutils.fileutils.pathutils import assemble_path
from polysaccharide2.genutils.logutils.IOHandlers import MSFHandlerFlex
from polysaccharide2.openfftools import TKREGS
from polysaccharide2.openfftools import topology
from polysaccharide2.monomers import MonomerGroup
from polysaccharide2.residues.partition import partition


# accept args from CLI
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-wdir', '--working_directory', help='Directory into which files generated should be saved'  , type=Path)
parser.add_argument('-pdb' , '--pdb_path'         , help='Path to PDB file from which to read molecule structure', type=Path)
parser.add_argument('-mono', '--monomer_path'     , help='Path to JSON file from which to read monomer chemistry', type=Path)
parser.add_argument('-n'   , '--mol_name'         , help='Name to assign to chemically-specified molecule'       , type=str)


# post-process args as needed
args = parser.parse_args()
args.working_directory.mkdir(exist_ok=True)
assert(args.pdb_path.suffix == '.pdb')
assert(args.monomer_path.suffix == '.json')


# main code
if __name__ == '__main__':
    # generate substructure cover and partition
    with MSFHandlerFlex(args.working_directory, proc_name=Path(__file__).stem, loggers='all') as log_handler:
        monogrp = MonomerGroup.from_file(args.monomer_path)
        offtop = Topology.from_pdb(args.pdb_path, _custom_substructures=monogrp.monomers, toolkit_registry=TKREGS['The RDKit'])
        was_partitioned = partition(offtop)
        assert(was_partitioned)

        # assign properties to Molecule
        offmol = topology.get_largest_offmol(offtop)
        offmol.name = args.mol_name
        offmol.properties['solvent'] = None
        offmol.properties['charge_method'] = None

        # save partitioned Topology
        mol_path = assemble_path(args.working_directory, args.mol_name, extension='sdf')
        topology.topology_to_sdf(mol_path, offtop=offtop, toolkit_registry=TKREGS['The RDKit'])
'''Assign charges to a chemically-defined Topology by a method of choice'''

import argparse
from pathlib import Path
import logging
LOGGER = logging.Logger(__name__)

from polysaccharide2.genutils.logutils.IOHandlers import MSFHandlerFlex
from polysaccharide2.genutils.fileutils.pathutils import assemble_path
from polysaccharide2.openfftools import topology
from polysaccharide2.openfftools.pcharge import MolCharger
from polysaccharide2.residues.rescharge.rctypes import ChargesByResidue


# accept args from CLI
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-wdir', '--working_directory', help='Directory into which files generated should be saved'  , type=Path)
parser.add_argument('-sdf' , '--sdf_path'         , help='Path to SDF file from which to read parameterized molecule structure', type=Path)

parser.add_argument('-cmet', '--charging_method'      , help='Charge assignment method to use to compute charges', type=str)
parser.add_argument('-lcp' , '--library_charge_path'  , help='Path to a JSON file containing library charges for the desired molecule', type=Path)

# post-process args as needed
args = parser.parse_args()
args.working_directory.mkdir(exist_ok=True)
assert(args.sdf_path.suffix == '.sdf')
assert(args.charging_method in MolCharger.subclass_registry)
assert(args.library_charge_path.suffix == '.json')


# main code
if __name__ == '__main__':
    with MSFHandlerFlex(args.working_dir, proc_name=__name__, loggers='all') as log_handler:
        offtop = topology.topology_from_sdf(args.sdf_path)
        offmol = topology.get_largest_offmol(offtop)
        mol_name = offmol.name

        ChargerType = MolCharger.subclass_registry[args.charging_method]
        if args.charging_method == 'RCT':
            if (args.library_charge_path is not None) and (args.library_charge_path.exists()):
                charger = ChargerType(ChargesByResidue.from_file(args.library_charge_path))
            else:
                raise FileExistsError
        else:
            charger = ChargerType()

        cmol = charger.charge_molecule(offmol, in_place=False)
        charged_top_path = assemble_path(args.working_dir, mol_name, extension='sdf', postfix=args.charging_method)
        topology.topology_to_sdf(charged_top_path, cmol.to_topology())
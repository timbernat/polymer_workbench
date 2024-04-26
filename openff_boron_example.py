# Import builtin/OpenFF
from pathlib import Path
from openff.interchange import Interchange
from openff.toolkit import Molecule, Topology, ForceField

# Import my custom tools which simplify parts of this process
from polymerist.openfftools import FFDIR, TKREGS, FF_DIR_REGISTRY, FF_PATH_REGISTRY
from polymerist.openfftools.topology import topology_to_sdf
from polymerist.openmmtools.serialization import serialize_openmm_pdb
from polymerist.openfftools.partialcharge import molchargers


# initialize molecule from SMILES string, generate conformation
smi = 'C=C-c1ccccc1-B(O)(O)'
offmol = Molecule.from_smiles(smi)
offmol.generate_conformers(n_conformers=1, toolkit_registry=TKREGS['OpenEye Toolkit'])

# assign partial charges
esp_charger = molchargers.EspalomaCharger()
esp_mol = esp_charger.charge_molecule(offmol)
esp_mol.partial_charges

# assign force field parameters
ffname = 'smirnoff99Frosst-1.0.0' # this was the only OpenFF forcefield I could find which could parameterize Boron 
ff = ForceField(f'{ffname}.offxml')
inc = ff.create_interchange(esp_mol.to_topology(), charge_from_molecules=[esp_mol])

# output PDB and (the preferable) SDF Topology files and GMX files
molname = '4-vinylphenyl-boronic_acid'
outdir = Path(molname)
outdir.mkdir(exist_ok=True)

topology_to_sdf(     outdir / f'{molname}.sdf', esp_mol.to_topology())
serialize_openmm_pdb(outdir / f'{molname}.pdb', inc.to_openmm_topology(), inc.positions.to_openmm())
inc.to_gromacs(prefix=f'{molname}/{molname}')
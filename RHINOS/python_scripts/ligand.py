import sys
from openff.toolkit import Molecule
from rdkit import Chem

openFF_mol = Molecule.from_file(str(sys.argv[1]))

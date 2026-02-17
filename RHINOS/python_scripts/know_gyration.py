from rdkit import Chem
from rdkit.Chem import Descriptors3D
import sys

molecule = Chem.MolFromPDBFile(str(sys.argv[1]))
gyration = Descriptors3D.RadiusOfGyration(molecule)

print(gyration * 2.9)

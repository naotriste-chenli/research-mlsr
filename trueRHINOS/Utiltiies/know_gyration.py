from rdkit import Chem
from rdkit.Chem import Descriptors3D
import sys

mol = str(sys.argv[1])
mol = Chem.MolFromMolFile(mol)
rog = Descriptors3D.RadiusOfGyration(mol) * 2.9

print(rog)
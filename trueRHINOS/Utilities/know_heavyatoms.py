from rdkit import Chem
from rdkit.Chem import Lipinski
import sys

mol = str(sys.argv[1])
mol = Chem.MolFromMolFile(mol)

heavy_atoms = Lipinski.HeavyAtomCount(mol)

print(str(heavy_atoms))
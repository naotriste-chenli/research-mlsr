from rdkit import Chem
from rdkit.Chem import AllChem
import sys
# 1. Load your stuff
mol = Chem.AddHs(Chem.MolFromSmiles("C#CC(OC)C(CO)CCC(=O)SCC(S)OO"))
core = Chem.MolFromMolFile(str(sys.argv[1])) # Your 3D reference

# 2. The "Easy" way to get the mapping
# This gets the indices of atoms in 'mol' that match 'core'
match = mol.GetSubstructMatch(core)

# 3. Create the 3D constraint (The "Force" step)
conf = core.GetConformer()
coord_map = {mol_idx: conf.GetAtomPosition(core_idx) 
             for core_idx, mol_idx in enumerate(match)}

# 4. Generate the 3D structure relative to that map
params = AllChem.ETKDGv3()
params.coordMap = coord_map
AllChem.EmbedMolecule(mol, params)

# 5. Clean up the new parts while keeping the scaffold fixed
AllChem.UFFOptimizeMolecule(mol, confId=0, vdwThresh=10.0)

Chem.MolToMolFile(mol, 'constrained_bioisostere.mol')

from rdkit import Chem
from rdkit.Chem import AllChem

suppl = Chem.SDMolSupplier('seeds.sdf', removeHs=False)
writer = Chem.SDWriter('seeds_cleaned.sdf')

for mol in suppl:
    if mol is None: continue
    # Add H's and generate a clean 3D conformation
    mol = Chem.AddHs(mol, addCoords=True)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)
    writer.write(mol)
writer.close()

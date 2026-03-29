from rdkit import Chem
from fragmenstein import Monster
import pyrosetta
import sys

pyrosetta.init(extra_options='-mute all -ignore_unrecognized_res true')


ref_hit = Chem.MolFromMolFile(sys.argv[1], sanitize=False)
ref_hit.SetProp("_Name", "seed_8") 

target_mol = Chem.MolFromMolFile(sys.argv[2], sanitize=False)

target_mol.SetProp("_Name", "placed_ligand")

monster = Monster(hits=[ref_hit])
monster.place(target_mol)

aligned_mol = monster.positioned_mol
Chem.MolToMolFile(aligned_mol, sys.argv[3])

print(f"Successfully placed {sys.argv[2]} onto {sys.argv[1]} as {sys.argv[3]}")

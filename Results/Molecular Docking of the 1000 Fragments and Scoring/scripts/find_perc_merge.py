from rdkit import Chem
from fragmenstein import Victor
import pyrosetta
import sys

pyrosetta.init(extra_options='-mute all -ignore_unrecognized_res true')

input_1 = Chem.MolFromMolFile(str(sys.argv[1]), sanitize=False)
input_2 = Chem.MolFromMolFile(str(sys.argv[2]), sanitize=False)
input_3 = str(sys.argv[4])

victor = Victor([input_1, input_2], pdb_filename=input_3)
victor.combine()
victor_pos = victor.monster.positioned_mol
victor_min = victor.minimized_mol

Chem.MolToMolFile(victor_min, str(sys.argv[3]))

from prody import *
import sys

pdb_file = str(sys.argv[1])
pdb_output = str(sys.argv[2])

protein = parsePDB(pdb_file)
select_prot = protein.select('protein')

writePDB(pdb_output, select_prot)
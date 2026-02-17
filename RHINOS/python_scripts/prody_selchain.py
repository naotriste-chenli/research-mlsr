from prody import *
import sys

file_name = str(sys.argv[1])
selchain = str(sys.argv[2])
output = str(sys.argv[3])

protein = parsePDB(file_name, chain=selchain)

writePDB(output, protein)
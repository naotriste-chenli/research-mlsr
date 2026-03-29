from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
import sys

input=str(sys.argv[1])
fpg = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048)

mols = Chem.SDMolSupplier(input, sanitize=False)
mols = [fpg.GetFingerprint(x) for x in mols]

for i in mols:
	tanimoto_similarity = DataStructs.TanimotoSimilarity(mols[0], i)
	print(tanimoto_similarity)

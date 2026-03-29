import pandas as pd
import sys

df = pd.read_csv(str(sys.argv[1]))
perc=df['VDss_Lombardo_drugbank_approved_percentile'].quantile(0.99)
top_smiles = df[df['VDss_Lombardo_drugbank_approved_percentile'] >= perc]['smiles']

print(top_smiles)


import pandas as pd
from numpy import int64

def transformTableAfterR(file, output):
    counts = pd.read_table(file)
    counts.index.name="Geneid"
    for col in counts:
        counts[col] = counts[col].astype(int64)
    counts.to_csv(output+'/count_table.tsv', sep='\t', index=True)
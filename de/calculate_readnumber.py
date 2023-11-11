import pandas as pd
df = pd.read_table("/home/ubuntu/splivardet/de/tables/result_table.tsv")
Total = df['baseMean'].sum()
print(Total)
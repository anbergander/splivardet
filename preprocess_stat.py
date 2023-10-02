
import sys
import pandas as pd
import numpy as np
import csv
import collections
from operator import add
import pickle

import warnings
warnings.filterwarnings('ignore')

def create_df_from_dict_entry(count_dict, grouplist, key):
    temp_dict = collections.defaultdict(list)
    for group in grouplist:
        for i in range(len(count_dict[key])):
            if count_dict["ids"][i] == group and count_dict[key][i] != group:
                temp_dict[group].append(count_dict[key][i])
        if len(temp_dict[group]) < max_length and key != "ids":
            temp_dict[group] = temp_dict[group] + ([np.nan] * (max_length - len(temp_dict[group])))

    if key != "ids":
        temp_df = pd.DataFrame.from_dict(temp_dict)

        names = list(temp_df.columns)
        temp_df_melt = pd.melt(temp_df.reset_index(), id_vars=['index'], value_vars=names)
        temp_df_melt.columns = ['index', 'group', 'value']
    else:
        return [], []
        
    return temp_df, temp_df_melt

def norm_by_exp(n, d):

    return 100*n/d if d else 0

def transform_z_o(pl, gl):
    n = len(gl)
    return list(map(lambda x: (x*(n-1)+0.5)/n, pl))
    return pl.apply(lambda x: (x*(n-1)+0.5)/n)
    

def identGroup(a):
    gr = a.split("_")[1]
    return gr


def identBarcode(a):
    bc = a.split("_")[0]
    return bc


def read_file_as_dict(path, delim):
    counts = {}
    with open(path, 'r') as file:
        csvreader = csv.reader(file, delimiter=delim)
        for row in csvreader:
            counts[row[0]] = row[1:]
    return counts


def read_tsv_as_dict(path):
    counts = read_file_as_dict(path, "\t")
    identifier_verbose = counts["ids"]
    counts["ids"] = list(map(identGroup, identifier_verbose))
    counts["barcodes"] = list(map(identBarcode, identifier_verbose))
    return counts

file = '/vol/fastq/hg38/barcode02.quantify.counts.tsv'
read_count_file = '/vol/fastq/hg38/read_counts.csv'


read_count_dict = read_file_as_dict(read_count_file, ",")
read_df = pd.DataFrame.from_dict(read_count_dict)
#print(read_df)
   
counts = read_tsv_as_dict(file)
raw_counts = dict(counts)
for key in counts:
    new_values = []
    raw_values = []
    if key == "ids" or key == "barcodes":
        continue
    for i in range(len(counts[key])):
        barcode = counts['barcodes'][i]
        read_count = int(read_count_dict[barcode][0])
        raw_count = float(counts[key][i])
        factor = 1000000 / read_count
        norm_count = raw_count * factor
        new_values.append(norm_count)
    counts[key] = new_values
    new_values = []



barcode_index= [int(x.split("barcode")[1])for x in counts["barcodes"]]
#res = dict(map(lambda i,j : (i,j) , barcode_index, counts["barcodes"]))
barcode_names = counts["barcodes"]
counts.pop("barcodes")
#sample_group = list(map(identGroup, df.index[1:]))

dataframes = {}
dataframes_raw = {}
dataframes_melt = {}
dataframes_melt_raw = {}
max_length = collections.Counter(counts["ids"]).most_common()[0][1]
anova_table_dict = collections.defaultdict(list)
p_value_list = []
gene_count_list = collections.defaultdict(list)


for key in counts.keys():
    if gene_count_list[key.split("_")[-1]] == []:
        gene_count_list[key.split("_")[-1]] = counts[key] 
        
    else:
        gene_count_list[key.split("_")[-1]] = list(map(add, gene_count_list[key.split("_")[-1]], counts[key]))

grouplist = set(counts["ids"])

gene_exp_df = pd.DataFrame.from_dict(gene_count_list, orient='index')
gene_exp_df.columns = barcode_names

for key in counts.keys():
    
    
    if key != "ids":
        counts[key] = list(map(norm_by_exp, counts[key], gene_count_list[key.split("_")[-1]]))

    
    
    group_list = set(counts["ids"])
    
    temp_df, temp_df_melt = create_df_from_dict_entry(counts, grouplist, key)
    temp_df_raw, temp_df_melt_raw = create_df_from_dict_entry(raw_counts, grouplist, key)

    if key != "ids":
 
        dataframes_melt[key] = temp_df_melt
        dataframes[key] = temp_df
        dataframes_melt_raw[key] = temp_df_melt_raw
        dataframes_raw[key] = temp_df_raw
        temp_df = dataframes_melt[key].dropna()
        join_df = False
        for group in grouplist:
            group_df = temp_df[temp_df['group'] == group]
            
            group_df['transf_value'] = transform_z_o(group_df['value'], counts['ids'])
            
            if type(join_df) != bool:
                join_df = pd.concat([join_df, group_df])
            else:
                join_df = group_df
            
        join_df.to_csv(key+".csv")
        join_df = False           
        with open("isoform_keys", 'a+') as f:
            f.write(str(key)+"\n")
            
count_df = pd.DataFrame.from_dict(counts, orient = 'index')
count_df.columns = barcode_names
count_df = count_df.drop(['ids'])
count_df = count_df.reset_index()
count_df = count_df.rename(columns={"index": "transcript"})
count_df.to_csv('norm_counts.tsv',sep = '\t', index = False)


with open('dataframes_melt_raw.pickle', 'wb') as f:
    pickle.dump(dataframes_melt_raw, f)
with open('dataframes_melt.pickle', 'wb') as f:
    pickle.dump(dataframes_melt, f)



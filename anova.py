
import sys
import os
import pandas as pd
import numpy as np
import csv
import collections
from operator import add, truediv

import scipy.stats as stats
import statsmodels.api as sm
from scipy.stats import shapiro
from scipy.stats import kruskal
import statsmodels
import plot_anova
from statsmodels.formula.api import ols

import warnings
warnings.filterwarnings("ignore")

def norm_by_exp(n, d):

    return n/d if d else 0

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


file = sys.argv[1]
read_count_file = sys.argv[2]
df = pd.read_csv(file, sep="\t").T

read_count_dict = read_file_as_dict(read_count_file, ",")
#read_df = pd.DataFrame.from_dict(read_count_dict)
#print(read_df)
    
counts = read_tsv_as_dict(file)
for key in counts:
    new_values = []
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
sample_group = list(map(identGroup, df.index[1:]))

dataframes = {}
dataframes_melt = {}
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

gene_exp_df.to_csv("normalized_gene_count.csv")
for key in counts.keys():
    
    
    if key != "ids":
        counts[key] = list(map(norm_by_exp, counts[key], gene_count_list[key.split("_")[-1]]))

    temp = collections.defaultdict(list)
    
    group_list = set(counts["ids"])
    
    for group in grouplist:
        for i in range(len(counts[key])):
            if counts["ids"][i] == group and counts[key][i] != group:
                temp[group].append(counts[key][i])
        if len(temp[group]) < max_length and key != "ids":
            temp[group] = temp[group] + ([np.nan] * (max_length - len(temp[group])))

    if key != "ids":
        temp_df = pd.DataFrame.from_dict(temp)

        names = list(temp_df.columns)
        temp_df_melt = pd.melt(temp_df.reset_index(), id_vars=['index'], value_vars=names)
        temp_df_melt.columns = ['index', 'group', 'value']
        #temp_df_melt['value'] = temp_df_melt['value']
        dataframes_melt[key] = temp_df_melt
        dataframes[key] = temp_df
        #dataframes[key] = temp_df.drop(['rbm20ttn'], axis=1)
        temp_df = dataframes_melt[key].dropna()
        join_df = False
        for group in grouplist:
            group_df = temp_df[temp_df['group'] == group]
            
            group_df['transf_value'] = transform_z_o(group_df['value'], counts['ids'])
            print(group_df)
            if type(join_df) != bool:
                join_df = pd.concat([join_df, group_df])
            else:
                join_df = group_df
            
        join_df.to_csv(key+".csv")
        join_df = False           
        
        dataset = [dataframes_melt[key][dataframes_melt[key]['group']==group]['value'].dropna().tolist() for group in grouplist]
        try:
            kruskal_table = kruskal(*dataset, nan_policy='omit')
            p_value = kruskal_table[1]
            p_value_list.append((key, p_value))
        except:
            pass
        
        #model = ols('value ~ C(group)', data=dataframes_melt[key]).fit()
        #anova_table = sm.stats.anova_lm(model, typ=2)
        #anova_table_dict[key] = anova_table
        #p_value = anova_table.at["C(group)", "PR(>F)"]
        
        
 
       # note: if the data is balanced (equal sample size for each group), Type 1, 2, and 3 sums of squares
        # (typ parameter) will produce similar results.
p_value = [y for (x, y) in p_value_list]
keys = np.array([x for (x, y) in p_value_list])
print(p_value)
a = statsmodels.stats.multitest.multipletests(p_value, method='fdr_bh')

candidate_list = keys[a[0]]
print(len(candidate_list))
candidate_file = open(os.path.basename(file).split(".")[0]+'.items.txt','w')
candidate_file.writelines([string + '\n' for string in candidate_list])
candidate_file.close()
for candidate in candidate_list:
    print (candidate)
    temp = dataframes_melt[candidate]
    temp.to_csv(path_or_buf= candidate+".csv", index=False)
    plot_anova.plot_posthoc(candidate+".csv")
    
print("end")
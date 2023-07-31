
import sys
import os
import pandas as pd
import numpy as np
import csv
import collections

import scipy.stats as stats
import statsmodels.api as sm
import statsmodels

from statsmodels.formula.api import ols



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

counts.pop("barcodes")
sample_group = list(map(identGroup, df.index[1:]))

dataframes = {}
dataframes_melt = {}
max_length = collections.Counter(counts["ids"]).most_common()[0][1]
anova_table_dict = collections.defaultdict(list)
p_value_list = []

for key in counts.keys():

    temp = collections.defaultdict(list)

    for group in set(counts["ids"]):
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
        temp_df_melt['value'] = temp_df_melt['value']
        dataframes_melt[key] = temp_df_melt
        dataframes[key] = temp_df.drop(['rbm20ttn'], axis=1)

        # dataframe_melt = dataframes_melt[key]
        # dataframe_melt.dropna(inplace=True)
        # ax = sns.boxplot(y='value', x='group', data=dataframe_melt, color='#99c2a2')
        # ax = sns.swarmplot(y='value', x='group', data=dataframe_melt, color='#7d0013')
        # plt.show()

        # stats f_oneway functions takes the groups as input and returns ANOVA F and p

        x = dataframes[key]
        values = [*x.to_numpy().T]
        values = [[x for x in value_list if ~np.isnan(x)] for value_list in values]
        fvalue, pvalue = stats.f_oneway(*values)

        # Ordinary Least Squares (OLS) model
        model = ols('value ~ C(group)', data=dataframes_melt[key]).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
        anova_table_dict[key] = anova_table
        p_value = anova_table.at["C(group)", "PR(>F)"]
        p_value_list.append((key, p_value))

        # note: if the data is balanced (equal sample size for each group), Type 1, 2, and 3 sums of squares
        # (typ parameter) will produce similar results.
p_value = [y for (x, y) in p_value_list]
keys = np.array([x for (x, y) in p_value_list])
a = statsmodels.stats.multitest.multipletests(p_value, method='fdr_bh')
candidate_list = keys[a[0]]
candidate_file = open(os.path.basename(file).split(".")[0]+'.items.txt','w')
candidate_file.writelines([string + '\n' for string in candidate_list])
candidate_file.close()
#for candidate in candidate_list:
#    temp = dataframes_melt[candidate]

#    temp.dropna(inplace=True)
#    temp = temp[temp.group != "rbm20ttn"]
#    print("bp")
#    res = stat()
#    res.tukey_hsd(df=temp, res_var='value', xfac_var='group', anova_model='value ~ C(group)')
#    res.tukey_summaryort
#    stat
print("end")

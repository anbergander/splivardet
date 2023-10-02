#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri Aug 25 17:22:47 2023

@author: denbi
"""
import pickle
import pandas as pd
from statsmodels.stats.multitest import multipletests
import numpy as np

file = open('dataframes_melt_raw.pickle','rb')
dataframes_melt_raw = pickle.load(file)
file.close()
file = open('dataframes_melt.pickle','rb')
dataframes_melt = pickle.load(file)
file.close()
posthoc_df = pd.read_csv('posthoc.csv')

posthoc_df[['p.value']] = posthoc_df[['p.value']].apply(pd.to_numeric)

#print(len(posthoc_df.loc[posthoc_df['p.value'] < 0.05]))
p_value = posthoc_df['p.value']
isoform_keys = np.array(posthoc_df['isoform_key'])

key_list = list(dataframes_melt.keys())
key = key_list[7737]
#print(key)
#print(dataframes_melt[key])
#print(dataframes_melt_raw[key])

a = multipletests(p_value, method='bonferroni', alpha=(0.05/3))
candidate_list = isoform_keys[a[0]]

temp_list = []
for candidate in candidate_list:
    read_list = list(map(float, (dataframes_melt_raw[candidate].dropna()['value'])))

    if sum(read_list)>=len(read_list):
        
        temp_list.append(candidate)
candidate_list = temp_list



ph_df = pd.read_csv("posthoc.csv")


#create unique list of names
isoform_keys = ph_df.isoform_key.unique()

#create a data frame dictionary to store your data frames
DataFrameDict = {elem : pd.DataFrame() for elem in isoform_keys}
grouplist = []
gl_does_not_exist = True
for key in candidate_list:
   
    DataFrameDict[key] = ph_df[:][ph_df.isoform_key == key].drop("Unnamed: 0",axis=1)
    DataFrameDict[key] = DataFrameDict[key].reset_index().drop(columns=['estimate', 'SE', 'df', 'z.ratio', 'index'])
    DataFrameDict[key] = DataFrameDict[key].drop_duplicates()
    if gl_does_not_exist:
        for i in DataFrameDict[key]['contrast'].tolist():
            grouplist.extend(i.split(' - '))
            gl_does_not_exist = False

        
grouplist = set(grouplist)
candidate_list =list(set(candidate_list))
for key in candidate_list:
    dummy_values = [1]*len(grouplist)
    temp_df = pd.DataFrame(data={i:dummy_values for i in grouplist})
    temp_df['gl'] = list(grouplist)
    temp_df.set_index('gl', inplace = True)
    temp_df.index.name = None

    for i in list(DataFrameDict[key][['contrast', 'p.value']].itertuples(index=False, name=None)):
        gg, pv = i
        groups = gg.split(" - ")
        temp_df.at[groups[0], groups[1]] = pv
        temp_df.at[groups[1], groups[0]] = pv
    DataFrameDict[key] = temp_df
keys = list(DataFrameDict.keys())
for key in keys:
    if key not in candidate_list:
        DataFrameDict.pop(key)
with open('posthoc.pickle', 'wb') as f:
    pickle.dump(DataFrameDict, f)

#print(DataFrameDict[isoform_keys[1]])

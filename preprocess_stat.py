
import sys
import pandas as pd
import numpy as np
import csv
import collections
from operator import add
import pickle
from gtfparse import read_gtf
from itertools import compress

import warnings
warnings.filterwarnings('ignore')


isofile = sys.argv[1]
wo = sys.argv[4]+'/'

#wo = '/vol/fastq/splivardet/'
#isofile = '/vol/fastq/hg38/barcode02.isoforms.gtf'

gtfp = read_gtf(isofile)

gl = list(gtfp['gene_id'])

uniprot = [ p for p in gl if not(p.startswith("chr") or p.startswith("ENST"))]
with open(wo+'resources/all_genes.txt', 'w+') as fp:
    fp.write('\n'.join(list(set(uniprot))))


uniprot = pd.read_csv(wo+"resources/uniprot.tsv", sep= "\t")

up = dict(zip(uniprot['Entry'],uniprot['Gene Names']))

file = open(wo+"resources/refseq_ensembl.txt")
annot = file.read()
file.close()

file = open(wo+"resources/transcript_to_gene.txt")
trans_to_gene = file.read()
file.close()

file = open(wo+"resources/ens_gensym.txt")
ens_gensym = file.read()
file.close()




a_d = {}
for n in annot.split("\n"):
    s = n.split("\t")
    if len(s)==2:
        a_d[s[1]] = s[0]
    
t_g = {}
for n in trans_to_gene.split("\n"):
    s = n.split("\t")
    if len(s)==2:
        t_g[s[0]] = s[1]

e_s = {}
s_e = {}
for n in ens_gensym.split("\n"):
    s = n.split("\t")
    if len(s)==2:
        e_s[s[1]] = s[0]
        s_e[s[0]] = s[1]


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

def transform_z_o(l, gl):

    x = float(l)
    x = 1 if x > 1 else x
    n = len(gl)
    return (x*(n-1)+0.5)/n

    

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

def gene_name(gene):

    try:
        
        if gene.startswith('ENST'):
            gn = t_g[gene.split('.')[0]]
        elif ':' in gene:
            gn = gene
        else:
            gn = e_s[str(up[gene]).split(' ')[0]]
    except: 
        gn = 'nan'

    return gn


file = sys.argv[3]
read_count_file = sys.argv[2]

#file='/vol/fastq/hg38/barcode02.quantify.counts.tsv'
#read_count_file = '/vol/fastq/hg38/read_counts.csv'

read_count_dict = read_file_as_dict(read_count_file, ",")
read_df = pd.DataFrame.from_dict(read_count_dict)

   
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
barcode_names = counts["barcodes"]
grouplist = set(counts["ids"])
ids = counts["ids"]

gr_bc_v = []
for i in ids:
    gr_bc_v.append(ids.count(i))

bc_counts = dict(zip(barcode_names, gr_bc_v))

counts.pop("barcodes")
counts.pop("ids")

count_df = pd.DataFrame.from_dict(counts, orient = 'index')
count_df.columns = barcode_names

count_df = count_df.reset_index()
count_df = count_df.rename(columns={"index": "transcript"})

known = count_df[count_df.transcript.str.contains('ENST')]
unknown = count_df[~count_df.transcript.str.contains('ENST')]
unknown = unknown[~unknown.transcript.str.contains('chr')]

a = known.transcript.str.split('_')
x = [i[0] for i in a]
y = [i[-1] for i in a]
f1 = [i.startswith('ENST') for i in x]
f2 = [i.startswith('ENST') for i in y]

known['ENST'] = x
known['ui'] = y


temp_1 = known[known.ui.str.startswith('ENST')]

ENST = temp_1.ui.tolist()
ENST = [i.split('.')[0]for i in ENST]
temp_1['ENST'] = ENST

temp_2 = known[~known.ui.str.startswith('ENST')]

ENST = temp_2.ENST.tolist()
ENST = [i.split('.')[0]for i in ENST]
temp_2['ENST'] = ENST

temp = pd.concat([temp_1, temp_2])
temp.drop(columns = 'ui', inplace = True)
ENST = temp.ENST.tolist()
ENSG = [t_g[i] if i in t_g else np.nan for i in ENST]
temp['ENSG'] = ENSG
temp.dropna(inplace=True)
temp['transcript'] = temp.ENST + '_' + temp.ENSG

known = temp

transc = unknown.transcript.tolist()
FLAIR = [i.split('_')[0]for i in transc]
ui = [i.split('_')[-1]for i in transc]
sym = [up[i] if i in up else np.nan for i in ui]
sym = [i.split(' ')[0] if type(i) != float else np.nan for i in sym]
ENSG = [e_s[i] if i in e_s else np.nan for i in sym]

unknown['ENSG'] = ENSG
unknown['ENST'] = FLAIR
unknown['transcript'] = unknown.ENST + '_' + unknown.ENSG
unknown.dropna(inplace=True)


count_df = pd.concat([known, unknown])

agg = dict.fromkeys(count_df.columns.difference(['ENST']), 'sum')
agg['ENSG'] = 'first'
agg['transcript'] = 'first'

aggr_count_df = count_df.groupby('ENST').agg(agg)

count_df.drop(columns=['ENST'], inplace = True)
aggr_count_df.reset_index(inplace=True)

temp = aggr_count_df.drop(columns=['transcript', 'ENST'])

groups = temp.groupby('ENSG')

agg = dict.fromkeys(temp.columns.difference(['ENSG']), 'sum')


norm_df = groups.agg(agg)


div_df = pd.DataFrame()

for e in aggr_count_df.ENSG.tolist():

    div_df = div_df.append(norm_df.loc[e], ignore_index=True)

aggr_count_df[aggr_count_df.columns[2:-1]] = aggr_count_df[aggr_count_df.columns[2:-1]]/div_df
aggr_count_df.fillna(0,inplace=True)
aggr_count_df[aggr_count_df.columns[2:-1]] = aggr_count_df[aggr_count_df.columns[2:-1]].applymap(transform_z_o, gl=ids).round(decimals=5)

aggr_count_df.set_index('transcript', inplace=True)
aggr_count_df.drop(columns=['ENSG','ENST'], inplace=True)

aggr_count_df.to_csv('norm_counts.tsv',sep = '\t', index = True)

count_df = aggr_count_df

tight = count_df.to_dict('tight')

counts = dict(zip(tight['index'], tight['data']))
counts['barcodes'] = barcode_names
counts['ids'] = ids

barcode_index= [int(x.split("barcode")[1])for x in barcode_names]

dataframes = {}

dataframes_melt = {}

max_length = collections.Counter(counts["ids"]).most_common()[0][1]

counts.pop('barcodes')

for key in counts.keys():

    
    group_list = set(counts["ids"])
    
    temp_df, temp_df_melt = create_df_from_dict_entry(counts, grouplist, key)
   
    if key != "ids":
 
        dataframes_melt[key] = temp_df_melt
        dataframes[key] = temp_df

        temp_df = dataframes_melt[key].dropna()
        join_df = False
        for group in grouplist:
            group_df = temp_df[temp_df['group'] == group]
            
            group_df['transf_value'] = group_df['value']
            
            if type(join_df) != bool:
                join_df = pd.concat([join_df, group_df])
            else:
                join_df = group_df
        if len(join_df.index) < 2:
            print(key)
            join_df = False
            continue
        join_df.to_csv(key+".csv")
        join_df = False           
        with open("isoform_keys", 'a+') as f:
            f.write(str(key)+"\n")
            

with open('dataframes_melt.pickle', 'wb') as f:
    pickle.dump(dataframes_melt, f)



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 16:30:00 2023

@author: denbi
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 10:50:50 2023

@author: denbi
"""
from swan_vis.utils import *
from swan_vis.talon_utils import *
from swan_vis.graph import *
import swan_vis as swan
import csv

from swan_vis.plottedgraph import PlottedGraph
from swan_vis.report import Report
import pandas as pd
import matplotlib as mpl
from scipy import sparse
import os
import sys

import warnings
warnings.filterwarnings("ignore")




# initialize a new SwanGraph
sg = swan.SwanGraph(sc=True)


bc = sys.argv[1]
candidates = 'candidates.txt'

cand_l = []
with open(candidates, 'r') as file:
  csvreader = csv.reader(file)
  for row in csvreader:
    if row[0].startswith('X') or row[0].startswith('NR'):
        continue
    entry = row[0].split(' ')[0]
    entry = entry.split('_')[-1]
    if entry.startswith('chr'):
        continue
    cand_l.append(entry)

cand_l = list(set(cand_l))

annot_gtf = 'swan_annot.gtf'
data_gtf = 'gtf.gtf'
ab_file = 'abundance.tsv'
meta = sys.argv[3]

meta = pd.read_csv(meta, sep = "\t",header=None)
meta.columns = ['dataset', 'cond', 'tissue', 'fastqpath'] 
meta.to_csv('temp_manifest.tab', index=False, sep='\t')
meta = 'temp_manifest.tab'

# add abundance for each dataset from the AnnData into the SwanGraph
sg.add_annotation(annot_gtf)
sg.add_transcriptome(data_gtf)
adata = sg.abundance_to_adata(ab_file, how='iso')
sg.add_abundance(ab_file)


sg.add_metadata(meta)

def plot_each_transcript(sg, tids, prefix,
                    indicate_dataset=False,
                    indicate_novel=False,
                    browser=False):
    """
    Plot each input transcript and automatically save figures

    Parameters:
        tids (list of str): List of transcript ids to plot
        prefix (str): Path and file prefix to automatically save
            the plotted figures
        indicate_dataset (str): Dataset name from SwanGraph to
            highlight with outlined nodes and dashed edges
            Incompatible with indicate_novel
            Default: False (no highlighting)
        indicate_novel (bool): Highlight novel nodes and edges by
            outlining them and dashing them respectively
            Incompatible with indicate_dataset
            Default: False
        browser (bool): Plot transcript models in genome browser-
            style format. Incompatible with indicate_dataset and
            indicate_novel
    """

    sg.check_plotting_args(indicate_dataset, indicate_novel, browser)

    # loop through each transcript in the SwanGraph object
    for tid in tids:
        sg.check_transcript(tid)

    for tid in tids:
        sg.pg.init_plot_settings(sg, tid=tid,
            indicate_dataset=indicate_dataset,
            indicate_novel=indicate_novel,
            browser=browser)
        fname = create_fname(prefix,
                             indicate_dataset,
                             indicate_novel,
                             browser,
                             ftype='path',
                             tid=tid)
        sg.pg.plot_graph()
        print('Saving transcript path graph for {} as {}'.format(tid, fname))
        save_fig(fname)

def plot_each_transcript_in_gene(sg, gid, prefix,
                         indicate_dataset=False,
                         indicate_novel=False,
                         browser=False):
    """
    Plot each transcript in a given gene and automatically save figures

    Parameters:
        gid (str): Gene id or gene name to plot transcripts from
        prefix (str): Path and file prefix to automatically save
            the plotted figures
        indicate_dataset (str): Dataset name from SwanGraph to
            highlight with outlined nodes and dashed edges
            Incompatible with indicate_novel
            Default: False (no highlighting)
        indicate_novel (bool): Highlight novel nodes and edges by
            outlining them and dashing them respectively
            Incompatible with indicate_dataset
            Default: False
        browser (bool): Plot transcript models in genome browser-
            style format. Incompatible with indicate_dataset and
            indicate_novel
    """

    if gid not in sg.t_df.gid.tolist():
        gid = sg.get_gid_from_gname(gid)
    sg.check_gene(gid)

    sg.check_plotting_args(indicate_dataset, indicate_novel, browser)

    # loop through each transcript in the SwanGraph object
    tids = sg.t_df.loc[sg.t_df.gid == gid, 'tid'].tolist()
    
    print()
    print('Plotting {} transcripts for {}'.format(len(tids), gid))
    for tid in tids:
        sg.pg.init_plot_settings(sg, tid=tid,
            indicate_dataset=indicate_dataset,
            indicate_novel=indicate_novel,
            browser=browser)
        fname = create_fname(prefix,
                             indicate_dataset,
                             indicate_novel,
                             browser,
                             ftype='path',
                             tid=tid)
        sg.pg.plot_graph()
        print('Saving transcript path graph for {} as {}'.format(tid, fname))
        save_fig(fname)

def subset_on_gene_sg(sg, gid=None, datasets=None):
        """
        Subset the swan Graph on a given gene and return the subset graph.

        Parameters:
            gid (str): Gene ID to subset on
            datasets (list of str): List of datasets to keep in the subset

        returns:
            subset_sg (swan Graph): Swan Graph subset on the input gene.
        """

        # didn't ask for either
        if not gid and not datasets:
            return sg

        # subset on gene
        if gid:
            # make sure this gid is even in the Graph
            sg.check_gene(gid)

            # get the strand
            strand = sg.get_strand_from_gid(gid)

            # subset t_df first, it's the easiest
            tids = sg.t_df.loc[sg.t_df.gid == gid].index.tolist()
            t_df = sg.t_df.loc[tids].copy(deep=True)
            t_df['path'] = sg.t_df.loc[tids].apply(
                    lambda x: copy.deepcopy(x.path), axis=1)
            t_df['loc_path'] = sg.t_df.loc[tids].apply(
                    lambda x: copy.deepcopy(x.loc_path), axis=1)

            # since we don't keep all transcripts in adata, make
            # sure to pare that down
            tids = list(set(tids)&set(sg.adata.var.index.tolist()))

            # subset loc_df based on all the locs that are in the paths from
            # the already-subset t_df
            paths = t_df['loc_path'].tolist()
            locs = [node for path in paths for node in path]
            locs = np.unique(locs)
            loc_df = sg.loc_df.loc[locs].copy(deep=True)

            # subset edge_df based on all the edges that are in the paths from
            # the alread-subset t_df
            paths = t_df['path'].tolist()
            edges = [node for path in paths for node in path]
            edges = np.unique(edges)
            edge_df = sg.edge_df.loc[edges].copy(deep=True)
        if not gid:
            t_df = sg.t_df.copy(deep=True)
            edge_df = sg.edge_df.copy(deep=True)
            loc_df = sg.loc_df.copy(deep=True)


        new_adatas = dict()
        # adatas = {'iso': sg.adata, 'edge': sg.edge_adata,
        #           'tss': sg.tss_adata, 'tes': sg.tes_adata}
        adatas = {'iso': sg.adata}
        for key, adata in adatas.items():

            if datasets and gid:
                new_adatas[key] = adata[datasets, tids]
            elif gid:
                new_adatas[key] = adata[:, tids]
            elif datasets:
                new_adatas[key] = adata[datasets, :]
            else:
                new_adatas[key] = adata

        # create a new graph that's been subset
        subset_sg = SwanGraph()
        subset_sg.loc_df = loc_df
        subset_sg.edge_df = edge_df
        subset_sg.t_df = t_df
        subset_sg.adata = new_adatas['iso']
        # subset_sg.edge_adata = new_adatas['edge']
        # subset_sg.tss_adata = new_adatas['tss']
        # subset_sg.tes_adata = new_adatas['tes']
        subset_sg.datasets = subset_sg.adata.obs.index.tolist()
        subset_sg.abundance = sg.abundance
        subset_sg.sc = sg.sc
        subset_sg.pg = sg.pg
        subset_sg.annotation = sg.annotation

        # renumber locs if using a gene
        if gid:
            if strand == '-':
                id_map = subset_sg.get_ordered_id_map(rev_strand=True)
                subset_sg.update_ids(id_map)
            else:
                subset_sg.update_ids()

            subset_sg.get_loc_types()

        # finally create the graph
        subset_sg.create_graph_from_dfs()

        return subset_sg

def get_tpm(sg, kind='iso'):
    """
    Retrieve TPM per dataset.

    Parameters:
        kind (str): {'iso', 'edge', 'tss', 'tes', 'ic'}
            Default: 'iso'

    Returns:
        df (pandas DataFrame): Pandas dataframe where rows are the different
            conditions from `dataset` and the columns are ids in the
            SwanGraph, and values represent the TPM value per
            isoform/edge/tss/tes/ic per dataset.
    """
    if kind == 'iso':
        adata = sg.adata
    elif kind == 'edge':
        adata = sg.edge_adata
    elif kind == 'tss':
        adata = sg.tss_adata
    elif kind == 'tes':
        adata = sg.tes_adata
    elif kind == 'ic':
        adata = sg.ic_adata

    adata.X = adata.layers['counts']
    df = pd.DataFrame.sparse.from_spmatrix(data=adata.X,
        index=adata.obs['dataset'].tolist(),
        columns=adata.var.index.tolist())
    return df

def calc_tpm(adata, obs_col='dataset', how='mean'):
	"""
	Calculate the TPM per condition given by `obs_col`.
	Default column to use is `adata.obs` index column, `dataset`.

	Parameters:
		adata (anndata AnnData): Annotated data object from the SwanGraph
		obs_col (str or list of str): Column name from adata.obs table to group on.
			Default: 'dataset'
		how (str): How to compute tpm across multiple datasets {'mean', 'max'}
		recalc (bool): Whether tpm data should be recalculated or tpm
			values should just be averaged
			Default: False

	Returns:
		df (pandas DataFrame): Pandas datafrom where rows are the different
			conditions from `obs_col` and the columns are transcript ids in the
			SwanGraph, and values represent the TPM value per isoform per
			condition.
	"""

	# only need to calculate tpm once when adding abundance


	# otherwise just grab the tpm from the adata

	data = adata.layers['counts'].toarray()

	# turn into a dataframe
	cols = adata.var.index.tolist()
	inds = adata.obs[obs_col].tolist()
	df = pd.DataFrame(data=data, columns=cols, index=inds)
	df.index.name = obs_col

	# average across tpm
	if obs_col != 'dataset':

		# keep track of original row order to sort by
		row_order = df.index.unique().tolist()
		row_map = {}
		for i, row in enumerate(row_order):
			row_map[row] = i

		df['row_order'] = df.index.map(row_map)
		df.reset_index(inplace=True)
		if how == 'mean':
			df = df.groupby(obs_col).mean()
		elif how == 'max':
			df = df.groupby(obs_col).max()
		df = df.sort_values(by='row_order', ascending=True)
		df.drop('row_order', inplace=True, axis=1)

	# make sparse
	data = sparse.csr_matrix(df.values)
	inds = df.index.tolist()
	cols = df.columns.tolist()
	df = pd.DataFrame.sparse.from_spmatrix(data, index=inds, columns=cols)

	return df

def gen_report(sg,
                   gid,
                   prefix,
                   datasets=None,
                   groupby=None,
                   metadata_cols=None,
                   novelty=False,
                      layer='tpm', # choose from tpm, pi
                   cmap='Spectral_r',
                   include_qvals=False,
                   q=0.05,
                   qval_obs_col=None,
                   qval_obs_conditions=None,
                   include_unexpressed=False,
                   indicate_novel=False,
                   display_numbers=False,
                   transcript_col='tid',
                   browser=False,
                   order='expression'):

        # check if groupby column is present
        multi_groupby = False
        indicate_dataset = False
        if groupby:
            # grouping by more than one column
            if type(groupby) == list and len(groupby) > 1:
                for g in groupby:
                    if g not in sg.adata.obs.columns.tolist():
                        raise Exception('Groupby column {} not found'.format(g))
                groupby = sg.add_multi_groupby(groupby)
                multi_groupby = True
            elif groupby not in sg.adata.obs.columns.tolist():
                raise Exception('Groupby column {} not found'.format(groupby))


        # check if metadata columns are present
        if metadata_cols:
            for c in metadata_cols:
                if c not in sg.adata.obs.columns.tolist():
                    raise Exception('Metadata column {} not found'.format(c))

                # if we're grouping by a certain variable, make sure
                # the other metadata cols we plan on plotting have unique
                # mappings to the other columns. if just grouping by dataset,
                # since each dataset is unique, that's ok
                if groupby and groupby != 'dataset':
                    if groupby == c:
                        continue

                    temp = sg.adata.obs[[groupby, c, 'dataset']].copy(deep=True)
                    temp = temp.groupby([groupby, c]).count().reset_index()
                    temp = temp.loc[temp.dataset!=0]
                    temp = temp.loc[~temp.dataset.isnull()]

                    # if there are duplicates from the metadata column, throw exception
                    if temp[groupby].duplicated().any():
                            raise Exception('Metadata column {} '.format(c)+\
                                'not compatible with groupby column {}. '.format(groupby)+\
                                'Groupby column has more than 1 unique possible '+\
                                'value from metadata column.')

        # check to see if input gene is in the graph
        if gid not in sg.t_df.gid.tolist():
            gid = sg.get_gid_from_gname(gid)
        sg.check_gene(gid)

        # check to see if these plotting settings will play together
        sg.check_plotting_args(indicate_dataset,
            indicate_novel, browser)

        # check to see if transcript column is present in sg.t_df
        if transcript_col not in sg.t_df.columns:
            raise Exception('Transcript identifier column {} '.format(transcript_col)+\
                            'not present in SwanGraph.')

        # get the list of columns to include from the input datasets dict
        if datasets:
            # get a df that is subset of metadata
            # also sort the datasets based on the order they appear in "datasets"
            i = 0
            sorters = []
            for meta_col, meta_cats in datasets.items():
                if meta_col not in sg.adata.obs.columns.tolist():
                    raise Exception('Metadata column {} not found'.format(meta_col))
                if type(meta_cats) == str:
                    meta_cats = [meta_cats]
                if i == 0:
                    temp = sg.adata.obs.loc[sg.adata.obs[meta_col].isin(meta_cats)]
                else:
                    temp = temp.loc[temp[meta_col].isin(meta_cats)]
                sort_ind = dict(zip(meta_cats, range(len(meta_cats))))
                sort_col = '{}_sort'.format(meta_col)
                temp[sort_col] = temp[meta_col].map(sort_ind).astype(int)
                sorters.append(sort_col)
                i += 1

            # sort the df based on the order that different categories appear in "datasets"
            temp.sort_values(by=sorters, inplace=True, ascending=True)
            temp.drop(sorters, axis=1, inplace=True)
            columns = temp.dataset.tolist()
            col_order = temp[groupby].unique().tolist()
            del temp
        else:
            columns = None
            col_order = None

        # if we've asked for novelty first check to make sure it's there
        if novelty:
            if not sg.has_novelty():
                raise Exception('No novelty information present in the graph. '
                    'Add it or do not use the "novelty" report option.')

        # abundance info to calculate TPM on - subset on datasets that will
        # be included
        #if columns or datasets:
            #subset_adata =     subset_on_gene_sg(sg, datasets=columns).adata
        #else:
        subset_adata = sg.adata

        # small SwanGraph with only this gene's data
        #sg =     subset_on_gene_sg(sg, gid=gid, datasets=columns)
        
        # if we're grouping data, calculate those new numbers
        # additionally order transcripts
        if groupby:
            if layer == 'tpm':
                # use whole adata to calc tpm
                t_df = tpm_df = calc_tpm(subset_adata, obs_col=groupby).transpose()

            elif layer == 'pi':
                # calc tpm just so we can order based on exp
                tpm_df = calc_tpm(subset_adata, obs_col=groupby).transpose()
                t_df, _ = calc_pi(sg.adata, sg.t_df, obs_col=groupby)
                t_df = t_df.transpose()

        else:
            if layer == 'tpm':
                # use whole adata to calc tpm
                t_df = tpm_df = get_tpm(sg).transpose()
                t_df = t_df[subset_adata.obs.dataset.tolist()]
                

        # order transcripts by user's preferences
        if order == 'pi' and not layer == 'pi':
            order = 'expression'

        if order == 'expression' and sg.abundance == False:
            order = 'tid'
        elif order == 'expression':
            order = 'log2tpm'
        tids = sg.t_df.loc[sg.t_df.gid == gid].index.tolist()
        
        tids = list(set(tids)&set(tpm_df.index.tolist()))
        tpm_df = tpm_df.loc[tids]
        t_df = t_df.loc[tids]

        if order == 'log2tpm':
            order_df = tpm_df
        elif order == 'pi':
            order_df = t_df
        else:
            order_df = t_df
        _, tids = sg.order_transcripts_subset(order_df, order=order)
        tpm_df = tpm_df.loc[tids]
        t_df = t_df.loc[tids]

    #     print('finished ordering transcripts')
        del tpm_df

        # remove unexpressed transcripts if desired
        if not include_unexpressed:
            t_df = t_df.loc[t_df.any(axis=1)]

        # make sure de has been run if needed
        qval_df = None

        # get tids in this report
        report_tids = t_df.index.tolist()

        # plot each transcript with these settings
        print()
        print('Plotting transcripts for {}'.format(gid))
        sg.plot_each_transcript(report_tids, prefix,
                                  browser=browser)
    #     print('finished plotting each transcript')

        # get a different prefix for saving colorbars and scales
        gid_prefix = prefix+'_{}'.format(gid)

        # if we're plotting tracks, we need a scale as well
        # also set what type of report this will be, 'swan' or 'browser'
        if browser:
            sg.pg.plot_browser_scale()
            save_fig(gid_prefix+'_browser_scale.png')
            report_type = 'browser'
        else:
            report_type = 'swan'

        # plot colorbar for either tpm or pi
        if layer == 'tpm':

            # take log2(tpm) (add pseudocounts)

            # min and max tpm vals
            g_max = t_df.max().max()
            g_min = t_df.min().min()

            # create a colorbar
            plt.rcParams.update({'font.size': 30})
            fig, ax = plt.subplots(figsize=(14, 1.5))
            fig.subplots_adjust(bottom=0.5)
            fig.patch.set_visible(False)
            ax.patch.set_visible(False)

            try:
                cmap = plt.get_cmap(cmap)
            except:
                raise ValueError('Colormap {} not found'.format(cmap))

            norm = mpl.colors.Normalize(vmin=g_min, vmax=g_max)

            cb = mpl.colorbar.ColorbarBase(ax,
                                cmap=cmap,
                                norm=norm,
                                orientation='horizontal')
            cb.set_label('transformed portion of total gene expression')
            plt.savefig(gid_prefix+'_colorbar_scale.png', format='png',
                bbox_inches='tight', dpi=200)
            plt.clf()
            plt.close()

        elif layer == 'pi':

            # min and max pi vals
            g_max = 100
            g_min = 0

            # create a colorbar between 0 and 1
            plt.rcParams.update({'font.size': 30})
            fig, ax = plt.subplots(figsize=(14, 1.5))
            fig.subplots_adjust(bottom=0.5)
            fig.patch.set_visible(False)
            ax.patch.set_visible(False)

            try:
                cmap = plt.get_cmap(cmap)
            except:
                raise ValueError('Colormap {} not found'.format(cmap))

            norm = mpl.colors.Normalize(vmin=0, vmax=100)

            cb = mpl.colorbar.ColorbarBase(ax,
                                cmap=cmap,
                                norm=norm,
                                orientation='horizontal')
            cb.set_label('Percent of isoform use (' +'$\pi$'+')')
            plt.savefig(gid_prefix+'_colorbar_scale.png', format='png',
                bbox_inches='tight', dpi=200)
            plt.clf()
            plt.close()

        # merge with sg.t_df to get additional columns
        if not col_order:
            datasets = t_df.columns
        cols = ['novelty', transcript_col]
        t_df = t_df.merge(sg.t_df[cols], how='left', left_index=True, right_index=True)

        #vc = t_df[t_df.columns[0:len(datasets)]]
        #try:
        #    nc = (vc/vc.sum())
        #except: 
        #    nc = vc.sum()
        #nc = nc.astype('float')

        #t_df[t_df.columns[0:len(datasets)]] = nc

        # order according to input specifications
        if col_order:
            col_order = [c for c in col_order if c in t_df.columns]
            datasets = col_order.copy()
            col_order+=list(set(t_df.columns)-set(col_order))
            t_df = t_df[col_order]

        # create report
        print('Generating report for {}'.format(gid))
        pdf_name = create_fname(prefix,
                     indicate_dataset,
                     indicate_novel,
                     browser,
                     ftype='report',
                     gid=gid)
        if transcript_col == 'tid':
            t_disp = 'Transcript ID'
        else:
            t_disp = 'Transcript Name'

        report = Report(gid_prefix,
                        report_type,
                        sg.adata.obs,
                        sg.adata.uns,
                        datasets=datasets,
                        groupby=groupby,
                        metadata_cols=metadata_cols,
                        novelty=novelty,
                        layer=layer,
                        cmap=cmap,
                        g_min=g_min,
                        g_max=g_max,
                        include_qvals=include_qvals,
                        qval_df=qval_df,
                        display_numbers=display_numbers,
                        t_disp=t_disp)
        report.add_page()

        # loop through each transcript and add it to the report


        
        for ind, entry in t_df.iterrows():
            tid = ind
            # display name for transcript
            if transcript_col == 'tid':
                t_disp = tid
            else:
                t_disp = entry[transcript_col]
                # t_disp = entry['tname']
            # else:
            #     t_disp = tid
            fname = create_fname(prefix,
                                 indicate_dataset,
                                 indicate_novel,
                                 browser,
                                 ftype='path',
                                 tid=tid)
            report.add_transcript(entry, fname, t_disp)
        report.write_pdf(pdf_name)

        # remove multi groupby column if necessary
        if multi_groupby:
            sg.rm_multi_groupby(groupby)

##########################################################################
############################# Data retrieval #############################
##########################################################################

sg.t_df['human_readable'] = sg.t_df['tname']

sg.t_df.loc[sg.t_df['tname'].str.contains('FLAIR'), 'human_readable'] = 'predicted'

sg.t_df = sg.t_df.drop_duplicates(subset='tid')
try:
	os.mkdir(bc+'/')
except:
	print()
for gene in cand_l:
    try:
        gen_report(sg, gene,
                   prefix=bc+'/'+gene,
                   cmap='viridis',
                   display_numbers=True,
                   novelty = False,
                   groupby='cond',
                   transcript_col='human_readable',
                   browser = True)
    except:
        print()


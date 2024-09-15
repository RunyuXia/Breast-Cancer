import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata as ad

import glob
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, fisher_exact, pearsonr
import matplotlib.gridspec as gridspec
import matplotlib as mpl

import matplotlib as mpl
import palettable
import matplotlib.patches as patches

def process(adata, ct):
    adata = adata[adata.obs['CellType_new'] == ct, :]
    hvg = open('/home/rx238/rds/hpc-work/aging/NMF_lym/counts/' + ct + '_HVGs.txt').read().split('\n')
    adata.X = adata.layers['raw'].copy()
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=10**4) ## TPT normalization
    adata1 = adata.copy()
    sc.pp.log1p(adata1)
    adata.raw = adata1
    hvg = hvg[:-1]
    adata = adata[:,hvg]
    sc.pp.scale(adata)
    sc.pp.neighbors(adata, use_rep='scvi', n_neighbors=15)
    sc.tl.umap(adata)
    
    return adata


def graphs(adata, ct, k, programs):
    
    gene_scores_file = pd.read_csv('/home/rx238/rds/hpc-work/aging/NMF_lym/cNMF/' + ct + '_cNMF/' + ct + '_cNMF.gene_spectra_score.k_' + k + '.dt_0_2.txt', delimiter='\t', index_col=0)
    gene_tpm_file = pd.read_csv('/home/rx238/rds/hpc-work/aging/NMF_lym/cNMF/' + ct + '_cNMF/tpm/' + ct + '_cNMF.gene_spectra_tpm.k_' + k + '.dt_0_2.txt', delimiter='\t', index_col=0)
    usage_file = pd.read_csv('/home/rx238/rds/hpc-work/aging/NMF_lym/cNMF/' + ct + '_cNMF/usage/' + ct + '_cNMF.usages.k_' + k + '.dt_0_2.consensus.txt', delimiter='\t', index_col=0)
    
    selected_usages = usage_file[programs]
    selected_usages.columns = ['Usage_%s' % i for i in selected_usages.columns]
    adata1 = adata.copy()
    adata1.obs = pd.merge(left=adata1.obs, right=selected_usages, how='left', left_index=True, right_index=True)
    sc.pl.umap(adata1, color=selected_usages.columns, use_raw=True, ncols=3, vmin=0, vmax='p99', save = '/' + ct + '_' + k + '_usage_all_new.png')
    df_transposed = gene_scores_file.T
    df_transposed.columns = [str(int(f)) for f in df_transposed.columns.to_list()]
    df_transposed = df_transposed[programs]
    df_transposed.to_csv('/home/rx238/rds/hpc-work/aging/NMF_lym_result/transposed/' + ct + k + '_gene_scores_per_program_all.csv', index=True, header=True)

    top_100 = df_transposed.apply(lambda x: x.nlargest(100).index)
    df_top_100_genes = pd.DataFrame(top_100)
    df_top_100_genes = df_top_100_genes[programs]
    df_top_100_genes.to_csv('/home/rx238/rds/hpc-work/aging/NMF_lym_result/top_100/' + ct + k + '_top_100_genes_per_program_all.csv', index=True, header=True)
    
    
def graph_batch(adata, ct, k, programs):
    print('Batch')
    base_dir = '/home/rx238/rds/hpc-work/aging/NMF_lym_result/results'  # Example cell type, replace with your variable ct
    directory = os.path.join(base_dir, ct)
# Check if the directory exists, and create it if it does not
    if not os.path.exists(directory):
        os.makedirs(directory)

    usage_matrix =pd.read_csv('/home/rx238/rds/hpc-work/aging/NMF_lym/cNMF/' + ct + '_cNMF/usage/' + ct + '_cNMF.usages.k_' + k + '.dt_0_2.consensus.txt', delimiter='\t', index_col=0)
    normalized_usage_matrix = usage_matrix.div(usage_matrix.sum(axis=1), axis=0)
    normalized_usage_matrix = normalized_usage_matrix[programs]


    clusters = adata.obs.copy()
    usage_by_cluster = normalized_usage_matrix.groupby(clusters['Batch'], axis=0).apply(lambda x: x.mean())

    batches = usage_by_cluster.index.tolist()
    n_batch = len(batches)

    programs_int = list(map(int, programs))
    gep_order = np.array(programs_int)
    num_programs = len(programs)
    
    width_per_program = 0.2
    fig_width = width_per_program * num_programs
    fig_height = 2  # Keep the height constant, adjust if necessary

    (fig, ax) = plt.subplots(1, 1, figsize=(fig_width, fig_height), dpi=500)
    
    sns.heatmap(usage_by_cluster, ax=ax)
    ax.set_title('GEP Usage by Batch')
    ax.set_ylabel('Batch')
    ax.set_xlabel('Inferred GEP')
    ax.set_yticks(np.arange(n_batch)+.5)
    
    ax.set_xticks(np.arange(num_programs)+.5)
    ax.set_xticklabels(programs, rotation=45, ha='right', fontsize=6)  # Set x-tick labels

    
    fig1_path = os.path.join(directory, f'{ct}_GEP_by_batch_{k}.dt_0_2.png')
    fig.savefig(fig1_path, dpi=600, bbox_inches='tight')

def graph_group(adata, ct, k, programs):
    
    base_dir = '/home/rx238/rds/hpc-work/aging/NMF_lym_result/results'  # Example cell type, replace with your variable ct
    directory = os.path.join(base_dir, ct)
# Check if the directory exists, and create it if it does not
    if not os.path.exists(directory):
        os.makedirs(directory)

    usage_matrix =pd.read_csv('/home/rx238/rds/hpc-work/aging/NMF_lym/cNMF/' + ct + '_cNMF/usage/' + ct + '_cNMF.usages.k_' + k + '.dt_0_2.consensus.txt', delimiter='\t', index_col=0)
    normalized_usage_matrix = usage_matrix.div(usage_matrix.sum(axis=1), axis=0)
    normalized_usage_matrix = normalized_usage_matrix[programs]


    clusters = adata.obs.copy()
    
    bins = [0, 30, 60, 90, 120]
    bin_labels = ['0-30', '30-60', '60-90', '90-120']
    clusters['Age_bin'] = pd.cut(clusters['Age'], bins, labels=bin_labels)
    
    clusters['Age_bin'] = clusters['Age_bin'].cat.add_categories(['Parous'])
    clusters['Age_bin'] = clusters['Age_bin'].cat.add_categories(['WKBR'])
    
    clusters.loc[clusters['Parity'] == 'Parous', 'Age_bin'] = 'Parous'
    clusters.loc[clusters['Genotype'] == 'WKBR', 'Age_bin'] = 'WKBR'
    
    
    usage_by_cluster = normalized_usage_matrix.groupby(clusters['Age_bin'], axis=0).apply(lambda x: x.mean())

    batches = usage_by_cluster.index.tolist()
    n_batch = len(batches)

    programs_int = list(map(int, programs))
    gep_order = np.array(programs_int)
    num_programs = len(programs)
    
    width_per_program = 0.2
    fig_width = width_per_program * num_programs
    fig_height = 2  # Keep the height constant, adjust if necessary

    (fig, ax) = plt.subplots(1, 1, figsize=(fig_width, fig_height), dpi=500)
    
    sns.heatmap(usage_by_cluster, ax=ax)
    ax.set_title('GEP Usage by Group')
    ax.set_ylabel('Group')
    ax.set_xlabel('Inferred GEP')
    ax.set_yticks(np.arange(n_batch)+.5)
    
    ax.set_xticks(np.arange(num_programs)+.5)
    ax.set_xticklabels(programs, rotation=45, ha='right', fontsize=6)  # Set x-tick labels
    
    
    fig1_path = os.path.join(directory, f'{ct}_GEP_by_group_{k}.dt_0_2.png')
    fig.savefig(fig1_path, dpi=600, bbox_inches='tight')


def gep_graph1(adata, ct, k, programs):
    print('Groups')
    base_dir = '/home/rx238/rds/hpc-work/aging/NMF_lym_result/results'  # Example cell type, replace with your variable ct
    directory = os.path.join(base_dir, ct)
# Check if the directory exists, and create it if it does not
    if not os.path.exists(directory):
        os.makedirs(directory)

    usage_matrix =pd.read_csv('/home/rx238/rds/hpc-work/aging/NMF_lym/cNMF/' + ct + '_cNMF/usage/' + ct + '_cNMF.usages.k_' + k + '.dt_0_2.consensus.txt', delimiter='\t', index_col=0)
    normalized_usage_matrix = usage_matrix.div(usage_matrix.sum(axis=1), axis=0)
    file_path = os.path.join(directory, f'{ct}_cNMF.usages.k_{k}.dt_0_2.consensus_normalized.csv')
    normalized_usage_matrix = normalized_usage_matrix[programs]
    normalized_usage_matrix.to_csv(file_path, index=True, header=True)

    clusters = adata.obs.copy()
    usage_by_cluster = normalized_usage_matrix.groupby(clusters['CellType_new'], axis=0).apply(lambda x: (x>.25).mean())

    cts = usage_by_cluster.index.tolist()
    n_ct = len(cts)

    programs_int = list(map(int, programs))
    gep_order = np.array(programs_int)

    (fig,ax) = plt.subplots(1,1,figsize=(3,2), dpi=500)
    sns.heatmap(usage_by_cluster, ax=ax)
    ax.set_title('% Cells With GEP Usage > 25%\nby Published Cluster')
    ax.set_ylabel('Clustered Cell-Type')
    ax.set_xlabel('Inferred GEP')
    ax.set_yticks(np.arange(n_ct)+.5)

    fig1_path = os.path.join(directory, f'{ct}_GEP_percent_{k}.dt_0_2.png')
    fig.savefig(fig1_path, dpi=600, bbox_inches='tight')

    bins = [0, 30, 60, 90, 120]
    bin_labels = ['0-30', '30-60', '60-90', '90-120']
    clusters['Age_bin'] = pd.cut(clusters['Age'], bins, labels=bin_labels)
    
    clusters['Age_bin'] = clusters['Age_bin'].cat.add_categories(['Parous'])
    clusters['Age_bin'] = clusters['Age_bin'].cat.add_categories(['WKBR'])
    
    clusters.loc[clusters['Parity'] == 'Parous', 'Age_bin'] = 'Parous'
    clusters.loc[clusters['Genotype'] == 'WKBR', 'Age_bin'] = 'WKBR'
    
    usage_info_ct = pd.merge(left=normalized_usage_matrix, right=clusters[['CellType_new', 'Age', 'Age_bin']], left_index=True, right_index=True)

    label_size = 8

    core_colors = type('CoreColors', (), {})
    cnames = ['red', 'blue', 'green', 'purple', 'orange', 'yellow', 'brown', 'pink', 'grey']

    def to_array_col(color):
        return np.array(color)/255.

    for cname,c in zip(cnames, palettable.colorbrewer.qualitative.Set1_9.colors):
        setattr(core_colors, cname, np.array(c)/255.)
    
    for cname, c in zip(['blue', 'green', 'red', 'orange', 'purple', 'brown', 'pink', 'grey'],
                        palettable.colorbrewer.qualitative.Paired_12.colors[::2]):
        setattr(core_colors, 'pale_'+cname, np.array(c)/255.)

    base_height_per_subplot = 0.5  # Height per subplot in inches
    fig_width = 3  # Width of the figure in inches
    total_fig_height = len(programs) * base_height_per_subplot  # Total height of the figure

# Create the figure with dynamic size based on the number of programs
    fig = plt.figure(figsize=(fig_width, total_fig_height), dpi=600)

    n_prog = len(programs)

    gs = gridspec.GridSpec(n_prog, 1, fig,
                    0.20, 0.11, 0.98, 0.92,
                    hspace=.3)
    
    col = 'Big GEP Group'
    used_groups = cts
    used_stim = ['0-30', '30-60', '60-90', '90-120', 'Parous', 'WKBR']

    clabels = cts
    nclusters = len(clabels)

    activity_programs = programs
        
    position = (np.arange(len(clabels))%6 - 1)*0.3
    n_groups = len(used_groups)*len(used_stim)
    positions = np.arange(n_groups) + (np.arange(n_groups)%6 - 1)*(-0.3)

    stim_facecolor = {'0-30': core_colors.pale_blue, '30-60': core_colors.pale_green, '60-90': core_colors.pale_orange, '90-120': core_colors.pale_red, 'Parous': core_colors.pale_purple, 'WKBR': core_colors.pale_brown}
    stim_edgecolor = {'0-30': core_colors.blue, '30-60': core_colors.green, '60-90': core_colors.orange,'90-120': core_colors.red, 'Parous': core_colors.purple, 'WKBR': core_colors.brown}
    subdat = usage_info_ct


    num_rows = k
    current_row = 0

    for (i,act_gep) in enumerate(activity_programs):
        current_row += 1
        
        ax = fig.add_subplot(gs[i,0], frameon=True)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        
        grp_def = []
        grouped_usages = []
        for grp in used_groups:
            for stim in used_stim:
                grp_def += [(grp, stim)]
                grouped_usages.append(subdat[(subdat['Age_bin']==stim) & (subdat['CellType_new']==grp)][act_gep].values)

        bp = ax.boxplot(grouped_usages, sym='o', notch=False, zorder=10, whis=[5,95], positions=positions,
                        patch_artist=True)
        for el in ['whiskers', 'caps']:
            for li,line in enumerate(bp[el]):
                line.set_color(stim_edgecolor[grp_def[int(np.floor(li/2))][1]])
                line.set_linewidth(0.5)
                line.set_zorder(-10)
                line.set_clip_on(False)
                
        for li,line in enumerate(bp['boxes']):
            line.set_linewidth(0.5)
            line.set_zorder(-10)
            line.set_facecolor(stim_facecolor[grp_def[li][1]])
            line.set_edgecolor(stim_edgecolor[grp_def[li][1]])
            line.set_clip_on(False)
        for li,line in enumerate(bp['medians']):
            line.set_zorder(-1)
            line.set_color(stim_edgecolor[grp_def[li][1]])
        for li,line in enumerate(bp['fliers']):
            line.set_markeredgecolor('none')
            line.set_markerfacecolor(stim_facecolor[grp_def[li][1]])
            line.set_alpha(0.3)
            line.set_markersize(1)
            line.set_rasterized(True)
        ax.set_xticks([])
        ax.set_yticks([0, .25, .5, .75, 1])
        ax.set_yticklabels([0, 25, 50, 75, 100], fontsize=3)
        ax.set_ylim([0, 1])
        ax.set_xlim(-0.5, len(grp_def)-0.5)
        for pi in range(1, len(used_groups), 2):
            ax.add_patch(
                patches.Rectangle(
                    (pi*6-0.5, 0.), 6.0, 1, 
                    facecolor="#f0f0f0", edgecolor='none', zorder=-20,
                    )
                )
            
        ax.set_ylabel(act_gep, fontsize=7)
        
        if current_row == num_rows:
                # ax.set_xticks(range(1, len(grp_def), 5))
                # ax.set_xticklabels(used_groups, rotation=45, ha='right',fontsize=5)
                ax.set_xlabel(ct, fontsize=7)
        if current_row == 1:
            ax.set_title(ct + ' GEP usage (%)', fontsize=7)
            
            stim_legend = [patches.Patch(facecolor=stim_facecolor['0-30'], linewidth=0.5,
                                        edgecolor=stim_edgecolor['0-30'], label='0-30 wks'),
                        patches.Patch(facecolor=stim_facecolor['30-60'], linewidth=0.5,
                                    edgecolor= stim_edgecolor['30-60'], label='30-60 wks'),
                        patches.Patch(facecolor=stim_facecolor['60-90'], linewidth=0.5,
                                    edgecolor= stim_edgecolor['60-90'], label='60-90 wks'),
                        patches.Patch(facecolor=stim_facecolor['90-120'], linewidth=0.5,
                                    edgecolor= stim_edgecolor['90-120'], label='90-120 wks'),
                        patches.Patch(facecolor=stim_facecolor['Parous'], linewidth=0.5,
                                    edgecolor= stim_edgecolor['Parous'], label='Parous'),
                        patches.Patch(facecolor=stim_facecolor['WKBR'], linewidth=0.5,
                                    edgecolor= stim_edgecolor['WKBR'], label='WKBR')]
            legend = ax.legend(handles=stim_legend, fontsize=3, loc=(.85,.95))
            legend.get_frame().set_linewidth(0.5)

    fig2_path = os.path.join(directory, f'{ct}_GEP_ct_{k}.png')
    fig.savefig(fig2_path, dpi=600, bbox_inches='tight')




sc.settings.figdir = '/home/rx238/rds/hpc-work/aging/NMF_lym_result/'

data = ad.read_h5ad('/home/rx238/rds/hpc-work/aging/final_data/adata_all_cells_renamed_add_lp_avd.h5ad')

data = data[(data.obs['Genotype'].isin(['WT', 'WKBR'])) & (data.obs['Parity'].isin(['NP', 'Parous'])) & (~data.obs["CellType_new"].isin(['Erythoid', 'Erythroid', 'Fb3 Uni2_2', 'Doublets', 'Caf', 'Tam Macrophage', 'Hs Nul2_2', 'Hs Uni2_1', 'Tumour', 'Myocytes'])) & (data.obs['Batch'] != 7) & (~data.obs['Sample'].isin(['ctrl2', 'ctrl1', 'NulControl'])),:]


k_values = pd.read_csv('/home/rx238/rds/hpc-work/aging/NMF_lym_result/k_selected.csv')
selected_p = pd.read_csv('/home/rx238/rds/hpc-work/aging/NMF_lym_result/all_p_selected.csv')

cts = k_values.columns.to_list()
cts = [name for name in cts if name not in ('Avd', 'Lp')]


k_dict = {}
program_dict = {}
object_dict = {}

for ct in cts:
    print(ct)
    adata_new = process(data, ct)
    k = k_values.at[0, ct]
    programs = selected_p[ct].dropna().to_list()
    object_dict[ct] = adata_new
    k_dict[ct] = str(k)
    filtered_programs_new = [str(int(item)) for item in programs]
    program_dict[ct] = filtered_programs_new


for ct in cts:
    adata = object_dict[ct]
    k = k_dict[ct]
    programs = program_dict[ct]
    print(ct)
    graphs(adata, ct, k, programs)
    gep_graph1(adata, ct, k, programs)
    graph_batch(adata, ct, k, programs)
    graph_group(adata, ct, k, programs)





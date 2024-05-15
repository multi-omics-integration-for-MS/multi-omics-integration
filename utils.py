import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.offline as pyo
import scanpy as sc
import seaborn as sns


def sankey_plot(
        labels,
        labels_titles=None,
        title=None,
        color_palette=sns.color_palette(),
        width=None,
        height=None
    ):
    '''
    This function plots a Sankey diagram of the sets of labels passed as arguments.

    :param labels1: list of labels list
    :param labels2: lables titles
    :param title: title of the plot
    :param color_palette: color palette to use
    '''

    sorted_pairs = sorted(zip(labels[0], labels[1]))
    sorted0, sorted1 = zip(*sorted_pairs)
    labels = [pd.Series(sorted0), pd.Series(sorted1)]

    n_types = [len(set(label_list)) for label_list in labels]

    plot_labels = []
    for i in labels:
        for j in i.unique():
            plot_labels.append(j)

    source = []
    target = []
    value = []
    for i in range(len(labels)-1):
        confusion_matrix = pd.crosstab(labels[i], labels[i+1])
        curr_source = []
        curr_target = []
        curr_value = []

        source_add = 0
        for j in range(0, i):
            source_add += n_types[j]
        target_add = source_add + n_types[i]

        for j in range(n_types[i]):
            for k in range(n_types[i+1]):
                if confusion_matrix.iloc[j, k] != 0:
                    curr_source.append(j+source_add)
                    curr_target.append(k+target_add)
                    curr_value.append(confusion_matrix.iloc[j, k])

        source += curr_source
        target += curr_target
        value += curr_value

    colors = []
    for i in range(len(labels)):
        colors += color_palette.as_hex()[:n_types[i]]

    fig = go.Figure(
        data=[
            go.Sankey(
                node = dict(
                    pad = 15,
                    thickness = 20,
                    line = dict(color = "black", width = 0.5),
                    color = colors
                ),
                link = dict(
                    source = source,
                    target = target,
                    value = value
                )
            )
        ]
    )

    for x_coordinate, column_name in enumerate(labels_titles):
        fig.add_annotation(
            x=x_coordinate,
            y=1.05,
            xref="x",
            yref="paper",
            text=column_name,
            showarrow=False
        )
    fig.update_layout(
        title_text=title, 
        xaxis={'showgrid': False, 'zeroline': False, 'visible': False},
        yaxis={'showgrid': False, 'zeroline': False, 'visible': False},
        plot_bgcolor='rgba(0,0,0,0)',
        font_size=8,
        width=width,
        height=height
    )

    """file_name = f'../sankey'
    if title is not None:
        camel_title = title.replace(' ', '_')
        file_name += f'_{camel_title}'
    file_name += '.html'
    pyo.plot(fig, filename=file_name, auto_open=False)"""
    fig.show()

def sankey_plot_with_labels(
        labels,
        labels_titles,
        title,
        path=None,
        colored_links=False,
        link_opacity=0.4,
        width=700,
        height=450,
    ):
    '''
    This function plots a Sankey diagram of the sets of labels passed as arguments.

    :param labels: list of labels list
    :param labels_titles: lables titles
    :param path: path to save the plot
    :param title: title of the plot
    '''

    n_clusters = [len(set(label_list)) for label_list in labels]

    plot_labels = []
    for i in range(len(labels)):
        plot_labels += np.unique(labels[i]).tolist()

    # Generate color palette for sankey nodes
    node_palette = sns.color_palette(None, len(plot_labels))
    link_palette = [f'rgba({int(r*255)}, {int(g*255)}, {int(b*255)}, {link_opacity})' for r, g, b in node_palette]

    source = []
    target = []
    value = []
    for i in range(len(labels)-1):
        confusion_matrix = pd.crosstab(labels[i], labels[i+1])
        curr_source = []
        curr_target = []
        curr_value = []

        source_add = 0
        for j in range(0, i):
            source_add += n_clusters[j]
        target_add = source_add + n_clusters[i]

        for j in range(n_clusters[i]):
            for k in range(n_clusters[i+1]):
                if confusion_matrix.iloc[j, k] != 0:
                    curr_source.append(j+source_add)
                    curr_target.append(k+target_add)
                    curr_value.append(confusion_matrix.iloc[j, k])

        source += curr_source
        target += curr_target
        value += curr_value

    fig = go.Figure(
        data=[
            go.Sankey(
                node = dict(
                    pad = 15,
                    thickness = 20,
                    line = dict(color = "black", width = 0.5),
                    label = plot_labels,
                    color = node_palette.as_hex()
                ),
                link = dict(
                    source = source,
                    target = target,
                    value = value,
                    color = [link_palette[i] for i in source] if colored_links else None
                )
            )
        ]
    )

    for x_coordinate, column_name in enumerate(labels_titles):
        fig.add_annotation(
            x=x_coordinate,
            y=1.05,
            xref="x",
            yref="paper",
            text=column_name,
            showarrow=False
        )
    fig.update_layout(
        title_text=title, 
        xaxis={'showgrid': False, 'zeroline': False, 'visible': False},
        yaxis={'showgrid': False, 'zeroline': False, 'visible': False},
        plot_bgcolor='rgba(0,0,0,0)',
        font_size=10,
        width=width,
        height=height
    )

    fig.show()
    if path is not None:
        pyo.plot(fig, filename=path, auto_open=False)
        fig.write_image(path.replace('.html', '.svg'))

def visualize_p_value_adj(adata_cell_type, cell_type_name, method='t-test'):
    fig, axs = plt.subplots(2, 1, figsize=(18, 6))

    axs[0].plot(adata_cell_type.uns[method]['pvals_adj']['MS'], label='MS cells')
    axs[0].plot([0, len(adata_cell_type.uns[method]['pvals_adj']['MS'])], [0.05, 0.05], 'r--')
    axs[0].set_xlabel('Genes ranked by score on MS cells')
    axs[0].set_ylabel('p-value adjusted')
    axs[0].set_xticks([])
    axs[0].legend()

    axs[1].plot(adata_cell_type.uns[method]['pvals_adj']['HC'], label='HC cells')
    axs[1].plot([0, len(adata_cell_type.uns[method]['pvals_adj']['HC'])], [0.05, 0.05], 'r--')
    axs[1].set_xlabel('Genes ranked by score on HC cells')
    axs[1].set_ylabel('p-value adjusted')
    axs[1].set_xticks([])

    plt.legend()
    plt.suptitle(f'{cell_type_name} p-values adjusted for multiple testing')
    plt.tight_layout()
    plt.show()

def dotplots_and_ranking_most_significant_genes(adata_cell_type, cell_type_name, n_genes=50, method='t-test'):
    sc.pl.rank_genes_groups(adata_cell_type, n_genes=n_genes, sharey=True, ncols=2, fontsize=6)

    sc.pl.dotplot(adata_cell_type, var_names=adata_cell_type.uns[method]['names']['MS'][:n_genes],
        groupby='MS/HC', expression_cutoff=0.1, title=f'{cell_type_name} most significant genes in MS cells (MS vs HC)')

    sc.pl.dotplot(adata_cell_type, var_names=adata_cell_type.uns[method]['names']['HC'][:n_genes],
        groupby='MS/HC', expression_cutoff=0.1, title=f'{cell_type_name} most significant genes in HC cells (HC vs MS)')

    sc.pl.rank_genes_groups_matrixplot(adata_cell_type, n_genes=20, key=method, groupby='MS/HC',
        title=f'{cell_type_name} most significant genes in MS cells')

def visualize_rank_genes_groups_violin(adata_cell_type, cell_type_name, method='t-test'):
    with warnings.catch_warnings():
        fig, axs = plt.subplots(1, 2, figsize=(15, 5))
        sc.pl.rank_genes_groups_violin(adata_cell_type, groups=['MS'], n_genes=10, key=method, show=False, ax=axs[0])
        axs[0].set_title('Genes most expressed in MS cells')

        sc.pl.rank_genes_groups_violin(adata_cell_type, groups=['HC'], n_genes=10, key=method, show=False, ax=axs[1])
        axs[1].set_title('Genes most expressed in HC cells')

        plt.suptitle(f'{cell_type_name} expression in MS and HC')
        plt.show()

        sc.pl.rank_genes_groups_stacked_violin(adata_cell_type, n_genes=20, groupby='MS/HC', key=method, figsize=(15, 3))

def visualize_venn_diagram_ttest_vs_wilcoxon(adata_cell_type, cell_type_name):
    wc_ms = sc.get.rank_genes_groups_df(adata_cell_type, group='MS', key='wilcoxon', pval_cutoff=0.01, log2fc_min=0)['names']
    wc_hc = sc.get.rank_genes_groups_df(adata_cell_type, group='HC', key='wilcoxon', pval_cutoff=0.01, log2fc_min=0)['names']
    tt_ms = sc.get.rank_genes_groups_df(adata_cell_type, group='MS', key='t-test', pval_cutoff=0.01, log2fc_min=0)['names']
    tt_hc = sc.get.rank_genes_groups_df(adata_cell_type, group='HC', key='t-test', pval_cutoff=0.01, log2fc_min=0)['names']

    from matplotlib_venn import venn2
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    venn2([set(tt_ms), set(wc_ms)], set_labels=('MS genes', 'HC genes'), ax=axs[0])
    axs[0].set_title('T-test')
    venn2([set(tt_hc), set(wc_hc)], set_labels=('MS genes', 'HC genes'), ax=axs[1])
    axs[1].set_title('Wilcoxon')

    plt.suptitle(f'Most expressed genes in {cell_type_name} (p-value < 0.01)')
    plt.tight_layout()
    plt.show()
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
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
                    #label = plot_labels, #FIXME: labels_titles
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
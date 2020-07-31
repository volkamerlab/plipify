"""
fp_visual.py

The visualization module for the calculated fingerprints.
This module takes the fingerprint in dataframe, processes the data further and creates mutliple different visualizations:
- Interactive stacked bar chart
- Heatmap
- Colour coded fingerprint table
"""

import matplotlib.pyplot as plt
import plotly.graph_objects as go
import seaborn as sns
import pandas as pd
from IPython.core.display import display, HTML
import ipywidgets as widgets
from ipywidgets import HBox

from .plip_fingerprints import divide_list, read_residues

interaction_colours = {
    "hydrophobic": "#33ccff",
    "hbond": "#ff6699",
    "waterbridge": "#4d4dff",
    "saltbridge": "#ff3300",
    "pistacking": "#00cc7a",
    "pication": "#cc66ff",
    "halogen": "#ace600",
    "metal": "#ff9933",
}


def fingerprint_barplot(fp_df):
    """
    Visualize fingerprint as barplot either based on count or fp_type.

    Parameters
    ----------
    fp_data = fingerpint in dataframe form

    """
    fig = go.Figure(
        data=[
            go.Bar(name=interaction, x=list(fp_df.index), y=list(fp_df[interaction]))
            for interaction in sorted(fp_df.columns, reverse=True)
        ],
    )
    # Change the bar mode
    fig.update_layout(
        barmode="stack",
        title_text="Residue Interactions",
        yaxis_title="Interactions",
        xaxis_title="Residues",
    )
    fig.show()


def fingerprint_heatmap(fp_df):
    """
    Visualize fingerprint as heatmap either based on count or fp_type.

    Parameters
    ----------
    fp_df = fingerpint in dataframe form

    """
    fig, ax = plt.subplots(figsize=(10, 7))  # plot size
    sns.heatmap(fp_df.T, annot=True, cmap="YlGnBu", ax=ax)
    plt.xlabel("Residues")
    plt.ylabel("Interaction Types")


def prepare_tabledata(fp_df):
    """
    Create interaction index dictionary for fingerprint table.
    The keys of this dictionary are the interaction types and the values for each interaction type are a list of the positions (indices) of the fingerprint array that represent this interaction type.

    Parameters
    ----------
    fp_df = fingerpint in dataframe form

    """
    residues = list(fp_df.index)
    interaction_types = list(fp_df.columns)
    res_fp = fp_df.values.tolist()
    fp_id = range(0, len(residues) * len(interaction_types))  # all fp indices
    interaction_list = interaction_types * len(residues)
    interaction_index = dict(zip(fp_id, interaction_list))
    return res_fp, interaction_index, residues


def cell_colour(fp_index, interaction_index):
    """
    Extract cell colour for specific residue interaction.

    Parameters
    ----------
    fp_index = the specific fingerprint position for which the colour is extracted
    interaction_index = dictionary of interaction types and their corresponding indices in the fingerprint

    """
    interaction_type = interaction_index[fp_index]
    interaction_colour = interaction_colours[interaction_type]
    return interaction_colour, interaction_type


def fingerprint_table(fp_df):
    """
    Create HTML and CSS table layout for the calculated fingerprint.

    Parameters
    ----------
    fp_data = fingerpint in dataframe form

    """
    res_fp, interaction_index, residues = prepare_tabledata(fp_df)

    html_legend = "<h3>Interactions in pre-defined binding site residues</h3><table><tr>"
    for key in interaction_colours:
        html_legend = (
            html_legend
            + '<td style="color:#fff; font-weight:bold; background-color:'
            + interaction_colours[key]
            + ';text-align:center">'
            + key
            + "</td>"
        )
    html_legend = html_legend + "</tr></table>"

    html_str = """<html>
    <style>
    .ttooltip {
    position: relative;
    display: inline-block;
    }

    .ttooltip .ttooltiptext {
    visibility: hidden;
    width: 120px;
    background-color: #000000;
    color: #fff;
    text-align: center;
    padding: 5px 0;
    border-radius: 6px;

    position: absolute;
    z-index: 1;
    top: -5px;
    left: 105%;

    opacity: 0.3;
    transition: opacity 0.3s;
    }

    .ttooltip .ttooltiptext::after {
    content: "";
    position: absolute;
    top: 100%;
    left: 50%;
    margin-left: -5px;
    border-width: 5px;
    border-style: solid;
    border-color:
    }


    .ttooltip:hover .ttooltiptext {
    visibility: visible;
    opacity: 0.7;
    }
    table{
        max-width: 600px;
    }
    </style>
    """

    html_str = html_str + '<table style="border:1px solid;border-color:#9698ed"><tr>'
    for res in residues:
        html_str = (
            html_str
            + '<th style="color:#fff;border:1px solid;border-color:#9698ed;background-color:#9698ed; text-align:center">'
            + str(res)
            + "</th>"
        )

    html_str = html_str + "</tr><tr>"

    fp_index = 0
    for res in res_fp:
        html_str = html_str + '<td style="border:1px solid;border-color:#9698ed;"><table><tr>'
        for i in res:
            if i == 0:
                html_str = html_str + '<td style="background-color:#ffffff"></td>'
            else:
                interaction_colour, interaction_type = cell_colour(fp_index, interaction_index)
                html_str = (
                    html_str
                    + '<td style="color:#fff; font-weight:bold; background-color:'
                    + interaction_colour
                    + '"><div class="ttooltip">'
                    + str(i)
                    + '<span class="ttooltiptext">'
                    + interaction_type
                    + "</span></div></td>"
                )
            fp_index += 1
        html_str = html_str + "</tr></table></td>"

    html_str = html_str + "</tr></table></html>"
    display(HTML(html_legend))
    display(HTML(html_str))

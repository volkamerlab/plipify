from plip_fingerprints import divide_list, read_residues
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import seaborn as sns
import pandas as pd
from IPython.core.display import display, HTML
import ipywidgets as widgets
from ipywidgets import HBox

interaction_colours = {
    "hydrophobic": "#4d4dff",
    "hbond": "#ff3300",
    "waterbridge": "#00cc7a",
    "saltbridge": "#cc66ff",
    "pistacking": "#ff9933",
    "pication": "#33ccff",
    "halogen": "#ff6699",
    "metal": "#ace600",
}


def prepare_plotdata(residue_file, interaction_types, fingerprint):
    """
    Prepare fingerprint data for plotting.
    """
    residues = read_residues(residue_file)
    res_fp = list(divide_list(fingerprint, len(interaction_types)))
    fp_df = pd.DataFrame(res_fp, columns=interaction_types, index=residues)
    fp_df = fp_df.T
    df_plot = fp_df.loc[
        :, (fp_df != 0).any(axis=0)
    ]  # change to eliminate redundant transpose
    plotdata = df_plot.T
    return plotdata


def fingerprint_barplot(residue_file, interaction_types, fingerprint):
    """
    Visualize fingerprint as barplot either based on count or fp_type
    """
    plotdata = prepare_plotdata(residue_file, interaction_types, fingerprint)

    # Plotly Bar Plot
    resis = map(str, list(plotdata.index))
    resis = ["res " + s for s in resis]
    fig = go.Figure(
        data=[
            go.Bar(name="hydrophobic", x=resis, y=list(plotdata["hydrophobic"])),
            go.Bar(name="hbond", x=resis, y=list(plotdata["hbond"])),
            go.Bar(name="waterbridge", x=resis, y=list(plotdata["waterbridge"])),
            go.Bar(name="saltbridge", x=resis, y=list(plotdata["saltbridge"])),
            go.Bar(name="pistacking", x=resis, y=list(plotdata["pistacking"])),
            go.Bar(name="pication", x=resis, y=list(plotdata["pication"])),
            go.Bar(name="halogen", x=resis, y=list(plotdata["halogen"])),
            go.Bar(name="metal", x=resis, y=list(plotdata["metal"])),
        ]
    )
    # Change the bar mode
    fig.update_layout(
        barmode="stack",
        title_text="Residue Interactions",
        yaxis_title="Interactions",
        xaxis_title="Residues",
    )
    fig.show()


def fingerprint_heatmap(residue_file, interaction_types, fingerprint):
    """
    Visualize fingerprint as heatmap either based on count or fp_type
    """
    plotdata = prepare_plotdata(residue_file, interaction_types, fingerprint)
    fig, ax = plt.subplots(figsize=(10, 7))  # plot size
    sns.heatmap(plotdata.T, annot=True, cmap="YlGnBu", ax=ax)
    plt.xlabel("Residues")
    plt.ylabel("Interaction Types")


def fingerprint_table(residue_file, interaction_types, count_fp, freq_fp, fp_type):
    """
    Create table visualization for fingerprint selection.
    """
    if fp_type == "frequency":
        fp_table_html(residue_file, interaction_types, freq_fp)
    else:
        fp_table_html(residue_file, interaction_types, count_fp)


def prepare_tabledata(residue_file, interaction_types, fingerprint):
    """
    Create interaction index dictionary for fingerprint table.
    """
    residues = read_residues(residue_file)
    res_fp = list(divide_list(fingerprint, len(interaction_types)))
    fp_id = range(0, len(residues) * len(interaction_types))  # all fp indices
    interaction_list = interaction_types * len(residues)
    interaction_index = dict(zip(fp_id, interaction_list))
    return res_fp, interaction_index


def cell_colour(fp_index, interaction_index):
    """
    Extract cell colour for specific residue interaction.
    """
    interaction_type = interaction_index[fp_index]
    interaction_colour = interaction_colours[interaction_type]
    return interaction_colour, interaction_type


def fp_table_html(residue_file, interaction_types, fingerprint):
    """
    Create HTML and CSS layout for the fingerprint table.
    """
    res_fp, interaction_index = prepare_tabledata(
        residue_file, interaction_types, fingerprint
    )
    html_legend = (
        "<h3>Interactions in pre-defined binding site residues</h3><table><tr>"
    )
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

    html_str = """<html><style>
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
    </style>
    """
    res_fp = list(
        divide_list(fingerprint, len(interaction_types))
    )  # divide fp for each redisue
    residues = read_residues(residue_file)
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
        html_str = (
            html_str + '<td style="border:1px solid;border-color:#9698ed;"><table><tr>'
        )
        for i in res:
            if i == 0:
                html_str = html_str + '<td style="background-color:#ffffff"></td>'
            else:
                interaction_colour, interaction_type = cell_colour(
                    fp_index, interaction_index
                )
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

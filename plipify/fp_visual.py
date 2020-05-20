from plip_fingerprints import divide_list, read_residues
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import seaborn as sns
import pandas as pd
from IPython.core.display import display, HTML
import ipywidgets as widgets
from ipywidgets import HBox


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

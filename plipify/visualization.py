"""
visualization.py

The visualization module for the calculated fingerprints.
This module takes the fingerprint in dataframe, processes the data
further and creates mutliple different visualizations:
- Interactive stacked bar chart
- Heatmap
- Colour coded fingerprint table
"""

import matplotlib.pyplot as plt
import plotly.graph_objects as go
import seaborn as sns


INTERACTION_PALETTE = {
    "hbond-don": "#22bbff",
    "hbond-acc": "#33ccff",
    "hydrophobic": "#ff6699",
    "waterbridge": "#4d4dff",
    "saltbridge": "#ff3300",
    "pistacking": "#00cc7a",
    "pication": "#cc66ff",
    "halogen": "#ace600",
    "metal": "#ff9933",
}


def fingerprint_barplot(fingerprint_df):
    """
    Visualize fingerprint as barplot either based on count or fp_type.

    Parameters
    ----------
    fp_data = fingerpint in dataframe form

    """
    fig = go.Figure(
        data=[
            go.Bar(
                name=interaction, x=list(fingerprint_df.index), y=list(fingerprint_df[interaction])
            )
            for interaction in sorted(fingerprint_df.columns, reverse=True)
        ],
    )
    # Change the bar mode
    fig.update_layout(
        barmode="stack",
        title_text="Residue Interactions",
        yaxis_title="Interactions",
        xaxis_title="Residues",
    )
    return fig


def fingerprint_heatmap(fingerprint_df, cmap="YlGnBu"):
    """
    Visualize fingerprint as heatmap either based on count or fp_type.

    Parameters
    ----------
    fingerprint_df = fingerpint in dataframe form

    """
    fig, ax = plt.subplots(figsize=(10, 7))  # plot size
    sns.heatmap(fingerprint_df.T, annot=True, cmap=cmap, ax=ax)
    ax.set_xlabel("Residues")
    ax.set_ylabel("Interaction Types")
    return fig


def _prepare_tabledata(fingerprint_df):
    """
    Create interaction index dictionary for fingerprint table.
    The keys of this dictionary are the interaction types and the values
    for each interaction type are a list of the positions (indices) of
    the fingerprint array that represent this interaction type.

    Parameters
    ----------
    fingerprint_df = fingerpint in dataframe form

    Returns
    -------
    fingerprint : list of int
    interaction_index : dict[int, str]
        Maps position in a fingerprint to the interaction type
    residues : list of str
        List of residue string representations
    """
    residues = list(fingerprint_df.index)
    interaction_types = list(fingerprint_df.columns)
    fingerprint = fingerprint_df.values.tolist()
    fp_id = range(len(residues) * len(interaction_types))  # all fp indices
    interaction_list = interaction_types * len(residues)
    interaction_index = dict(zip(fp_id, interaction_list))
    return fingerprint, interaction_index, residues


_TABLE_CSS = """
    table.plipify-legend {
        text-align: center;
        color: #fff;
    }

    table.plipify-legend td {
        padding: 3px;
    }
    table.plipify-interactions {
        border: 1px solid #9698ed;
        text-align: center;
        color: #fff;
        font-weight: bold;
    }

    table.plipify-interactions th {
        background-color: #9698ed;
        padding: 3px;
    }

    .plipify-ttooltip {
        position: relative;
        display: inline-block;
    }

    .plipify-ttooltip .plipify-ttooltiptext {
        visibility: hidden;
        width: 120px;
        background-color: #000000;
        color: #fff;
        text-align: center;
        padding: 5px 0;
        border-radius: 5px;

        position: absolute;
        z-index: 1;
        top: -5px;
        left: 105%;

        opacity: 0.3;
        transition: opacity 0.3s;
    }

    .plipify-ttooltip .plipify-ttooltiptext::after {
        content: "";
        position: absolute;
        top: 100%;
        left: 50%;
        margin-left: -5px;
        border: 5px solid;
    }


    .plipify-ttooltip:hover .plipify-ttooltiptext {
        visibility: visible;
        opacity: 0.7;
    }

    """


def fingerprint_table(fingerprint_df, as_widget=True):
    """
    Create HTML and CSS table layout for the calculated fingerprint.

    Parameters
    ----------
    fp_data = fingerpint in dataframe form
    as_widget = whether to build a IPyWidget object or just return the HTML string
    """
    fingerprint, interaction_index, residues = _prepare_tabledata(fingerprint_df)

    html = f"""
    <style>
    {_TABLE_CSS}
    </style>
    <h3>Interactions in pre-defined binding site residues</h3>
    <table class="plipify-legend">
        <tr>
    """
    for key in INTERACTION_PALETTE:
        html += f'<td style="background-color:{INTERACTION_PALETTE[key]}">{key}</td>'

    html += """
        </tr>
    </table>

    <table class="plipify-interactions">
        <tr>
    """

    for residue in residues:
        html += f"<th>{residue}</th>"
    html += "</tr><tr>"

    fp_index = 0
    for residue_fp in fingerprint:
        html += """<td>
            <table>
                <tr>
        """
        for bit in residue_fp:
            if bit == 0:  # no interactions reported at this position
                html += "<td></td>"
            else:
                interaction_type = interaction_index[fp_index]
                interaction_colour = INTERACTION_PALETTE[interaction_type]
                html += f"""
                    <td style="background-color: {interaction_colour}">
                        <div class="plipify-ttooltip">{bit}
                            <span class="plipify-ttooltiptext">{interaction_type}</span>
                        </div>
                    </td>
                    """
            fp_index += 1
        html += """
                </tr>
            </table>
        </td>
        """

    html += """
        </tr>
    </table>
    """

    if as_widget:
        from ipywidgets import HTML

        return HTML(html)
    return html


def fingerprint_nglview(fingerprint_df, protein_path):
    """
    Depict interaction hotspots in NGLView

    Parameters
    ----------
    fingerprint_df: pandas.Dataframe
    protein_path: str or pathlib.Path

    Returns
    -------
    view : nglview.NGLWidget
    """
    import nglview as nv

    view = nv.NGLWidget(height="600px")
    structure = nv.adaptor.FileStructure(str(protein_path))
    protein_component = view.add_component(structure)

    selection, selection_ons = [], []
    tooltips = {}
    for residue, row in fingerprint_df.iterrows():
        resid = residue[4:].split(".")[0]
        values = sorted(zip(fingerprint_df.columns, row), key=lambda kv: kv[1], reverse=True)
        tooltips[resid] = ", ".join([f"{val}x{col.title()}" for (col, val) in values if val])

        # Display interacting residues
        selection.append(f"({resid} and not _H)")
        selection_ons.append(f"({resid} and ((_O) or (_N) or (_S)))")

    protein_component.add_ball_and_stick(
        sele=" or ".join(selection), colorScheme="chainindex", aspectRatio=1.5
    )
    protein_component.add_ball_and_stick(
        sele=" or ".join(selection_ons), colorScheme="element", aspectRatio=1.5
    )

    # Patch tooltips
    # Adapted from https://github.com/umesh-timalsina/mbuild/blob/master/mbuild/utils/jsutils.py

    tooltip_js = (
        f"""
        const extraTooltips = {tooltips};
        """
        # We split the fstring in two so the curly braces below do not need to escaped
        """
        this.stage.mouseControls.add('hoverPick', (stage, pickingProxy) => {
            let tooltip = this.stage.tooltip;
            if(pickingProxy && pickingProxy.atom && !pickingProxy.bond){
                let atom = pickingProxy.atom;
                tooltip.innerText = "ATOM: " + atom.qualifiedName();
                if(atom.resno in extraTooltips) {
                    tooltip.innerText += ". Residue interactions: " + extraTooltips[atom.resno];
                }
            }
        });
        """
    )

    infotext_js = (
        f"""
        const extraTooltips = {tooltips};
        """
        """
        this.stage.signals.clicked.removeAll();
        this.stage.signals.clicked.add((pickingProxy) => {
                if(pickingProxy){
                   let pickingText = null;
                   this.model.set('picked', {});
                   this.touch();
                   let currentPick = {};
                   if(pickingProxy.atom){
                        currentPick.atom1 = pickingProxy.atom.toObject();
                        currentPick.atom1.name = pickingProxy.atom.qualifiedName();
                        pickingText = "Atom: " + currentPick.atom1.name;
                        if(currentPick.atom1.resno in extraTooltips) {
                            pickingText += ". Residue interactions: " + extraTooltips[currentPick.atom1.resno];
                        }
                   }

                   if(pickingProxy.instance){
                        currentPick.instance = pickingProxy.instance;
                   }
                   var nComponents = this.stage.compList.length;
                   for(let i = 0; i < nComponents; i++){
                        let comp = this.stage.compList[i];
                        if(comp.uuid == pickingProxy.component.uuid){
                            currentPick.component = i;
                        }
                   }
                   this.model.set('picked', currentPick);
                   this.touch();
                   this.$pickingInfo.text(pickingText);
                }
        });
        """
    )
    view._js(tooltip_js)
    view._js(infotext_js)

    return view
"""
visualization.py

The visualization module for the calculated fingerprints.
This module takes the fingerprint in dataframe, processes the data
further and creates mutliple different visualizations:
- Interactive stacked bar chart
- Heatmap
- Colour coded fingerprint table
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import MDAnalysis as mda
import plotly.graph_objects as go
import pymol
import seaborn as sns
from pymol import cmd, util

from plipify.core import Structure

INTERACTION_PALETTE = {
    "hbond-don": "#0173b2",
    "hbond-acc": "#de8f05",
    "hbond-don-bb": "#0173b2",
    "hbond-don-sc": "#08C1C9",
    "hbond-acc-bb": "#de8f05",
    "hbond-acc-sc": "#F57B07",
    "waterbridge": "#d55e00",
    "saltbridge": "#cc78bc",
    "pistacking": "#ca9161",
    "pication": "#fbafe4",
    "halogen": "#949494",
    "metal": "#ece133",
    "hydrophobic": "#029e73",
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
                name=interaction,
                x=list(fingerprint_df.index),
                y=list(fingerprint_df[interaction]),
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
    fig.update_xaxes(type="category")

    return fig


def fingerprint_heatmap(fingerprint_df, cmap="YlGnBu"):
    """
    Visualize fingerprint as heatmap either based on count or fp_type.

    Parameters
    ----------
    fingerprint_df = fingerpint in dataframe form

    """
    fig, ax = plt.subplots(figsize=(10, 7))  # plot size
    sns.heatmap(fingerprint_df, annot=True, cmap=cmap, ax=ax, fmt="g")
    ax.set_xlabel("Interaction Types")
    ax.set_ylabel("Residues")
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

    ## Reorder fingerprint_df columns to match the INTERACTION PALETTE
    selected_interaction_palette = {int_type: values for int_type, values in INTERACTION_PALETTE.items() if
                                    int_type in list(fingerprint_df.columns)}
    fingerprint_df = fingerprint_df.reindex(columns=selected_interaction_palette.keys())

    ## so that the resulting table will be in the same order as the INTERACTION PALETTE
    interaction_types = list(fingerprint_df.columns)
    fingerprint = fingerprint_df.values.tolist()
    fp_id = range(len(residues) * len(interaction_types))  # all fp indices
    interaction_list = interaction_types * len(residues)
    interaction_index = dict(zip(fp_id, interaction_list))
    return fingerprint_df, fingerprint, interaction_index, residues


_TABLE_CSS = """

    table.plipify-legend {
        text-align: center;
        color: #fff;
    }

    table.plipify-legend td {
        padding: 3px;
        border: 3px solid white;
    }

    table.plipify-interactions {
        border: 1px solid black;
        text-align: center;
        color: white;
        font-weight: bold;
    }

    table.plipify-interactions th {
        text-align: center;
        background-color: #005473;
        padding: 3px;
        border-right: 1px solid black;
        border-left: 1px solid black;
    }

    table.plipify-interactions tr {
        text-align: center;
        background-color: white;
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


def fingerprint_table(fingerprint_df, as_widget=True, structure=None):
    """
    Create HTML and CSS table layout for the calculated fingerprint.

    Parameters
    ----------
    fingerprint_df : pd.DataFrame
        Fingerpint in dataframe form
    as_widget: bool
        Whether to build a IPyWidget object or just return the HTML string
    structure : plipify.core.Structure
        If supplied, will use this to generate residue labels in the form A123 etc

    """
    fingerprint_df, fingerprint, interaction_index, residues = _prepare_tabledata(fingerprint_df)

    # if a structure is supplied, get one-letter residues e.g. H163
    ol_residues = None

    if structure:
        try:
            fingerprint_df.attrs["residues"] = {
                r.seq_index: f"{structure.residues[i].one_letter_code}{r.seq_index}"
                for i, r in enumerate(structure.residues)
            }

            ol_residues = fingerprint_df.attrs["residues"]

        except AttributeError:
            print("Supplied structure argument is not a plipify.core.Structure object!")



    html = f"""
    <style>
    {_TABLE_CSS}
    </style>
    <h3>Interactions in pre-defined binding site residues</h3>
    <table class="plipify-legend">
        <tr>
    """
    for key in INTERACTION_PALETTE:
        if key in fingerprint_df.columns.to_list():
            html += f'<td style="background-color:{INTERACTION_PALETTE[key]}">{key}</td>'

    html += """
        </tr>
    </table>

    <table class="plipify-interactions">
        <tr>
    """

    for residue in residues:
        if ol_residues and residue in ol_residues:
            residue = ol_residues[residue]

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
                fp_index += 1
                continue
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


def fingerprint_nglview(fingerprint_df, structure, fp_index_to_residue_id=None):
    """
    Depict interaction hotspots in NGLView

    Parameters
    ----------
    fingerprint_df: pandas.Dataframe
    structure: plipify.core.Structure
    fp_key_to_residue_id : dict
        Maps fingerprint index to the actual residue ids in `structure`.
        Same as expected by InteractionFingerprint.calculate_fingerprint
        `residue_indices`.

    Returns
    -------
    view : nglview.NGLWidget
    """
    view = structure.view(solvent_selection_query="NOT all")

    selection, selection_ons = [], []
    tooltips = {}
    for resid, row in fingerprint_df.iterrows():
        values = sorted(zip(fingerprint_df.columns, row), key=lambda kv: kv[1], reverse=True)
        if fp_index_to_residue_id is not None:
            residue = structure.get_residue_by(**fp_index_to_residue_id[resid])
            if residue:
                resid = residue.seq_index
            else:
                print(
                    resid,
                    "->",
                    fp_index_to_residue_id[resid],
                    "not found in this structure",
                )
                continue

        tooltips[resid] = ", ".join([f"{val}x{col.title()}" for (col, val) in values if val])

        # Display interacting residues
        selection.append(f"({resid} and not _H)")
        selection_ons.append(f"({resid} and ((_O) or (_N) or (_S)))")

    view.add_ball_and_stick(sele=" or ".join(selection), colorScheme="chainindex", aspectRatio=1.5)
    view.add_ball_and_stick(
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


def nglview_color_side_chains_by_frequency(fp_focused, selected_structure_pdb, ligand="LIG"):
    """
    Adapted from this [gist](https://gist.github.com/dominiquesydow/1aa19cc5583a31bfbd58bc8a066a7705)
    by Dominique Sydow.
    """
    import nglview as nv
    from matplotlib import cm, colors

    # add sum column
    fp_color = fp_focused.copy()
    fp_color["sum"] = fp_color.sum(axis=1)
    # prepare color maps
    max_n_interactions = fp_color["sum"].max()
    max_n_interactions
    cmap = cm.get_cmap("Reds", max_n_interactions)
    # get residue-color map
    # NGLview needs color as hex code
    fp_color["color"] = fp_color["sum"].apply(lambda x: colors.rgb2hex(cmap(x)))
    # NGLview needs residue IDs as str
    fp_color["residue_id_str"] = fp_color.index.astype(str)
    # convert to list as needed by nglview colorscheme
    color_residue_map = fp_color[["color", "residue_id_str"]].to_numpy().tolist()
    color_scheme = nv.color._ColorScheme(color_residue_map, label="#interactions")

    # nglview
    view = nv.NGLWidget()
    structure = nv.adaptor.FileStructure(str(selected_structure_pdb))
    view.add_component(structure, default_representation=False)

    view.add_representation(
        "cartoon",
        selection="protein",
        color=color_scheme,
    )

    # add selected side chains
    residue_list = fp_color["residue_id_str"].tolist()
    view.add_representation(
        "licorice",
        selection=f"protein and ({' '.join([str(x) for x in residue_list])})",
        color=color_scheme,
    )

    view.add_ball_and_stick(ligand)
    view.center(ligand)
    return view


def fingerprint_writepdb(
    fingerprint_df,
    structure,
    output_path,
    ligand=False,
    ligand_name="LIG",
    summed=False,
    verbose=True,
) -> dict:
    """
    Write interaction hotspots to a PDB file

    Parameters
    ----------
    fingerprint_df: pandas.Dataframe
    structure: plipify.core.Structure
    output_path: pathlib.PosixPath
        A path to store the output files
    ligand: bool
        Set to True to write out the bound ligand structure
    ligand_name: str
        The name of the bound ligand
    summed : bool
        Specify whether to write out summed interaction PDB, default = False
    verbose : bool
        Print some information as we process files or not, default = True
    Returns
    -------
    systems: dict[str: plipify.core.Structure]
    """

    def _load_universe(input_structure):
        u = mda.Universe(input_structure)

        return u

    # create output directory to store generated PDBs
    OUTDIR = Path(output_path) / "interaction_pdbs"
    OUTDIR.mkdir(exist_ok=True, parents=True)

    # set selection strings
    if ligand:
        sel_string = f"protein or resname {ligand_name}"
    else:
        sel_string = "protein"

    systems = {interaction: None for interaction in list(fingerprint_df)}

    for interaction_col in fingerprint_df:  # loop over interaction types

        if verbose:
            print(f"Analysing {interaction_col} interactions...")

        # create MDAnalysis.Universe
        u_int = _load_universe(structure._path)

        # set temperature factors, default = 0 for all residues
        u_int.add_TopologyAttr("tempfactors")

        sys_int = u_int.select_atoms(sel_string)

        for resid, value in fingerprint_df[interaction_col].iteritems():  # loop over residues

            # assign temperature facture value based on interaction value
            sel = sys_int.select_atoms(f"resid {str(resid)}")
            sel.atoms.tempfactors = value

        # write out new pdb for interaction type
        try:
            with mda.Writer(OUTDIR / f"sys_int_{str(interaction_col)}.pdb", sys_int.n_atoms) as W:
                W.write(sys_int)
        except TypeError:
            print(f"Warning! Couldn't write sys_int_{str(interaction_col)}.pdb")
            continue
        else:
            if verbose:
                print(f"Written sys_int_{str(interaction_col)}.pdb to {OUTDIR}")

        systems[interaction_col] = Structure.from_pdbfile(
            str(OUTDIR / f"sys_int_{str(interaction_col)}.pdb")
        )

    if summed:

        # sum interactions to get total number per residue
        fingerprint_df_sum = fingerprint_df.T.sum()

        u_summed = _load_universe(structure._path)
        u_summed.add_TopologyAttr("tempfactors")
        sys_summed = u_summed.select_atoms(sel_string)

        for resid, summed_ints in fingerprint_df_sum.iteritems():

            sel_summed = sys_summed.select_atoms(f"resid {str(resid)}")
            sel_summed.atoms.tempfactors = summed_ints

        # write out new pdb for total interactions
        try:
            with mda.Writer(OUTDIR / "sys_summed_interactions.pdb", sys_summed.n_atoms) as W:
                W.write(sys_summed)
        except TypeError:
            print("Warning! Couldn't write sys_summed_interactions.pdb")
        else:
            if verbose:
                print(f"Written sys_summed_interactions.pdb to {OUTDIR}")

        systems["summed_interactions"] = Structure.from_pdbfile(
            str(OUTDIR / "sys_summed_interactions.pdb")
        )

    return systems


class PymolVisualizer(object):
    """
    A series of tools to create publication ready images using PyMol

    Example usage:

    v = PymolVisualizer(pdb=pdb_file.pdb)
    v.set_style()
    v.create_image(surface=True, ligand_col="green", spectrum_col="white_pink")
    v.render(name="example")

    Parameters
    ----------
    pdb : str
       Path to PDB file to load e.g. plipify.core.Structure._path
    ligand_interactions: bool
        Specify whether to highlight residues in contact with the ligand, default = True.
        This requires a pdb file with bfactor values > 0
        e.g. from plipify.visualization.fingerprint_writepdb
    """

    def __init__(self, pdb: str, ligand_interactions: bool = True, verbose: bool = True):

        self._pdb = pdb
        self._ligand_interactions = ligand_interactions
        self._verbose = verbose

        # Launch pymol session
        if self._verbose:
            print("Launching PyMol session...")
        pymol.pymol_argv = ["pymol", "-qc"]
        pymol.finish_launching()

        # Load file
        if self._verbose:
            print("Loading supplied PDB file...")
        if os.path.exists(self._pdb):
            cmd.load(self._pdb)
        else:
            raise FileNotFoundError(f"{self._pdb} not found")

        # Rename objects
        # NOTE not sure this step is needed
        name = os.path.splitext(os.path.basename(self._pdb))[0]
        cmd.set_name(name, "system")

        # Remove solvent
        # TODO check for more names
        cmd.remove("resn HOH")
        cmd.remove("resn SOL")
        cmd.remove("resn NA")
        cmd.remove("resn CL")

    def set_style(
        self,
        ambient: float = 0.4,
        ambient_occlusion_mode: int = 1,
        ambient_occlusion_scale: float = 15,
        background_col: str = "white",
        antialias: int = 1,
        ortho: int = 1,
        ray_trace_mode: int = 0,
        flat_sheets: int = 1,
    ):

        """
        Set style using PyMol parameters

        Parameters
        ----------
        ambient : float
            Set the amount of ambient light (0-1), default = 0.4
        ambient_occlusion_mode : int
            Darken areas of the protein that are occluded (0 = off, 1 = on), default = 1
        ambient_occlusion_scale : float
            The scale by which ambient occlusion values are modified, default = 15
        background_col : str
            The background colour, default = "white"
        antialias : int
            Control the amount of antialiasing PyMol does using ray (0-4), default = 1
        ortho : int
            Specify whether to use an orthographic view (0 = off, 1 = on), default = 1
        ray_trace_mode : int
            Control in-built PyMol ray trace mode (0-3), default = 0
        flat_sheets: int
            Specify whether beta sheets should be artistically flattened (set to 0 to follow backbone), default = 1
        """

        # Add filter (ambient) and other misc settings
        # TODO add other options

        cmd.set("ambient", ambient)
        cmd.set("ambient_occlusion_mode", ambient_occlusion_mode)
        cmd.set("ambient_occlusion_scale", ambient_occlusion_scale)
        cmd.bg_colour(background_col)
        cmd.set("antialias", antialias)
        cmd.set("ortho", ortho)
        cmd.set("ray_trace_mode", ray_trace_mode)

        # Set the beta sheet style
        cmd.set("cartoon_flat_sheets", flat_sheets)

    def create_image(
        self,
        highlight_style: str = "sticks",
        show_backbone: bool = False,
        spectrum_col: str = "white_green",
        surface: bool = False,
        surface_mode: int = 3,
        surface_col: str = "white",
        transparency: float = 0.75,
        protein_style: str = "cartoon",
        protein_col: str = "white",
        show_ligand: bool = True,
        ligand_resname: str = "LIG",
        ligand_style: str = "sticks",
        ligand_col: str = "cyan",
        cnc_protein: bool = False,
        cnc_ligand: bool = True,
        view: str = "ligand",
        viewport_x: int = 720,
        viewport_y: int = 720,
    ):

        """
        Create a PyMol image

        Parameters
        ----------
        highlight_style : str
            Set the style for residues to be highlighted as, default = "sticks"
        show_backbone : bool
            Specify whether to show the backbone of highlighted residues, default = False
        spectrum_col : str
            Set the spectrum colour to use for highlighting interaction sites, default = "white_green"
        surface : bool
            Specify whether to show the protein surface, default = False
        surface_mode : int
            Set what atoms are used to specify the surface (0-2), default = 3
        surface_col : str
            Set the colour of the protein surface, default = "white"
        transparency : float
            Set the transparency of the surface if specified, default = 0.75
        protein_style : str
            Set the protein style, default = "cartoon"
        protein_col : str
            Set the protein colour, default = "white"
        show_ligand : bool
            Specify whether to show a bound ligand based on ligand_resname, default = True
        ligand_resname : str
            The name of the bound ligand, default = "LIG"
        ligand_style : str
            Set the ligand style, default = "sticks"
        ligand_col : str
            Set the ligand colour, default = "cyan"
        cnc_protein : bool
            Specify whether to only colour Carbon atoms in the protein, default = False
        cnc_ligand : bool
            Specify whether to only colour Carbon atoms in the ligand, default = True
        view: str
            Specify where to centre the camera on, can take "ligand" to centre on supplied ligand or the output from PyMol "get_view" command, default = "ligand"
        viewport_x : int
            Set the viewport dimensions in x, default = 720
        viewport_y : int
            Set the viewport dimensions in y, default = 720
        """

        self._viewport_x = viewport_x
        self._viewport_y = viewport_y

        # Set colours
        cmd.set_color("cborange", [0.8706, 0.5608, 0.0196])
        cmd.set_color("cbpurple", [0.66, 0.35, 0.63])
        cmd.set_color("cblight_blue", [0.52, 0.75, 0.98])
        cmd.set_color("cbgreen", [0.0078, 0.6196, 0.4510])
        cmd.set_color("cbblue", [0.0039, 0.4510, 0.6980])
        cmd.set_color("cbred", [0.8353, 0.3686, 0])
        cmd.set_color("cbbrown", [0.7921, 0.5686, 0.3804])
        cmd.set_color("cbgrey", [0.5804, 0.5804, 0.5804])
        cmd.set_color("cbgray", [0.5804, 0.5804, 0.5804])

        # Set colour dictionary
        c_dict = {
            "yellow": cmd.util.cbay,
            "green": cmd.util.cbag,
            "cyan": cmd.util.cbac,
            "light magenta": cmd.util.cbam,
            "salmon": cmd.util.cbam,
            "white": cmd.util.cbaw,
            "slate": cmd.util.cbab,
            "orange": cmd.util.cbao,
            "purple": cmd.util.cbap,
            "pink": cmd.util.cbak,
        }

        # Select protein residues (backbone and sidechain)
        cmd.select("prot", "bb. or sc.")

        if show_ligand:
            cmd.select("ligand", f"resn {ligand_resname}")

        # Sort out protein and ligand colours
        if show_ligand:
            c_dict[ligand_col]("ligand")
        c_dict[protein_col]("prot")

        # Hide everything
        cmd.hide("all")
        # Ensure secondary structure is show correctly
        cmd.dss("prot")

        # Colour protein residues based on bfactor value, if specified
        if self._ligand_interactions:
            if self._verbose:
                print("Highlighting ligand interaction sites...")
            cmd.spectrum("b", spectrum_col, "prot")
            if show_backbone:
                cmd.show(highlight_style, "prot and not hydrogen and b > 0")
            else:
                cmd.select(
                    "hotspots",
                    "not resn PRO and prot and not name N and not name C and not name O and not hydrogen and b > 0",
                )
                cmd.select(
                    "hotspots_pro",
                    "resn PRO and not name C and not name O and not hydrogen and b > 0",
                )
                cmd.show(
                    highlight_style,
                    "hotspots",
                )
                # Handle Proline
                cmd.show(
                    highlight_style,
                    "hotspots_pro",
                )

        # Show the whole protein
        cmd.show(protein_style, "prot and not hydrogen")
        if surface:
            cmd.show("surface", "prot and not hydrogen")
        cmd.enable("prot")

        # Show the ligand, if specified
        if show_ligand:
            cmd.enable("ligand")
            cmd.show(ligand_style, "ligand and not hydrogen")
            if cnc_ligand:
                util.cnc("ligand")

        # Colour only Carbon atoms, if specified
        if cnc_protein:
            util.cnc("prot")

        # Add a protein surface, if specified
        if surface:
            cmd.set("surface_mode", surface_mode)
            cmd.set("transparency", transparency)
            cmd.set("surface_color", surface_col)

            cmd.set("transparency", 0.75, "b>0")

        # Set the viewport and view
        if self._verbose:
            print("Setting PyMol view...")
        cmd.viewport(self._viewport_x, self._viewport_y)

        if view == "ligand":
            if self._verbose:
                print("Focussing view on ligand and binding site")
            cmd.center("ligand or hotspots or hotspots_pro")
            cmd.zoom("ligand or hotspots or hotspots_pro", 5)

        else:
            if len(view.split()) == 18:
                cmd.set_view(f"{view}")
                if self._verbose:
                    print("Focussing view on user-supplied coordinates")
            else:
                if self._verbose:
                    print("Couldn't parse supplied view!")

    def render(self, name, save_path="./", dpi=300):

        """
        Render PyMol image

        Parameters
        ----------

        name : str
            The name of the image
        save_path : str
            The path of where to save the PyMol image, default = "./"
        dpi : float
            The dots per inch of the image, default = 300
        """

        d = Path(save_path)
        d.mkdir(exist_ok=True, parents=True)

        # Create image
        if self._verbose:
            print("Rendering PyMol image...")
        cmd.set("ray_transparency_contrast", 3.0)
        cmd.ray(self._viewport_x, self._viewport_y)
        filename = str(d / f"{name}.png")
        cmd.png(filename, dpi=dpi)
        if self._verbose:
            print(f"Image created! Saving to {filename}")

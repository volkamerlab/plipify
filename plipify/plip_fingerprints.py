import csv
import os
import pandas as pd
from tqdm.notebook import tqdm
from plip.modules.preparation import PDBComplex
from plip.modules.report import BindingSiteReport
from IPython.display import display, Markdown

##### PREPARATION ######
# fixed for now
interaction_types = [
    "hydrophobic",
    "hbond",
    "waterbridge",
    "saltbridge",
    "pistacking",
    "pication",
    "halogen",
    "metal",
]


def read_residues(path):
    """
    Reads residue file.
    """
    with open(path) as file:
        residue_reader = csv.reader(file, delimiter=",")
        residues = residue_reader.next()
    residues = map(int, residues)
    return residues


def divide_list(list_name, n):
    """
    Split list into lists of given size.
    """
    for i in range(0, len(list_name), n):
        yield list_name[i : i + n]


def residue_dictionary(residues):
    """
    Create index dictionary for residues.
    """
    residue_i = range(0, len(residues) * len(interaction_types))
    residue_i = list(divide_list(residue_i, len(interaction_types)))
    residue_dict = dict(zip(residues, residue_i))
    return residue_dict


def interaction_dictionary(interaction_types):
    """
    Create index dictionary for interactions.
    """
    interaction_i = range(0, len(interaction_types))
    interaction_dict = dict(zip(interaction_types, interaction_i))
    return interaction_dict


###### PLIP DATA ######
def analyze_interactions(pdbfile):
    protlig = PDBComplex()
    protlig.load_pdb(pdbfile)  # load the pdb
    for ligand in protlig.ligands:
        protlig.characterize_complex(ligand)  # find ligands and analyze interactions
    sites = {}
    for key, site in sorted(protlig.interaction_sets.items()):
        binding_site = BindingSiteReport(site)  # collect data about interactions
        # tuples of *_features and *_info will be converted to pandas df
        keys = (
            "hydrophobic",
            "hbond",
            "waterbridge",
            "saltbridge",
            "pistacking",
            "pication",
            "halogen",
            "metal",
        )
        interactions = {
            k: [getattr(binding_site, k + "_features")] + getattr(binding_site, k + "_info")
            for k in keys
        }
        sites[key] = interactions
    return sites


# for displaying data
def site_to_dataframes(site):
    keys = [
        "hydrophobic",
        "hbond",
        "waterbridge",
        "saltbridge",
        "pistacking",
        "pication",
        "halogen",
        "metal",
    ]
    dfs = {}
    for key in keys:
        records = site[key][1:]
        if not records:
            dfs[key] = None
        else:
            dfs[key] = pd.DataFrame.from_records(records, columns=site[key][0])
    return dfs


def get_plip_data(data, name_file):
    """
    Get data from plip and store it in a dictionary.
    """
    with open(os.path.join(data, name_file)) as f:
        non_covalent_filenames = [line.strip() for line in f if line.strip()]
    interactions = {}
    for filename in tqdm(non_covalent_filenames):
        full_filename = os.path.join(data, filename)
        if not os.path.isfile(full_filename):
            print("File", full_filename, "does not exist?")
            continue
        interactions[filename] = analyze_interactions(full_filename)
    return interactions


def show_plip_data(interactions):
    for structure, sites in interactions.items():
        display(Markdown("# Structure {}".format(structure)))
        for site_name, site_interactions in sites.items():
            if not site_name.startswith("LIG"):
                continue  # fragments are labeled as LIG; other "sites" detected by PLIP are XRC artefacts
            display(Markdown("## Site {}".format(site_name)))
            for interaction_type, dataframe in site_to_dataframes(site_interactions).items():
                if dataframe is not None:
                    display(Markdown("### {}".format(interaction_type)))
                    display(dataframe)


###### FINGERPRINTS ######
def interaction_fingerprint(residue_dictionary, interaction_dict, residue_list, interaction_type):
    """
    Create one indivual fingerprint for one structure.
    """
    fp_length = len(residue_dictionary) * len(interaction_dict)
    fp = [0] * fp_length
    for residue in residue_list:
        fp_index = residue_dictionary[residue][interaction_dict[interaction_type]]
        fp[fp_index] = fp[fp_index] + 1
    return fp


def interaction_fingerprint_list(interactions, residue_dict, interaction_dict):
    """
    Create list of fingerprints for all given structures.
    """
    fp_list = []
    for sites in interactions.items():
        for site_name, site_interactions in sites.items():
            if not site_name.startswith("LIG"):
                continue  # fragments are labeled as LIG; other "sites" detected by PLIP are XRC artefacts
            for interaction_type, dataframe in site_to_dataframes(site_interactions).items():
                if dataframe is not None:
                    residue_nos = dataframe["RESNR"].tolist()
                    fp = interaction_fingerprint(
                        residue_dict, interaction_dict, residue_nos, interaction_type
                    )
                    fp_list.append(fp)
    return fp_list


def frequency_interaction_fingerprint(fp_list):
    """
    Create frequency fingerprint from list of fingerprints.
    """
    count_fp = [sum(i) for i in zip(*fp_list)]
    frequency_fp = [float(i) / len(fp_list) for i in count_fp]
    return frequency_fp

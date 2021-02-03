"""
core.py

This module collects the base objects that will be
used across the different stages of the pipeline.

Namely:

- Interaction, which "connects" two residues
- Residue, the most basic interacting component
- Structure, a list of binding sites
- BindingSite, a list of residues, within a Structure
- Dataset, a collection of structures
"""

from collections import defaultdict, Counter
from pathlib import Path

from plip.structure.preparation import PDBComplex
from plip.exchange.report import BindingSiteReport
from Bio.Data import IUPACData

###
# Interaction Types
###


class BaseInteraction:
    """
    This the base class for all interaction types

    Parameters
    ----------
    interaction : dict
        Unspecified schema. Might be different per interaction type!
        For example, HBonds might have donor and acceptors,
        but hydrophobic are not aware of such a concept.
    """

    shorthand = ""

    def __init__(self, interaction):
        self.interaction = interaction

    def __getitem__(self, value):
        return self.interaction[value]

    def __repr__(self):
        return f"<{self.__class__.__name__} with {self.interaction}>"

    def to_dataframe(self):
        import pandas as pd

        return pd.DataFrame.from_dict(self.interaction, orient="index").T

    def _ipython_display_(self):
        if self.interaction:
            from IPython.display import display

            display(self.to_dataframe())


class HydrophobicInteraction(BaseInteraction):
    """
    HydrophobicInteraction, a subclass of BaseInteraction
    """

    shorthand = "hydrophobic"


class HbondInteraction(BaseInteraction):
    """
    HbondInteraction, a subclass of BaseInteraction
    """

    shorthand = "hbond"


class HbondDonorInteraction(HbondInteraction):
    """
    HbondInteraction, a subclass of BaseInteraction
    """

    shorthand = "hbond-don"


class HbondAcceptorInteraction(HbondInteraction):
    """
    HbondInteraction, a subclass of BaseInteraction
    """

    shorthand = "hbond-acc"


class WaterbridgeInteraction(BaseInteraction):
    """
    WaterbridgeInteraction, a subclass of BaseInteraction
    """

    shorthand = "waterbridge"


class SaltbridgeInteraction(BaseInteraction):
    """
    SaltbridgeInteraction, a subclass of BaseInteraction
    """

    shorthand = "saltbridge"


class PistackingInteraction(BaseInteraction):
    """
    PistackingInteraction, a subclass of BaseInteraction
    """

    shorthand = "pistacking"


class PicationInteraction(BaseInteraction):
    """
    PicationInteraction, a subclass of BaseInteraction
    """

    shorthand = "pication"


class HalogenInteraction(BaseInteraction):
    """
    HalogenInteraction, a subclass of BaseInteraction
    """

    shorthand = "halogen"


class MetalInteraction(BaseInteraction):
    """
    MetalInteraction, a subclass of BaseInteraction
    """

    shorthand = "metal"


###
#  Structural objects
###


class BaseResidue:
    """
    A collection of covalently bonded atoms
    """

    _ALLOWED_RESIDUE_NAMES = []

    def __init__(self, name):
        self.name = self._check_valid_name(name)

    def _check_valid_name(self, name):
        if not self._ALLOWED_RESIDUE_NAMES:  # no checks defined!
            return name
        if name in self._ALLOWED_RESIDUE_NAMES:
            return name
        raise ValueError(
            "Residue name {} is not valid! Must be one of ({})".format(
                name, ", ".join(self._ALLOWED_RESIDUE_NAMES)
            )
        )


class ProteinResidue(BaseResidue):
    """
    A residue belonging to a protein sequence
    """

    # TODO: Fill list in!
    _ALLOWED_RESIDUE_NAMES = []

    def __init__(self, name, seq_index, chain, interactions=None, structure=None):
        self.seq_index = seq_index
        self.name = name
        self.chain = chain
        self.interactions = interactions or []
        self.structure = structure

    def count_interactions(self):
        interaction_types = [interaction.shorthand for interaction in self.interactions]
        counter = Counter(interaction_types)
        return counter

    def __repr__(self):
        if self.interactions:
            return "<ProteinResidue {}, and {} interactions>".format(
                self.identifier, len(self.interactions)
            )
        return "<ProteinResidue {}>".format(self.identifier)

    @property
    def identifier(self):
        return "{}:{}.{}".format(self.name, self.seq_index, self.chain)

    def is_protein(self):
        return self.name.title() in IUPACData.protein_letters_3to1

    @property
    def one_letter_code(self):
        return IUPACData.protein_letters_3to1[self.name.title()]

    @property
    def three_letter_code(self):
        return self.name.title()


class LigandResidue(BaseResidue):
    """
    A small molecule in the vicinity of a binding site
    """

    # TODO: Fill list in!
    _ALLOWED_RESIDUE_NAMES = []


class BindingSite:
    """
    A collection of protein residues and ligands
    """

    def __init__(self, interactions, name=None):
        self.interactions = interactions
        self.name = name

    def __repr__(self):
        return "<BindingSite with name='{}' and {} interactions>".format(
            self.name,
            len([interaction for interaction in self.interactions.values() if interaction]),
        )


class Structure:
    """
    A collection of ProteinResidue objects, and potentially LigandResidue
    objects, which can define one or more BindingSite objects.

    Parameters
    ----------
    residues : list of ProteinResidue
        Order is relevant, it's assumed that the order of the residues
        dictate the sequence of the protein.
    ligands: list of PLIP Ligand objects
    binding_sites: dict of {str: }

    """

    INTERACTIONS = (
        HydrophobicInteraction,
        HbondInteraction,
        WaterbridgeInteraction,
        SaltbridgeInteraction,
        PistackingInteraction,
        PicationInteraction,
        HalogenInteraction,
        MetalInteraction,
    )

    def __init__(self, residues=None, ligands=None, binding_sites=None):
        self.residues = residues or []
        self.ligands = ligands or []
        self.binding_sites = binding_sites or []
        self._path = None

    @classmethod
    def from_pdbfile(cls, path, only_ligands=None, ligand_identifier=None):
        """
        Read pdb file and then collect and process PLIP data.

        Parameters
        ----------
        path : str
            Path to an existing PDB file
        only_ligands : list of str
            Ligand names to characterize. Other ligands will be ignored.
        ligand_identifier : str
            A String that the binding site names start with, when they are ligand binding sites.
            If filled, only binding sites with this identifier will be considered.
        """
        pdbcomplex = PDBComplex()
        pdbcomplex.load_pdb(path)

        structure = cls()
        structure._path = path

        residues = []
        for r in pdbcomplex.resis:
            residue = ProteinResidue(
                name=r.GetName(),
                seq_index=r.GetNum(),
                chain=r.GetChain(),
                structure=structure,
            )
            residues.append(residue)

        structure.residues = residues

        ligands = pdbcomplex.ligands
        for ligand in ligands:
            if only_ligands is not None and ligand.longname not in only_ligands:
                continue
            pdbcomplex.characterize_complex(ligand)

        binding_sites = []
        for key, site in sorted(pdbcomplex.interaction_sets.items()):
            report = BindingSiteReport(site)
            interactions = []
            if ligand_identifier is None or key.startswith(ligand_identifier):
                for InteractionType in cls.INTERACTIONS:
                    shorthand = InteractionType.shorthand
                    if shorthand == "hbond-acc":
                        shorthand == "hbond"
                    elif shorthand == "hbond-don":
                        continue  # skip, we processed it already
                    features = getattr(report, shorthand + "_features")
                    # list of BaseInteraction Subclasses (depending on type)
                    for interaction_data in getattr(report, shorthand + "_info"):
                        interaction_dict = dict(zip(features, interaction_data))
                        seq_index, chain = (
                            interaction_dict["RESNR"],
                            interaction_dict["RESCHAIN"],
                        )
                        residue = structure.get_residue_by(seq_index=seq_index, chain=chain)
                        if shorthand == "hbond":
                            if interaction_dict["PROTISDON"]:
                                interaction_obj = HbondDonorInteraction(
                                    interaction=interaction_dict
                                )
                            else:
                                interaction_obj = HbondAcceptorInteraction(
                                    interaction=interaction_dict
                                )
                        else:
                            interaction_obj = InteractionType(interaction=interaction_dict)
                        interactions.append(interaction_obj)
                        residue.interactions.append(interaction_obj)

                interactions_by_type = defaultdict(list)
                for i in interactions:
                    interactions_by_type[i.shorthand].append(i)

                binding_site = BindingSite(interactions_by_type, name=key)
                binding_sites.append(binding_site)

        structure.binding_sites = binding_sites
        return structure

    def get_residue_by(self, index=None, seq_index=None, chain=None):
        """
        Access residues by either position in the Python list (index)
        or their sequence identifier in the PDB records (seq_index)

        Parameters
        ----------
        index, seq_index : int
        chain : str or None
            If chain is None or "any", behaviour is not specified
            shall there exist more than one residue with that
            seq_index. If None, a warning will be emitted; this
            warning is omitted if chain == "any".
        """

        if index is None and seq_index is None:
            return None
        if index is not None and seq_index is not None:
            raise ValueError("Can only specify index OR seq_index, not both")
        if index is not None and chain not in (None, "any"):
            raise ValueError("`index` does not support `chain` spec")

        if index is not None:
            return self.residues[index]
        if seq_index is not None:
            for residue in self.residues:
                if residue.seq_index == seq_index:
                    if chain is None:  # TODO: Check there are no other residues with same index!
                        print(
                            "! Warning, you didn't select a chain. "
                            "First match will be returned but there could be more."
                        )
                        break
                    elif chain in ("any", residue.chain):
                        break
            else:  # break not reached!
                raise ValueError(
                    "No residue with such sequence index: {}, {}!".format(seq_index, chain)
                )
            return residue

    def sequence(self, with_gaps=True):
        if not with_gaps:
            return "".join([r.one_letter_code for r in self.residues])

        aabypos = {r.seq_index: r.one_letter_code for r in self.residues}
        max_residue = max(r.seq_index for r in self.residues)
        return "".join([aabypos.get(i, "-") for i in range(1, max_residue + 1)])

    @property
    def identifier(self):
        if self._path:
            return Path(self._path).stem

    @property
    def description(self):
        s = (
            f"{self.__class__.__name__} with {len(self.residues)} residues, "
            f"{len(self.ligands)} ligands "
            f"and {len(self.binding_sites)} binding sites"
        )
        if self._path:
            s += f" (loaded from file `{self._path}`)"
        return s

    def __repr__(self) -> str:
        return f"<{self.description}>"

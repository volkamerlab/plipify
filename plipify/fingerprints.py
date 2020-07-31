"""
fingerprints.py
---------------

Factories that take a Structure or multiple structures and produce
an interaction fingerprint.

"""

from collections import Counter, defaultdict
import numpy as np
import pandas as pd
from .core import Structure, ProteinResidue


class InteractionFingerprint:
    """
    This class will take a protein-ligand structure,
    analyze its interactions and report a fingerprint
    per residue.

    Parameters
    ----------
    residue_indices : list of indices
    """

    def __init__(
        self,
        structures,
        residue_indices,
        interaction_types=(
            "hydrophobic",
            "hbond",
            "waterbridge",
            "saltbridge",
            "pistacking",
            "pication",
            "halogen",
            "metal",
        ),
    ):
        self.structures = structures
        self.residue_indices = residue_indices
        self.interaction_types = interaction_types

    def calculate_fingerprint(
        self,
        labeled=True,
        cumulative=True,
        as_dataframe=False,
        remove_non_interacting_residues=False,
    ):
        """
        Calulative interaction fingerprint for one or multiple strcutures.

        Parameters
        ----------
        labeled = boolean deciding whether to make each fingerprint bit a labeled value or simple integer
        cumulative = defines if the fp is a summed up fp or multiple structures
        as_dataframe = if true return fp as data_frame, else as array
        remove_non_interacting_residues = if true, remove all fp bits that belong to residues for which there are no interactions
        """

        # TODO: Some boolean paths are not covered here! Provide errors or implement missing path.
        fingerprints = []
        for structure in self.structures:
            fingerprint = self._calculate_fingerprint_one_structure(structure, labeled=labeled)
            fingerprints.append(fingerprint)

        if cumulative:
            cumul_fp = self._acumulate_fingerprints(fingerprints)
            if labeled and as_dataframe:
                plotdata = defaultdict(list)
                for entry in cumul_fp:
                    plotdata[entry.label["type"]].append(entry)
                labels = [
                    labeled_value.label["residue"].identifier
                    for labeled_value in plotdata["hbond"]
                ]
                df = pd.DataFrame.from_dict(
                    {k: [x.value for x in v] for (k, v) in plotdata.items()}
                )
                df.index = labels
                # change to eliminate redundant transpose
                if remove_non_interacting_residues:
                    df = df.T
                    df = df.loc[:, (df != 0).any(axis=0)]
                    return df.T
                else:
                    return df

    def _acumulate_fingerprints(self, fingerprints):
        """
        Calculate the cumulative fingerprint from fingerprints of multiple structures.

        Parameters
        ----------
        fingerprints = list of fingperprints to sum up
        """
        summed_fp = []
        # Iterate over the positions in the finger print
        # [ 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 ]
        # [ 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 ]
        # [ 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 ]
        # [ 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 ]
        #   ^ -->
        for position in zip(*fingerprints):
            total = sum([getattr(structure, "value", structure) for structure in position])
            if hasattr(position[0], "label"):  #  this is the labeled fingerprint!
                labels = [structure.label for structure in position]
                # Check all residues are equivalent!
                for attr in ("name", "seq_index", "chain"):
                    attrs = [getattr(label["residue"], attr) for label in labels]
                    assert all([attrs[0] == attr_ for attr_ in attrs[1:]])
                assert all([labels[0]["type"] == label["type"] for label in labels[1:]])
                old_res = labels[0]["residue"]
                residue = ProteinResidue(old_res.name, old_res.seq_index, old_res.chain)
                new_label = {"residue": residue, "type": labels[0]["type"]}
                summed_fp.append(_LabeledValue(value=total, label=new_label))
            else:
                summed_fp.append(total)
        return summed_fp

    def _calculate_fingerprint_one_structure(self, structure, labeled=False):
        """
        Calculate the interaction fingerprint for a single structure.

        Parameters
        ----------
        structure = structure object based on pdb file
        labeled = boolean deciding whether to make each fingerprint bit a labeled value or simple integer
        """
        fp_length = len(self.residue_indices) * len(self.interaction_types)
        fingerprint = []
        for index in self.residue_indices:
            residue = structure.get_residue_by(seq_index=index)
            counter = residue.count_interactions()
            for interaction in self.interaction_types:
                if labeled:
                    label = {"residue": residue, "type": interaction}
                    n_interactions = _LabeledValue(counter[interaction], label=label)
                else:
                    n_interactions = counter[interaction]
                fingerprint.append(n_interactions)
        assert len(fingerprint) == fp_length, "Expected length not matched"
        if not labeled:
            return np.asarray(fingerprint)
        return fingerprint

    def clear_fingerprint(self):
        self._fingerprint = None


class _LabeledValue:
    """
    This class is used to assign additional information to a value in the fingerprint.
    """

    def __init__(self, value, label):
        self.value = value
        self.label = label

    def __repr__(self):
        return "<LabeledValue {} labeled with object {}>".format(self.value, self.label)

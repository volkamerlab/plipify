"""
fingerprints.py
---------------

Factories that take a Structure and produce
an interaction fingerprint
"""

from collections import Counter

from core import Structure


class BaseFingerprint:
    """
    Base class for Fingerprint factories
    """

    def calculate(self, structure):
        """
        Parameters
        ----------
        structure: plipify.core.Structure
        """
        raise NotImplementedError("Implement in subclass!")


class MappedFingerprint:
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
        residue_indices,
        interactions=(
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
        self.residue_indices = residue_indices
        self.interactions = interactions

    def calculate(self, structure):
        fp_length = len(self.residue_indices) * len(self.interactions)
        fingerprint = []
        for index in self.residue_indices:
            residue = structure.get_residue_by(seq_index=index)
            interactions = residue.interactions  # TODO
            interaction_types = [interaction.__class__ for interaction in interactions]
            counter = Counter(interaction_types)
            fingerprint.extend(
                [
                    counter[Structure.INTERACTION_KEYS[interaction]]
                    for interaction in self.interactions
                ]
            )
        assert len(fingerprint) == fp_length, "Expected length not matched"
        return fingerprint

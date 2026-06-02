"""
fingerprints.py
---------------

Factories that take a Structure or multiple structures and produce
an interaction fingerprint.

"""
from collections import defaultdict, Counter
from tempfile import TemporaryDirectory
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd
from Bio.AlignIO.FastaIO import MultipleSeqAlignment, Seq, SeqRecord
from Bio.AlignIO import write as write_alignment, read as read_alignment

from .core import ProteinResidue


class InteractionFingerprint:
    """
    This class will take a protein-ligand structure,
    analyze its interactions and report a fingerprint
    per residue.

    """

    def __init__(
        self,
        interaction_types=(
            "hydrophobic",
            "hbond-don",
            "hbond-acc",
            "waterbridge",
            "saltbridge",
            "pistacking",
            "pication",
            "halogen",
            "metal",
            "covalent",
        ),
            split_backbone_sidechain_hbonds=False
    ):
        self.indices = None
        self.split_backbone_sidechain_hbonds = split_backbone_sidechain_hbonds
        self.interaction_types = interaction_types

    def count_interactions_with_hbond_split(self, residue):
        """
        The purpose of this function is to enable the split of sidechain and backbone hydrogen bonds
        """
        interaction_types = []
        for interaction in residue.interactions:
            int_type = interaction.shorthand
            if int_type == 'hbond-don' or int_type == 'hbond-acc':
                sidechain = interaction.interaction["SIDECHAIN"]
                if sidechain:
                    int_type += '-sc'
                else:
                    int_type += '-bb'
            interaction_types.append(int_type)
        counter = Counter(interaction_types)
        return counter

    def calculate_fingerprint(
        self,
        structures,
        residue_indices=None,
        labeled=True,
        cumulative=True,
        as_dataframe=False,
        remove_non_interacting_residues=False,
        remove_empty_interaction_types=False,
        ensure_same_sequence=True,
    ):
        """
        Cumulative interaction fingerprint for one or multiple structures.

        Parameters
        ----------
        structures : list of core.Structure objects
        residue_indices :  list of dict[int, <int or None>], or None
            list of dictionaries (one per structure) that maps
            unaligned position in sequence vs aligned position (after
            running Muscle on all the sequences). If not provided,
            it will be auto computed with `self.calculate_indices_mapping`
        labeled : bool
            decide whether to make each fingerprint bit a labeled value
            or simple integer
        cumulative : bool
            defines if the fp is a summed up fp or multiple structures
        as_dataframe : bool
            if true return fp as data_frame, else as array
        remove_non_interacting_residues : bool
            remove all fp bits that belong to residues for which
            there are no interactions
        remove_empty_interaction_types : bool
            remove interaction types that do not report any residues
        ensure_same_sequence : bool
            if true, check that all residues are identical for each position
            across structures. Only meaningful if cumulative=True
        """
        if residue_indices is None:
            residue_indices = self.calculate_indices_mapping(structures)
        if len(structures) != len(residue_indices):
            raise ValueError(
                f"Number of residue indices mappings ({len(residue_indices)}) "
                f"does not match number of structures ({len(structures)})"
            )
        # TODO: Some boolean paths are not covered here! Provide errors or implement missing path.
        fingerprints = []
        for structure, indices in zip(structures, residue_indices):
            try:
                fingerprints.append(
                    self._calculate_fingerprint_one_structure(
                        structure, indices.values(), labeled=labeled
                    )
                )
            except Exception as e:
                print(
                    f"! Warning, could not process structure {structure} "
                    f"due to error `{type(e).__name__}`: {e}"
                )

        if cumulative:
            cumul_fp = self._acumulate_fingerprints(
                fingerprints, ensure_same_sequence=ensure_same_sequence
            )
            if labeled and as_dataframe:
                plotdata = defaultdict(list)
                for entry in cumul_fp:
                    plotdata[entry.label["type"]].append(entry)
                df = pd.DataFrame.from_dict(
                    {k: [x.value for x in v] for (k, v) in plotdata.items()}
                )
                df.index = residue_indices[0].keys()
                # change to eliminate redundant transpose
                if remove_non_interacting_residues:
                    # remove all zero rows
                    df = df.loc[(df != 0).any(axis=1)]
                if remove_empty_interaction_types:
                    # remove all zero columns
                    df = df.loc[:, (df != 0).any(axis=0)]
                return df

        return fingerprints

    def _acumulate_fingerprints(self, fingerprints, ensure_same_sequence=True):
        """
        Calculate the cumulative fingerprint from fingerprints of multiple structures.

        Parameters
        ----------
        fingerprints = list of fingperprints to sum up
        ensure_same_sequence = if true, check that all residues are identical
            for each position across structures.
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
            if hasattr(position[0], "label"):  # this is the labeled fingerprint!
                labels = [structure.label for structure in position]
                if ensure_same_sequence:
                    # Check all residues are equivalent!
                    for attr in ("name", "seq_index", "chain"):
                        attrs = set(getattr(label["residue"], attr) for label in labels)
                        if len(attrs) > 1:
                            raise ValueError(
                                f"Residue at position {position[0].value} should be the same "
                                f"one across structures. Too many seen values for `{attr}`: {attrs}! "
                                f"Your structures might not be sequence-aligned."
                            )
                types = set(label["type"] for label in labels)
                if len(types) > 1:
                    raise ValueError(
                        f"Position {position[0].value} contains more than one type: {types}."
                    )
                old_res = labels[0]["residue"]
                residue = ProteinResidue(old_res.name, old_res.seq_index, old_res.chain)
                new_label = {"residue": residue, "type": labels[0]["type"]}
                summed_fp.append(_LabeledValue(value=total, label=new_label))
            else:
                summed_fp.append(total)
        return summed_fp

    def _calculate_fingerprint_one_structure(self, structure, indices, labeled=False):
        """
        Calculate the interaction fingerprint for a single structure.

        Parameters
        ----------
        structure = structure object based on pdb file
        indices = list of dict
            each dict contains kwargs that match Structure.get_residue_by
            so it can return a Residue object. For example:
            {"seq_index": 1, "chain": "A"}
        """
        empty_counter = Counter()
        fp_length = len(indices) * len(self.interaction_types)
        fingerprint = []
        for index_kwargs in indices:
            residue = structure.get_residue_by(**index_kwargs)
            if residue:
                if self.split_backbone_sidechain_hbonds:
                    counter = self.count_interactions_with_hbond_split(residue)
                else:
                    counter = residue.count_interactions()
            else:
                # FIXME: This is a bit hacky. Let's see if we can
                # come up with something more elegant.
                residue = ProteinResidue("GAP", 0, None)
                counter = empty_counter
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

    @staticmethod
    def calculate_indices_mapping(structures):
        """
        Align sequences of `structures` and provide a mapping of
        sequence positions to alignment positions (accounting for gaps).

        Only matching residue types are reported.

        Parameters
        ----------
        structures : list of plipify.core.Structure

        Returns
        -------
        indices : list of dict[int, int]
        """
        sequences = [s.sequence() for s in structures]
        maxlen = max(len(s) for s in sequences)
        # pad with ending -
        sequences = [s if len(s) == maxlen else (s + "-" * (maxlen - len(s))) for s in sequences]

        identifiers = [s.identifier for s in structures]
        records = [SeqRecord(Seq(s), id=i) for s, i in zip(sequences, identifiers)]
        unaligned = MultipleSeqAlignment(records)
        with TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            infile = str(tmp / "in.fasta")
            outfile = str(tmp / "out.fasta")
            logfile = str(tmp / "log.txt")
            write_alignment(unaligned, infile, "fasta")
            subprocess.run(['muscle', '-align', infile, '-output', outfile, '-log', logfile])
            aligned = read_alignment(outfile, "fasta")

        offset = unaligned.get_alignment_length() - aligned.get_alignment_length()

        old2new = []
        for old in unaligned:
            new = next(n for n in aligned if n.id == old.id)
            new = "-" * offset + new
            old2new.append({})
            gaps = 0
            keep = None
            for i in range(unaligned.get_alignment_length()):
                oldchar = old[i]
                newchar = new[i]
                if oldchar == newchar:
                    if oldchar == "-":
                        continue
                    old2new[-1][i + 1] = {"seq_index": i + 1, "chain": "any"}
                    if keep is not None:
                        old2new[-1][keep] = {"seq_index": keep + gaps, "chain": "any"}
                        gaps = 0
                        keep = None
                else:
                    keep = i + 1
                    gaps += 1
        return old2new


class _LabeledValue:
    """
    This class is used to assign additional information to a value in the fingerprint.
    """

    def __init__(self, value, label):
        self.value = value
        self.label = label

    def __repr__(self):
        return "<LabeledValue {} labeled with object {}>".format(self.value, self.label)

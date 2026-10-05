"""
fingerprints.py
---------------

Factories that take a Structure or multiple structures and produce
an interaction fingerprint.

"""
import subprocess
from collections import defaultdict, Counter
from tempfile import TemporaryDirectory
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.AlignIO.FastaIO import Seq, SeqRecord
from Bio.AlignIO import read as read_alignment
from Bio.SeqIO import write as write_sequences

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
        Align sequences of `structures` with MUSCLE and map every alignment
        column to the residue each structure has at that column.

        Only columns where all structures with a residue agree on the
        residue type are reported. Structures with a gap at a reported
        column get ``seq_index=None``, which is fingerprinted as a GAP.

        Parameters
        ----------
        structures : list of plipify.core.Structure

        Returns
        -------
        indices : list of dict[int, dict]
            One dict per structure, all with the same keys (1-based alignment
            columns). Values are kwargs for ``Structure.get_residue_by``.
        """
        # sequence() has "-" for unresolved residues, so the n-th letter
        # of a sequence sits at seq_index (position in the string + 1)
        sequences = [s.sequence() for s in structures]
        seq_indices = [[i + 1 for i, c in enumerate(s) if c != "-"] for s in sequences]
        ungapped = [s.replace("-", "") for s in sequences]

        # positional identifiers: structure identifiers may be missing or duplicated
        identifiers = [f"s{i}" for i in range(len(structures))]
        records = [SeqRecord(Seq(s), id=i) for s, i in zip(ungapped, identifiers)]

        with TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            infile = str(tmp / "in.fasta")
            outfile = str(tmp / "out.fasta")
            logfile = str(tmp / "log.txt")
            write_sequences(records, infile, "fasta")
            subprocess.run(['muscle', '-align', infile, '-output', outfile, '-log', logfile])
            aligned = read_alignment(outfile, "fasta")

        # muscle reorders sequences in its output
        aligned_by_id = {record.id: str(record.seq) for record in aligned}
        rows = [aligned_by_id[i] for i in identifiers]
        for row, seq in zip(rows, ungapped):
            assert row.replace("-", "") == seq, "MUSCLE altered an input sequence"

        # for each structure, the seq_index at each alignment column (None for gaps)
        columns = []
        for row, indices in zip(rows, seq_indices):
            residues = iter(indices)
            columns.append([None if c == "-" else next(residues) for c in row])

        old2new = [{} for _ in structures]
        for col in range(aligned.get_alignment_length()):
            residue_types = {row[col] for row in rows} - {"-"}
            if len(residue_types) != 1:
                continue
            for mapping, structure_columns in zip(old2new, columns):
                mapping[col + 1] = {"seq_index": structure_columns[col], "chain": "any"}
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

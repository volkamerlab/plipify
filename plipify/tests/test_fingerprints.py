"""
Tests for the sequence alignment mapping in plipify.fingerprints.

Run with `pytest -s` to see the printed alignments.
"""

import shutil

import pytest

from plipify.fingerprints import InteractionFingerprint

pytestmark = pytest.mark.skipif(shutil.which("muscle") is None, reason="muscle not installed")


class FakeStructure:
    """Minimal stand-in for plipify.core.Structure: only what the mapping needs."""

    def __init__(self, identifier, sequence):
        self.identifier = identifier
        self._sequence = sequence

    def sequence(self):
        return self._sequence

    def __repr__(self):
        return f"FakeStructure({self.identifier!r}, {self._sequence!r})"


def map_and_show(*structures):
    """Run the mapping, print it as a small table and return seq_index lists per structure."""
    mappings = InteractionFingerprint.calculate_indices_mapping(list(structures))
    columns = list(mappings[0])
    print()
    print("column    " + "".join(f"{c:>3}" for c in columns))
    for structure, mapping in zip(structures, mappings):
        indices = [mapping.get(c, {}).get("seq_index") for c in columns]
        letters = ["-" if i is None else structure.sequence()[i - 1] for i in indices]
        print(f"{structure.identifier:<10}" + "".join(f"{letter:>3}" for letter in letters))
        print("seq_index " + "".join(f"{'.' if i is None else i:>3}" for i in indices))
    return [list(mapping) for mapping in mappings], [
        [v["seq_index"] for v in mapping.values()] for mapping in mappings
    ]


def test_identical_sequences():
    """
    Same sequence twice: no gaps, column i is residue i in both.

    FAIL if any column maps to a different residue number.
    """
    keys, (a, b) = map_and_show(FakeStructure("a", "MKTAYW"), FakeStructure("b", "MKTAYW"))

    assert a == b == [1, 2, 3, 4, 5, 6]


def test_insertion():
    """
    b has an extra G, so muscle puts a gap in a at column 3:

        a  MK-TAYW
        b  MKGTAYW

    FAIL if a and b get different keys, or if a loses residues after the gap.
    The old code did both: it compared a before and after alignment letter by
    letter, so after the gap nothing matched and a kept only residues 1-2.
    """
    keys, (a, b) = map_and_show(FakeStructure("a", "MKTAYW"), FakeStructure("b", "MKGTAYW"))

    assert keys[0] == keys[1], "a and b must have the same keys"
    assert a == [1, 2, None, 3, 4, 5, 6]
    assert b == [1, 2, 3, 4, 5, 6, 7]


def test_missing_residues():
    """
    b has residues 1, 2 and 6 missing in the crystal structure ("-" in sequence()):

        a  MKTAYWQR
        b  --TAY-QR

    FAIL if b's residues lose their real numbers (3, 4, 5, 7, 8).
    """
    keys, (a, b) = map_and_show(FakeStructure("a", "MKTAYWQR"), FakeStructure("b", "--TAY-QR"))

    assert keys[0] == keys[1], "a and b must have the same keys"
    assert a == [1, 2, 3, 4, 5, 6, 7, 8]
    assert b == [None, None, 3, 4, 5, None, 7, 8]


def _seq_indices(structures, **kwargs):
    """Mapping as {column: [seq_index per structure]} for the given options."""
    mappings = InteractionFingerprint.calculate_indices_mapping(list(structures), **kwargs)
    return {col: [m[col]["seq_index"] for m in mappings] for col in mappings[0]}


# nine identical sequences and one with T instead of A at position 4
_MISMATCH = [FakeStructure(f"wt{i}", "MKTAYW") for i in range(9)] + [FakeStructure("odd", "MKTTYW")]


def test_strict_consensus_drops_mismatch_column():
    """
    min_consensus=1.0 (the original rule): one odd residue removes column 4 for everyone.

    FAIL if column 4 is kept.
    """
    cols = _seq_indices(_MISMATCH, min_consensus=1.0)
    assert sorted(cols) == [1, 2, 3, 5, 6]


def test_consensus_masks_mismatch():
    """
    min_consensus=0.9: column 4 is 9/10 = 90% A, so it is kept, and the odd
    structure gets None there (contributes nothing) while the others keep residue 4.

    FAIL if column 4 is dropped or the odd structure is not masked.
    """
    cols = _seq_indices(_MISMATCH, min_consensus=0.9)
    assert sorted(cols) == [1, 2, 3, 4, 5, 6]
    assert cols[4] == [4] * 9 + [None]


def test_no_masking_keeps_every_residue():
    """
    min_consensus=0.0 and mask_mismatches=False: every column is kept and every
    structure maps to its own residue, including the odd T at position 4.

    FAIL if any structure gets None (there are no gaps in this example).
    """
    cols = _seq_indices(_MISMATCH, min_consensus=0.0, mask_mismatches=False)
    assert sorted(cols) == [1, 2, 3, 4, 5, 6]
    assert all(v == [col] * 10 for col, v in cols.items())

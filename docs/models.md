# Plipify's object model

We have several important modules:

* `core`: this module is responsible for I/O, PLIP execution
    * This is the conversion from `plip_fingerprints`. Other content is in `fingerprints`
* `fingerprints`: this module is responsible for PLIP->bitstring conversion
* `reports`: this module takes a PLIP report and outputs several representations or plots (stacked bars, heatmaps, interaction table).
    * Right now the functionality for this module is on fp_visual.


## CORE

* BaseInteraction as a base class for all kinds of nonbonded interactions, as defined by PLIP. This behaves as dictionary with a type, basically.
* Structural objects
    * BaseResidue, with subclasses: ProteinResidue, LigandResidue
        * ProteinResidue is part of a larger, linear sequence (as in a protein sequence)
        * LigandResidue is an isolated entity (no sequence)
    * BindingSite, is a collection of ProteinResidues (those relevant for the interaction) and LigandResidues
    * Structure: A collection of ProteinResidue objects, and potentially LigandResidue objects, which can define one or more BindingSite objects.

## FINGERPRINTS

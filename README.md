PLIPify: Protein-Ligand Interaction Frequencies across Multiple Structures
==============================
[![Actions Status](https://github.com/volkamerlab/plipify/workflows/CI/badge.svg)](https://github.com/volkamerlab/plipify/actions)

_Powered by: [Volkamer lab](http://volkamerlab.org/)_

:construction_worker: **Note: This repo is still work-in-progress.**

### Background

Protein-ligand interactions are an essential part of research in structural bioinformatics and drug discovery.
Tools such as the Protein-Ligand Interaction Profiler ([PLIP](https://plip-tool.biotec.tu-dresden.de/plip-web/plip/index)) enable us to get detailed interaction profiles for single structures. However, combining this data for multiple structures of a protein to identify possible interaction hotspots across them, e.g. when bound to different ligands, remains difficult.
The aim of `plipify` is to create and visualize a fingerprint that represents the protein-ligand interaction frequencies over multiple structures of the same protein.

Note that full credits for protein-ligand profile computation go to [PLIP](https://plip-tool.biotec.tu-dresden.de/plip-web/plip/index) [1]. `plipify` provides a wrapper around PLIP, which allows to digest multiple structures at once, performs the mapping of the individual profiles to fingerprints and reports protein-ligand interaction frequencies.

> [1] Salentin, S. et al. [PLIP: fully automated protein-ligand interaction profiler](https://academic.oup.com/nar/article/43/W1/W443/2467865). Nucl. Acids Res., 2015, 43 (W1): W443-W447.

### Project 01

#### Exploring SARS-CoV-2 Main Protease Interaction Hotspots Using `plipify`

This is part of a community effort to rapidly find new hits to target the virus main protease.

The COVID-19 (coronavirus disease 2019) pandemic, caused by severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2) has become a global health emergency and there is still an urgent need for effective anti-COVID drugs (see [COVID-19 Science Report: Therapeutics](https://sph.nus.edu.sg/covid-19/research/)). A promising target is the main protease Mpro of SARS-CoV-2 (first crystal structure, 01-2020, [6LU7](https://www.rcsb.org/structure/6LU7)). [UK’s Diamond Light Source](https://www.diamond.ac.uk/covid-19/for-scientists/Main-protease-structure-and-XChem/Downloads.html) has performed a large crystal-based fragment screen on Mpro, yielding by now over 400 complex structures.
A huge crowed sourcing campaign - the [COVID moonshot](https://postera.ai/moonshot) project - was invoked by [PostEra](https://covid.postera.ai/covid) and partners which encourage researchers from around the world to use the fragment hits as a starting point and contribute, amongst others, by suggesting potential inhibitors (effective and easy-to-make).

Besides other structures-based attempts (see our [repo](https://github.com/volkamerlab/covid19-SBapproach)), we applied `plipify` to the set of available structures, to generate more insides about the common bindig modes.

For more details, please see the [fragalysis.ipynb](https://github.com/volkamerlab/plipify/blob/master/projects/01/fragalysis.ipynb).

### Project 02/03

Coming soon ...

### Installation using conda

#### Prerequisite
Anaconda and Git should be pre-installed. See [Anaconda's website](https://www.anaconda.com/products/individual) and [Git's website](https://git-scm.com/downloads) for download.

#### How to

1. Clone the github repository:

```console
git clone https://github.com/volkamerlab/plipify.git
```

2. Change directory:

```console
cd plipify
```

3. Create the conda environment:

```console
conda env create -f devtools/conda-envs/test_env.yaml
```

4. Activate the environment:

```console
conda activate plipify
```

5. Install plipify package:

```console
pip install -e .
```

### Contributors

* Methodology: Jaime Rodríguez-Guerra, Franziska Fritz, Andrea Volkamer
    * plipify nb: Franziska Fritz, Jaime Rodríguez-Guerra
* Projects:
    * 01: One protein against many ligands: Jaime Rodríguez-Guerra, William Glass, Andrea Volkamer
    * WIP: 02: One ligand against several targets: Jaime Rodríguez-Guerra, David Schaller, Andrea Volkamer
    * WIP: 03: Automated interaction statistics for any protein in the PDB: David Schaller, Jaime Rodríguez-Guerra, Andrea Volkamer

### Repository structure and important files

```
|-- LICENSE
|-- README.md
|-- devtools    <- environment file
|-- plipify     <- plipify code
|-- projects
|   |-- 01      <- One protein against many ligands
|   |-- 02      <- One ligand against several targets
```

### Copyright

Copyright (c) 2021, Volkamer Lab

#### Acknowledgements

Project based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.2.

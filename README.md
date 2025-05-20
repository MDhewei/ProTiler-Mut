[![](https://img.shields.io/badge/Pypi-v0.1.0-519dd9.svg)](https://pypi.org/project/protiler-mut/)
[![License: GUN](https://img.shields.io/badge/License-GUN-yellow.svg)](https://github.com/MDhewei/MOFF/blob/master/LICENSE)
![](https://img.shields.io/badge/language-python-orange.svg)

# ProTiler-Mut

## Introduction of ProTiler-Mut 
<div align="left"><img src="StaticFiles/ProTiler-Mut_logo.png"  height="140" width="1000"></div>
                                     

    Hi, this is ProTiler-Mut, a computational pipeline designed for comprehensive analysis of tiling mutagenesis screen data
    
    I have three major functions: 1). Clustering and Categorization of functional mutations from tiling mutagenesis screens
                                  2). "3D-RRA" module for robust identification of functional substructures from identified mutation clusters
                                  3). PPI-mapping for specific mutation or substructure to identify mutaiton-associated PPIs
    
    
    Hope you enjoy playing with me ^o^!
                                     
    Any questions or bugs, please contact hwkobe.1027@gmail.com or whe3@mdanderson.org
                                     

## Installation

#### If Anaconda (or miniconda) is not installed with Python 3, it is highly recommended to download and install Python3 Anaconda from here: https://www.anaconda.com/download/

### Dependencies

> **Python Packages**:
> Following are the specific versions used when developing the tool, other versions should also be OK,
> if it is not,please install the corresponding version instead

     biopython==1.79, matplotlib==3.5.3, mygene==3.2.2, numpy==1.21.6, 
     pandas==1.3.5, Requests==2.32.3, rich==14.0.0, scikit_learn==0.20.0, scipy==1.7.3, 
     seaborn==0.13.2, setuptools==68.0.0, statsmodels==0.13.5, umap_learn==0.5.3

> **Pymol is required for ProTiler-Mut, install it using following command:**
```console   
conda install -c conda-forge pymol-open-source
```
> **On macOS, you need this command in addition:**
```console  
pip install PyQt5
```

### There are three ways to install ProTiler-Mut
#### Install ProTiler-Mut through pip
 ```console     
 pip install ProTiler-Mut
 ```
    
#### OR you can install ProTiler-Mut through git clone
```console   
git clone https://github.com/MDhewei/ProTiler-Mut.git
cd ProTiler-Mut
pip install -r requirements.txt .
```

#### OR you can install ProTiler-Mut through Docker (Still under development)
With Docker no installation is required, the only dependence is Docker itself. Users will completely get rid of all the installation and configuration issues. Docker will do all the dirty work for you!

Docker can be downloaded freely from here: https://store.docker.com/search?offering=community&type=edition

To get an image of ProTiler-Mut, simply execute the following command:
```console   
$ docker pull MDhewei/ProTiler-Mut
 ```

## How to use ProTiler-Mut
### There are two ways to use ProTiler-Mut

### 1. Integrate ProTiler-Mut Functional Modules into Your Python Code (suitable for users with basic python programing skill)
ProTiler-Mut can be seamlessly imported into your Python environment and integrated directly into custom scripts or pipelines.
A demo notebook demonstrating how to use ProTiler-Mut as a Python package is available here: [Demo.ipynb](https://github.com/MDhewei/ProTiler-Mut/blob/main/Demo.ipynb)

### 2. Use ProTiler-Mut via Command Line as a Standalone Tool (Other users, can be excecuted in any terminal in both Linux or MacOS system)
#### (1) ProTiler-Mut cluster: Perform the clustering and categorization of functional mutations

       usage: protiler-mut cluster [-h] -i INPUTFILE -g GENE_ID -s SAMPLES -c CONTROL [-p PDB] [-n N_CLUSTERS] [-m METHOD]
                                      [-d METRIC] [--pdf-report PDF_REPORT] [-o OUTPUT_FOLDER]

       optional arguments:
       -h, --help            show this help message and exit

       Required arguments for clustering.:

       -i INPUTFILE, --inputfile INPUTFILE
                        The inputfile contains information of tiling mutagenesis screens including symbol of target
                        gene(s),targeted residue position, mutation types and phenotypic scores. Accept .txt, .cvs or
                        .xlsx fileformats. 
       -g GENE_ID, --gene_id GENE_ID
                        The symbol of targeted protein-coding gene, for example: ERCC2
       -s SAMPLES, --samples SAMPLES
                        Comma-separated sample column names.eg., "CISP,OLAP,DOX,CPT"
       -c CONTROL, --control CONTROL
                        Comma-separated control column names.eg., T0

       Optional arguments for clustering.:

       -p PDB, --pdb PDB     File path to the PDB of targeted protein structure.
       -n N_CLUSTERS, --n-clusters N_CLUSTERS
                        Number of clusters for clustering analysis.
       -m METHOD, --method METHOD
                        Clustering linkage method (default: average).
       -d METRIC, --metric METRIC
                        Clustering metric (default: euclidean).
       --pdf-report PDF_REPORT
                        Generate pdf report of clustering, visualization and annotation.
        -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                        Output folder for saving the results.

#### (2) ProTiler-Mut 3d-rra: Perform "3D-RRA" to call significant substructures in specific mutation clusters

       usage: protiler-mut 3d-rra [-h] -g GENE_ID -i INPUTFILE -p PDB -n N [-r NUM_PERMUTATIONS] [-t1 DISTANCE_THRESHOLD1]
                              [-t2 DISTANCE_THRESHOLD2] [-o OUTPUT_FOLDER]

       optional arguments:
       -h, --help            show this help message and exit

       Required arguments for 3D-RRA.:

       -g GENE_ID, --gene_id GENE_ID
                        The symbol of targeted protein-coding gene, for example: ERCC2
       -i INPUTFILE, --inputfile INPUTFILE
                        Path output tables file generated in cluster module which annotat the significant mutations, their
                        cluster assignment and residue position
       -p PDB, --pdb PDB     File path to the PDB of targeted protein structure
       -n N, --n N           Number of mutation samples for RRA analysis

       Optional arguments for 3D-RRA.:

       -r NUM_PERMUTATIONS, --num-permutations NUM_PERMUTATIONS
                        Number of permutations (default: 10000).
       -t1 DISTANCE_THRESHOLD1, --distance-threshold1 DISTANCE_THRESHOLD1
                        Distance threshold to identify clusters of seed mutations on 3D structure(default: 10.0 Å).
       -t2 DISTANCE_THRESHOLD2, --distance-threshold2 DISTANCE_THRESHOLD2
                        Distance threshold to identify surrounding signals near identified seed mutations(default: 5.0 Å).
       -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                         Output folder for results.

 #### (3). ProTiler-Mut ppi-mapping: Perform PPI-mapping for specific mutation or substructure to identify mutaiton-associated PPIs

       usage: protiler-mut.py ppi-mapping [-h] -g GENE_ID -i INPUTFILE -f PDB_FILES -b CHAINS [-t DISTANCE_THRESHOLD]
                                   [-o OUTPUT_FOLDER]

       optional arguments:
       -h, --help            show this help message and exit

       Required arguments for PPI-mapping.:

       -g GENE_ID, --gene_id GENE_ID
                        The symbol of targeted protein-coding gene, for example: ERCC2
       -i INPUTFILE, --inputfile INPUTFILE
                        Path output tables file generated in cluster module which annotat the significant mutations, their
                        cluster assignment and residue position, See example file
       -f PDB_FILES, --pdb-files PDB_FILES
                        Comma-separated list of paths of protein complex PDB files involving the target protein.
       -b CHAINS, --chains CHAINS
                        Comma-separated list of corresponding chain IDs of the target protein(e.g., A,B,A).

       Optional arguments for PPI mapping.:

       -t DISTANCE_THRESHOLD, --distance-threshold DISTANCE_THRESHOLD
                        Distance threshold to determine whether two residues interact between among different
                        chains(default: 5.0 Å).
       -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                        Output folder for results.

 ### For a detailed tutorial to run ProTiler-Mut, please refer to the video at Youtube https://www.youtube.com/@bioinforbricker  

 


[![](https://img.shields.io/badge/Pypi-v1.2.3-519dd9.svg)](https://pypi.org/project/MOFF/)
[![License: GUN](https://img.shields.io/badge/License-GUN-yellow.svg)](https://github.com/MDhewei/MOFF/blob/master/LICENSE)
![](https://img.shields.io/badge/language-python-orange.svg)

# ProTiler-Mut
Computational pipeline for comprehensive analysis of tiling mutagenesis screen data 

### Installation 
 Step1: Install Anaconda (highly recomended)
    
 ```console
 wget https://repo.continuum.io/archive/Anaconda2-2018.12-Linux-x86_64.sh 
 bash Anaconda2-2018.12-Linux-x86_64.sh 
 ```

 Step2: Install ProTiler-Mut through pip
 ```console     
 pip install ProTiler-Mut
 ```
    
 Step3: **OR** you can install ProTiler-Mut through git clone
 ```console   
 git clone https://github.com/MDhewei/ProTiler-Mut.git
 cd ProTiler-Mut
 python setup.py install
 ```

Step4: **OR** you can install ProTiler-Mut through Docker

## How to use ProTiler-Mut

### 1. ProTiler-Mut cluster: Perform the clustering and categorization of functional mutations

       usage: protiler-mut.py cluster [-h] -i INPUTFILE -g GENE_ID -s SAMPLES -c CONTROL [-p PDB] [-n N_CLUSTERS] [-m METHOD]
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

### 2. ProTiler-Mut 3d-rra: Perform "3D-RRA" to call significant substructures in specific mutation clusters

       usage: protiler-mut.py 3d-rra [-h] -g GENE_ID -i INPUTFILE -p PDB -n N [-r NUM_PERMUTATIONS] [-t1 DISTANCE_THRESHOLD1]
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

 ### 3. ProTiler-Mut ppi-mapping: Perform PPI-mapping for specific mutation or substructure to identify mutaiton-associated PPIs

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

   

 


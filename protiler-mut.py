#!/home/whe/miniconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on March 2 Sunday 13:32:33 2025
@author: Wei He
@email: whe3@mdanderson.org
ProTiler-Mut: A tool for analysis for tiling mutagenesis screen data including clustering, 3D-RRA, and PPI mapping.
"""

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # Suppress TensorFlow INFO and WARNING messages

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

__version__ = "1.0.0"

import argparse, logging, sys, pkg_resources
import pandas as pd

# Import all functions from your three modules.
from Clustering import *
from ThreeD_RRA import *
from PPI_mapping import *


from rich.console import Console
from rich.panel import Panel
from rich.text import Text


def print_banner():
    """
    Print a colorful ASCII art banner using the rich library.
    """
    console = Console()
    banner_text = r"""
      ██████╗ ██████╗  ██████╗ ████████╗██╗██╗     ███████╗██████╗     ███╗   ███╗██╗   ██╗████████╗
      ██╔══██╗██╔══██╗██╔═══██╗╚══██╔══╝██║██║     ██╔════╝██╔══██╗    ████╗ ████║██║   ██║╚══██╔══╝
    ██████╔╝██████╔╝██║   ██║   ██║   ██║██║     █████╗  ██████╔╝    ██╔████╔██║██║   ██║   ██║    
    ██╔═══║ ██╔══██╗██║   ██║   ██║   ██║██║     ██╔══╝  ██╔══██╗    ██║╚██╔╝██║██║   ██║   ██║  
    ██║     ██║  ██║╚██████╔╝   ██║   ██║███████╗███████╗██║  ██║    ██║ ╚═╝ ██║╚██████╔╝   ██║    
    ╚═╝     ╚═╝  ╚═╝ ╚═════╝    ╚═╝   ╚═╝╚══════╝╚══════╝╚═╝  ╚═╝    ╚═╝     ╚═╝ ╚═════╝    ╚═╝       
    """ 
   
    panel = Panel(Text(banner_text, justify="center", style="bold green"), title="ProTiler-Mut", 
                  subtitle="Advanced Tool for Tiling Mutagenesis Screen Data Analysis", expand=False)
    console.print(panel)


def ProTilerMain():
    print_banner()
    import pkg_resources
    import os

    # Get the full path to the 'StaticFiles' directory within your package.
    #static_files_path = pkg_resources.resource_filename(__name__, 'StaticFiles')
    #print("Static files are located at:", static_files_path)
#     logging.info('Loading and reading genome-wide exons information files...')
#     exons_dic_file = os.path.join(static_files_path,'Exons_dic_hg38.json')
#     exons_dic = json.loads( open(exons_dic_file).read())

#     logging.info('Loading and reading genome-wide domain annotation...')
#     domain_dic_file =  os.path.join(static_files_path,'Pfam_domain_proteome.json') 
#     domain_dic = json.loads( open(domain_dic_file).read())

#     logging.info('Loading and reading clinvar dataset...')
#     clinvar_dic_file =  os.path.join(static_files_path,'Clinvar_genome.json') 
#     clinvar_dic = json.loads( open(domain_dic_file).read())

#     logging.info('Loading and reading genome-wide post-translational modification file...')
#     po_file = os.path.join(static_files_path,'Phosphosite_genome.json')
#     po_dic = json.loads( open(po_file).read())
    
    ## Set logging format
    logging.basicConfig(level=logging.DEBUG,  
                        format='%(levelname)s:%(asctime)s @%(message)s',  
                        datefmt='%a, %d %b %Y %H:%M:%S',
                        filemode='a')
    
    ## Add arguments for user input with command line.
    parser = argparse.ArgumentParser(
        description='ProTiler-Mut: An advanced tool for comprehensive analysis of tiling mutagenesis screen data.'
    )
    
    ## Add subparsers for the three sub-functions.
    subparsers = parser.add_subparsers(help='Commands to run ProTiler-Mut', dest='subcmd')
    
    # ----- Cluster sub-command -----
    
    subp_cluster = subparsers.add_parser('cluster', help='Perform clustering of mutations, 1D&3D visualization and amino acid annotation')
    req_cluster = subp_cluster.add_argument_group(title='Required arguments for clustering.', description='')
    
    req_cluster.add_argument('-i','--inputfile',type=str,help='The inputfile contains information of tiling mutagenesis screens including symbol of target gene(s),targeted residue position, mutation types and phenotypic scores. Accept .txt, .cvs or .xlsx fileformats',required=True)
    
    req_cluster.add_argument('-g','--gene_id',type=str,help='The symbol of targeted protein-coding gene, for example: ERCC2',required=True)
    req_cluster.add_argument('-s','--samples', type=str, required=True, help='Comma-separated sample column names.eg., "CISP,OLAP,DOX,CPT"')
    req_cluster.add_argument('-c','--control', type=str, required=True, help='Comma-separated control column names.eg., T0')
    
    opt_cluster = subp_cluster.add_argument_group(title='Optional arguments for clustering.', description='')
    opt_cluster.add_argument('-p','--pdb', type=str, default=None, help='File path to the PDB of targeted protein structure.')
    opt_cluster.add_argument('-n','--n-clusters', type=int, default=3, help='Number of clusters for clustering analysis.')
    opt_cluster.add_argument('-m','--method', type=str, default='average', help='Clustering linkage method (default: average).')
    opt_cluster.add_argument('-d','--metric', type=str, default='euclidean', help='Clustering metric (default: euclidean).')
    opt_cluster.add_argument('--pdf-report', help='Generate pdf report of clustering, visualization and annotation.')
    opt_cluster.add_argument('-o','--output-folder', type=str, default='results', help='Output folder for saving the results.')
    
    # ----- 3D-RRA sub-command -----
    subp_rra = subparsers.add_parser('3d-rra', help='Perform 3D-RRA analysis to identify "hotspot" funcitonal substructure')
    req_rra = subp_rra.add_argument_group(title='Required arguments for 3D-RRA.', description='')
    req_rra.add_argument('-g','--gene_id',type=str,help='The symbol of targeted protein-coding gene, for example: ERCC2',required=True)
    req_rra.add_argument('-i','--inputfile', type=str, required=True, help='Path output tables file generated in cluster module which annotat the significant mutations, their cluster assignment and residue position')
    req_rra.add_argument('-p','--pdb', type=str, required=True, help='File path to the PDB of targeted protein structure')
    req_rra.add_argument('-n','--n', type=int, required=True, help='Number of mutation samples for RRA analysis')
    
    opt_rra = subp_rra.add_argument_group(title='Optional arguments for 3D-RRA.', description='')
    opt_rra.add_argument('-r','--num-permutations', type=int, default=10000, help='Number of permutations (default: 10000).')
    opt_rra.add_argument('-t1','--distance-threshold1', type=float, default=10.0, help='Distance threshold to identify clusters of seed mutations on 3D structure(default: 10.0 Å).')
    opt_rra.add_argument('-t2','--distance-threshold2', type=float, default=5.0, help='Distance threshold to identify surrounding signals near identified seed mutations(default: 5.0 Å).')
    opt_rra.add_argument('-o','--output-folder', type=str, default='results', help='Output folder for results.')
    
    # ----- PPI-mapping sub-command -----
    subp_ppi = subparsers.add_parser('ppi-mapping', help='Perform PPI mapping to identify mutation affected PPI interfaces.')
    req_ppi = subp_ppi.add_argument_group(title='Required arguments for PPI-mapping.', description='')
    
    req_ppi.add_argument('-g','--gene_id',type=str,help='The symbol of targeted protein-coding gene, for example: ERCC2',required=True)
    req_ppi.add_argument('-i','--inputfile', type=str, required=True, help='Path output tables file generated in cluster module which annotat the significant mutations, their cluster assignment and residue position')
    req_ppi.add_argument('-f','--pdb-files', type=str, required=True, help='Comma-separated list of paths of protein complex PDB files involving the target protein.')
    req_ppi.add_argument('-b','--chains', type=str, required=True, help='Comma-separated list of corresponding chain IDs of the target protein(e.g., A,B,A).')
        
    opt_ppi = subp_ppi.add_argument_group(title='Optional arguments for PPI mapping.', description='')
    opt_ppi.add_argument('-t','--distance-threshold', type=float, default=5.0, help='Distance threshold to determine whether two residues interact between among different chains(default: 5.0 Å).')
    opt_ppi.add_argument('-o','--output-folder', type=str, default='results', help='Output folder for results.')
    
    args = parser.parse_args()
    
    if args.subcmd == "cluster":
        
        df_gene = pd.read_csv(args.gene_data)
        sample_list = args.samples.split(',')
        control_list = args.control.split(',')
        df_clust = clustering(df_gene, args.gene, sample_list, control_list, args.n_clusters, args.method, args.metric, args.output_folder)
        Annotation(df_clust,gene,output_folder,pdb_path=None)
        
        logging.info("Clustering analysis completed. Results saved in %s", args.output_folder)
    
    elif args.subcmd == "3d-rra":
        perform_3d_rra(args.gene, args.gene_data, args.pdb, args.n, args.num_permutations, args.output_folder)
        logging.info("3D-RRA analysis completed. Results saved in %s", args.output_folder)
    
    elif args.subcmd == "ppi-mapping":
        mut_numbers = [int(x) for x in args.mutations.split(',')]
        chains = args.chains.split(',')
        pdb_files = args.pdb_files.split(',')
        table = build_ppi_interface_table(mut_numbers, chains, pdb_files, args.distance_threshold, args.output_folder)
        logging.info("PPI mapping analysis completed. Results saved in %s", args.output_folder)
        print(table)
    
    else:
        parser.print_help()

if __name__ == "__main__":
    ProTilerMain()

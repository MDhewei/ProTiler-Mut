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
from .cluster import *
from .threeD_rra import *
from .ppi_mapping import *


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


def main():
    print_banner()
    import pkg_resources
    import os
    
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
    opt_cluster.add_argument('-t','--significance_threshold', type=int, default=2, help='Threshold to determine significant mutation in the screen')
    opt_cluster.add_argument('-p','--pdb', type=str, default=None, help='File path to the PDB of targeted protein structure.')
    opt_cluster.add_argument('-n','--n-clusters', type=int, default=3, help='Number of clusters for clustering analysis.')
    opt_cluster.add_argument('-m','--method', type=str, default='average', help='Clustering linkage method (default: average).')
    opt_cluster.add_argument('-d','--metric', type=str, default='euclidean', help='Clustering metric (default: euclidean).')
    opt_cluster.add_argument('--pdf-report', help='Generate pdf report of clustering, visualization and annotation.')
    opt_cluster.add_argument('-o','--outputfolder', type=str, default='Clustering_results', help='Output folder for saving the results.')
    
    
    # ----- 3D-RRA sub-command -----
    subp_rra = subparsers.add_parser('3d-rra', help='Perform 3D-RRA analysis to identify "hotspot" funcitonal substructure')
    req_rra = subp_rra.add_argument_group(title='Required arguments for 3D-RRA.', description='')
    req_rra.add_argument('-i','--inputfile',type=str,help='The inputfile contains information of tiling mutagenesis screens including symbol of target gene(s),targeted residue position, mutation types and phenotypic scores. Accept .txt, .cvs or .xlsx fileformats',required=True)
    req_rra.add_argument('-g','--gene_id',type=str,help='The symbol of targeted protein-coding gene, for example: ERCC2',required=True)
    req_rra.add_argument('-c','--clustertable', type=str, required=True, help='Path output tables file generated in cluster module which annotat the significant mutations, their cluster assignment and residue position')
    req_rra.add_argument('-p','--pdb', type=str, required=True, help='File path to the PDB of targeted protein structure')
    #req_rra.add_argument('-n','--n', type=int, required=True, help='Number of mutation samples for RRA analysis')
    
    opt_rra = subp_rra.add_argument_group(title='Optional arguments for 3D-RRA.', description='')
    opt_rra.add_argument('-t','--significance_threshold', type=int, default=2, help='Threshold to determine significant mutation in the screen')
    opt_rra.add_argument('-r','--num_permutations', type=int, default=10000, help='Number of permutations (default: 10000).')
    opt_rra.add_argument('-t3','--top_cutoff', type=float, default=0.25, help='Top percentile signals used for RRA analysis (default: 0.25.')
    opt_rra.add_argument('-t1','--distance_threshold1', type=float, default=15.0, help='Distance threshold to identify clusters of seed mutations on 3D structure(default: 15.0 Å).')
    opt_rra.add_argument('-t2','--distance_threshold2', type=float, default=5.0, help='Distance threshold to identify surrounding signals near identified seed mutations(default: 5.0 Å).')
    opt_rra.add_argument('-o','--outputfolder', type=str, default='3D_RRA_results', help='Output folder for results.')
    
    # ----- PPI-mapping sub-command -----
    subp_ppi = subparsers.add_parser('ppi-mapping', help='Perform PPI mapping to identify mutation affected PPI interfaces.')
    req_ppi = subp_ppi.add_argument_group(title='Required arguments for PPI-mapping.', description='')
    
    req_ppi.add_argument('-g','--gene_id',type=str,help='The symbol of targeted protein-coding gene, for example: ERCC2',required=True)
    req_ppi.add_argument('-m','--mut_number', type=str, required=True, help='Comma-separated list of amino acids in specific substructure')
    req_ppi.add_argument('-f','--pdb_path', type=str, required=True, help='Path to the folder containing all the protein complex PDB files involving the target protein.')
    req_ppi.add_argument('-b','--chains', type=str, required=True, help='Comma-separated list of corresponding chain IDs of the target protein(e.g., A,B,A).')
        
    opt_ppi = subp_ppi.add_argument_group(title='Optional arguments for PPI mapping.', description='')
    opt_ppi.add_argument('-t','--distance_threshold', type=float, default=5.0, help='Distance threshold to determine whether two residues interact between among different chains(default: 5.0 Å).')
    opt_ppi.add_argument('-o','--outputfolder', type=str, default='PPI_Mapping_results', help='Output folder for results.')
    
    args = parser.parse_args()
    
    if args.subcmd == "cluster":
        df_all = pd.read_csv(args.inputfile)
        df_all.index = df_all['Function']
        sample_list = args.samples.split(',')
        control = args.control
        gene = args.gene_id; n_cluster = args.n_clusters
        method = args.method; metric = args.metric
        outputfolder = args.outputfolder
        th = args.significance_threshold
        pdb = args.pdb
        df_sig = df_all[abs(df_all['score'])>th]
        df_gene = df_sig[df_sig['Gene']==gene]
        
        df_clust = clustering(df_gene, gene, sample_list, control, 
                              n_cluster, method, metric, outputfolder)
        
        dicts = download_and_load_jsons()
        #df_clust['AA'] = df_clust['AA'].astype['Int64']
        df_clust = df_clust.dropna(subset=['AA'])
        df_anno = annotation(dicts,df_clust,gene,outputfolder,pdb)
        visualization_1d(dicts,df_all,df_clust,gene,outputfolder)
        
        if pdb != None:
            cluster_ls = []
            df_clust = df_clust[df_clust['Function']=='missense']
            for c in sorted(set(df_clust['Cluster'])):
                df_c = df_clust[df_clust['Cluster']==c]
                aa_ls = [aa for aa in df_c['AA'] if aa!='NA']
                cluster_ls.append(aa_ls)
            visualization_3d(cluster_ls,gene,outputfolder,pdb)
        logging.info("Clustering analysis completed. Results saved in %s", args.outputfolder)
     
    elif args.subcmd == "3d-rra":
        th = args.significance_threshold
        gene = args.gene_id; pdb = args.pdb; th1 = args.distance_threshold1; th2 = args.distance_threshold2
        df_all = pd.read_csv(args.inputfile)
        outputfolder = args.outputfolder
        df_all['Rank'] = abs(df_all['score']).rank(ascending=False)/df_all.shape[0]
        
        df_sig = df_all[abs(df_all['score'])>th]
        df_gene = df_sig[df_sig['Gene']==gene]
        df_gene = df_gene.dropna(subset=['AA'])
        df_clust = pd.read_csv(args.clustertable)
        df_clust = df_clust.dropna(subset=['AA'])
        df_substr = RRA_3D(gene, df_all, df_gene, df_clust, pdb, th1, th2, 
                           outputfolder,args.num_permutations,args.top_cutoff)
        visualize_substructure(df_substr,pdb,gene,outputfolder)
        logging.info("3D-RRA analysis completed. Results saved in %s", args.outputfolder)
    
    elif args.subcmd == "ppi-mapping":
        mut_numbers = args.mut_number.split(',')
        print(mut_numbers)
        gene = args.gene_id
        outputfolder = args.outputfolder
        #[int(x) for x in args.mutations.split(',')]
        chains = args.chains.split(',')
        pdb_path = args.pdb_path
        pdb_files = [os.path.join(pdb_path,s) for s in os.listdir(pdb_path) if '.pdb' in s]
        print(chains,pdb_files)
        df_ppi = build_ppi_interface_table(gene, mut_numbers,  pdb_files, chains, args.distance_threshold)
        visualize_interfaces(df_ppi,outputfolder)
        if outputfolder is not None:
            csv_path = os.path.join(os.getcwd(), outputfolder, gene+'_'+args.mut_number+"_ppi_interface_mutations.csv")
            df_ppi.to_csv(csv_path, index=False)
            print(f"Interface table saved to: {csv_path}")
        
        logging.info("PPI mapping analysis completed. Results saved in %s", outputfolder)
        
    else:
        parser.print_help()

if __name__ == "__main__":
    main()

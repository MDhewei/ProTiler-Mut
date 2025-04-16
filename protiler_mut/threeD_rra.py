#!/usr/bin/env python3
"""
3D RRA Method for Detecting Important 3D Substructures in Protein Structures

This script combines multiple functions:
  - get_uniprot_id_from_gene_symbol: Retrieve UniProt IDs from gene symbols.
  - get_seed_matrix: Compute a distance matrix between significant residues.
  - cluster_merge: Iteratively merge residue clusters based on a distance threshold.
  - structure_motif: Expand a residue cluster into a structural motif.
  - get_aa_cluster_3d: Identify 3D clusters of amino acids from gene data.
  - beta_calc: Compute the aggregated beta statistic.
  - efficient_rra_permutation: Perform a fast, vectorized permutation test for RRA.
  
A main() function demonstrates how these functions might be used in a 3D RRA workflow.
"""
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

import os
import sys
import random
import bisect
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
from scipy.stats import beta
from statsmodels.stats.multitest import fdrcorrection
from Bio.PDB.PDBParser import PDBParser
import mygene
from pymol import cmd
import pandas as pd


def get_uniprot_id_from_gene_symbol(gene_symbol: str) -> Optional[str]:
    """
    Retrieve the Swiss-Prot UniProt ID given a gene symbol using MyGene.info.

    Parameters:
        gene_symbol (str): Gene symbol to query.

    Returns:
        Optional[str]: The UniProt ID if found, otherwise None.
    """
    mg = mygene.MyGeneInfo()
    result = mg.query(gene_symbol, fields='uniprot.Swiss-Prot', species='9606')
    if 'hits' in result and len(result['hits']) > 0:
        for hit in result['hits']:
            if 'uniprot' in hit:
                uniprot_id = hit['uniprot'].get('Swiss-Prot')
                if uniprot_id:
                    return uniprot_id
    return None


def get_seed_matrix(df_cluster: pd.DataFrame, pdb_file: str) -> pd.DataFrame:
    """
    Build a seed matrix of pairwise distances between significant residues based on a z-score threshold.

    Parameters:
        df_gene (pd.DataFrame): DataFrame containing gene information with columns 'AA' and 'Zscore'.
        pdb_file (str): Path to the PDB file.
        zscore_threshold (float): Threshold for selecting significant residues.

    Returns:
        pd.DataFrame: A symmetric DataFrame of minimum pairwise distances (in angstroms) between residues.
    """
    # Filter significant residues
    #df_sig = df_gene[df_gene['Zscore'] >= zscore_threshold]
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = parser.get_structure('', pdb_file)
    chain = structure[0]['A']

    # Get chain residues and their positions
    residues = list(chain.get_residues())
    chain_res_indices = [res.id[1] for res in residues]

    # Collect unique significant residue indices present in the chain
    significant_indices = []
    for aa in df_cluster['AA']:
        print(aa)
        if aa in chain_res_indices and aa not in significant_indices:
            significant_indices.append(int(aa))

    # Calculate pairwise minimum distances between atoms of significant residues
    distance_matrix = []
    for aa1 in significant_indices:
        distances_row = []
        res1 = next(res for res in residues if res.id[1] == aa1)
        for aa2 in significant_indices:
            res2 = next(res for res in residues if res.id[1] == aa2)
            atom_distances = [atom1 - atom2 for atom1 in res1.get_atoms() for atom2 in res2.get_atoms()]
            distances_row.append(min(atom_distances))
        distance_matrix.append(distances_row)

    # Create DataFrame with labels (e.g., 'res23')
    labels = [f"res{idx}" for idx in significant_indices]
    seed_df = pd.DataFrame(distance_matrix, columns=labels, index=labels)
    return seed_df


def cluster_merge(seed_df: pd.DataFrame, merge_threshold: float) -> List[List[str]]:
    """
    Merge clusters iteratively based on the distance matrix.
    
    Parameters:
        distance_data (pd.DataFrame): Symmetric DataFrame containing pairwise distances.
        col_names (List[str]): List of residue labels (e.g., ['res23', 'res45']).
        merge_threshold (float): Distance threshold for merging clusters.
        
    Returns:
        List[List[str]]: A list of clusters (each a list of residue labels).
    """
    link_set = []
    merged_elements = []

    while True:
        min_val = sys.maxsize
        merge_pair = (None, None)

        # Find the pair with the smallest distance (excluding self comparisons)
        for t1 in seed_df.columns:
            for t2 in seed_df.columns:
                if t1 != t2 and seed_df[t1][t2] < min_val:
                    min_val = seed_df[t1][t2]
                    merge_pair = (t1, t2)
        col, row = merge_pair

        # Update distances: average distances for the merged residue and set the partner's distances high
        for c in seed_df.columns:
            if c not in (col, row):
                avg_distance = np.mean([seed_df[col][c], seed_df[row][c]])
                seed_df[col][c] = avg_distance
                seed_df[c][col] = avg_distance
            
            seed_df[row][c] = 100
            seed_df[c][row] = 100

        if min_val < merge_threshold:
            merged = False
            merge_indices = []
            # Check if either residue is already in an existing cluster
            for i, cluster in enumerate(link_set):
                if set(cluster) & {row, col}:
                    link_set[i] = list(set(cluster) | {row, col})
                    merged_elements.extend([row, col])
                    merge_indices.append(i)
                    merged = True
            # Merge multiple clusters if necessary
            if len(merge_indices) > 1:
                combined = []
                for idx in sorted(merge_indices, reverse=True):
                    combined.extend(link_set.pop(idx))
                link_set.append(list(set(combined)))
            elif not merged:
                link_set.append([row, col])
                merged_elements.extend([row, col])
        else:
            break

    # Residues not merged become singleton clusters
    singleton_clusters = [[c] for c in seed_df.columns if c not in merged_elements]
    cluster_list = singleton_clusters + link_set
    return cluster_list


def structure_motif(cluster: List[str], pdb_file: str, distance_cutoff: float) -> List[int]:
    """
    Identify a structural motif by expanding a given cluster of residues. For each residue in the cluster,
    find all residues that are within the specified distance cutoff.
    
    Parameters:
        cluster (List[str]): List of residue labels (e.g., ['res23', 'res45']).
        pdb_file (str): Path to the PDB file.
        distance_cutoff (float): Distance cutoff (in angstroms) to include a residue in the motif.
    
    Returns:
        List[int]: List of residue numbers (as integers) forming the structural motif.
    """
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = parser.get_structure('', pdb_file)
    chain = structure[0]['A']
    residues = list(chain.get_residues())
    chain_res_indices = [res.id[1] for res in residues]

    motif = []
    for res_label in cluster:
        try:
            residue_number = int(res_label[3:])  # Assuming label like "res23"
            print(residue_number)
        except ValueError as e:
            print(e)
            continue
        if residue_number in chain_res_indices:
            res1 = next(res for res in residues if res.id[1] == residue_number)
            for res2 in residues:
                atom_distances = [atom1 - atom2 for atom1 in res1.get_atoms() for atom2 in res2.get_atoms()]
                if min(atom_distances) <= distance_cutoff and res2.id[1] not in motif:
                    motif.append(res2.id[1])
    return motif


def get_aa_cluster_3d(motif, df_gene: pd.DataFrame) -> Tuple[List[int], List]:
    """
    Identify amino acid clusters in 3D space by comparing candidate residues to significant residues.
    
    Parameters:
        aa_sig_list (List): List of significant amino acid positions.
        df_gene (pd.DataFrame): DataFrame containing gene data with columns 'AA' and 'Rank'.
        
    Returns:
        Tuple[List[int], List]: A tuple with a list of amino acid positions in the cluster and their corresponding ranks.
    """

    aa_cluster = []
    rank_list = []
    for _, row in df_gene.iterrows():
        try:
            aa_candidate = row['AA']
            rank = row['Rank']
            #print(aa_candidate)
            if pd.isna(aa_candidate):
                continue
            if int(aa_candidate) in motif:
                #print('Yes')
                aa_cluster.append(int(aa_candidate))
                rank_list.append(rank)
        except Exception:
            continue   
    return aa_cluster, rank_list


def beta_calc(values: List[float]) -> float:
    """
    Calculate the minimum beta CDF value over a list of values.
    
    Parameters:
        values (List[float]): List of rank values.
    
    Returns:
        float: The minimum beta CDF value; returns 1.0 if the list is empty.
    """
    n = len(values)
    if n > 0:
        beta_values = [beta.cdf(values[k], k + 1, n - k) for k in range(n)]
        return min(beta_values)
    else:
        return 1.0


def efficient_rra_permutation(rank_list, df_be, n, num_permutations=100000, significance_cutoff=0.25) -> Tuple[float, float, float]:
    """
    Perform an efficient, vectorized Robust Rank Aggregation (RRA) permutation test.
    
    Parameters:
        rank_list (array-like): Observed ranks.
        df_be (pd.DataFrame): DataFrame containing a 'Rank' column with the overall rank distribution.
        n (int): Number of ranks to sample for aggregation.
        num_permutations (int): Number of permutations.
        significance_cutoff (float): Only ranks below this cutoff are considered.
    
    Returns:
        Tuple[float, float, float]: (Ro_obs, p_value, fdr)
            - Ro_obs: Aggregated beta statistic for the observed sample.
            - p_value: Fraction of permutation statistics as extreme as Ro_obs.
            - fdr: Estimated false discovery rate.
    """
    # Overall rank distribution as numpy array
    rank_all = np.array(df_be['Rank'])
    
    rank_list = np.array(rank_list)
    if len(rank_list) < n:
        additional = np.random.choice(rank_all, size=n - len(rank_list), replace=False)
        observed_sample = np.concatenate([rank_list, additional])
    else:
        observed_sample = np.random.choice(rank_list, size=n, replace=False)
    observed_filtered = np.sort(observed_sample[observed_sample < significance_cutoff])
    
    def compute_beta_statistic(ranks):
        n_val = len(ranks)
        if n_val == 0:
            return 1.0
        a_params = np.arange(1, n_val + 1)
        b_params = n_val - np.arange(0, n_val)
        beta_values = beta.cdf(ranks, a_params, b_params)
        return np.min(beta_values)
    
    Ro_obs = compute_beta_statistic(observed_filtered)
    
    # Vectorized permutation sampling: shape (num_permutations, n)
    #print(len(rank_all),n)
    #perm_samples = np.random.choice(rank_all, size=(num_permutations, n), replace=False)
    perm_samples = np.array([np.random.choice(rank_all, size=n, replace=False)
                             for _ in range(num_permutations)])

    perm_samples.sort(axis=1)
    perm_samples_masked = np.where(perm_samples < significance_cutoff, perm_samples, 1.0)
    
    a_params = np.arange(1, n + 1)
    b_params = n - np.arange(0, n)
    beta_vals = beta.cdf(perm_samples_masked, a_params, b_params)
    perm_beta_stats = np.min(beta_vals, axis=1)
    
    #p_value = np.mean(perm_beta_stats <= Ro_obs)
    sorted_perm = np.sort(perm_beta_stats)
    p_value = (np.searchsorted(sorted_perm, Ro_obs) + 1) / (len(sorted_perm) + 1)
    perm_sig = sorted_perm[0:int(len(sorted_perm)*0.05)]
    fdr = (np.searchsorted(perm_sig, Ro_obs) + 1) / (len(perm_sig) + 1)
    
    return Ro_obs, p_value, fdr


def RRA_3D(gene,df_all,df_gene,df_clust,pdbfile,th1,th2,output_folder,num_permutations=10000,significance_cutoff=0.25):
    
    if output_folder is not None and not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    cluster_assign = []; substr_ls = []; sigaa_ls = []; 
    Ro_ls = []; p_ls = []; fdr_ls = []
    
    for c in set(df_clust['Cluster']):
        df_c = df_clust[df_clust['Cluster']==c]
        seed_df = get_seed_matrix(df_c, pdbfile)
        
        if seed_df.shape[0]>1:
            cluster_ls = cluster_merge(seed_df,th1)

            for clust in cluster_ls:
                substr = structure_motif(clust,pdbfile,th2)
                print(clust,substr)
                aa_clust_ls,rank_ls = get_aa_cluster_3d(substr, df_gene)
                print(aa_clust_ls,rank_ls)
                if len(aa_clust_ls) >= 3:
                    Ro_obs, p_value, fdr = efficient_rra_permutation(rank_ls,df_all,len(rank_ls),
                                                                    num_permutations,significance_cutoff)
                    print(Ro_obs, p_value, fdr)
                    cluster_assign.append('cluster_'+str(int(c)+1))
                    sigaa_ls.append(aa_clust_ls)
                    substr_ls.append(substr)
                    Ro_ls.append(Ro_obs)
                    p_ls.append(p_value)
                    fdr_ls.append(fdr)
                    
    df_substr = pd.DataFrame({'Gene':[gene]*len(p_ls),'Cluster':cluster_assign,'Seed_aa':sigaa_ls,'Substructure':substr_ls,
                              'Ro_value':Ro_ls,'Pvalue':p_ls,'FDR':fdr_ls})
    print(df_substr)
    if output_folder is not None:
        table_path = os.path.join(os.getcwd(),output_folder, f"{gene}_Substructures_list.csv")
        df_substr.to_csv(table_path, index=False)
        print(f"Substructure table saved to: {table_path}")
                          
    return df_substr
                             
                             
from pymol import cmd, finish_launching
import os
import pandas as pd

def visualize_substructure(df_substr: pd.DataFrame,
                          pdb_file: str,
                          gene,outputfolder,
                          object_name: str = 'structure'
                          ):
    """
    For each cluster in df_substr, load the PDB, highlight its substructures,
    and save the entire PyMOL session as a .pse file.

    Parameters:
        df_substr (pd.DataFrame): must have columns
            - 'Cluster' (e.g. 'cluster_1')
            - 'Substructure' (list of residue numbers)
        pdb_file (str): path to the PDB file
        object_name (str): name to give the loaded object in PyMOL
        output_dir (str): directory where .pse files will be saved
    """
    #finish_launching()
    # Ensure output directory exists
    #os.makedirs(output_dir, exist_ok=True)

    # Load structure once
    cmd.reinitialize()
    cmd.load(pdb_file, object_name)
    cmd.show('cartoon', object_name)
    cmd.color('gray70', object_name)

    # A simple palette of PyMOL colors
    palette = ['red', 'blue', 'green', 'yellow', 'magenta',
               'cyan', 'orange', 'purple', 'salmon', 'lime']

    clusters = sorted(df_substr['Cluster'].unique())
    for ci, cluster in enumerate(clusters):
        # Hide any previous sticks so only this cluster shows
        cmd.hide('sticks', object_name)

        # Subset to just this cluster’s rows
        df_clust = df_substr[df_substr['Cluster'] == cluster].reset_index(drop=True)
        if df_clust.empty:
            continue

        # Highlight each sub‐motif in the cluster
        for pi, row in df_clust.iterrows():
            residues = row['Substructure']
            if not residues:
                continue
            resi_str = '+'.join(str(r) for r in residues)
            sel_name = gene+f"_{cluster}_part{pi+1}"
            cmd.select(sel_name, f"{object_name} and resi {resi_str}")
            color = palette[(ci + pi) % len(palette)]
            cmd.color(color, sel_name)
            cmd.show('sticks', sel_name)
            cmd.set('stick_radius', 0.2, sel_name)

        # Frame the view on the current highlights
        cmd.zoom(object_name)

        # Save the session
        session_path = os.path.join(outputfolder, gene+f"_{cluster}.pse")
        cmd.save(session_path)
        print(f"Saved session for {cluster} → {session_path}")

__all__ = [
    "RRA_3D","visualize_substructure","get_uniprot_id_from_gene_symbol"
]
# Example usage:
if __name__ == '__main__':
    df_substr = threeD_rra(df_all, df_gene, df_clust, pdbfile, th1, th2)
    save_cluster_sessions(df_substr, pdbfile, output_dir='cluster_sessions')
            
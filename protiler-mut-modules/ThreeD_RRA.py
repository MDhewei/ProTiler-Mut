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


def get_seed_matrix(df_gene: pd.DataFrame, pdb_file: str, zscore_threshold: float) -> pd.DataFrame:
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
    df_sig = df_gene[df_gene['Zscore'] >= zscore_threshold]
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = parser.get_structure('', pdb_file)
    chain = structure[0]['A']

    # Get chain residues and their positions
    residues = list(chain.get_residues())
    chain_res_indices = [res.id[1] for res in residues]

    # Collect unique significant residue indices present in the chain
    significant_indices = []
    for aa in df_sig['AA']:
        if aa in chain_res_indices and aa not in significant_indices:
            significant_indices.append(aa)

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


def cluster_merge(distance_data: pd.DataFrame, col_names: List[str], merge_threshold: float) -> List[List[str]]:
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
        for t1 in col_names:
            for t2 in col_names:
                if t1 != t2 and distance_data.at[t1, t2] < min_val:
                    min_val = distance_data.at[t1, t2]
                    merge_pair = (t1, t2)
        col, row = merge_pair

        # Update distances: average distances for the merged residue and set the partner's distances high
        for c in col_names:
            if c not in (col, row):
                avg_distance = np.mean([distance_data.at[col, c], distance_data.at[row, c]])
                distance_data.at[col, c] = avg_distance
                distance_data.at[c, col] = avg_distance
            distance_data.at[row, c] = 100
            distance_data.at[c, row] = 100

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
    singleton_clusters = [[c] for c in col_names if c not in merged_elements]
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
        except ValueError:
            continue
        if residue_number in chain_res_indices:
            res1 = next(res for res in residues if res.id[1] == residue_number)
            for res2 in residues:
                atom_distances = [atom1 - atom2 for atom1 in res1.get_atoms() for atom2 in res2.get_atoms()]
                if min(atom_distances) <= distance_cutoff:
                    motif.append(res2.id[1])
    return motif


def get_aa_cluster_3d(aa_sig_list: List, df_gene: pd.DataFrame, pdb_file: str) -> Tuple[List[int], List]:
    """
    Identify amino acid clusters in 3D space by comparing candidate residues to significant residues.
    
    Parameters:
        aa_sig_list (List): List of significant amino acid positions.
        df_gene (pd.DataFrame): DataFrame containing gene data with columns 'AA' and 'Rank'.
        pdb_file (str): Path to the PDB file.
        
    Returns:
        Tuple[List[int], List]: A tuple with a list of amino acid positions in the cluster and their corresponding ranks.
    """
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = parser.get_structure('', pdb_file)
    chain = structure[0]['A']

    aa_cluster = []
    rank_list = []
    for _, row in df_gene.iterrows():
        try:
            aa_candidate = row['AA']
            rank = row['Rank']
            if pd.isna(aa_candidate):
                continue
            for aa_sig in aa_sig_list:
                if pd.isna(aa_sig):
                    continue
                atom_distances = [
                    atom1 - atom2
                    for atom1 in chain[int(aa_candidate)]
                    for atom2 in chain[int(aa_sig)]
                ]
                if min(atom_distances) <= 6:
                    aa_cluster.append(int(aa_candidate))
                    rank_list.append(rank)
                    break
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


def efficient_rra_permutation(rank_list, df_be, n, num_permutations=1000000, significance_cutoff=0.25) -> Tuple[float, float, float]:
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
    perm_samples = np.random.choice(rank_all, size=(num_permutations, n), replace=False)
    perm_samples.sort(axis=1)
    perm_samples_masked = np.where(perm_samples < significance_cutoff, perm_samples, 1.0)
    
    a_params = np.arange(1, n + 1)
    b_params = n - np.arange(0, n)
    beta_vals = beta.cdf(perm_samples_masked, a_params, b_params)
    perm_beta_stats = np.min(beta_vals, axis=1)
    
    p_value = np.mean(perm_beta_stats <= Ro_obs)
    sorted_perm = np.sort(perm_beta_stats)
    fdr = (np.searchsorted(sorted_perm, Ro_obs) + 1) / (len(sorted_perm) + 1)
    
    return Ro_obs, p_value, fdr


def main():
    """
    Demonstration of the 3D RRA workflow.
    
    This example assumes you have:
      - A tiling mutagenesis screen data file (CSV) with columns 'AA', 'Zscore', and 'Rank'
      - A PDB file for the protein structure.
    Adjust file paths, thresholds, and parameters as needed.
    """
    # Example inputs (replace these with your actual data paths and parameters)
    gene_symbol = "TP53"
    pdb_file = "example.pdb"        # Path to your PDB file
    gene_data_file = "gene_data.csv"  # CSV file with gene data
    
    zscore_threshold = 2.0
    merge_threshold = 6.0
    distance_cutoff = 6.0
    n_permutation_sample = 10       # Number of ranks to sample for RRA (example)
    num_permutations = 10000        # Use higher (e.g., 1e6) for final analysis
    
    # Load gene expression data (assumes CSV with columns: 'AA', 'Zscore', 'Rank')
    df_gene = pd.read_csv(gene_data_file)
    
    # (Optional) Retrieve UniProt ID for the gene
    uniprot_id = get_uniprot_id_from_gene_symbol(gene_symbol)
    print("UniProt ID:", uniprot_id)
    
    # Compute the seed matrix for significant residues
    seed_matrix = get_seed_matrix(df_gene, pdb_file, zscore_threshold)
    print("Seed matrix:")
    print(seed_matrix)
    
    # Merge clusters based on the seed matrix
    clusters = cluster_merge(seed_matrix.copy(), list(seed_matrix.columns), merge_threshold)
    print("Clusters:")
    print(clusters)
    
    # Identify structural motifs from each cluster
    motifs = []
    for cluster in clusters:
        motif = structure_motif(cluster, pdb_file, distance_cutoff)
        motifs.append(motif)
    print("Structural motifs:")
    print(motifs)
    
    # Identify 3D amino acid clusters from the gene data
    aa_sig_list = list(df_gene[df_gene['Zscore'] >= zscore_threshold]['AA'])
    aa_cluster, rank_list = get_aa_cluster_3d(aa_sig_list, df_gene, pdb_file)
    print("Amino acid cluster positions:", aa_cluster)
    print("Corresponding ranks:", rank_list)
    
    # Perform the permutation-based Robust Rank Aggregation (RRA)
    Ro_obs, p_value, fdr = efficient_rra_permutation(
        rank_list, df_gene, n_permutation_sample, num_permutations=num_permutations
    )
    print("RRA Results:")
    print("Aggregated beta statistic (Ro_obs):", Ro_obs)
    print("P-value:", p_value)
    print("FDR:", fdr)


if __name__ == "__main__":
    main()

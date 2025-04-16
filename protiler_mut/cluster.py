#!/usr/bin/env python3
from __future__ import division

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

import os
import json
import math
import re
import random
import requests
from typing import Dict

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize
from matplotlib import cm
from matplotlib.patches import Rectangle
from textwrap import wrap

from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.neighbors import KernelDensity

from scipy.stats import pearsonr
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.ResidueDepth import get_surface
from Bio.PDB.ResidueDepth import residue_depth
from Bio.PDB.DSSP import DSSP


import umap
import pymol

#=============================================================================

def clustering(
    df_gene: pd.DataFrame,
    gene: str,
    sample_list: list,
    sample_control: list,
    n_clusters: int,
    clustering_method: str,
    clustering_metric: str,
    output_folder: str = None,
) -> pd.DataFrame:
    """
    Perform hierarchical clustering on screen data, plot a heatmap with dendrograms, and display UMAP, PCA,
    and TSNE visualizations. Save the results if an output folder is provided.

    Parameters:
        df_gene (pd.DataFrame): DataFrame containing gene data with at least columns for expression.
        gene (str): Gene name.
        sample_list (list): List of sample column names to include in the clustering.
        sample_control (list): List of control sample column names.
        n_clusters (int): Number of clusters to generate.
        clustering_method (str): Linkage method (e.g., 'average') for clustering.
        clustering_metric (str): Distance metric (e.g., 'euclidean') for clustering.
        output_folder (str, optional): Directory path to save outputs (table and figures).

    Returns:
        pd.DataFrame: Updated DataFrame with added 'Cluster' and 'Correlation' columns.
        
    Note:
        This function expects that a global dictionary 'corr_dic' exists to provide correlation values
        for each sgID.
    """
    # Create output folder if provided
    if output_folder is not None and not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    #df_gene.index = df_gene['Function']+'_'+df_gene['Mutation']
    
    # Subtract control samples from the screen data
    df_x = df_gene.loc[:, sample_list]
    if sample_control!='':
        df_x = df_x.subtract(df_x[sample_control], axis=0)

    # Perform hierarchical clustering
    cluster_model = AgglomerativeClustering(
        n_clusters=n_clusters, linkage=clustering_method, affinity=clustering_metric
    )
    cluster_model.fit(df_x)
    df_gene["Cluster"] = cluster_model.labels_

    # Set up the heatmap figure
    plt.figure(figsize=(8, 6), dpi=300)
    plt.rcParams["font.size"] = 28
    plt.rcParams["font.family"] = "Helvetica"

    # Generate colors for clusters
    cluster_labels = cluster_model.labels_
    color_options = [ "blue", "red","green", "orange", "purple",
                     "yellow", "cyan", "pink", "teal", "brown"]
    top_colors = color_options[: len(set(cluster_labels))]
    top_colors_dict = dict(zip(set(cluster_labels), top_colors))
    top_colors_list = [top_colors_dict[label] for label in cluster_labels]

    # Create the clustermap
    g = sns.clustermap(
        df_x.T,
        metric=clustering_metric,
        method=clustering_method,
        cmap="coolwarm",
        dendrogram_ratio=(0.1, 0.1),
        vmin=-1,
        vmax=1,
        row_cluster=False,
        col_cluster=True,
        cbar_pos=(0.02, 0.35, 0.03, 0.2),
        figsize=(20, 10),
        xticklabels=False,
        col_colors=top_colors_list,
    )

    cbar = g.ax_cbar
    cbar.set_title("LFC", pad=10)
    
    # Annotate mutation categories based on sgID naming
    mutation_categories = []
    for idx in g.dendrogram_col.reordered_ind:
        sgid = str(df_x.index[idx])
        if "nonsense" in sgid or "splice" in sgid:
            mutation_categories.append("Nonsense &Splicing")
        else:
            mutation_categories.append("Missense")
    
    bottom_colors_dict = {"Nonsense &Splicing": "purple", "Missense": "grey"}
    bottom_colors_list = [bottom_colors_dict[cat] for cat in mutation_categories]
    
    title_font = {"size": 40, "weight": "bold"}
    g.ax_heatmap.set_title(gene, fontdict=title_font, pad=100)
    
    # Add a bottom annotation bar for mutation categories
    ax = g.ax_heatmap
    for idx, color in enumerate(bottom_colors_list):
        ax.add_patch(
            Rectangle(
                (idx, len(sample_list)+0.2),  # x, y
                1,
                0.2,
                facecolor=color,
                transform=ax.transData,
                clip_on=False,
            )
        )
    
    # Create legends for clusters and mutation categories
    for label, color in top_colors_dict.items():
        g.ax_col_dendrogram.bar(0, 0, color=color, label=f"Cluster {label+1}", linewidth=0)
    for label, color in bottom_colors_dict.items():
        ax.bar(0, 0, color=color, label="\n".join(wrap(label, width=9)), linewidth=0)
    
    g.ax_col_dendrogram.legend(
        loc="center", title="Clusters", bbox_to_anchor=(-0.1, 0), fontsize=25
    )
    ax.legend(
        loc="center", title="Muation Category", bbox_to_anchor=(-0.12, 0.02), fontsize=25
    )
    
    # Add correlation annotation bars using the global 'corr_dic'
    norm = Normalize(vmin=-1, vmax=1)
    cmap_obj = cm.get_cmap("bwr")
    
    df_lof = df_gene[(df_gene['Function']=='nonsense')|(df_gene['Function'].str.contains('splice')==True)]
    df_lof_data = df_lof.loc[:,sample_list]
    df_lof_data.reset_index(drop=True, inplace=True)
    df_gene_data = df_gene.loc[:,sample_list]
    
    corr_vals = []
    
    if df_lof.shape[0] > 0:
        for ix1 in g.dendrogram_col.reordered_ind:
            #print(df_gene_data.iloc[ix1,:])
            ls1 = list(df_gene_data.iloc[ix1,:])
            corr_ls = []
            for ix2 in df_lof_data.index:
                ls2 = list(df_lof_data.iloc[ix2,:])
                #print(ls1,ls2)
                r = pearsonr(ls1,ls2)[0]
                corr_ls.append(r)
            corr_vals.append(np.mean(corr_ls))
                
    else:
        corr_vals = [0] * df_gene.shape[0]
            
    for i, cval in enumerate(corr_vals):
        rect_color = cmap_obj(norm(cval))
        ax.add_patch(
            Rectangle(
                (i, len(sample_list)+0.7),  # x, y
                1,
                0.5,
                facecolor=rect_color,
                transform=ax.transData,
                clip_on=False,
            )
        )
        anno = "NA" if cval == 0 else ""
        ax.text(
            i + 0.5,
            len(sample_list)+0.95,
            anno,
            ha="center",
            va="center",
            color="black",
            fontsize=20,
            fontweight="bold",
            transform=ax.transData,
        )
    
    sm = cm.ScalarMappable(norm=norm, cmap=cmap_obj)
    sm.set_array([])
    cbar_ax = g.fig.add_axes([-0.02, -0.12, 0.1, 0.05])
    cbar2 = plt.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar2.set_label("Correlation", fontsize=25)
    cbar2.ax.tick_params(labelsize=26)
    
    df_gene["Correlation"] = corr_vals
    
    if output_folder is not None:
        clustermap_path = os.path.join(os.getcwd(),output_folder, f"{gene}_ClusterMap.png")
        g.fig.savefig(clustermap_path, dpi=300, facecolor="white", bbox_inches="tight")
        print(f"Clustermap figure saved to: {clustermap_path}")
        table_path = os.path.join(os.getcwd(),output_folder, f"{gene}_Cluster_Table.csv")
        df_gene.to_csv(table_path, index=False)
        print(f"Clustering table saved to: {table_path}")
        
    plt.show()
    
    # Dimensionality reduction
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df_x)
    
    # UMAP
    reducer = umap.UMAP(random_state=47, n_neighbors=6, min_dist=0.3, metric="euclidean")
    embedding = reducer.fit_transform(scaled_data)
    
    # PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_data)
    
    # TSNE (highlighting potential outliers with 'nonsense' or 'splice' in index)
    lof_idx = [
        ix 
        for ix in range(df_x.shape[0])
        if "splice" in str(df_x.index[ix]) or "nonsense" in str(df_x.index[ix])
    ]
    tsne = TSNE(n_components=2, random_state=47, perplexity=5)
    tsne_result = tsne.fit_transform(scaled_data)
    
    # Plot UMAP, PCA, and TSNE results
    plt.figure(figsize=(12, 4), dpi=300)
    #sns.set_style("ticks")
    plt.rcParams["font.size"] = 16
    plt.rcParams["font.family"] = "Helvetica"
    
    plt.subplot(1, 3, 1)
    plt.scatter(embedding[:, 0], embedding[:, 1], c=top_colors_list, s=80)
    
    plt.scatter(
        embedding[lof_idx, 0],
        embedding[lof_idx, 1],
        c="purple",
        s=85,
        label="Nonsense&Splicing",
    )
    plt.xlabel("UMAP Dimension 1")
    plt.ylabel("UMAP Dimension 2")
    
    plt.subplot(1, 3, 2)
    plt.scatter(pca_result[:, 0], pca_result[:, 1], c=top_colors_list, s=80)
    plt.scatter(
        pca_result[lof_idx, 0],
        pca_result[lof_idx, 1],
        c="purple",
        s=85,
        label="Nonsense&Splicing",
    )
    plt.xlabel("PCA Dimension 1")
    plt.ylabel("PCA Dimension 2")
    
    plt.subplot(1, 3, 3)
    plt.scatter(tsne_result[:, 0], tsne_result[:, 1], c=top_colors_list, s=80)
    plt.scatter(
        tsne_result[lof_idx, 0],
        tsne_result[lof_idx, 1],
        c="purple",
        s=85,
        label="Nonsense&Splicing",
    )
    plt.xlabel("TSNE Dimension 1")
    plt.ylabel("TSNE Dimension 2")
    
    plt.tight_layout()
    
    if output_folder is not None:
        dr_path = os.path.join(os.getcwd(),output_folder, f"{gene}_DimReduction.png")
        plt.savefig(dr_path, dpi=300, facecolor="white", bbox_inches="tight")
        print(f"Dimensionality reduction plot saved to: {dr_path}")
    
    plt.show()
    return df_gene
    
def download_and_load_jsons() -> Dict[str, dict]:
    """
    Download and load all of our annotation JSONs from Figshare.

    Returns:
        A dict with keys:
            - 'exons_dic'
            - 'domain_dic'
            - 'alpha_dic'
            - 'po_dic'
            - 'clinvar_dic'
        each mapping to the loaded JSON content.
    """
    # 1) Define your Figshare download URLs here:
    url_mapping = {
        "exons_dic": "https://figshare.com/ndownloader/files/53596592",
        "domain_dic": "https://figshare.com/ndownloader/files/53596595",
        "alpha_dic": "https://figshare.com/ndownloader/files/53596607",
        "po_dic": "https://figshare.com/ndownloader/files/53596601",
        "clinvar_dic": "https://figshare.com/ndownloader/files/53596589"
    }

    # 2) (Optional) Where to cache the files locally
    cache_dir = os.path.expanduser("~/.cache/protein_annotations")
    os.makedirs(cache_dir, exist_ok=True)

    results: Dict[str, dict] = {}

    for key, url in url_mapping.items():
        local_path = os.path.join(cache_dir, f"{key}.json")

        # Download if not already cached
        if not os.path.exists(local_path):
            resp = requests.get(url, stream=True)
            resp.raise_for_status()
            with open(local_path, "wb") as fh:
                for chunk in resp.iter_content(1 << 20):
                    fh.write(chunk)
            print(f"[downloaded] {key} → {local_path}")
        else:
            print(f"[cached]    {key} → {local_path}")

        # Load JSON into memory
        with open(local_path, "r") as fh:
            results[key] = json.load(fh)

    return results


def annotation(dicts,df_clust,gene,output_folder,pdb_path=None):
    exons_dic   = dicts["exons_dic"]
    domain_dic  = dicts["domain_dic"]
    alpha_dic   = dicts["alpha_dic"]
    po_dic      = dicts["po_dic"]
    clinvar_dic = dicts["clinvar_dic"]
    # Create output folder if provided
    if output_folder is not None and not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    df_clust = df_clust.dropna(subset=['AA'])
    ## Annotate alphafold score
    if gene in alpha_dic:
        alpha_ls = []
        for aa in df_clust['AA']:
            score = 'NA'
            if aa != 'NA' and str(int(aa)) in alpha_dic[gene]:
                score = round(np.mean(alpha_dic[gene][str(int(aa))]),3)
            alpha_ls.append(score)
        df_clust['Alphafold.score'] = alpha_ls
    else:
        df_clust['Alphafold.score'] = ['NA'] * df_clust.shape[0]
             
    ## Domain annotation
    if gene in domain_dic:
        dom_ls = []   
        for aa in df_clust['AA']:
            domain = 'NA'
            if aa != 'NA':
                for dom in domain_dic[gene]:
                    name = dom[0]; start = dom[1]; end = dom[2]
                    if start <= aa <= end:
                        domain = name
                        break
            dom_ls.append(domain)
        df_clust['Pfam.domain'] = dom_ls
        
    else:
        df_clust['Pfam.domain'] = ['NA'] * df_clust.shape[0]
        
    ## Phosphorylation annotation
    if gene in po_dic:
        po_list = [s[0] for s in po_dic[gene] if s[1]>0]
        po_ls = []
        for aa in df_clust['AA']:
            po = 0
            if aa != 'NA':
                if int(aa) in po_list:
                    po = 1
            po_ls.append(po)
        df_clust['Phospho.site'] = po_ls
    else:
        df_clust['Phospho.site'] = ['NA'] * df_clust.shape[0]
            
    ## Clivar annotation:
    if gene in clinvar_dic:
        cli_ls = []
        for aa in df_clust['AA']:
            cli = 'NA'
            if aa != 'NA':
                for mut in clinvar_dic[gene]:
                    mut_aa = int(mut[3:-3])
                    if abs(aa - mut_aa)<=5:
                        cli = mut
                        break
            cli_ls.append(cli)
        df_clust['ClinVar'] = cli_ls
            
    else:
        df_clust['ClinVar'] = ['NA'] * df_clust.shape[0]
    
    if pdb_path is None:
        print(f"No pdb path provided for gene {gene}. Skipping structure score calculation.")
        # Fill with 'NA' if you wish to maintain columns
        df_clust['ResidueDepth'] = ['NA'] * df_clust.shape[0]
        df_clust['ASA'] = ['NA'] * df_clust.shape[0]
    
    else:
        print('Attempt to parse the PDB structure')
        try:
            parser = PDBParser(PERMISSIVE=1, QUIET=True)
            structure = parser.get_structure('', pdb_path)
            model = structure[0]
        except Exception as e:
            print(f"Error parsing pdb file for gene {gene}: {e}")
            model = None
        
        if model is not None:
            
            try:
                surface = get_surface(model)
            except Exception as e:
                print(f"Error calculating surface for gene {gene}: {e}")
                surface = None
            
            try:
                chain = model['A']
            except Exception as e:
                print(f"Error accessing chain A for gene {gene}: {e}")
                chain = None
            
            try:
                dssp = DSSP(model, pdb_path, dssp='mkdssp')
            except Exception as e:
                print(f"Error calculating DSSP for gene {gene}: {e}")
                dssp = None

            rd_ls = []; asa_ls = []
            # Assuming df_clust['AA'] holds residue numbers (integers)
            for aa in df_clust['AA']:
                aa = int(aa)
                try:
                    if chain is not None and surface is not None and dssp is not None:
                        # Get the residue from chain using the residue number directly
                        residue = chain[aa]
                        rd = residue_depth(residue, surface)
                        # Construct the DSSP key; adjust based on your DSSP output format
                        key = ('A', (' ', aa, ' '))
                        asa = dssp[key][3]
                        rd_ls.append(rd)
                        asa_ls.append(asa)
                    else:
                        rd_ls.append('NA')
                        asa_ls.append('NA')
                except Exception as e:
                    rd_ls.append('NA')
                    asa_ls.append('NA')
            
            df_clust['ResidueDepth'] = rd_ls
            df_clust['ASA'] = asa_ls
        
        else:
            df_clust['ResidueDepth'] = ['NA'] * df_clust.shape[0]
            df_clust['ASA'] = ['NA'] * df_clust.shape[0]
    
    # Save the annotated table if an output folder is provided
    if output_folder is not None:
        table_path = os.path.join(os.getcwd(),output_folder, f"{gene}_Missense_Annotation_Table.csv")
        df_clust.to_csv(table_path, index=False)
        print(f"Annotation table saved to: {table_path}")
    
    return df_clust


def kde_sklearn(x: np.ndarray, x_grid: np.ndarray, bandwidth: float) -> np.ndarray:
    """
    Compute Kernel Density Estimation using Scikit-learn.

    Parameters:
        x (np.ndarray): Input data points of shape (n_samples, 1).
        x_grid (np.ndarray): Points to evaluate the density on.
        bandwidth (float): Bandwidth for the kernel.

    Returns:
        np.ndarray: Density estimates evaluated on x_grid.
    """
    kde_model = KernelDensity(kernel="gaussian", bandwidth=bandwidth)
    kde_model.fit(x)
    log_pdf = np.exp(kde_model.score_samples(x_grid))
    return log_pdf


def get_sift_kde(con_list: list, aa_list: list, bw: float) -> np.ndarray:
    """
    Generate a kernel density estimate (KDE) for SIFT scores.

    Parameters:
        con_list (list): List of conservation scores.
        aa_list (list): List of amino acid positions.
        bw (float): Bandwidth for KDE.

    Returns:
        np.ndarray: KDE values evaluated on the amino acid grid.
    """
    binned = list(pd.cut(con_list, 100, labels=False))
    x_grid = np.array(aa_list).reshape(-1, 1)
    weight_list = []
    for i, count in enumerate(binned):
        weight_list += [(i + 1)] * int(count)
    x = np.array(weight_list).reshape(-1, 1)
    kde_values = kde_sklearn(x, x_grid, bw)
    return kde_values


def plot_exons(res_max: float, pos: float, exon_list: list) -> None:
    """
    Plot exon regions as bars.

    Parameters:
        res_max (float): Maximum residue value.
        pos (float): Vertical position for plotting.
        exon_list (list): List of exon boundaries.
    """
    plt.text(-(res_max / 7), pos + 0.2, "Exons list")
    for i in range(len(exon_list) - 1):
        exon_width = exon_list[i + 1] - exon_list[i]
        color = "silver" if (i + 1) % 2 == 0 else "black"
        plt.bar(
            exon_list[i] + exon_width / 2,
            0.6,
            width=exon_width,
            bottom=pos,
            facecolor=color,
            alpha=0.2,
        )


def plot_domain(res_max: float, pos: float, domains: list) -> None:
    """
    Plot Pfam domain annotations.

    Parameters:
        res_max (float): Maximum residue value.
        pos (float): Vertical position for plotting.
        domains (list): List of domains, each represented as (name, start, end).
    """
    plt.bar(res_max / 2, 0.6, width=res_max, bottom=pos, color="silver", alpha=0.2)
    plt.text(-(res_max / 7), pos + 0.2, "Pfam domain")
    for domain in domains:
        dom_name, dom_start, dom_end = domain
        loc = int((dom_start + dom_end) / 2)
        length = dom_end - dom_start
        plt.bar(loc, 0.6, width=length, bottom=pos, facecolor="lightcoral", alpha=0.8)
        label_y = pos + 0.2 if len(dom_name) <= 15 else pos + 0.1
        label_text = dom_name if len(dom_name) <= 15 else dom_name[:15]
        plt.text(loc, label_y, label_text, fontsize=8, horizontalalignment="center")


def plot_contlist(
    res_max: float, cont_list: list, pos: float, label: str, palette: str, reverse: bool
) -> None:
    """
    Plot a continuous data list as a bar plot.

    Parameters:
        res_max (float): Maximum residue value.
        cont_list (list): List of continuous values.
        pos (float): Vertical position for plotting.
        label (str): Label text to display.
        palette (str): Seaborn color palette to use.
        reverse (bool): If True, sort in descending order.
    """
    plt.text(-(res_max / 7), pos + 0.2, label)
    pair_list = [(i + 1, cont_list[i]) for i in range(len(cont_list))]
    pair_list.sort(key=lambda x: x[1], reverse=reverse)
    color_list = sns.color_palette(palette, len(cont_list))
    min_val = min(cont_list)
    zero_list = [pair for pair in pair_list if pair[1] == min_val]
    non_zero = [pair for pair in pair_list if pair[1] != min_val]
    plt.bar(
        [pair[0] for pair in zero_list],
        height=[0.6] * len(zero_list),
        width=1,
        bottom=pos,
        color=color_list[0],
    )
    plt.bar(
        [pair[0] for pair in non_zero],
        height=[0.6] * len(non_zero),
        width=1,
        bottom=pos,
        color=color_list[len(zero_list):],
    )


def plot_bilist(res_max: float, loc_list: list, pos: float, label: str, color: str) -> None:
    """
    Plot a binary list as a bar plot.

    Parameters:
        res_max (float): Maximum residue value.
        loc_list (list): List of locations to mark.
        pos (float): Vertical position for plotting.
        label (str): Label text to display.
        color (str): Color for the bars.
    """
    plt.text(-(res_max / 7), pos + 0.2, label)
    plt.bar(res_max / 2, 0.6, width=res_max, bottom=pos, color="silver", alpha=0.2)
    plt.bar(loc_list, [0.6] * len(loc_list), width=res_max / 1000, bottom=pos, color=color)


def visualization_1d(dicts,df_all,
     df_clust: pd.DataFrame, gene: str, output_folder: str
) -> None:
    """
    Generate a one-dimensional visualization of tiling mutagenesis screen data for missense mutations.

    Parameters:
        df_all (pd.DataFrame): DataFrame containing all sgRNA data.
        df_clust (pd.DataFrame): DataFrame with sgRNA clustering results.
        gene (str): Gene name.
        output_folder (str): Directory path to save the resulting figure.
        
    Note:
        This function expects global dictionaries 'exons_dic', 'alpha_dic', 'po_dic',
        'clinvar_dic', and 'domain_dic' to be defined.
    """
    exons_dic   = dicts["exons_dic"]
    domain_dic  = dicts["domain_dic"]
    alpha_dic   = dicts["alpha_dic"]
    po_dic      = dicts["po_dic"]
    clinvar_dic = dicts["clinvar_dic"]
    # Retrieve exon information
    res_max = exons_dic[gene][0]["AA.num"]
    exon_list = exons_dic[gene][0]["AA.list"]

    # Compute Alphafold score list if available
    if gene in alpha_dic:
        score_list = [np.mean(alpha_dic[gene][str(n)]) for n in range(1, len(alpha_dic[gene]) + 1)]
    else:
        score_list = [0] * res_max

    kde_score_list = get_sift_kde(score_list, list(range(1, res_max)), 10)

    # Retrieve phosphosite and ClinVar data if available
    if gene in po_dic:
        po_list = [s[0] for s in po_dic[gene] if s[1] > 0]
    else:
        po_list = [0] * res_max

    if gene in clinvar_dic:
        clinvar_list = [int(s[3:-3]) for s in clinvar_dic[gene] if s[-3:] != "Ter"]
    else:
        clinvar_list = [0] * res_max

    # Retrieve domain data if available
    dom_list = domain_dic.get(gene, [])

    plt.figure(figsize=(8, 5), dpi=300)
    plt.rcParams["font.size"] = 10
    plt.rcParams["font.family"] = "Helvetica"
    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Plot various tracks
    plot_contlist(res_max, kde_score_list, 6, "Alphafold", "Blues", reverse=False)
    plot_exons(res_max, 5, exon_list)
    plot_domain(res_max, 4, dom_list)
    plot_bilist(res_max, po_list, 3, "Phosphosite", "black")
    plot_bilist(res_max, clinvar_list, 2, "Pathogenicity", "red")

    # Plot all missense sgRNA locations
    df_gene_missense = df_all[(df_all["Gene"] == gene) & (df_all["Function"] == "missense")]
    aa_list_all = [aa for aa in list(df_gene_missense["AA"]) if aa != "NA"]
    plot_bilist(res_max, aa_list_all, 1, "sgRNAs all", "gray")

    # Plot significant sgRNAs per cluster for missense mutations
    df_clust_missense = df_clust[df_clust["Function"] == "missense"]
    color_list = ["red", "blue", "green", "orange", "purple", "yellow", "cyan", "pink", "teal", "brown"]
    for cluster in sorted(set(df_clust_missense["Cluster"])):
        df_cluster = df_clust_missense[df_clust_missense["Cluster"] == cluster]
        aa_list_res = []
        for aa in df_cluster["AA"]:
            if str(aa) != "NA":
                offset = random.randint(1, int(res_max / 100)) if aa not in aa_list_res else random.randint(2, int(res_max / 100))
                aa_list_res.append(aa + offset)
        plot_bilist(res_max, aa_list_res, -cluster, f"Cluster_{cluster+1}", color_list[cluster])

    plt.title(f"Tiling mutagenesis screen data of Missense sgRNAs targeting {gene}")
    plt.yticks([])
    plt.xticks([i - i % 10 for i in range(0, int(res_max), int(res_max / 100) * 10)])
    plt.xlabel("AA_Location")
    plt.xlim(0, res_max + 5)

    fig_save_path = os.path.join(os.getcwd(), output_folder, f"{gene}_ClusterDistribution_1D.png")
    plt.savefig(fig_save_path, dpi=300, facecolor="white", bbox_inches="tight")
    print(f"1D clusering map saved to: {fig_save_path}")
    plt.show()
    

def visualization_3d(gene, clust_ls, output_folder, pdb_path=None):
    
    if output_folder is not None and not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    if pdb_path is None:
        print(f"No pdb path provided for gene {gene}. Skip.")
        # Fill with 'NA' if you wish to maintain columns
        return None
    
    pymol.cmd.load(pdb_path)
    #pymol.cmd.hide('surface')
    pymol.cmd.hide('everything')
    pymol.cmd.color('gray')
    pymol.cmd.bg_color('white')
    
    pymol.cmd.show('cartoon')
    pymol.cmd.set('cartoon_transparency',0.7)
    #pymol.cmd.spectrum('b','blue_white_red',minimum=50,maximum=100)
    
    color_list = ["red", "blue", "green", "orange", "purple", 
                  "yellow", "cyan", "pink", "teal", "brown"]
    
    one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
                 'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y', \
                 'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A', \
                 'GLY':'G', 'PRO':'P', 'CYS':'C'}
    
    n = 0
    for ix_ls in clust_ls:
        sel = '+'.join([str(num) for num in ix_ls])
        print(sel)
        
        color = color_list[n]

        pymol.cmd.select('Cluster'+str(n+1), "resi "+sel)
        pymol.cmd.show('sticks', 'Cluster'+str(n+1))
        pymol.cmd.color(color, 'Cluster'+str(n+1))
        

        pymol.cmd.label('n. CA and i. '+sel, 'resi')
        pymol.cmd.set('label_size',10)
        pymol.cmd.set('stick_radius',0.5)
        
        n += 1
    
    pymol.cmd.set('stick_transparency',0)  
    
    print('Save pymol session and images ...')
    
    #figsavepath = os.path.join(os.getcwd(),outputfolder,gene+'_BE_Cluster_3D.png')
    pymolsavepath = os.path.join(os.getcwd(),output_folder,gene+'_ClustersDistribution_3D.pse')
    pymol.cmd.save(pymolsavepath)
    pymol.cmd.reinitialize()
    #pymol.cmd.png(figsavepath,width=1200, height=1200, dpi=300, ray=0)
    print('Finished') 

__all__ = [
    "clustering","annotation","download_and_load_jsons",
    "visualization_1d","visualization_3d",
]

if __name__ == "__main__":
    main()

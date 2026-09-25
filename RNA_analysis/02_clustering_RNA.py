from __future__ import annotations

import anndata as ad
import scanpy as sc

sc.set_figure_params(dpi=100, facecolor="white")

adata_RNA = sc.read_h5ad("after_qc.h5ad")
print(adata_RNA)


#6 PCA → DIMENSIONALITY REDUCTION
#PCA reveals main axes of variation + does denoising

sc.tl.pca(adata_RNA, n_comps=30, use_highly_variable=True)
sc.pl.pca_variance_ratio(adata_RNA, n_pcs=30, log=True, save="_variance.pdf") #variance plot
sc.pl.pca( # principal components plot
    adata_RNA,
    color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
    save="_qc_metrics_plot.pdf",
)

#Reading the PCA Plots
#Top 2 (Cells by gene expression): Each color is a sample (iPSC, classical, etc.)
    #Cells close together = similar gene expression
    #Long stretch right = different population of cells
#Bottom 2 (By Mt genes): Purple = low mt gene %, yellow/green = high mt gene %, tail right is higher mt %

#Reading PCA variance plot
    #Looking for a curve which flats out (elbow shape)
    #Where elbow first appears (so PC # whatever) means 0-# captures most variation
    #So in neighbors step youd set n_pcs=15 or something close


#7 NEAREST NEIGHBOR GRAPHS (UMAP)
#If you inspect batch effects in your UMAP it can be beneficial
    #to integrate across samples and perform batch correction/integration.

sc.pp.neighbors(adata_RNA, n_pcs=15)
sc.tl.umap(adata_RNA)
sc.pl.umap(
    adata_RNA,
    color=["sample", "total_counts", "n_genes_by_counts", "pct_counts_mt"],
    size=2,
    ncols=2,
    save="_neighbor_umap.pdf",
)


#8 CLUSTERING (LEIDEN)

# Using the igraph implementation and a fixed number of iterations can be significantly faster,
# especially for larger datasets
sc.tl.leiden(adata_RNA, flavor="igraph", n_iterations=2, directed=False)
sc.pl.umap(adata_RNA, color=["leiden"], save="_leiden_cluster.pdf")

#Reading Leiden UMAP
    #Each dot = 1 cell
    #Cluster of cells = similar gene expression
    #Each color = one cluster IDd by Leiden alg
    #Distinct "blobs" are cell clusters with similar transcriptomes


#9 RE-ACCESS QUALITY CONTROL/FILTERING
#Visualize QC metrics in UMAP

#BELOW UMAP DOES NOT APPLY TO ncMo and cMo
sc.pl.umap(
    adata_RNA,
    color=["leiden", "predicted_doublet", "doublet_score"],
    #increase horizontal space between panels
    wspace=0.5,
    size=3,
    save="_doublet_quality_umap.pdf",
)

sc.pl.umap(
    adata_RNA,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
    save=_"gene_count_umap.pdf",
)

#Reading QC UMAPS
    #Doublet UMAP
        #Middle Graph (Predicted): Blue means singlets, orange means likely doublets
        #Right (Doublet Score): Purple = low doublet score, Yellow/Green = High score
    #Mt and gene and total UMAP
        #Total_Count: Purple = low total counts, yellow = high total counts (should be uniform across UMAP)
        #Mt count: Purple = low mito %, yellow = high mito %
        #Genes by count: Purple = few genes detected, yellow = many genes detected (very high count could indicate doublet)


#Batch Correction (If needed) NOT APPLICABLE HERE

#Picking Ideal Resolution

sc.tl.leiden(adata_RNA, resolution=0.02, key_added='leiden_res_0.02')
sc.tl.leiden(adata_RNA, resolution=0.50, key_added='leiden_res_0.50')
sc.tl.leiden(adata_RNA, resolution=2.00, key_added='leiden_res_2.00')

sc.pl.umap(
    adata_RNA,
    color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
    legend_loc="on data",
)

#10 MARKER GENE SET DEFINITION
#use Buchrieser and calderon for iPSC mono and primary mono respectively
#can also use CellMarker and PanglaoDB

marker_genes = {
    "iMono": ["CD14", "FCGR3A", "CD163", "POU5F1", "SOX2", "NANOG", "KLF4", "DNMT3B", "EPCAM"], #FCGR3A is CD16, POU5F1 is OCT4
    "cMo": ["CD14", "S100A8", "S100A9", "CCR2", "LYZ"], #CellMarker 2.0 and Ziegler (for all)
    #"intermediate mono": ["FCGR3A", "CD14", "HLA-DR", "MARCO"], #CellMarker 2.0 #not used
    "ncMo": ["FCGR3A", "CD14", "HLA-DRA", "ITGAM", "CX3CR1"], #ITGAM is CD11B

}

all_genes = [g for gs in marker_genes.values() for g in gs]
missing = [g for g in all_genes if g not in adata_RNA.var_names]
print("Missing:", missing)
marker_genes = {k: [g for g in v if g in adata_RNA.var_names] for k, v in marker_genes.items()}

sc.pl.dotplot(adata_RNA, marker_genes, groupby="leiden_res_0.50", standard_scale="var") #can change res

sc.pl.dotplot(adata_RNA, marker_genes, groupby="sample", standard_scale="var", save="_markers_by_sample.pdf")

imono_sub = adata_RNA[adata_RNA.obs["sample"] == "iMono"].copy()

imono_panel = {
    "monocyte": ["CD14", "FCGR3A", "LYZ", "S100A8", "CD163", "ITGAM"],
    "pluripotency": ["POU5F1", "NANOG", "SOX2", "LIN28A", "DNMT3B"],
    "proliferation": ["MKI67", "TOP2A"],
}

missing_panel = [g for gs in imono_panel.values() for g in gs if g not in imono_sub.var_names]
print("Missing from iMono panel:", missing_panel)
imono_panel = {k: [g for g in v if g in imono_sub.var_names] for k, v in imono_panel.items()}

sc.pl.dotplot(imono_sub, imono_panel, groupby="leiden_res_0.50", standard_scale="var", save="_markers_imono.pdf")

#11 DIFFERENTLY EXPRESSED GENES AS MARKERS

# Obtain cluster-specific differentially expressed genes
sc.tl.rank_genes_groups(imono_sub, groupby="leiden_res_0.50", method="wilcoxon")

#Visualize top 5 DEGs on dot plot
sc.pl.rank_genes_groups_dotplot(imono_sub, groupby="leiden_res_0.50", standard_scale="var", n_genes=5)

#Can also get top DEGs as a table
for cl in imono_sub.obs["leiden_res_0.50"].cat.categories:
    print(f"--- cluster {cl} ---")
    print(sc.get.rank_genes_groups_df(imono_sub, group=cl).head(5))

#And Cluster DEGs #May be used after initial run
#dc_cluster_genes = sc.get.rank_genes_groups_df(adata_RNA, group="7").head(5)["names"]
#sc.pl.umap(
    #adata_RNA,
    #color=[*dc_cluster_genes, "leiden_res_0.50"],
    #legend_loc="on data",
    #frameon=False,
    #ncols=3,
#)

adata_RNA.write_h5ad("after_clustering.h5ad")
imono_sub.write_h5ad("imono_clustered.h5ad")


#then move onto DESeq2



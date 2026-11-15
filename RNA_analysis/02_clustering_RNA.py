from __future__ import annotationsa

import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")

sc.set_figure_params(dpi=100, facecolor="white")
sc.settings.figdir = "figures/clustering"

adata_RNA = sc.read_h5ad("after_qc.h5ad")
print(adata_RNA)

group_of = {"MV1": "hiPSC", "MV2": "hiPSC", "MV3": "hiPSC", "MV5": "hiPSC", "MV6": "hiPSC", "Control": "hiPSC",
            "iMono": "iMono", "cMo": "cMo", "intMo": "intMo", "ncMo": "ncMo"}
adata_RNA.obs["group"] = pd.Categorical(adata_RNA.obs["sample"].map(group_of),
                                        categories=["hiPSC", "iMono", "cMo", "intMo", "ncMo"])
print(adata_RNA.obs["group"].value_counts())

#5 Highly Variable Genes 
sc.pp.highly_variable_genes(adata_RNA, n_top_genes=2000) #can change this value if needed (2000 is default)

#6 PCA → DIMENSIONALITY REDUCTION
#PCA reveals main axes of variation + does denoising

sc.tl.pca(adata_RNA, n_comps=30, use_highly_variable=True)
sc.pl.pca_variance_ratio(adata_RNA, n_pcs=30, log=True, save="_variance.pdf") #variance plot
sc.pl.pca( # principal components plot
    adata_RNA,
    color=["group", "group", "pct_counts_mt", "pct_counts_mt"],
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

N_pcs=15 #can edit

sc.pp.neighbors(adata_RNA, n_pcs=15)
sc.tl.umap(adata_RNA)
sc.pl.umap(adata_RNA, color=["group", "sample"], size=2, wspace=0.4, save="_by_sample.pdf")
sc.pl.umap(
    adata_RNA,
    color=["total_counts", "n_genes_by_counts", "pct_counts_mt"],
    size=2,
    ncols=3,
    save="_qc_metrics.pdf",
)


#8 CLUSTERING (LEIDEN)

# Using the igraph implementation and a fixed number of iterations can be significantly faster,
# especially for larger datasets
for res in [0.1, 0.5, 1.0]:
    sc.tl.leiden(adata_RNA, resolution=res, key_added=f"leiden_res_{res}",
                 flavor="igraph", n_iterations=2, directed=False)
    
sc.pl.umap(
    adata_RNA,
    color=["leiden_res_0.1", "leiden_res_0.5", "leiden_res_1.0"],
    legend_loc="on data",
    size=2,
    save="_leiden_resolutions.pdf",
)
Cluster_key = "leiden_res_0.5"

#Checking if the iMonos are closer to hiPSCs or primaries 
ct = pd.crosstab(adata_RNA.obs[Cluster_key], adata_RNA.obs["sample"])
ct.to_csv("cluster_by_sample.csv")
print(ct)

#Reading Leiden UMAP
    #Each dot = 1 cell
    #Cluster of cells = similar gene expression
    #Each color = one cluster IDd by Leiden alg
    #Distinct "blobs" are cell clusters with similar transcriptomes


#9 RE-ACCESS QUALITY CONTROL/FILTERING
#Visualize QC metrics in UMAP

sc.pl.umap(
    adata_RNA,
    color=[Cluster_key, "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
    size=2,
    save="_qc_by_cluster.pdf",
)

sc.pl.violin(adata_RNA, ["n_genes_by_counts", "pct_counts_mt"], groupby=Cluster_key,
             stripplot=False, save="_qc_by_cluster.pdf")
#Reading QC UMAPS
    #Doublet UMAP
        #Middle Graph (Predicted): Blue means singlets, orange means likely doublets
        #Right (Doublet Score): Purple = low doublet score, Yellow/Green = High score
    #Mt and gene and total UMAP
        #Total_Count: Purple = low total counts, yellow = high total counts (should be uniform across UMAP)
        #Mt count: Purple = low mito %, yellow = high mito %
        #Genes by count: Purple = few genes detected, yellow = many genes detected (very high count could indicate doublet)


#10 MARKER GENE SET DEFINITION
#use Buchrieser and calderon for iPSC mono and primary mono respectively
#can also use CellMarker and PanglaoDB

marker_genes = {
    "monocyte": ["CD14", "LYZ", "S100A8", "S100A9", "CSF1R", "ITGAM", "CD300E", "HLA-DRA"],   #ITGAM = CD11B
    "cMo": ["CCR2", "VCAN", "FCN1"],
    "ncMo": ["FCGR3A", "CX3CR1"],                              #FCGR3A = CD16
    "macrophage": ["CD163", "MRC1", "APOE", "SPP1"],
    "pluripotency": ["POU5F1", "SOX2", "NANOG", "LIN28A", "DNMT3B"],      #POU5F1 is OCT4
    "proliferation": ["MKI67", "TOP2A"],
}

missing = [g for gs in marker_genes.values() for g in gs if g not in adata_RNA.var_names]
print("Missing:", missing)
marker_genes = {k: [g for g in v if g in adata_RNA.var_names] for k, v in marker_genes.items()}
marker_genes = {k: v for k, v in marker_genes.items() if v}

sc.pl.dotplot(adata_RNA, marker_genes, groupby="sample", standard_scale="var", save="_markers_by_sample.pdf")
sc.pl.dotplot(adata_RNA, marker_genes, groupby=Cluster_key, standard_scale="var", save="_markers_by_cluster.pdf")

#11 DIFFERENTLY EXPRESSED GENES AS MARKERS

sc.tl.rank_genes_groups(adata_RNA, groupby=Cluster_key, method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(adata_RNA, groupby=Cluster_key, standard_scale="var", n_genes=5,
                                save="_top_DEGs_by_cluster.pdf")
deg = sc.get.rank_genes_groups_df(adata_RNA, group=None)
deg.to_csv("DEGs_by_cluster.csv", index=False)
for cl in adata_RNA.obs[Cluster_key].cat.categories:
    print(f"--- cluster {cl} ---")
    print(deg[deg["group"] == cl].head(5)[["names", "logfoldchanges", "pvals_adj"]])

#And Cluster DEGs #May be used after initial run
#dc_cluster_genes = sc.get.rank_genes_groups_df(adata_RNA, group="7").head(5)["names"]
#sc.pl.umap(
    #adata_RNA,
    #color=[*dc_cluster_genes, "leiden_res_0.50"],
    #legend_loc="on data",
    #frameon=False,
    #ncols=3,
#)

#12 Making iMONO subclusters

imono_sub = adata_RNA[adata_RNA.obs["sample"] == "iMono"].copy()
sc.pp.neighbors(imono_sub, n_pcs=N_pcs)
sc.tl.umap(imono_sub)
sc.tl.leiden(imono_sub, resolution=0.5, key_added="leiden_imono",
             flavor="igraph", n_iterations=2, directed=False)
sc.pl.umap(imono_sub, color=["leiden_imono", "pct_counts_mt", "n_genes_by_counts"],
           size=5, wspace=0.4, save="_imono_subclusters.pdf")

imono_panel = {
    "monocyte": ["CD14", "FCGR3A", "LYZ", "S100A8", "CD163", "ITGAM"],
    "macrophage": ["APOE", "SPP1", "MRC1"],
    "pluripotency": ["POU5F1", "NANOG", "SOX2", "LIN28A", "DNMT3B"],
    "proliferation": ["MKI67", "TOP2A"],
}
imono_panel = {k: [g for g in v if g in imono_sub.var_names] for k, v in imono_panel.items()}
imono_panel = {k: v for k, v in imono_panel.items() if v}
sc.pl.dotplot(imono_sub, imono_panel, groupby="leiden_imono", standard_scale="var", save="_markers_imono.pdf")

sc.tl.rank_genes_groups(imono_sub, groupby="leiden_imono", method="wilcoxon")
sc.get.rank_genes_groups_df(imono_sub, group=None).to_csv("DEGs_imono_subclusters.csv", index=False)

#SAVING
adata_RNA.write_h5ad("after_clustering.h5ad")
imono_sub.write_h5ad("imono_clustered.h5ad")
print("Finished with no errors. Figures in figures/clustering, tables: cluster_by_sample.csv, DEGs_by_cluster.csv, DEGs_imono_subclusters.csv")


#then move onto DESeq2



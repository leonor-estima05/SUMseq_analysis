#Installed scanpy anndata
#Installed leidenalg igraph
#Installed pooch (background manager for scanpy to import datasets)-not needed in actual data
#Installed scikit-image
#Installed scanpy[scrublet]
#Installed harmonypy 0.0.6


from __future__ import annotations

import anndata as ad
import pooch
import scanpy as sc

sc.set_figure_params(dpi=100, facecolor="white")

#1. DATA LOADING
#Real Data Formats:
#RNA: barcodes.tsv(cell IDs), features.tsv(gene data), matrix.mtx and a merged seurat object .qs
#ATAC: fragments as .bed.gz, ArchR .arrow (not scanpy!)

adata_RNA = sc.read_10x_mtx("file path to folder with the 3 RNA files") #or other filetype

#samples = {
    #"iMono: "file path to folder with the 3 RNA files for iPSCs",
    #"cMo": "file path to folder with the 3 RNA files for monocytes",
    #"intMo": "",
    #"ncMo": "",
#}

#RNA_datasets = {}

#for sample_id, filename in samples.items():
    #sample_adata = sc.read_10x_mtx(path) #or _h5 depending on file format
    #sample_adata.var_names_make_unique()
    #RNA_datasets[sample_id] = sample_adata

#adata_RNA = ad.concat(RNA_datasets, label = "sample")
#adata_RNA.obs_names_make_unique()
#print(adata_RNA.obs["sample"].value_counts()) #shows how many cells per cell type
#print(adata_RNA) #will show number of cells x genes


#2 QUALITY CONTROL
#mostly uses sc.pp.calculate_qc_metrics()
#measures mitochondrial, ribosomal and hemoglobin (important bc it would indicate RBC contamination)

# mitochondrial genes, "MT-" for human, "Mt-" for mouse
#adata_RNA.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
#adata_RNA.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
#adata_RNA.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

#sc.pp.calculate_qc_metrics(adata_RNA, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

#sc.pl.violin(
    #adata_RNA,
    #["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb""],
    #groupby="sample",
    #jitter=0.4,
    #multi_panel=True,
#)

#For my violin and scatter plots im seperating by cell type
#Violin plots tell you # of genes expressed in count matrix, total counts per cell, % of counts in mitochondrial genes
#Reading violin plots
#Wide part: many cells have this value
#Narrow: few cells have value
#Mt: Most cells should cluster low, cells on long tail upwards are cells with high mt % (cut off tail)
#Gene Count: Cells at flat bottom, few genes. Outliers at top (may be doublets). Want normal distribution
#Total Count: Want normal distribution. Outliers remove on either side.

# mt vs total counts
#sc.pl.scatter(adata_RNA, x='total_counts', y='pct_counts_mt', color='sample')
# hb vs total counts
#sc.pl.scatter(adata_RNA, x='total_counts', y='pct_counts_hb', color='sample')
# ribo vs total counts
#sc.pl.scatter(adata_RNA, x='total_counts', y='pct_counts_ribo', color='sample')
# genes vs total counts (for doublet detection)
#sc.pl.scatter(adata_RNA, x='total_counts', y='n_genes_by_counts', color='sample')

#Based on the QC Metric plots filter out low quality cells/genes-change based on findings
#They filtered cells with less than 100 genes expressed and genes detected in less than 3 cells
#Set cell # by 0.1% of total cells

#print(f"Before filtering: {adata_RNA.shape[0]} cells, {adata_RNA.shape[1]} genes")


# filter by mitochondrial %
#adata_RNA = adata_RNA[adata_RNA.obs['pct_counts_mt'] < 20] #cells with more than 20% mt genes
# filter by hemoglobin %
#adata_RNA = adata_RNA[adata_RNA.obs['pct_counts_hb'] < 1] #cells with more than 1% hb genes
#sc.pp.filter_cells(adata_RNA, min_genes=100) #cells with less than 100 genes expressed
#sc.pp.filter_genes(adata_RNA, min_cells=3) #genes in less than 3 cells

#print(f"After filtering: {adata_RNA.shape[0]} cells, {adata_RNA.shape[1]} genes")
#print("Filtering complete")


#3 DOUBLET REMOVAL (Warning: Scrublet takes a while)
#Run doublet detection algorithm called Scrublet (Wolock et. al., 2019)

#sc.pp.scrublet(adata_RNA, batch_key="sample")
# filter out predicted doublets
#adata_RNA = adata_RNA[adata_RNA.obs['predicted_doublet'] == False]
#print(f"After doublet removal: {adata_RNA.shape[0]} cells")


#Saving count data before normalization for DESeq2

#adata_forDES = adata.copy()
#adata_forDES.write_h5ad("adata_forDES.h5ad")

#4 NORMALIZATION
#Using count depth scaling + log + 1 transformation
#Via target_sum in pp.normalize_total

# Normalizing to median total counts
#sc.pp.normalize_total(adata_RNA)
# Logarithmize the data
#sc.pp.log1p(adata_RNA)


#5 FEATURE SELECTION
#Used to include only most informative genes for dimensionality reduction
#Uses scanpy function pp.highly_variable_genes

#sc.pp.highly_variable_genes(adata_RNA, n_top_genes=2000, batch_key="sample")
#sc.pl.highly_variable_genes(adata_RNA)
# then filter to keep only highly variable genes
#adata_RNA = adata_RNA[:, adata_RNA.var['highly_variable']]


#6 PCA → DIMENSIONALITY REDUCTION
#PCA reveals main axes of variation + does denoising

#sc.tl.pca(adata_RNA)
#sc.pl.pca_variance_ratio(adata_RNA, n_pcs=50, log=True) #variance plot
#sc.pl.pca( # principal components plot
    #adata_RNA,
    #color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
    #dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    #ncols=2,
    #size=2,
#)

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
    #We recommend checking out scanorama and scvi-tools for batch integration.

#sc.pp.neighbors(adata_RNA)
#sc.tl.umap(adata_RNA)
#sc.pl.umap(
    #adata_RNA,
    #color="sample",
    #size=2,
#)


#8 CLUSTERING (LEIDEN)

# Using the igraph implementation and a fixed number of iterations can be significantly faster,
# especially for larger datasets
#sc.tl.leiden(adata_RNA, flavor="igraph", n_iterations=2)
#sc.pl.umap(adata_RNA, color=["leiden"])

#Reading Leiden UMAP
    #Each dot = 1 cell
    #Cluster of cells = similar gene expression
    #Each color = one cluster IDd by Leiden alg
    #Distinct "blobs" are cell clusters with similar transcriptomes


#9 RE-ACCESS QUALITY CONTROL/FILTERING
#Visualize QC metrics in UMAP

#sc.pl.umap(
    #adata,
    #color=["leiden", "predicted_doublet", "doublet_score"],
    # increase horizontal space between panels
    #wspace=0.5,
    #size=3,
#)

#sc.pl.umap(
    #adata,
    #color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    #wspace=0.5,
    #ncols=2,
#)

#Reading QC UMAPS
    #Doublet UMAP
        #Middle Graph (Predicted): Blue means singlets, orange means likely doublets
        #Right (Doublet Score): Purple = low doublet score, Yellow/Green = High score
    #Mt and gene and total UMAP
        #Total_Count: Purple = low total counts, yellow = high total counts (should be uniform across UMAP)
        #Mt count: Purple = low mito %, yellow = high mito %
        #Genes by count: Purple = few genes detected, yellow = many genes detected (very high count could indicate doublet)


#Batch Correction (If needed)

#sc.external.pp.harmony_integrate(adata_RNA, key='sample')
#sc.pp.neighbors(adata_RNA, use_rep='X_pca_harmony')
#sc.tl.umap(adata_RNA)
#sc.tl.leiden(adata_RNA, resolution=0.5)
#sc.pl.umap(adata_RNA, color=['leiden', 'sample'])


#Picking Ideal Resolution

#sc.tl.leiden(adata_RNA, resolution=0.02, key_added='leiden_res_0.02')
#sc.tl.leiden(adata_RNA, resolution=0.50, key_added='leiden_res_0.50')
#sc.tl.leiden(adata_RNA, resolution=2.00, key_added='leiden_res_2.00')

#sc.pl.umap(
    #adata_RNA,
    #color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
    #legend_loc="on data",
#)

#10 MARKER GENE SET DEFINITION
#use Buchrieser and calderon for iPSC mono and primary mono respectively
#can also use CellMarker and PanglaoDB

#marker_genes = {
    #"iPSC mono": ["CD14", "FCGR3A", "CD163", "POU5F1", "SOX2", "NANOG", "KLF4", "DNMT38", "EPCAM"], #FCGR3A is CD16, POU5F1 is OCT4
    #"classical mono": ["CD14", "S100A8", "S100A9", "CCR2", "LYZ"], #CellMarker 2.0 and Ziegler (for all)
    #"intermediate mono": ["FCGR3A", "CD14", "HLA-DR", "MARCO"], #CellMarker 2.0
    #"nonclassical mono": ["FCGR3A", "CD14", "HLA-DR", "ITGAM", "CX3CR1"], #ITGAM is CD11B

#}

#sc.pl.dotplot(adata_RNA, marker_genes, groupby="leiden_res_0.5", standard_scale="var") #can change res

#Then label according to cell type

#adata.obs["cell_type_lvl1"] = adata.obs["leiden_res_0.02"].map( #ONLY KNOWN AFTER DOTPLOTS
    #{
        #"0": "iPSC mono",
        #"1": "classical mono",
        #"2": "intermediate mono",
        #"3": "nonclassical mono",
    #}
#)

#sc.pl.dotplot(adata_RNA, marker_genes, groupby="leiden_res_0.50", standard_scale="var")

#11 DIFFERENTLY EXPRESSED GENES AS MARKERS

# Obtain cluster-specific differentially expressed genes
#sc.tl.rank_genes_groups(adata_RNA, groupby="leiden_res_0.50", method="wilcoxon")

#Visualize top 5 DEGs on dot plot
#sc.pl.rank_genes_groups_dotplot(adata_RNA, groupby="leiden_res_0.50", standard_scale="var", n_genes=5)

#Can also get top DEGs as a table
#sc.get.rank_genes_groups_df(adata_RNA, group="7").head(5)

#And Cluster DEGs
#dc_cluster_genes = sc.get.rank_genes_groups_df(adata_RNA, group="7").head(5)["names"]
#sc.pl.umap(
    #adata_RNA,
    #color=[*dc_cluster_genes, "leiden_res_0.50"],
    #legend_loc="on data",
    #frameon=False,
    #ncols=3,
#)

#then move onto DESeq2







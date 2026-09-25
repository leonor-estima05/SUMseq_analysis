#Installed scanpy anndata
#Installed leidenalg igraph
#Installed pooch (background manager for scanpy to import datasets)-not needed in actual data
#Installed scikit-image
#Installed scanpy[scrublet]
#Installed harmonypy 0.0.6 #will not be used because it makes little biological sense


from __future__ import annotations

import anndata as ad
import pooch
import scanpy as sc

sc.set_figure_params(dpi=100, facecolor="white")

#1. DATA LOADING
#Real Data Formats:
#RNA: barcodes.tsv(cell IDs), features.tsv(gene data), matrix.mtx and a merged seurat object .qs
#ATAC: fragments as .bed.gz, ArchR .arrow (not scanpy!)

base = "/home/l.estima/scanpy_input" #dir with individual sample mtx files (wrapped into .gz so read_10x_mtx will work)

samples = {
    "iMono": f"{base}/iMono",
    "cMo": f"{base}/cMo",
    #"intMo": f"{base}/intMo", #intMo will not be run, at 3 median genes per cell it would produce too much noise.
    "ncMo": f"{base}/ncMo",
}

RNA_datasets = {}

for sample_id, path in samples.items():
    sample_adata = sc.read_10x_mtx(path) #or _h5 depending on file format
    sample_adata.var_names_make_unique()
    RNA_datasets[sample_id] = sample_adata

adata_RNA = ad.concat(RNA_datasets, label = "sample", index_unique = "_")
print(adata_RNA.obs["sample"].value_counts())
print(adata_RNA) #will show number of cells x genes


#2 QUALITY CONTROL
#mostly uses sc.pp.calculate_qc_metrics()
#measures mitochondrial, ribosomal and hemoglobin (important bc it would indicate RBC contamination)

#mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata_RNA.var["mt"] = adata_RNA.var_names.str.startswith("MT-")
#ribosomal genes
adata_RNA.var["ribo"] = adata_RNA.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata_RNA.var["hb"] = adata_RNA.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(adata_RNA, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

sc.pl.violin(
    adata_RNA,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
    groupby="sample",
    jitter=0.4,
    multi_panel=True,
)

#For my violin and scatter plots im seperating by cell type
#Violin plots tell you # of genes expressed in count matrix, total counts per cell, % of counts in mitochondrial genes
#Reading violin plots
#Wide part: many cells have this value
#Narrow: few cells have value
#Mt: Most cells should cluster low, cells on long tail upwards are cells with high mt % (cut off tail)
#Gene Count: Cells at flat bottom, few genes. Outliers at top (may be doublets). Want normal distribution
#Total Count: Want normal distribution. Outliers remove on either side.

# mt vs total counts
sc.pl.scatter(adata_RNA, x='total_counts', y='pct_counts_mt', color='sample')
# hb vs total counts
sc.pl.scatter(adata_RNA, x='total_counts', y='pct_counts_hb', color='sample')
# ribo vs total counts
sc.pl.scatter(adata_RNA, x='total_counts', y='pct_counts_ribo', color='sample')
# genes vs total counts (for doublet detection)
sc.pl.scatter(adata_RNA, x='total_counts', y='n_genes_by_counts', color='sample')

#Based on the QC Metric plots filter out low quality cells/genes-change based on findings
#They filtered cells with less than 20 genes expressed and genes detected in less than 3 cells, this is low due to low sequencing depth
#Set cell # by 0.1% of total cells

print(f"Before filtering: {adata_RNA.shape[0]} cells, {adata_RNA.shape[1]} genes")


# filter by mitochondrial %
adata_RNA = adata_RNA[adata_RNA.obs['pct_counts_mt'] < 20].copy() #cells with more than 20% mt genes
# filter by hemoglobin %
adata_RNA = adata_RNA[adata_RNA.obs['pct_counts_hb'] < 1].copy() #cells with more than 1% hb genes
sc.pp.filter_cells(adata_RNA, min_genes=20) #cells with less than 20 genes expressed, quite low due to low depth reads in this data
sc.pp.filter_genes(adata_RNA, min_cells=3) #genes in less than 3 cells

print(f"After filtering: {adata_RNA.shape[0]} cells, {adata_RNA.shape[1]} genes")
print("Filtering complete")

#Save full unfiltered counts (all genes, all cells) as a DESeq2 fallback
adata_RNA_full = adata_RNA.copy()
adata_RNA_full.write_h5ad("adata_full_counts.h5ad")


#3 DOUBLET REMOVAL (Warning: Scrublet takes a while)
#Run doublet detection algorithm called Scrublet (Wolock et. al., 2019)
#This will not be run on ncMo and cMo due to low median/cell gene counts

imono = adata_RNA[adata_RNA.obs["sample"] == "iMono"].copy()
sc.pp.scrublet(imono)
#will output predicted doublet scores
adata_RNA.obs["predicted_doublet"] = False
adata_RNA.obs.loc[imono.obs_names, "predicted_doublet"] = imono.obs["predicted_doublet"]
adata_RNA.obs["doublet_score"] = float("nan")
adata_RNA.obs.loc[imono.obs_names, "doublet_score"] = imono.obs["doublet_score"]
print(f"Predicted doublets in iMono: {imono.obs['predicted_doublet'].sum()} of {imono.n_obs} cells")


#Saving filtered count data before normalization for DESeq2

adata_RNA_forDES = adata_RNA.copy()
adata_RNA_forDES.write_h5ad("adata_filtered_counts.h5ad")

adata_RNA.layers["counts"] = adata_RNA.X.copy()

#4 NORMALIZATION
#Using count depth scaling + log + 1 transformation
#Via target_sum in pp.normalize_total

#Normalizing to median total counts
sc.pp.normalize_total(adata_RNA, target_sum=1e4)
# Logarithmize the data
sc.pp.log1p(adata_RNA)


#5 FEATURE SELECTION
#Used to include only most informative genes for dimensionality reduction
#Uses scanpy function pp.highly_variable_genes

sc.pp.highly_variable_genes(adata_RNA, n_top_genes=500) #low top gene number used due to low count data
sc.pl.highly_variable_genes(adata_RNA)
# then filter to keep only highly variable genes, not doing this due to later gene expression analysis combined with low count data
#adata_RNA = adata_RNA[:, adata_RNA.var['highly_variable']].copy()


#6 PCA → DIMENSIONALITY REDUCTION
#PCA reveals main axes of variation + does denoising

sc.tl.pca(adata_RNA, n_comps=30, use_highly_variable=True)
sc.pl.pca_variance_ratio(adata_RNA, n_pcs=30, log=True) #variance plot
sc.pl.pca( # principal components plot
    adata_RNA,
    color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
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
)


#8 CLUSTERING (LEIDEN)

# Using the igraph implementation and a fixed number of iterations can be significantly faster,
# especially for larger datasets
sc.tl.leiden(adata_RNA, flavor="igraph", n_iterations=2, directed=False)
sc.pl.umap(adata_RNA, color=["leiden"])

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
)

sc.pl.umap(
    adata_RNA,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
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

sc.pl.dotplot(adata_RNA, marker_genes, groupby="sample", standard_scale="var")

imono_sub = adata_RNA[adata_RNA.obs["sample"] == "iMono"].copy()

imono_panel = {
    "monocyte": ["CD14", "FCGR3A", "LYZ", "S100A8", "CD163", "ITGAM"],
    "pluripotency": ["POU5F1", "NANOG", "SOX2", "LIN28A", "DNMT3B"],
    "proliferation": ["MKI67", "TOP2A"],
}

missing_panel = [g for gs in imono_panel.values() for g in gs if g not in imono_sub.var_names]
print("Missing from iMono panel:", missing_panel)
imono_panel = {k: [g for g in v if g in imono_sub.var_names] for k, v in imono_panel.items()}

sc.pl.dotplot(imono_sub, imono_panel, groupby="leiden_res_0.50", standard_scale="var")

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


#then move onto DESeq2







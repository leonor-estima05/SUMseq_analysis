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
    save="_qc_metrics.pdf",
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
sc.pl.scatter(adata_RNA, x="total_counts", y="pct_counts_mt", color="sample", save="_mt_vs_counts.pdf")
# hb vs total counts
sc.pl.scatter(adata_RNA, x="total_counts", y="pct_counts_hb", color="sample", save="_hb_vs_counts.pdf")
# ribo vs total counts
sc.pl.scatter(adata_RNA, x="total_counts", y="pct_counts_ribo", color="sample", save="_ribo_vs_counts.pdf")
# genes vs total counts (for doublet detection)
sc.pl.scatter(adata_RNA, x="total_counts", y="n_genes_by_counts", color="sample", save="_genes_vs_counts.pdf")

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
sc.pl.highly_variable_genes(adata_RNA, save="_hvg.pdf")
# then filter to keep only highly variable genes, not doing this due to later gene expression analysis combined with low count data
#adata_RNA = adata_RNA[:, adata_RNA.var['highly_variable']].copy()

adata_RNA.write_h5ad("after_qc.h5ad")
print(f"Saved: {adata_RNA.shape[0]} cells, {adata_RNA.shape[1]} genes")

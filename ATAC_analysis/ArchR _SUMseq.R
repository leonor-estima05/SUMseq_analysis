
library(ArchR)
set.seed(1)

addArchRThreads(threads=8)
addArchRGenome("hg38")

#1 DATA LOADING

in_rds <- "/home/projects/icell/sumseq_pilot/ATAC/ATAC_outdir/ATAC/4.ArchR_projectFolder/ArchRProject.minFrags0.rds"
out_dir <- "/home/projects/icell/sumseq_pilot/ATAC/archr_analysis/01_qc"
plot_dir <- "/home/projects/icell/sumseq_pilot/ATAC/archr_analysis/plots/01_qc"


MIN_FRAGS <- 1000 #minimum unique fragments per cell 
MIN_TSS <- 4 #minimum TSS enrichment score
FILTER_RATIO <- 1 #for filterdoublets

DROP_SAMPLES <- c(0) #kept to 0 for now but can change to c("cMo_s8") later

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

proj <- readRDS(in_rds)
cat("\nBarcodes per sample (unfiltered):\n")
print(table(proj$Sample)


#2 QUALITY CONTROL

#Quality Control Per Cell-likely already done by snakemake
#1. Number of unique nuclear fragments (not MT DNA), cells with too little fragments are not useful
#2. Signal to background ratio, low ratio is often dead/dying cells
#3. Fragment size distribution. Expect depletion of fragments that are length of DNA wrapped around a nucleosome (~147bp)

qc_df <- data.frame(
  Sample = proj$Sample,
  logFrags = log10(proj$nFrags),
  TSS = proj$TSSEnrichment
)

qc_df <- qc_df <- qc_df[qc_df$logFrags >= log10(300), ]   # drop near-empty barcodes to keep the plot light

qc_plot <- ggplot(qc_df, aes(x = logFrags, y = TSS)) +
  geom_bin2d(bins = 100) +
  scale_fill_viridis_c(trans = "log10", name = "Barcodes") +
  geom_hline(yintercept = MIN_TSS, lty = "dashed") +
  geom_vline(xintercept = log10(MIN_FRAGS), lty = "dashed") +
  facet_wrap(~ Sample) +
  coord_cartesian(ylim = c(0, 30)) +
  labs(x = "Unique fragments (log10)", y = "TSS enrichment") +
  theme_bw()
ggsave(file.path(plot_dir, "QC_TSS_vs_frags_per_sample.pdf"), qc_plot, width = 12, height = 9)

#Filter out cells

keep <- proj$cellNames[
  proj$nFrags >= MIN_FRAGS &
  proj$TSSEnrichment >= MIN_TSS &
  !(proj$Sample %in% DROP_SAMPLES)
]
proj <- proj[keep, ]

proj@sampleColData <- proj@sampleColData[unique(proj$Sample), , drop = FALSE]
 
cat("\nCells per sample after fragment/TSS filtering:\n")
print(table(proj$Sample))


#3 DOUBLET REMOVAL

proj <- addDoubletScores(
  input = proj,
  k = 10, #nearest neighbor number
  knnMethod = "UMAP",
  LSIMethod = 1
)

n_before <- table(proj$Sample)
proj <- filterDoublets(proj, filterRatio = FILTER_RATIO)
n_after <- table(proj$Sample)[names(n_before)]

cat("\nDoublet filtering per sample:\n")
print(data.frame(
  before      = as.vector(n_before),
  after       = as.vector(n_after),
  removed_pct = round(100 * (1 - as.vector(n_after) / as.vector(n_before)), 1),
  row.names   = names(n_before)
))

#SAVING

write.csv(as.data.frame(getCellColData(proj)),
          file.path(plot_dir, "cells_passing_QC.csv"))
 
proj <- saveArchRProject(
  ArchRProj       = proj,
  outputDirectory = out_dir,
  dropCells       = TRUE,   # shrink arrow copies to the kept cells
  load            = TRUE
)

#Seperate script here

-------------------------------------------------------------------------------------------------

#4 DIMENSIONALITY REDUCTION WITH LATENT SEMANTIC INDEXING (LSI) 

proj <- addIterativeLSI(
  ArchRproj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30
)


#5 BATCH CORRECTION (May not be used)

proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)


#6 CLUSTERING

proj <- addClusters(
  input = proj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

head(proj$Clusters) #shows all cluster names
table(proj$Clusters) #table with number of cells in each cluster 
cM <- confusionMatrix(paste0(proj$Clusters), paste0(proj$Sample))
cM #checking cluster composition by cell typw

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "Harmony",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP") #Cluster UMAP
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "Sample", embedding = "UMAP") #Sample Distribution UMAP (clusters should be well mixed within sample)

p1
p2

#Saving UMAPs
plotPDF(p1, name = "Clusters-UMAP.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
plotPDF(p2, name = "Sample-UMAP.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj, outputDirectory = "Save-Proj", load = FALSE)

#7 GENE MARKER UMAPS

markergenes <- list(                                   
  "monocyte general": c("PTPRC", "CD33", "HLA-DRA"),
  "classical monocyte": c("CD14"),
  "intermediate monocyte": c("CD14", "FCGR3A"),
  "nonclassical monocyte": c("CD14", "FCGR3A", "CD300E")
)
  

markergenes_vector <- list(
  "PTPRC", "CD33", "HLA-DRA", #general
  "CD14",                     #classical
  "CD14", "FCGR3A",          #intermediate
  "CD14", "FCGR3A", "CD300E" #nonclassical
)

marker_UMAP <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneScoreMatrix",
  name = markergenes_vector,
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = getImputeWeights(proj),
  plotAs = "points",
  highlightCells = getCellNames(ArchRProj = proj)[proj$Clusters == "C1"] #highlighting C1 cluster, adjust as needed
)

#plot_track <- plotBrowserTrack(
#  ArchRProj = proj,
#  groupBy = "Clusters",
#  geneSymbol = "CD14", #Adjust as needed
#  upstream = 50000,
 # downstream = 50000
#)


#8 ANNOTATING CLUSTERS

markers_annotate <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C1  # check top marker genes for cluster 1, etc.

#Now cross-reference with gene marker UMAPs to determine which cluster is what

#Cluster Labelling
proj$cell_type_lvl1 <- mapLabels(
  proj$Clusters,
  newLabels = c(
    "C1" = "iPSC Mono",
    "C2" = "Classical Mono",
    "C3" = "Intermediate Mono",
    "C4" = "Nonclassical Mono"
  ),
  oldLabels = c("C1", "C2", "C3", "C4")
)

#Visualizing Labels

plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "cell_type_lvl1",
  embedding = "UMAP",
)

#9 CHECKING FOR PLURIPOTENCY-ASSOCIATED LOCI (May add proliferation-associated)

pluripotencyGenes <- c("POU5F1", "NANOG", "SOX2", "KLF4", "MYC", "DNMT3B", "EPCAM") #SOURCE

pluripotency_plot <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneScoreMatrix",
  name = pluripotencyGenes,
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = getImputeWeights(proj)
)

pluripotency_trackplot <- plotBrowserTrack(
  ArchRProj = proj,
  groupBy = "cell_type_lvl1",
  geneSymbol = c("POU5F1", "NANOG", "SOX2"),
  upstream = 50000,
  downstream = 50000
)

grid::grid.draw(pluripotency_trackplot$POU5F1) #Adjust as needed


#10 MAC2 CALLING PEAKS

# Install MACS2 externally by running either of the code below on a terminal or WSL
# pip install MACS2
# (or via conda: conda install -c bioconda macs2)

# Find MACS2 to see if installation worked
# pathToMacs2 <- findMacs2()

#Create pseudo-bulk replicates per cell type

proj <- addGroupCoverage(
  ArchRProj = proj,
  groupBy = "cell_type_lvl1"
)

proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "cell_type_lvl1",
  pathToMacs2 = "/home/projects/icell/sumseq_pilot/envs/macs2/bin/macs2"
)

getPeakSet(proj)

proj <- addPeakMatrix(proj)

#Save Again
saveArchRProject(ArchRProj = proj, outputDirectory = "Save-ProjPeaks", load=FALSE)

peaks_heatmap <- plotMarkerHeatMap(               #THIS IS IMPORTANT!!!!!
  seMarker = markers_annotate,
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  transpose = TRUE
)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")


#11 MOTIF ENRICHMENT AND CHROME VAR (identify which TF binding motifs are enriched in accessible domain regions of each cell type)

proj <- addMotifAnnotations(
  ArchRProj = proj,
  motifSet = "cisbp",
  namae = "Motif"
) 

proj <- addBgdPeaks(proj) #adding background peaks for later chromVAR analysis

proj <- addDeviationsMatrix(  #This compares per cell how much/less accessible each TF motifs binding sites are compared to background (chromVAR deviation scores)
  ArchRProj = proj,
  peakAnnotation = "Motif",
  force = TRUE
) 

motifsUp <- peakAnnoEnrichment(
  seMarker = markers_annotate,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

#Ranked Motif Enrichment Plot

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE)]
df$rank <- seq_len(nrow(df))

ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj > 2)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(data = df[1:10,], aes(label = TF), size = 2) +
  theme_minimal()

#Visualize Motif Deviation Scores (UMAP)
 
motifs <- c("topmotifsfromdata") #Adjust as needed, based on the top TF enriched sites from the enrichment plot
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse = "|"), useMatrix = "MotifMatrix")

p <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "MotifMatrix",
  name = markerMotifs,
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)

#12 PEAK TO GENE (for Scenic+)

proj <- addPeak2GeneLinks(
  ArchRProj = proj,
  reducedDims = "Harmony"
)

p2g <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0.45,
  resolution = 1000,
  returnLoops = FALSE
)

p <- plotBrowserTrack(  #shows the peak-gene links visually
  ArchRProj = proj,
  groupBy = "cell_type_lvl1",
  geneSymbol = "CD14", #Adjust gene of interest
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(proj)
)

grid:grid.draw(p$CD14) #Adjust gene of interest 

#13 HEATMAP OF PEAK TO GENE LINKS

p2g_heatmap <- plotPeak2GeneHeatmap( #shows peak-gene pair clusters for each cell type
  ArchRProj = proj,
  groupBy = "cell_type_lvl1"
)

saveArchRProject(ArchRProj = proj, outputDirectory = "Save-ProjFinal", load = FALSE)








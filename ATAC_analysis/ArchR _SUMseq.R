

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
install.packages("renv") #installing renv enviorment
renv::init() #initialize renv enviorment

#Installing & Initializing ArchR
devtools::install_github("GreenleafLab/ArchR", ref="dev", repos = BiocManager::repositories())
library(ArchR)
ArchR::installExtraPackages()
set.seed(1)

#1 DATA LOADING

input_files <- list.files(path = "Path to arrow files", 
                          pattern = ".arrow$",
                          full.names = TRUE)

input_files #check if they have been loaded properly 

proj <- ArchRProject(
  ArrowFiles = input_files,
  outputDirectory = "ArchROutput",
  copyArrows = TRUE
)

proj #check project creation worked, also gives you # of sample cells

addArchRThreads(threads=1)
addArchRGenome("hg38")

#2 QUALITY CONTROL

#Quality Control Per Cell-likely already done by snakemake
#1. Number of unique nuclear fragments (not MT DNA), cells with too little fragments are not useful
#2. Signal to background ratio, low ratio is often dead/dying cells
#3. Fragment size distribution. Expect depletion of fragments that are length of DNA wrapped around a nucleosome (~147bp)

plotFragmentSizes(ArchRproj = proj) #fragment size distribution plot, to check data quality. Reference included

proj <- proj[proj$Frags > 1000, ] #set to 1000 unique fragments per cell. Adjust as needed

proj <- proj[proj$TSSEnrichment > 4, ] #Adjust as needed, based on next plot

df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))
df

qc_plot <- ggPoint( #combined QC plot
  x=df[,1],
  y=df[,2],
  colorDensity = TRUE,
  xlabel = "Unique Fragments Log10",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + 
  geom_vline(xintercept = 3, lty = "dashed")

qc_plot


#3 DOUBLET REMOVAL

proj <- addDoubletScores(
  input = proj,
  k = 10, #nearest neighbor number
  knnMethod = "UMAP",
  LSIMethod = 1
)

proj <- filterDoublets(proj, filterRatio = 1.5) #Adjust ratio as needed


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
  pathToMacs2 = pathToMacs2
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








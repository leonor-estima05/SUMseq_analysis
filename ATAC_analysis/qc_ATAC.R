
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

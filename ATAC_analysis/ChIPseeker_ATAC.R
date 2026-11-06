#I will be using this dataset's peaks as the background for GO pathway enrichment comparison 



#Need to convert current data to a GRange object (which ChIP-seq expects) first

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GenomicRanges)

#Data loading + genome (hg38)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
dir <- "/home/projects/icell/sumseq_pilot/ATAC/snapatac2_out/peaks_motifs/diff_motifs"

#Convering the current text files into the needed Grange object

read_peaks <- function(file) {
  p <- readLines(file)
  p <- p[p !=""]
  GRanges(sub(":.*", "", p),
          IRanges(as.numeric(sub(".*:(\\d+)-.*", "\\1", p)),
                  as.numeric(sub(".*-", "", p))))
  }

#my background for peak reference
bg_gr <- readPeakFile("/home/projects/icell/sumseq_pilot/ATAC/snapatac2_out/peaks_motifs/merged_peaks.bed")
bg_anno <- annotatePeak(bg_gr, TxDb = txdb, tssRegion = c(-3000, 3000))
universe <- unique(as.data.frame(bg_anno)$geneId)


#starting with compring iMono and hiPSCs
gr <- read_peaks(file.path(dir, "peaks_iMono_up_vs_hiPSC.txt"))

#Now to annotate each peak with its nearest gene and feature
anno <- annotatePeak(gr, TxDb = txdb, annoDb = "org.Hs.eg.db", tssRegion = c(-3000, 3000))
plotAnnoPie(anno)
write.csv(as.data.frame(anno), file.path(dir, "annotated_iMono_up_vs_hiPSC.csv"), row.names = FALSE)

#Pathway Enrichment
genes <- unique(as.data.frame(anno)$geneId)
ego <- enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", readable = TRUE, universe = universe)
dotplot(ego, showCategory = 15)

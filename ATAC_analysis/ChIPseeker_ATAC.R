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

#starting with compring iMono and hiPSCs
gr <- read_peaks(file.path(dir, "peaks_iMono_up_vs_hiPSC.txt"))

#Now to annotate each peak with its nearest gene and feature
anno <- annotatePeak(gr, TxDb = txdb, annoDB = "org.Hs.eg.db", tssRegion = c(-3000, 3000))
plotAnnoPie(anno)
write.csv(as.data.frame(anno), file.path(dir, "annotated_iMono_up_vs_hiPSC.csv"), row.names = FALSE)

#Pathway Enrichment
genes <- unique(as.data.frame(anna)$geneId)
ego <- enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", readable = TRUE)
dotplot(ego, showCategory = 15)

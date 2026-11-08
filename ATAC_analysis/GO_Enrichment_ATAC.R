library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

dir <- "/home/projects/icell/sumseq_pilot/ATAC/snapatac2_out/peaks_motifs/diff_motifs/annotating_peaks/" 
out <- file.path(dir, "GO")
dir.create(out, showWarnings = FALSE)

read_genes <- function(f) {
  g <- readLines(f)
  g[g !=""]
  }

universe <- read_genes(file.path(dir, "universe_genes.txt"))

files <- list.files(dir, pattern = "^genes_.*\\.txt$", full.names = TRUE)

for (f in files) {
  key <- sub("^genes_(.*)\\.txt$", "\\1", basename(f))
  genes <- read_genes(f)
  cat(key, ":", length(genes), "genes\n")

  ego <- enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  universe = universe, pvalueCutoff = 0.05, qvalueCutoff = 0.2) 

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat("  no significant GO terms\n")
    next
  }

  write.csv(as.data.frame(ego), file.path(out, paste0("GO_", key, ".csv")), row.names = FALSE)
  p <- dotplot(ego, showCategory = 15) + ggtitle(key)
  ggsave(file.path(out, paste0("GO_", key, ".png")), p, width = 7, height = 6, dpi = 300)
}


gene_lists <- lapply(files, read_genes)
names(gene_lists) <- sub("^genes_(.*)\\.txt$", "\\1", basename(files))

cc <- compareCluster(gene_lists, fun = "enrichGO", OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL", ont = "BP", universe = universe)

p <- dotplot(cc, showCategory = 5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out, "GO_all_sets_compared.png"), p, width = 10, height = 9, dpi = 300)

cat("Done\n")

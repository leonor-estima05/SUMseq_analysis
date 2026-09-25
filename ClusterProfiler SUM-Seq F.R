BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#https://www.youtube.com/watch?v=KimtzOf3Qz8


#Opening files from DESeq2
res_iPSC_vs_classical <- readRDS("res_iPSC_vs_classical.rds")
res_iPSC_vs_intermediate <- readRDS("res_iPSC_vs_intermediate.rds")
res_iPSC_vs_nonclassical <- readRDS("res_iPSC_vs_nonclassical.rds")



#iPSC v Classical 

#Filter Significant Genes
sig_genes_classical <- res[which(res$padj < 0.05)]
sig_genes_classical_ids <- rownames(sig_genes_classical)
length(sig_genes_classical_ids)

#Run GO Analysis
go_results_classical <- enrichGO(gene=sig_genes_classical_ids,
                                 OrgDb=org.Hs.eg.db,
                                 keyType="ENSEMBL",
                                 pAdjustMethod="BH",
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.05)
#Results Data Frame
head(as.data.frame(go_results_classical))

#Dotplot Visualization
dotplot(go_results_classical, showCategory = 15)

#Barplot (Gene Count)
barplot(go_results_classical, showCategory = 15)

#Gene Concept Network (Which genes drive which GOs)
cnetplot(go_results_classical, showCategory = 5)

#Enrichment Map (Overlap between GO Terms)
go_results_classical_2 <- pairwise_termsim(go_results_classical)
emapplot(go_results_classical_2, showCategory = 15)





#iPSC v Intermediate

sig_genes_intermediate <- res[which(res$padj < 0.05)]
sig_genes_intermediate_ids <- rownames(sig_genes_intermediate)
length(sig_genes_intermediate_ids)

go_results_intermediate <- enrichGO(gene=sig_genes_intermediate_ids,
                                 OrgDb=org.Hs.eg.db,
                                 keyType="ENSEMBL",
                                 pAdjustMethod="BH",
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.05)

head(as.data.frame(go_results_intermediate))

dotplot(go_results_intermediate, showCategory = 15)


barplot(go_results_intermediate, showCategory = 15)


cnetplot(go_results_intermediate, showCategory = 5)


go_results_intermediate_2 <- pairwise_termsim(go_results_intermediate)
emapplot(go_results_intermediate_2, showCategory = 15)





#iPSC v Nonclassical

sig_genes_nonclassical <- res[which(res$padj < 0.05)]
sig_genes_nonclassical_ids <- rownames(sig_genes_nonclassical)
length(sig_genes_nonclassical_ids)

go_results_nonclassical <- enrichGO(gene=sig_genes_nonclassical_ids,
                                 OrgDb=org.Hs.eg.db,
                                 keyType="ENSEMBL",
                                 pAdjustMethod="BH",
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.05)

head(as.data.frame(go_results_nonclassical))

dotplot(go_results_nonclassical, showCategory = 15)

barplot(go_results_nonclassical, showCategory = 15)


cnetplot(go_results_nonclassical, showCategory = 5)


go_results_nonclassical_2 <- pairwise_termsim(go_results_nonclassical)
emapplot(go_results_nonclassical_2, showCategory = 15)



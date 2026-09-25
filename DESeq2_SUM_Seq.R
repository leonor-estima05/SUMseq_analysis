

#Load Libraries
library(DESeq2)
library(tidyverse)
library(ComplexHeatMaps)
library(circlize)
library(RColorBrewer)

countData <- read.csv("counts_data.csv", row.names=1)
metaData <- read.csv("sample_info.csv", row.names=1)
#Genes must be rownames not columns!! if they are you dont need the countData$ensgene

#1. Correctly Format the Data

#Checking for correct structure of both
head(countData) 
head(metaData)
dim(countData)
dim(metaData)

#Setting row names and cleaning up countdata, only need this if genes arent rows and are columns!!
#rownames(countData) <- countData$ensgene
#countData$ensgene <- NULL

#Checking it worked
head(countData)
dim(countData)

#2 Create dds Object

dds <- DESeqDataSetFromMatrix(countData=countData,
                              colData=metaData,
                              design = ~cell_type_lvl1)

dds


keep <- rowSums(counts(dds)) >= 10 #keeping rows that have at least 10 reads
dds <- dds[keep,]

dds

#set factor (which condition should DESeq compare against)
dds$cell_type_lvl1 <- relevel(dds$cell_type_lvl1, ref = "classical mono")
levels(dds$cell_type_lvl1) #Checking if levels are correct for iPSC_mono, classical, etc.


#3 Run DESeq
dds <- DESeq(dds)

#Comparison Information
resultsNames(dds)

#4  contrasts
res_iPSC_vs_classical <- results(dds, contrast = c("cell_type_lvl1", "iPSC mono", "classical mono"))
res_iPSC_vs_intermediate <- results(dds, contrast = c("cell_type_lvl1", "iPSC mono", "intermediate mono"))
res_iPSC_vs_nonclassical <- results(dds, contrast = c("cell_type_lvl1", "iPSC mono", "nonclassical mono"))

summary(res_iPSC_vs_classical)
summary(res_iPSC_vs_intermediate)
summary(res_iPSC_vs_nonclassical)

#MA Plot for one comparison
plotMA(res_iPSC_vs_classical)
plotMA(res_iPSC_vs_intermediate)
plotMA(res_iPSC_vs_nonclassical)

#PCA Plot - shows ALL groups at once, only needs to be done once
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "cell_type_lvl1")

#Volcano Plot - for one comparison
res_df <- as.data.frame(res_iPSC_vs_classical)
res_df$significant <- res_df$padj < 0.05
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point() +
  theme_minimal()

res_df2 <- as.data.frame(res_iPSC_vs_intermediate)
res_df2$significant <- res_df2$padj < 0.05
ggplot(res_df2, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point() +
  theme_minimal()

res_df3 <- as.data.frame(res_iPSC_vs_nonclassical)
res_df3$significant <- res_df3$padj < 0.05
ggplot(res_df3, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point() +
  theme_minimal()




#Heat Map iPSC v Classical

#Getting top 30 expressed genes
top_genes_classical <- head(order(res_iPSC_vs_classical$padj, na.last = TRUE), 30)
heatmapData_classical <- assay(vsd)[top_genes_classical, ]
heatmapData_classical_scaled <- t(scale(t(heatmapData_classical))) #scaling per gene (row) to normalize gene expression across sample


#Creating HeatMap 
max(heatmapData_classical_scaled)
min(heatmapData_classical_scaled) #checking min and max values to then choose color range (-6,0,6)

col_fun1 <- colorRamp2(c(-6,0,6), c("blue", "lightgray", "red")) #color palette options
col_fun2 <- RColorBrewer::brewer.pal(name = "Spectral", n = 11)

annotation <- HeatmapAnnotation(
  cell_type = metaData$cell_type_lvl1,
  col = list(cell_type = c(
    "iPSC mono" = "red",
    "classical mono" = "lightblue",
    "intermediate mono" = "lightgreen",
    "nonclassical mono" = "purple"
  )),
  annotation_legend_param = list(
    cell_type = list(title = "Cell Type")
  )
)

gene_label <- ifelse(res_iPSC_vs_classical$log2FoldChange[top_genes_classical] > 0, 
                            "up in iPSC", "down in iPSC")

row_annotation <- rowAnnotation(
  direction = gene_label,
  col = list(direction = c("up in iPSC" = "red", "down in iPSC" = "lightblue"))
)

custom_color_heatmap_classical <- Heatmap(
  heatmapData_classical_scaled,
  name = "Expression",
  col = col_fun2, #color palette 
  top_annotation = annotation,
  right_annotation = row_annotation, #row annotation
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_split = metaData$cell_type_lvl1,
  column_title = "iPSC vs Classical Monocytes Top 30 DEGs",
  heatmap_legend_param = list(
    title = "Scaled Expression",
    legend_direction = "horizontal",
    legend_width = unit(6, "cm")
  )
)

draw(custom_color_heatmap_classical, heatmap_legend_side = "bottom")



#HeatMap iPSC v Intermediate

#Getting top 30 expressed genes
top_genes_intermediate <- head(order(res_iPSC_vs_intermediate$padj, na.last = TRUE), 30)
heatmapData_intermediate <- assay(vsd)[top_genes_intermediate, ]
heatmapData_intermediate_scaled <- t(scale(t(heatmapData_intermediate))) #scaling per gene (row) to normalize gene expression across sample


#Creating HeatMap 


custom_color_heatmap_intermediate <- Heatmap(
  heatmapData_intermediate_scaled,
  name = "Expression",
  col = col_fun2, #color palette 
  top_annotation = annotation,
  right_annotation = row_annotation, #row annotation
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_split = metaData$cell_type_lvl1,
  column_title = "iPSC vs Intermediate Monocytes Top 30 DEGs",
  heatmap_legend_param = list(
    title = "Scaled Expression",
    legend_direction = "horizontal",
    legend_width = unit(6, "cm")
  )
)

draw(custom_color_heatmap_intermediate, heatmap_legend_side = "bottom")






#HeatMap iPSC v Nonclassical


#Getting top 30 expressed genes
top_genes_nonclassical <- head(order(res_iPSC_vs_nonclassical$padj, na.last = TRUE), 30)
heatmapData_nonclassical <- assay(vsd)[top_genes_nonclassical, ]
heatmapData_nonclassical_scaled <- t(scale(t(heatmapData_nonclassical))) #scaling per gene (row) to normalize gene expression across sample


#Creating HeatMap 


custom_color_heatmap_nonclassical <- Heatmap(
  heatmapData_nonclassical_scaled,
  name = "Expression",
  col = col_fun2, #color palette 
  top_annotation = annotation,
  right_annotation = row_annotation, #row annotation
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_split = metaData$cell_type_lvl1,
  column_title = "iPSC vs Non-Classical Monocytes Top 30 DEGs",
  heatmap_legend_param = list(
    title = "Scaled Expression",
    legend_direction = "horizontal",
    legend_width = unit(6, "cm")
  )
)

draw(custom_color_heatmap_nonclassical, heatmap_legend_side = "bottom")

  



#Saving for ClusterProfiler
saveRDS(res_iPSC_vs_classical, "res_iPSC_vs_classical.rds")
saveRDS(res_iPSC_vs_intermediate, "res_iPSC_vs_intermediate.rds")
saveRDS(res_iPSC_vs_nonclassical, "res_iPSC_vs_nonclassical.rds")

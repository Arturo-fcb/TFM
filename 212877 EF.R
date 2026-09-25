# El enriquecimiento funcional se calculó sobre la versión sin voom de las listas de DEGs, tal como se indica explícitamente en el apartado 6.2.3 del TFM.
# ==================================================
# 0. LIBRERÍAS
# ==================================================
library(edgeR)
library(limma)
library(pheatmap)
library(ggplot2)
 
# Enriquecimiento
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
 
# ==================================================
# 1. RUTA DE TRABAJO
# ==================================================
setwd("C:/ARTURO/Arturo universidad/Master/UE/TFM/dATOS/Nueva carpeta/212877")
 
# ==================================================
# 2. LEER DATOS 
# ==================================================
files <- list.files(pattern = "GSM.*_C[1-3]_RNASeq.txt.gz$", full.names = TRUE)
 
data_list <- lapply(files, function(f) {
  read.delim(f, header = TRUE, comment.char = "#", stringsAsFactors = FALSE, fill = TRUE)
})
 
genes <- gsub("\\s.*", "", data_list[[1]][, 1])
 
counts <- do.call(cbind, lapply(data_list, function(x) {
  tmp <- x[, ncol(x)]
  as.numeric(gsub(".*\\s", "", tmp))
}))
 
rownames(counts) <- genes
colnames(counts) <- gsub(".txt.gz", "", basename(files))
 
# ==================================================
# 3. METADATA
# ==================================================
cell_line <- sub(".*_(.*)_C[0-9]_RNASeq.*", "\\1", colnames(counts))
cell_line <- factor(cell_line)
 
# ==================================================
# 4. NORMALIZACIÓN
# ==================================================
dge <- DGEList(counts = counts, group = cell_line)
keep <- filterByExpr(dge, group = cell_line)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE)
 
# ==================================================
# 5. PCA
# ==================================================
pca <- prcomp(t(logCPM), scale. = TRUE)
 
pca_df <- data.frame(PC1 = pca$x[,1],
                     PC2 = pca$x[,2],
                     cell_line = cell_line)
 
png("PCA.png", 1200, 800)
ggplot(pca_df, aes(PC1, PC2, color = cell_line)) +
  geom_point(size = 8) + 
  theme_minimal() +
  theme(
    legend.text = element_text(size = 20),   
    legend.title = element_text(size = 22),  
    axis.title = element_text(size = 18)     
  ) +
  guides(color = guide_legend(override.aes = list(size = 8))) 
dev.off()
 
 
 
# ==================================================
# 6. DEFINICIÓN ADRN vs MES
# ==================================================
ADRN_lines <- c("BE2", "SY5Y", "55N", "WN")
MES_lines  <- c("SHEP")
 
phenotype <- ifelse(cell_line %in% ADRN_lines, "ADRN",
                    ifelse(cell_line %in% MES_lines, "MES", "OTHER"))
phenotype <- factor(phenotype)
 
sel <- phenotype %in% c("ADRN", "MES")
 
logCPM_sub <- logCPM[, sel]
phenotype_sub <- droplevels(phenotype[sel])
 
# ==================================================
# 7. EXPRESIÓN DIFERENCIAL (limma)
# ==================================================
design <- model.matrix(~ phenotype_sub)
 
fit <- lmFit(logCPM_sub, design)
fit <- eBayes(fit)
 
results <- topTable(fit, coef = "phenotype_subMES", number = Inf)
write.csv(results, "DEGs_ADRN_vs_MES.csv")
 
# ==================================================
# VOLCANO MEJORADO ADRN vs MES
# ==================================================
 
results$threshold <- "NS"
results$threshold[results$logFC > 1 & results$adj.P.Val < 0.05] <- "Up in MES"
results$threshold[results$logFC < -1 & results$adj.P.Val < 0.05] <- "Up in ADRN"
 
results$threshold <- factor(results$threshold, 
                            levels = c("Up in ADRN", "NS", "Up in MES"))
 
png("Volcano_ADRN_vs_MES.png", 1400, 900)
 
ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  
  geom_point(alpha = 0.7, size = 2) +
  
  # Líneas de corte
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  # Colores claros
  scale_color_manual(values = c(
    "Up in ADRN" = "#2C7BB6",   # azul
    "NS" = "grey70",
    "Up in MES" = "#D7191C"     # rojo
  )) +
  
  theme_minimal() +
  
  labs(
    title = " MES vs ADRN",
    subtitle = "Genes sobreexpresados en MES (derecha) vs ADRN (izquierda)",
    x = "log2 Fold Change (MES / ADRN)",
    y = "-log10 adjusted p-value",
    color = "Categoria"
  ) +
  
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
  )
 
dev.off()
 
# ==================================================
# 9. GENES SIGNIFICATIVOS
# ==================================================
sig <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
 
genes_up_MES  <- rownames(sig[sig$logFC > 0, ])
genes_up_ADRN <- rownames(sig[sig$logFC < 0, ])
 
writeLines(genes_up_MES, "genes_up_MES.txt")
writeLines(genes_up_ADRN, "genes_up_ADRN.txt")
 
 
#####
top20_MES <- head(sig[order(-sig$logFC), ], 20)
 
top20_ADRN <- head(sig[order(sig$logFC), ], 20)
###
df_MES <- top20_MES
df_MES$gene <- rownames(df_MES)
 
ggplot(df_MES, aes(x = reorder(gene, logFC), y = logFC)) +
  geom_col(fill = "#D7191C") +
  coord_flip() +
  labs(title = "Top 20 genes up en MES",
       x = "Gen",
       y = "log2 Fold Change")
 
 
#
df_ADRN <- top20_ADRN
df_ADRN$gene <- rownames(df_ADRN)
 
ggplot(df_ADRN, aes(x = reorder(gene, logFC), y = logFC)) +
  geom_col(fill = "#2C7BB6") +
  coord_flip() +
  labs(title = "Top 20 genes up en ADRN",
       x = "Gen",
       y = "log2 Fold Change")
 
 
# ==================================================
# 10. CONVERSIÓN A ENTREZ
# ==================================================
convert_to_entrez <- function(symbols) {
  res <- bitr(symbols, fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)
  return(res$ENTREZID)
}
 
entrez_MES  <- convert_to_entrez(genes_up_MES)
entrez_ADRN <- convert_to_entrez(genes_up_ADRN)
 
# ==================================================
# 11. GO ENRICHMENT
# ==================================================
ego_MES <- enrichGO(entrez_MES, OrgDb=org.Hs.eg.db,
                    ont="BP", pAdjustMethod="BH")
 
ego_ADRN <- enrichGO(entrez_ADRN, OrgDb=org.Hs.eg.db,
                     ont="BP", pAdjustMethod="BH")
 
write.csv(as.data.frame(ego_MES), "GO_MES.csv")
write.csv(as.data.frame(ego_ADRN), "GO_ADRN.csv")
 
# ==================================================
# 12. KEGG
# ==================================================
kegg_MES <- enrichKEGG(entrez_MES, organism="hsa")
kegg_ADRN <- enrichKEGG(entrez_ADRN, organism="hsa")
 
write.csv(as.data.frame(kegg_MES), "KEGG_MES.csv")
write.csv(as.data.frame(kegg_ADRN), "KEGG_ADRN.csv")
 
# ==================================================
# 13. GSEA
# ==================================================
gene_list <- results$logFC
names(gene_list) <- rownames(results)
gene_list <- sort(gene_list, decreasing = TRUE)
 
gsea <- gseGO(geneList = gene_list,
              OrgDb = org.Hs.eg.db,
              ont = "BP",
              keyType = "SYMBOL")
 
write.csv(as.data.frame(gsea), "GSEA_GO.csv")
 
# ==================================================
# FIN
# ==================================================
cat("Pipeline completo ejecutado\n")
###### El codigo es elaboracion propia pero estructurado con IA

5. GSE62564 — Carga de datos, expresión diferencial HR vs. non-HR y enriquecimiento

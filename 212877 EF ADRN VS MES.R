# Versión final reportada en el TFM (incluye el ajuste con voom() sobre el DGEList filtrado, que sustituye a una primera versión sobre logCPM directamente).
# =============================
# 0. LIBRERÍAS
# =============================
library(edgeR)
library(limma)
library(pheatmap)
library(ggplot2)
 
# =============================
# 1. RUTA DE TRABAJO
# =============================
setwd("C:/ARTURO/Arturo universidad/Master/UE/TFM/dATOS/Nueva carpeta/212877")
 
# =============================
# 2. LEER ARCHIVOS RNA-SEQ (SOLO CONTROLES)
# =============================
files <- list.files(pattern = "GSM.*_C[1-3]_RNASeq.txt.gz$", full.names = TRUE)
length(files)   # 24
 
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
 
dim(counts)
head(counts[, 1:4])
 
# =============================
# 3. METADATA: LÍNEA CELULAR
# =============================
cell_line <- sub(".*_(.*)_C[0-9]_RNASeq.*", "\\1", colnames(counts))
cell_line <- factor(cell_line)
table(cell_line)
 
# =============================
# 4. OBJETO edgeR + FILTRADO + NORMALIZACIÓN
# =============================
dge <- DGEList(counts = counts, group = cell_line)
keep <- filterByExpr(dge, group = cell_line)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE)
 
# =============================
# 5. PCA (TODAS LAS LÍNEAS CELULARES)
# =============================
pca <- prcomp(t(logCPM), scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[, 1],
                     PC2 = pca$x[, 2],
                     cell_line = cell_line)
 
png("PCA_8_cell_lines.png", 1200, 800)
ggplot(pca_df, aes(PC1, PC2, color = cell_line)) +
  geom_point(size = 6) +
  theme_minimal() +
  ggtitle("PCA - 8 líneas celulares (controles)") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 22))
dev.off()
 
# =============================
# 6. MDS PLOT
# =============================
png("MDS_8_cell_lines.png", 1200, 800)
plotMDS(dge, labels = cell_line, col = as.numeric(cell_line))
title("MDS - 8 líneas celulares (controles)")
dev.off()
 
# =============================
# 7. HEATMAP DE CORRELACIÓN ENTRE MUESTRAS
# =============================
cor_mat <- cor(logCPM)
annotation <- data.frame(cell_line = cell_line)
rownames(annotation) <- colnames(logCPM)
 
png("Correlation_heatmap_8_cell_lines.png", 1200, 1000)
pheatmap(cor_mat,
         main = "Correlación entre muestras (logCPM)",
         annotation_col = annotation)
dev.off()
 
# =============================
# 8. CLUSTERING JERÁRQUICO
# =============================
hc <- hclust(dist(t(logCPM)))
png("Dendrogram_8_cell_lines.png", 1200, 800)
plot(hc, main = "Clustering jerárquico - 8 líneas celulares")
dev.off()
 
# =============================
# 9. HEATMAP TOP 500 GENES MÁS VARIABLES
# =============================
var_genes <- apply(logCPM, 1, var)
top_var <- names(sort(var_genes, decreasing = TRUE))[1:500]
 
annotation <- data.frame(cell_line = cell_line)
rownames(annotation) <- colnames(logCPM)
 
png("Heatmap_top500_var_8_cell_lines.png", 1200, 1000)
pheatmap(logCPM[top_var, ],
         scale = "row",
         annotation_col = annotation,
         main = "Top 500 genes más variables")
dev.off()
 
# =============================
# 10. DISEÑO PARA COMPARACIÓN ADRN vs MES 
# =============================
# ADRN: BE2, SY5Y, 55N
# MES:  SHEP 
ADRN_lines <- c("BE2", "SY5Y", "55N", "WN")
 
MES_lines  <- c("SHEP")
 
phenotype <- ifelse(cell_line %in% ADRN_lines, "ADRN",
                    ifelse(cell_line %in% MES_lines, "MES", "OTHER"))
phenotype <- factor(phenotype, levels = c("ADRN", "MES", "OTHER"))
table(phenotype)
 
# Seleccionar solo las muestras ADRN y MES
sel <- phenotype %in% c("ADRN", "MES")
logCPM_sub <- logCPM[, sel]
phenotype_sub <- droplevels(phenotype[sel])
 
# Comprobar cuántas muestras de cada grupo
table(phenotype_sub)
 
# =============================
# 11. LIMMA PARA ADRN vs MES
# =============================
design <- model.matrix(~ phenotype_sub)   # Intercept + coeficiente MES vs ADRN
fit <- lmFit(logCPM_sub, design)
fit <- eBayes(fit)
 
results_ADRN_MES <- topTable(fit, coef = "phenotype_subMES", number = Inf)
write.csv(results_ADRN_MES, "DEGs_ADRN_vs_MES.csv")
 
# =============================
# 12. VOLCANO PLOT
# =============================
results_ADRN_MES$threshold <- as.factor(abs(results_ADRN_MES$logFC) > 1 &
                                          results_ADRN_MES$adj.P.Val < 0.05)
 
png("Volcano_ADRN_vs_MES.png", 1200, 800)
ggplot(results_ADRN_MES, aes(logFC, -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  ggtitle("Volcano - MES vs ADRN") +
  xlab("log2 Fold Change (MES/ADRN)") +
  ylab("-log10(adj.P.Value)")
dev.off()
 
# =============================
# 13. HEATMAP TOP 50 DEGs
# =============================
top_genes_ADRN_MES <- rownames(results_ADRN_MES)[1:50]
mat_ADRN_MES <- logCPM_sub[top_genes_ADRN_MES, ]
 
annotation <- data.frame(phenotype = phenotype_sub)
rownames(annotation) <- colnames(mat_ADRN_MES)
 
png("Heatmap_top50_ADRN_vs_MES.png", 1000, 1000)
pheatmap(mat_ADRN_MES,
         scale = "row",
         annotation_col = annotation,
         main = "Top 50 DEGs - MES vs ADRN")
dev.off()
 
# =============================
# 14. EXTRACCIÓN DE LISTAS DE GENES SIGNIFICATIVOS
# =============================
sig <- results_ADRN_MES[results_ADRN_MES$adj.P.Val < 0.05 & abs(results_ADRN_MES$logFC) > 1, ]
genes_up_MES  <- rownames(sig[sig$logFC > 0, ])
genes_up_ADRN <- rownames(sig[sig$logFC < 0, ])
 
cat("Número de genes sobreexpresados en MES:", length(genes_up_MES), "\n")
cat("Número de genes sobreexpresados en ADRN:", length(genes_up_ADRN), "\n")
 
# Guardar listas 
writeLines(genes_up_MES, "genes_up_MES.txt")
writeLines(genes_up_ADRN, "genes_up_ADRN.txt")
 
 
 
# =============================
# 10. DISEÑO PARA COMPARACIÓN ADRN vs MES (CON VOOM)
# =============================
# ADRN: BE2, SY5Y, 55N, WN
# MES:  SHEP
ADRN_lines <- c("BE2", "SY5Y", "55N", "WN")
MES_lines  <- c("SHEP")
 
phenotype <- ifelse(cell_line %in% ADRN_lines, "ADRN",
                    ifelse(cell_line %in% MES_lines, "MES", "OTHER"))
phenotype <- factor(phenotype, levels = c("ADRN", "MES", "OTHER"))
table(phenotype)
 
# Seleccionar solo las muestras ADRN y MES
sel <- phenotype %in% c("ADRN", "MES")
phenotype_sub <- droplevels(phenotype[sel])
table(phenotype_sub)
 
# ── CAMBIO CLAVE ─────────────────────────────────────────────────────────
# voom() trabaja sobre CONTEOS (el DGEList ya filtrado con filterByExpr),
# no sobre logCPM. Por eso subseteamos el objeto dge original, no logCPM.
dge_sub <- dge[, sel]
 
# Recalcular los factores de normalización tras el subsetting es buena
# práctica: al quitar muestras, los factores TMM originales ya no son
# óptimos para este subconjunto de 15 muestras.
dge_sub <- calcNormFactors(dge_sub)
 
design <- model.matrix(~ phenotype_sub)
 
# =============================
# 11. VOOM + LIMMA PARA ADRN vs MES
# =============================
# voom() estima la relación media-varianza propia de RNA-seq y genera un
# peso de precisión por observación (gen x muestra), que lmFit() usa para
# ponderar el ajuste. Sin esto, limma trata los logCPM como si tuvieran
# varianza homogénea, lo cual no es cierto en RNA-seq (los genes de baja
# expresión son sistemáticamente más ruidosos).
png("voom_mean_variance_ADRN_vs_MES.png", 1000, 800)
v <- voom(dge_sub, design, plot = TRUE)   # el plot muestra la curva mean-variance
dev.off()
 
fit <- lmFit(v, design)
fit <- eBayes(fit)
 
results_ADRN_MES <- topTable(fit, coef = "phenotype_subMES", number = Inf)
write.csv(results_ADRN_MES, "DEGs_ADRN_vs_MES_voom.csv")
 
# =============================
# 12. VOLCANO PLOT (igual que antes, sobre los nuevos resultados)
# =============================
results_ADRN_MES$threshold <- "NS"
results_ADRN_MES$threshold[results_ADRN_MES$logFC >  1 & results_ADRN_MES$adj.P.Val < 0.05] <- "Up in MES"
results_ADRN_MES$threshold[results_ADRN_MES$logFC < -1 & results_ADRN_MES$adj.P.Val < 0.05] <- "Up in ADRN"
results_ADRN_MES$threshold <- factor(results_ADRN_MES$threshold,
                                     levels = c("Up in ADRN", "NS", "Up in MES"))
 
png("Volcano_ADRN_vs_MES_voom.png", 1400, 900)
ggplot(results_ADRN_MES, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("Up in ADRN" = "#2C7BB6", "NS" = "grey70", "Up in MES" = "#D7191C")) +
  theme_minimal() +
  labs(title = "MES vs ADRN (voom)",
       subtitle = "Genes sobreexpresados en MES (derecha) vs ADRN (izquierda)",
       x = "log2 Fold Change (MES / ADRN)", y = "-log10 adjusted p-value", color = "Categoria")
dev.off()
 
# =============================
# 13. GENES SIGNIFICATIVOS
# =============================
sig <- results_ADRN_MES[results_ADRN_MES$adj.P.Val < 0.05 & abs(results_ADRN_MES$logFC) > 1, ]
genes_up_MES  <- rownames(sig[sig$logFC > 0, ])
genes_up_ADRN <- rownames(sig[sig$logFC < 0, ])
 
cat("Número de genes sobreexpresados en MES:", length(genes_up_MES), "\n")
cat("Número de genes sobreexpresados en ADRN:", length(genes_up_ADRN), "\n")
 
writeLines(genes_up_MES,  "genes_up_MES_voom.txt")
writeLines(genes_up_ADRN, "genes_up_ADRN_voom.txt")
 
# =============================
# 14. COMPROBACIÓN DIRECTA DE LOS 4 GENES CANDIDATOS
# =============================
# Esto es lo que más me interesa que me pegues de vuelta: si ARPP21 o DLG2
# cruzan el umbral de significación con voom, es un resultado que vale la
# pena remarcar en el TFM.
candidatos <- c("ARPP21", "DLG2", "AXL", "MYCN")
results_ADRN_MES[candidatos, c("logFC", "P.Value", "adj.P.Val")]

4. GSE212877 — Enriquecimiento funcional

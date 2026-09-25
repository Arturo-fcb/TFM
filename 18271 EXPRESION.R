 
#  GSE18271 
#  Arturo Fernández
 
 
# -----------------------------
# -----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
 
BiocManager::install(c("affy", "limma", "hgu133plus2.db", "annotate"), ask=FALSE)
 
install.packages(c("ggplot2", "pheatmap", "ggrepel"))
 
library(affy)
library(limma)
library(hgu133plus2.db)
library(annotate)
library(ggplot2)
library(pheatmap)
library(ggrepel)
 
 
# 1. CARGAR DATOS 
# -----------------------------
setwd("C:/ARTURO/Arturo universidad/Master/UE/TFM/dATOS/Nueva carpeta/18271")
 
data <- ReadAffy()
 
 
# 2. CONTROL DE CALIDAD
 
png("boxplot_raw.png", width=2000, height=1400, res=300)
boxplot(data,
        main="Distribución de intensidades (datos crudos)",
        las=2,
        col="lightblue")
dev.off()
 
 
# 3. NORMALIZACIÓN (RMA)
 
eset <- rma(data)
expr <- exprs(eset)
 
png("boxplot_normalized.png", width=2000, height=1400, res=300)
boxplot(expr,
        main="Datos normalizados (RMA)",
        las=2,
        col="lightgreen")
dev.off()
 
 
# 4. GRUPOS
 
group <- factor(c(
  "gn","gn","gnb","gn","gn","gn","gnb","gnb","gn","gnb","gnb",
  "gn","gnb","gnb","gn","gn","gn","gnb","gnb","gn","gnb","gn"
))
 
table(group)
 
 
# 5. DISEÑO 
 
design <- model.matrix(~ group)
colnames(design) <- c("Intercept","GNB_vs_GN")
 
 
# 6. EXPRESIÓN DIFERENCIAL
 
fit <- lmFit(expr, design)
fit <- eBayes(fit)
 
results <- topTable(fit,
                    coef="GNB_vs_GN",
                    number=Inf,
                    adjust.method="fdr")
 
 
# 7. ANOTACIÓN
 
results$Symbol <- getSYMBOL(rownames(results), "hgu133plus2.db")
 
write.csv(results, "DEG_results_GSE18271.csv")
 
 
# 8. VOLCANO PLOT 
 
results$threshold <- as.factor(abs(results$logFC) > 1 & results$adj.P.Val < 0.05)
 
png("volcano.png", width=2000, height=1400, res=300)
 
ggplot(results, aes(x=logFC, y=-log10(adj.P.Val), color=threshold)) +
  geom_point(alpha=0.7, size=2) +
  scale_color_manual(values=c("grey","red"),
                     labels=c("No significativo","Significativo")) +
  theme_minimal(base_size = 16) +
  labs(
    title="Volcano Plot: GN vs GNB",
    x="Log2 Fold Change",
    y="-log10(p ajustado)",
    color="Clasificación"
  ) +
  theme(plot.title = element_text(hjust=0.5, face="bold"))
 
dev.off()
 
 
# 9. PCA 
 
pca <- prcomp(t(expr), scale.=TRUE)
pca_df <- data.frame(pca$x, group=group)
 
var_explained <- summary(pca)$importance[2,]
 
png("PCA.png", width=2000, height=1400, res=300)
 
ggplot(pca_df, aes(PC1, PC2, color=group)) +
  geom_point(size=5) +
  theme_minimal(base_size = 18) +
  labs(
    title = "Análisis de Componentes Principales (PCA)",
    x = paste0("PC1 (", round(var_explained[1]*100,1), "% varianza)"),
    y = paste0("PC2 (", round(var_explained[2]*100,1), "% varianza)"),
    color = "Grupo"
  ) +
  theme(
    plot.title = element_text(hjust=0.5, face="bold"),
    legend.position = "right"
  )
 
dev.off()
 
 
# 10. HEATMAP 
 
annotation_col <- data.frame(Grupo = group)
rownames(annotation_col) <- colnames(expr)
 
top_genes <- head(rownames(results), 50)
 
png("heatmap.png", width=2000, height=1600, res=300)
 
pheatmap(expr[top_genes, ],
         scale="row",
         annotation_col = annotation_col,
         show_rownames = FALSE,
         fontsize_col = 10,
         main = "Top 50 genes diferencialmente expresados")
 
dev.off()
 
 
# 11. EXPORTAR DATOS
 
write.csv(expr, "RMA_normalized_expression.csv")
 
# TOP 20
top20 <- results[order(results$adj.P.Val), ][1:20,
                                             c("Symbol","logFC","adj.P.Val","AveExpr","B")]
 
colnames(top20) <- c("Gen","logFC","p_adj","Exp_media","B")
 
write.csv(top20, "Top20_DEGs.csv")
 
 
 

2. GSE18271 — Enriquecimiento funcional

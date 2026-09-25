library(GEOquery)
setwd("C:/ARTURO/Arturo universidad/Master/UE/TFM/dATOS/Nueva carpeta/62564")
# ── 1. Cargar metadatos clínicos desde el series matrix ──────────────────────
gse <- getGEO(filename = 'GSE62564_series_matrix.txt')
meta <- pData(gse)
 
# ── 2. Extraer variables clínicas limpias ────────────────────────────────────
# Las características están en columnas characteristics_ch1.X
# Hay que parsearlas porque cada fila mezcla variables distintas por muestra
 
# Función para extraer valor de una variable por nombre
extract_char <- function(meta, pattern) {
  chars <- meta[, grep("characteristics_ch1", colnames(meta))]
  apply(chars, 1, function(row) {
    hit <- grep(pattern, row, value = TRUE)
    if (length(hit) == 0) return(NA)
    as.numeric(gsub(paste0(".*", pattern, ": "), "", hit[1]))
  })
}
 
meta$os_days  <- extract_char(meta, "os day")
meta$os_event <- extract_char(meta, "os bin")
meta$efs_days <- extract_char(meta, "efs day")
meta$efs_event<- extract_char(meta, "efs bin")
meta$age_days <- extract_char(meta, "age")
meta$high_risk <- ifelse(grepl("HR", apply(
  meta[, grep("characteristics_ch1", colnames(meta))], 1, paste, collapse=" "
)), 1, 0)
 
# ── 3. Cargar matriz de expresión ────────────────────────────────────────────
expr <- read.delim('GSE62564_SEQC_NB_RNA-Seq_log2RPM.txt',
                   header=TRUE, row.names=1, check.names=FALSE)
 
# ── 4. Verificar alineación muestras ─────────────────────────────────────────
# Los IDs de expresión son SEQC_NB001...NB498
# Los IDs de meta son GSM...  → usar títulos
sample_titles <- gsub(" \\[2\\]", "", meta$title)  # "SEQC_NB001"
rownames(meta) <- sample_titles
 
# Intersección de muestras comunes
common <- intersect(rownames(meta), colnames(expr))
expr_filt <- expr[, common]
meta_filt  <- meta[common, ]
 
cat("Muestras con expresión Y metadatos:", length(common), "\n")
cat("OS disponible:", sum(!is.na(meta_filt$os_days)), "\n")
cat("EFS disponible:", sum(!is.na(meta_filt$efs_days)), "\n")
 
 
 
 
library(survival)
library(survminer)
library(AnnotationDbi)
library(org.Hs.eg.db)
 
# ── 1. Convertir símbolos génicos a RefSeq IDs ───────────────────────────────
genes_interes <- c("ARPP21", "DLG2", "AXL", "MYCN")
 
refseq_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = genes_interes,
  keytype  = "SYMBOL",
  columns  = c("SYMBOL", "REFSEQ")
)
 
# Ver qué RefSeq IDs corresponden a cada gen
print(refseq_map)
 
# Quedarnos solo con RefSeq de ARNm (NM_)
refseq_map <- refseq_map[grepl("^NM_", refseq_map$REFSEQ), ]
 
# Ver cuáles están presentes en nuestra matriz
refseq_map$en_expr <- refseq_map$REFSEQ %in% rownames(expr_filt)
print(refseq_map)
 
 
 
 
# ── 2. Seleccionar transcrito canónico de cada gen ───────────────────────────
# Usamos el transcrito principal (NM_ más corto / canónico) presente en la matriz
 
transcrito_seleccionado <- c(
  ARPP21 = "NM_016300",    # transcrito canónico principal
  DLG2   = "NM_001364",    # transcrito canónico principal
  AXL    = "NM_001699",    # transcrito canónico principal
  MYCN   = "NM_005378"     # único disponible
)
 
# ── 3. Extraer expresión y unir con metadatos ────────────────────────────────
surv_data <- meta_filt[, c("os_days", "os_event", "efs_days", "efs_event",
                           "age_days", "high_risk")]
 
for (gen in names(transcrito_seleccionado)) {
  refseq_id <- transcrito_seleccionado[gen]
  surv_data[[gen]] <- as.numeric(expr_filt[refseq_id, rownames(surv_data)])
}
 
# Verificar que todo está bien
cat("Dimensiones:", nrow(surv_data), "muestras x", ncol(surv_data), "variables\n")
summary(surv_data[, names(transcrito_seleccionado)])
 
# ── 4. Función para hacer KM de un gen ───────────────────────────────────────
km_gen <- function(gen, tiempo, evento, titulo) {
  
  df         <- surv_data[!is.na(surv_data[[tiempo]]) & !is.na(surv_data[[gen]]), ]
  mediana    <- median(df[[gen]], na.rm = TRUE)
  df$grupo   <- factor(ifelse(df[[gen]] >= mediana, "Alto", "Bajo"), levels = c("Bajo", "Alto"))
  df$t_years <- df[[tiempo]] / 365.25
  df$ev      <- df[[evento]]
  
  km_fit <- survfit(Surv(t_years, ev) ~ grupo, data = df)
  
  ggsurvplot(
    km_fit,
    data              = df,
    pval              = TRUE,
    risk.table        = TRUE,
    conf.int          = TRUE,
    palette           = c("#2C7BB6", "#D7191C"),
    legend.labs       = c(paste(gen, "bajo"), paste(gen, "alto")),
    xlab              = "Tiempo (años)",
    ylab              = "Probabilidad de supervivencia",
    title             = paste0(titulo, " — ", gen),
    ggtheme           = theme_minimal(base_size = 14),
    risk.table.height = 0.25
  )
}
 
km_gen("ARPP21", "os_days", "os_event", "Overall Survival")
km_gen("DLG2",   "os_days", "os_event", "Overall Survival")
km_gen("AXL",    "os_days", "os_event", "Overall Survival")
km_gen("MYCN",   "os_days", "os_event", "Overall Survival")
 
 
# Guardar todos en PNG con buena resolución
for (gen in c("ARPP21", "DLG2", "AXL", "MYCN")) {
  png(paste0("KM_OS_", gen, ".png"), width = 2400, height = 2000, res = 300)
  print(km_gen(gen, "os_days", "os_event", "Overall Survival"))
  dev.off()
}
library(survival)
 
# Modelo univariante para cada gen
for (gen in c("ARPP21", "DLG2", "AXL", "MYCN")) {
  cox_uni <- coxph(Surv(os_days, os_event) ~ surv_data[[gen]], 
                   data = surv_data)
  cat("\n=====", gen, "=====\n")
  print(summary(cox_uni)$coefficients)
}
 
# Modelo multivariante — ajustando por edad, high_risk, y MYCN
cox_multi <- coxph(
  Surv(os_days, os_event) ~ ARPP21 + DLG2 + AXL + MYCN + age_days + high_risk,
  data = surv_data
)
summary(cox_multi)
 
 
 
library(limma)
 
# ── 1. Verificar distribución de grupos ──────────────────────────────────────
table(surv_data$high_risk)
# Debería dar algo como: 0 = no-HR, 1 = HR
 
# ── 2. Preparar grupos y matriz ──────────────────────────────────────────────
group_hr <- factor(surv_data$high_risk, levels = c(0, 1),
                   labels = c("nonHR", "HR"))
 
# La matriz expr_filt ya está en log2RPM — usarla directamente
# Asegurarse de que el orden de columnas coincide con surv_data
all(colnames(expr_filt) == rownames(surv_data))  # debe dar TRUE
 
# ── 3. Modelo lineal limma ───────────────────────────────────────────────────
design <- model.matrix(~ group_hr)
colnames(design) <- c("Intercept", "HR_vs_nonHR")
 
fit  <- lmFit(expr_filt, design)
fit  <- eBayes(fit)
 
results_hr <- topTable(
  fit,
  coef          = "HR_vs_nonHR",
  number        = Inf,
  adjust.method = "fdr",
  sort.by       = "p"
)
 
# Resumen
DEGs_hr <- subset(results_hr, adj.P.Val < 0.05 & abs(logFC) > 1)
cat("Total DEGs HR vs nonHR:", nrow(DEGs_hr), "\n")
cat("Sobreexpresados en HR: ", sum(DEGs_hr$logFC > 0), "\n")
cat("Infraexpresados en HR: ", sum(DEGs_hr$logFC < 0), "\n")
 
head(results_hr[, c("logFC", "adj.P.Val", "AveExpr")], 20)
 
 
 
 
 
 
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
 
# ── 1. Anotar RefSeq → símbolo génico ───────────────────────────────────────
results_hr$refseq <- rownames(results_hr)
 
conv <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = results_hr$refseq,
  keytype = "REFSEQ",
  columns = c("REFSEQ", "SYMBOL")
)
 
# Quedarse con la primera coincidencia por RefSeq
conv <- conv[!duplicated(conv$REFSEQ), ]
results_hr$Symbol <- conv$SYMBOL[match(results_hr$refseq, conv$REFSEQ)]
 
# Cuántos anotados
cat("Genes anotados con símbolo:", sum(!is.na(results_hr$Symbol)), "\n")
 
# ── 2. DEGs anotados ─────────────────────────────────────────────────────────
DEGs_hr$refseq  <- rownames(DEGs_hr)
DEGs_hr$Symbol  <- conv$SYMBOL[match(DEGs_hr$refseq, conv$REFSEQ)]
DEGs_hr_annot   <- DEGs_hr[!is.na(DEGs_hr$Symbol), ]
 
cat("DEGs con símbolo génico:", nrow(DEGs_hr_annot), "\n")
 
# Top 10 sobreexpresados y top 10 infraexpresados
top_up   <- head(DEGs_hr_annot[order(-DEGs_hr_annot$logFC), 
                               c("Symbol","logFC","adj.P.Val")], 10)
top_down <- head(DEGs_hr_annot[order(DEGs_hr_annot$logFC),  
                               c("Symbol","logFC","adj.P.Val")], 10)
 
cat("\nTop 10 sobreexpresados en HR:\n"); print(top_up)
cat("\nTop 10 infraexpresados en HR (más altos en nonHR):\n"); print(top_down)
 
# ── 3. Volcano plot ──────────────────────────────────────────────────────────
results_hr$categoria <- "No significativo"
results_hr$categoria[results_hr$logFC >  1 & results_hr$adj.P.Val < 0.05] <- "Sobreexpresado en HR"
results_hr$categoria[results_hr$logFC < -1 & results_hr$adj.P.Val < 0.05] <- "Infraexpresado en HR"
 
ggplot(results_hr, aes(x=logFC, y=-log10(adj.P.Val), color=categoria)) +
  geom_point(alpha=0.6, size=1.5) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="gray50") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="gray50") +
  scale_color_manual(values=c("Sobreexpresado en HR"="#D7191C",
                              "No significativo"="grey70",
                              "Infraexpresado en HR"="#2C7BB6")) +
  labs(title="Volcano Plot — HR vs non-HR (GSE62564)",
       x="log2 Fold Change", y="-log10 (p-valor ajustado)", color="") +
  theme_minimal(base_size=14)
 
# ── 4. ORA — GO Biological Process ──────────────────────────────────────────
genes_up   <- DEGs_hr_annot$Symbol[DEGs_hr_annot$logFC > 0]
genes_down <- DEGs_hr_annot$Symbol[DEGs_hr_annot$logFC < 0]
 
ego_up <- enrichGO(gene=genes_up, OrgDb=org.Hs.eg.db, keyType="SYMBOL",
                   ont="BP", pAdjustMethod="BH", qvalueCutoff=0.05, readable=TRUE)
 
ego_down <- enrichGO(gene=genes_down, OrgDb=org.Hs.eg.db, keyType="SYMBOL",
                     ont="BP", pAdjustMethod="BH", qvalueCutoff=0.05, readable=TRUE)
 
dotplot(ego_up,   showCategory=15) + ggtitle("GO BP — Sobreexpresados en HR")
dotplot(ego_down, showCategory=15) + ggtitle("GO BP — Infraexpresados en HR (↑ nonHR)")
 
# ── 5. Verificar genes candidatos en los resultados ──────────────────────────
candidatos <- c("ARPP21","DLG2","AXL","MYCN")
results_hr_cand <- results_hr[results_hr$Symbol %in% candidatos, 
                              c("Symbol","logFC","adj.P.Val")]
print(results_hr_cand[order(results_hr_cand$adj.P.Val), ])
cox.zph(cox_multi)
 

6. GSE62564 — Supervivencia: Kaplan-Meier y regresión de Cox por gen

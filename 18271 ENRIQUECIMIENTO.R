library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
 
 
# Filtrar DEGs significativos
DEGs <- subset(results, adj.P.Val < 0.05 & abs(logFC) > 1)
 
# Extraer símbolos válidos
genes <- DEGs$Symbol
genes <- genes[!is.na(genes)]
 
 
ego <- enrichGO(
  gene          = genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
 
# Visualización
dotplot(ego, showCategory = 20) + ggtitle("GO Biological Process")
 
 
# Convertir símbolos a ENTREZ (KEGG usa ENTREZ)
entrez_genes <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
 
ekegg <- enrichKEGG(
  gene         = entrez_genes,
  organism     = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)
 
dotplot(ekegg, showCategory = 20) + ggtitle("KEGG Pathways")
 
 
 
gene_list <- results$logFC
names(gene_list) <- results$Symbol
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!is.na(names(gene_list))]
 
 
 
gsea_go <- gseGO(
  geneList     = gene_list,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  keyType      = "SYMBOL",
  pAdjustMethod = "BH",
  verbose      = FALSE
)
 
ridgeplot(gsea_go) + ggtitle("GSEA - GO Biological Process")
 
 
 
# Necesita ENTREZ
gene_list_entrez <- bitr(names(gene_list), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
gene_list_kegg <- gene_list[gene_list_entrez$SYMBOL]
names(gene_list_kegg) <- gene_list_entrez$ENTREZID
 
gsea_kegg <- gseKEGG(
  geneList = gene_list_kegg,
  organism = "hsa",
  pAdjustMethod = "BH",
  verbose = FALSE
)
 
gseaplot2(gsea_kegg, geneSetID = 1)
 
 
write.csv(as.data.frame(ego), "GO_BP_results.csv")
write.csv(as.data.frame(ekegg), "KEGG_results.csv")
write.csv(as.data.frame(gsea_go), "GSEA_GO_results.csv")
write.csv(as.data.frame(gsea_kegg), "GSEA_KEGG_results.csv")
 
 

3. GSE212877 — Expresión diferencial ADRN vs. MES (con voom)

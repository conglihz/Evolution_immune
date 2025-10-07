library("DESeq2")
library(ggplot2)
library(dplyr)
library("RColorBrewer")
library(pheatmap)
library(ggrepel)
library(tibble)
library(DEGreport)
library(cowplot)

dmel_dsim <- read.csv("dvir_up_genes_not_dmel_up_pr_0812.txt", header = FALSE)
dmel_no_dsim <- read.csv("dmel_up_genes_not_dsim_up_pr_0812.txt", header = FALSE)
## GO enrichment analysis
library(clusterProfiler)
gene = as.list(dmel_dsim[, 1]) 
gene_pr = as.list(dmel_no_dsim[, 1]) 
head(gene)

library(org.Dm.eg.db)

ego <- enrichGO(gene          = gene,
                universe      = keys(org.Dm.eg.db, keytype = "FLYBASE"),
                OrgDb         = org.Dm.eg.db,
                keyType       = "FLYBASE",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
ego_pr <- enrichGO(gene       = gene_pr,
                universe      = keys(org.Dm.eg.db, keytype = "FLYBASE"),
                OrgDb         = org.Dm.eg.db,
                keyType       = "FLYBASE",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# Visualize enriched GO terms as a directed acyclic graph
goplot(ego)

# Bar plot
library(enrichplot)

go_bar_ef <- barplot(ego, showCategory=5, order=T) +
              #coord_cartesian(xlim = c(0, 80)) + 
              ggtitle(expression(italic("Enterococcus faecalis")))+
              theme(
                plot.title = element_text(hjust = 0, size = 8), 
                axis.title = element_text(size = 8), 
                axis.text.x = element_text(size = 5),  
                axis.text.y = element_text(size = 5), 
                legend.title = element_text(size = 5),
                legend.text = element_text(size = 5),
                legend.key.size = unit(0.4, "cm"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()
              )
print(go_bar_ef)
ggsave("dmel_dsim.pdf", go_bar_ef, width = 3, height = 3, dpi = 300)


go_bar_pr <- barplot(ego_pr, showCategory=10, order=T) + 
              #coord_cartesian(xlim = c(0, 80)) + 
              ggtitle(expression(italic("Providencia rettgeri")))+
              theme(
                plot.title = element_text(hjust = 0, size = 8), 
                axis.title = element_text(size = 8), 
                axis.text.x = element_text(size = 5),
                axis.text.y = element_text(size = 5),
                legend.title = element_text(size = 5),
                legend.text = element_text(size = 5),
                legend.key.size = unit(0.4, "cm"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()
              )
print(go_bar_pr)
ggsave("dmel_no_dsim_pr0812.pdf", go_bar_pr, width = 3, height = 3, dpi = 300)


bar_plot <- plot_grid(go_bar_ef+ theme(legend.position="none"), 
                      go_bar_pr+ theme(legend.position="none"))#, labels = c('C', 'D'))
print(bar_plot)
legend_ef <- get_legend(go_bar_ef)
ggsave("go_bar_legend1_0220.pdf", legend_ef, width = 2, height = 2, dpi = 300)
legend_pr <- get_legend(go_bar_pr)
ggsave("go_bar_legend2_0220.pdf", legend_pr, width = 2, height = 2, dpi = 300)

ggsave("go_bar_0220_3.pdf", bar_plot, width = 4, height = 3, dpi = 300)

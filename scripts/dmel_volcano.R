library("DESeq2")
library(ggplot2)
library(dplyr)
library("RColorBrewer")
library(pheatmap)
library(ggrepel)
library(tibble)
library(DEGreport)
library(cowplot)


## volcano plot
plot_pr <- ggplot(res_table_pr_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = gene_category), size = 0.3) +
  ggtitle(expression(italic("Providencia rettgeri"))) +
  xlab("Log2 fold change") + 
  ylab("Adjusted p-value (-log10)") +
  geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "dashed", size = 0.3) +
  geom_vline(xintercept = c(-1,1), color = "grey", linetype = "dashed", size = 0.3) +
  scale_x_continuous(limits = c(-5,12)) +
  scale_y_continuous(limits = c(0,50)) +
  scale_color_manual(values = c("Upregulated" = "#F8766D", "Downregulated" = "#00BA38", "Non-significant" = "grey")) +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.line = element_line(color = "black", linewidth = 0.4),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )  
print(plot_pr)

legend <- get_legend(plot_ef)
volcano <- plot_grid(plot_ef+ theme(legend.position="none"),
                     plot_pr+ ylab(NULL))#, labels = c('A', 'B'))
volcano <- plot_grid(volcano, legend, rel_widths = c(4, 1))
print(volcano)
ggsave("volcano_0220.pdf", plot = volcano, width = 4, height = 2) 
ggsave("legend_0220.pdf", plot = legend, width = 2, height = 2) 




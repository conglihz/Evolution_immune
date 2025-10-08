
# Rapid evolution of innate immunity revealed by species-biased responses to bacteria in Drosophila


This repository contains scripts and data to support the paper.

## Scripts

bacterial_load.R
- plot bacterial load at day 1.

survival_curve.R
- plot a survival curve with survival data using Kaplan-Meier estimate.

bacterial_load_correlation_with_survival.ipynb
- do a correlaton between bacterial load and survival considering phylogenetic relationships.

dmel_volcano.R
- plot a volcano plot using deseq2 result table.

plotGO.R
- plot GO enrichment bar plot.

pivot_table.ipynb
- get a pivot table with result from deseq2 comparing untreat vs bacteria_infected groups.

EDS&spearman.ipynb
- calculate spearman's correlation and EDS score based on log2 fold changes between untreat vs bacteria_infected groups across five species.

pagels_lambda_calculation.R
- calculate pagel's lambda values.

pagel_lambda_comparison.ipynb
- compare pagel's lambda values between upregulated genes and other genes.

check_gene_conservation&expression.ipynb
- based on orthomcl results and expression data, check the fold changes of upregulated genes across species (if one gene is upregulated in another species).

find_pattern_genes.ipynb
- find genes with specific expression patterns.

plot_individual_gene_expression.ipynb
- for individual genes, plot its expression across species.

check_allgene_expression.ipynb
- check fold changes of conserved/nonconserved genes, if a gene is regulated in the same direction across species.

get_new_transcripts.ipynb
- get new transcripts from RNAseq results that are upregulated.

get_orf.ipynb
- get orfs from upregulated candidate transcripts.

all_orfs.ipynb
- get all possible orfs from genome/transcriptome and calculate their physicochemical properties.


## Data
groups_g.txt
- orthomcl results in orthologous gene clusters.

combined_res_singlecopy_ef.txt
- combined result table from Deseq2 comparing untreated vs bacteria-infected, converting gene ID in non-dmel to corresponding gene ID in dmel from E. faecalis dataset.

combined_res_singlecopy_pr.txt
- combined result table from Deseq2 comparing untreated vs bacteria-infected, converting gene ID in non-dmel to corresponding gene ID in dmel from P. rettgeri dataset.

setwd("~/Documents/Zhao_Lab/manuscript/figures/RNAseq_updated/5species/pagels_lambda")
library(ape)
library(phytools)

# load data
df <- read.csv("pivot_df_pr0708.csv")

# load and prepare the tree
tree <- read.tree("tree.nwk")
species_cols <- c("dmel", "dsim", "dana", "dvir", "sleb")

# ensure tree tips match species columns
common_species <- intersect(tree$tip.label, species_cols)
tree <- drop.tip(tree, setdiff(tree$tip.label, common_species))

# initialize results list
lambda_results <- data.frame(
  gene = df$gene,
  lambda = NA,
  logL = NA,
  p_value = NA
)

# loop over genes
for (i in 1:nrow(df)) {
  gene_id <- df$gene[i]
  trait_vector <- as.numeric(df[i, common_species])
  names(trait_vector) <- common_species
  
  # skip if all NA or no variation
  if (all(is.na(trait_vector)) || var(trait_vector, na.rm = TRUE) == 0) {
    next
  }
  
  # try-catch to handle potential fitting errors
  try({
    result <- phylosig(tree, trait_vector, method = "lambda", test = TRUE)
    lambda_results$lambda[i] <- result$lambda
    lambda_results$logL[i] <- result$logL
    lambda_results$p_value[i] <- result$P
  }, silent = TRUE)
}

# save results
write.csv(lambda_results, "pagel_lambda_by_gene_pr.csv", row.names = FALSE)

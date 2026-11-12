# Example: microscopy eukaryotic phytoplankton at the genus level.
# Keep all inputs in the working directory. Abundance tables have taxa in rows and samples in columns.

library(ape)
library(NST)
library(iCAMP)
library(vegan)

# 1. Inputs and parameters
taxonomy_file <- "IdenPlantGenu_Asm_Taxonomy.csv"
abundance_file <- "IdenPlantGenu_Asm.csv"
group_file <- "treatment.csv"
group_order <- c("C", "E", "P", "EP", "W", "WE", "WP", "WEP")
seed <- 123
rand_number <- 1000
n_worker <- 1
dataset <- sub("\\.csv$", "", abundance_file)

# 2. Build a taxonomy-derived proxy tree with unit branch lengths
taxonomy <- read.csv(taxonomy_file, check.names = FALSE)
# For genus-level analysis, use Phylum/Class/Order/Family/Genus; for species-level analysis, also read the Species column and add /Species to the formula below.
taxonomy <- as.data.frame(lapply(taxonomy[, c("Phylum", "Class", "Order", "Family", "Genus")], as.factor))
tree <- as.phylo(~Phylum/Class/Order/Family/Genus, data = taxonomy, collapse = FALSE)
tree$edge.length <- rep(1, nrow(tree$edge))
write.tree(tree, file = paste0(dataset, "_ProxyTree.tre"))

# 3. Match taxa and samples, then convert to relative abundance
abundance <- read.csv(abundance_file, row.names = 1, check.names = FALSE)
comm <- t(as.matrix(abundance))
comm <- comm[rowSums(comm) > 0, , drop = FALSE]
comm <- comm[, colSums(comm) > 0, drop = FALSE]
stopifnot(!anyDuplicated(tree$tip.label), setequal(colnames(comm), tree$tip.label))
comm <- comm[, tree$tip.label, drop = FALSE]
comm_rel <- decostand(comm, method = "total")
pd <- cophenetic.phylo(tree)
treatment <- read.csv(group_file, row.names = 1, check.names = FALSE)
groups <- treatment[rownames(comm), "Treatment", drop = FALSE]
names(groups) <- "group"
groups$group <- droplevels(factor(groups$group, levels = group_order))

# 4. pNST: obtain pairwise pMST and betaNTI from the same null-model run
set.seed(seed)
phylogenetic_result <- pNST(
  comm = comm_rel,
  pd = pd,
  group = groups,
  abundance.weighted = TRUE,
  rand = rand_number,
  nworker = n_worker,
  phylo.shuffle = TRUE,
  taxo.null.model = NULL,
  dirichlet = TRUE,
  SES = TRUE,
  RC = FALSE,
  output.rand = FALSE
)
write.csv(phylogenetic_result$index.pair.grp, paste0(dataset, "_ProxyTree_pMST_pairwise.csv"), row.names = FALSE)

# 5. Calculate RCbray independently, using all retained samples as one regional pool
set.seed(seed)
rc_result <- RC.cm(
  comm = comm_rel,
  rand = rand_number,
  nworker = n_worker,
  weighted = TRUE,
  meta.group = NULL,
  sig.index = "RC",
  taxo.metric = "bray",
  dirichlet = TRUE,
  na.zero = FALSE,
  output.bray = FALSE,
  detail.null = FALSE
)

# 6. Save betaNTI and RCbray for within-treatment sample pairs
beta_pairs <- phylogenetic_result$index.pair
pairwise <- data.frame(name1 = beta_pairs$name1, name2 = beta_pairs$name2, betaNTI = beta_pairs$bNTI.wt)
row_index <- match(pairwise$name1, rownames(rc_result$index))
col_index <- match(pairwise$name2, colnames(rc_result$index))
pairwise$RCbray <- rc_result$index[cbind(row_index, col_index)]
group1 <- as.character(groups[pairwise$name1, "group"])
group2 <- as.character(groups[pairwise$name2, "group"])
pairwise$group <- group1
pairwise <- pairwise[group1 == group2, c("group", "name1", "name2", "betaNTI", "RCbray")]
stopifnot(all(is.finite(pairwise$betaNTI)), all(is.finite(pairwise$RCbray)))
write.csv(pairwise, paste0(dataset, "_ProxyTree_betaNTI_RCbray_pairwise.csv"), row.names = FALSE)

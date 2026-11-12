# Calculate betaNTI and RCbray for within-treatment sample pairs.
# Keep abundance tables, trees, the grouping table and this script in the working directory.

# 1. Inputs and parameters
input_files <- rbind(
  c("IdenPlantGenu_Asm.csv", "IdenPlantGenu_Asm_ProxyTree.tre"),
  c("IdenPlantSpec_Asm.csv", "IdenPlantSpec_Asm_ProxyTree.tre"),
  c("IdenZooGenu_Asm.csv", "IdenZooGenu_Asm_ProxyTree.tre"),
  c("IdenZooSpec_Asm.csv", "IdenZooSpec_Asm_ProxyTree.tre"),
  c("SeqPlantOtu_Asm.csv", "SeqPlantOtu_Asm_SeqTree.nwk"),
  c("SeqPlantSpec_Asm.csv", "SeqPlantSpec_Asm_SeqTree.nwk"),
  c("SeqPlantGenu_Asm.csv", "SeqPlantGenu_Asm_SeqTree.nwk"),
  c("SeqZooOtu_Asm.csv", "SeqZooOtu_Asm_SeqTree.nwk"),
  c("SeqZooSpec_Asm.csv", "SeqZooSpec_Asm_SeqTree.nwk"),
  c("SeqZooGenu_Asm.csv", "SeqZooGenu_Asm_SeqTree.nwk"),
  c("SeqPlantGenu_Asm.csv", "SeqPlantGenu_Asm_ProxyTree.tre"),
  c("SeqPlantSpec_Asm.csv", "SeqPlantSpec_Asm_ProxyTree.tre"),
  c("SeqZooGenu_Asm.csv", "SeqZooGenu_Asm_ProxyTree.tre"),
  c("SeqZooSpec_Asm.csv", "SeqZooSpec_Asm_ProxyTree.tre"),
  c("IdenPlantGenu_AsmCy.csv", "IdenPlantGenu_AsmCy_ProxyTree.tre"),
  c("IdenPlantSpec_AsmCy.csv", "IdenPlantSpec_AsmCy_ProxyTree.tre")
)
group_file <- "treatment.csv" # Sample names are in the first column; Treatment identifies the treatment group.
group_order <- c("C", "E", "P", "EP", "W", "WE", "WP", "WEP") 
seed <- 123 
rand_number <- 1000 
n_worker <- 1 

# 2. Load packages and sample groups
library(NST) 
library(iCAMP) 
library(vegan) 
library(ape) 
treatment <- read.csv(group_file, row.names = 1, check.names = FALSE)

# 3. Analyze each abundance-table and tree combination
# Abundance tables: taxa in rows, samples in columns, unique taxon names in the first column; finite, nonnegative values without missing data.
for (i in seq_len(nrow(input_files))) { 
  abundance_file <- input_files[i, 1] 
  tree_file <- input_files[i, 2] 
  dataset <- sub("\\.(tre|nwk)$", "", tree_file) 
  abundance <- read.csv(abundance_file, row.names = 1, check.names = FALSE) 
  comm <- t(as.matrix(abundance)) 
  comm <- comm[rowSums(comm) > 0, , drop = FALSE] 
  comm <- comm[, colSums(comm) > 0, drop = FALSE] 

  tree <- read.tree(tree_file) 
  colnames(comm) <- gsub("[[:space:]]+", "_", sub("^[a-z]__", "", trimws(colnames(comm)))) 
  tree$tip.label <- gsub("[[:space:]]+", "_", sub("^[a-z]__", "", trimws(tree$tip.label)))
  tree <- keep.tip(tree, which(tree$tip.label %in% colnames(comm) & !duplicated(tree$tip.label))) 
  stopifnot(setequal(colnames(comm), tree$tip.label)) 
  comm <- comm[, tree$tip.label, drop = FALSE] 
  comm_rel <- decostand(comm, method = "total") 
  pd <- cophenetic.phylo(tree) 
  groups <- treatment[rownames(comm), "Treatment", drop = FALSE] 
  names(groups) <- "group" 
  groups$group <- droplevels(factor(groups$group, levels = group_order)) 

  # 4. betaNTI
  set.seed(seed) 
  beta_result <- pNST( 
    comm = comm_rel,
    pd = pd, 
    group = groups, 
    abundance.weighted = TRUE,
    phylo.shuffle = TRUE, 
    taxo.null.model = NULL, 
    rand = rand_number, 
    nworker = n_worker, 
    SES = TRUE, 
    RC = FALSE, 
    output.rand = FALSE 
  )
  beta_pairs <- beta_result$index.pair 

  # 5. RCbray
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

  # 6. Save within-treatment pairwise metrics
  pairwise <- data.frame(name1 = beta_pairs$name1, name2 = beta_pairs$name2, betaNTI = beta_pairs$bNTI.wt) 
  row_index <- match(pairwise$name1, rownames(rc_result$index)) 
  col_index <- match(pairwise$name2, colnames(rc_result$index)) 
  pairwise$RCbray <- rc_result$index[cbind(row_index, col_index)] 
  group1 <- as.character(groups[pairwise$name1, "group"]) 
  group2 <- as.character(groups[pairwise$name2, "group"]) 
  pairwise$group <- group1
  pairwise <- pairwise[group1 == group2, c("group", "name1", "name2", "betaNTI", "RCbray")]
  stopifnot(all(is.finite(pairwise$betaNTI)), all(is.finite(pairwise$RCbray)))
  write.csv(pairwise, paste0(dataset, "_betaNTI_RCbray_pairwise.csv"), row.names = FALSE)
}

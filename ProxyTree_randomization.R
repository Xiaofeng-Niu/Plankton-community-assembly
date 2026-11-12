# Proxy-tree randomization: compare pairwise pMST, betaNTI and RCbray for the reference tree and 999 randomized trees.
# Keep input files and this script in the working directory.

# 1. Inputs and parameters
input_files <- rbind(
  c("IdenPlantGenu_Asm.csv", "IdenPlantGenu_Asm_ProxyTree.tre"),
  c("IdenPlantSpec_Asm.csv", "IdenPlantSpec_Asm_ProxyTree.tre"),
  c("IdenZooGenu_Asm.csv", "IdenZooGenu_Asm_ProxyTree.tre"),
  c("IdenZooSpec_Asm.csv", "IdenZooSpec_Asm_ProxyTree.tre"),
  c("IdenPlantGenu_AsmCy.csv", "IdenPlantGenu_AsmCy_ProxyTree.tre"),
  c("IdenPlantSpec_AsmCy.csv", "IdenPlantSpec_AsmCy_ProxyTree.tre")
)
group_file <- "treatment.csv" # Sample names are in the first column; Treatment identifies the treatment group.
group_order <- c("C", "E", "P", "EP", "W", "WE", "WP", "WEP")
seed <- 123
tree_number <- 999
rand_number <- 1000
n_worker <- 1

# 2. Load packages and sample groups
library(NST)
library(iCAMP)
library(vegan)
library(ape)
treatment <- read.csv(group_file, row.names = 1, check.names = FALSE)

# 3. Prepare the six datasets
for (i in seq_len(nrow(input_files))) {
  abundance_file <- input_files[i, 1]
  tree_file <- input_files[i, 2]
  dataset <- sub("\\.(tre|nwk)$", "", tree_file)
  # Abundance tables: taxa in rows, samples in columns, unique taxon names in the first column; finite, nonnegative values without missing data.
  abundance <- read.csv(abundance_file, row.names = 1, check.names = FALSE)
  comm <- t(as.matrix(abundance))
  comm <- comm[rowSums(comm) > 0, , drop = FALSE]
  comm <- comm[, colSums(comm) > 0, drop = FALSE]
  tree <- read.tree(tree_file)
  colnames(comm) <- gsub("[[:space:]]+", "_", sub("^[a-z]__", "", trimws(colnames(comm))))
  tree$tip.label <- gsub("[[:space:]]+", "_", sub("^[a-z]__", "", trimws(tree$tip.label)))
  stopifnot(!anyDuplicated(colnames(comm)), !anyDuplicated(tree$tip.label), setequal(colnames(comm), tree$tip.label))
  stopifnot(all(is.finite(tree$edge.length)), all(tree$edge.length > 0), is.ultrametric(tree))
  comm <- comm[, tree$tip.label, drop = FALSE]
  comm_rel <- decostand(comm, method = "total")
  groups <- treatment[rownames(comm), "Treatment", drop = FALSE]
  names(groups) <- "group"
  groups$group <- droplevels(factor(groups$group, levels = group_order))

  # 4. Calculate RCbray once per dataset, using all retained samples as one regional pool
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

  # 5. Reference and randomized proxy trees
  pair_results <- vector("list", tree_number + 1)
  for (r in 0:tree_number) { # Replicate 0 is the reference tree; replicates 1-999 are randomized trees.
    current_tree <- tree
    if (r > 0) {
      set.seed(seed + i * 1000000 + r)
      current_tree <- multi2di(tree, random = TRUE)
      new_edges <- which(current_tree$edge.length == 0)
      current_tree$edge.length[new_edges] <- runif(length(new_edges), 0.01, 0.10)
      tip_depth <- node.depth.edgelength(current_tree)[seq_len(Ntip(current_tree))]
      terminal_edges <- match(seq_len(Ntip(current_tree)), current_tree$edge[, 2])
      current_tree$edge.length[terminal_edges] <- current_tree$edge.length[terminal_edges] + max(tip_depth) - tip_depth
      current_tree <- reorder.phylo(current_tree, order = "cladewise")
    }
    pd <- cophenetic.phylo(current_tree) # Recalculate taxon distances from the current tree in every iteration.
    set.seed(seed + i * 2000000 + r)
    result <- pNST(
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

    # 6. Match within-treatment pairwise metrics by sample names
    mst <- result$index.pair.grp
    beta <- result$index.pair
    mst_key <- paste(pmin(mst$name1, mst$name2), pmax(mst$name1, mst$name2), sep = "|")
    beta_key <- paste(pmin(beta$name1, beta$name2), pmax(beta$name1, beta$name2), sep = "|")
    beta_index <- match(mst_key, beta_key)
    rc_row <- match(mst$name1, rownames(rc_result$index))
    rc_col <- match(mst$name2, colnames(rc_result$index))
    pairwise <- data.frame(replicate = r, group = mst$group, name1 = mst$name1, name2 = mst$name2,
                           pMST = mst$MST.ij.bMNTD, betaNTI = beta$bNTI.wt[beta_index],
                           RCbray = rc_result$index[cbind(rc_row, rc_col)])
    stopifnot(!anyDuplicated(mst_key), all(is.finite(as.matrix(pairwise[, c("pMST", "betaNTI", "RCbray")]))))
    pair_results[[r + 1]] <- pairwise
  }
  pairwise <- do.call(rbind, pair_results)
  write.csv(pairwise, paste0(dataset, "_randomization_pairwise.csv"), row.names = FALSE)
}

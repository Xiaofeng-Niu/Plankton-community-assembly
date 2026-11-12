# tNST: calculate pairwise modified stochasticity ratios within treatments.
# Keep abundance tables, the grouping table and this script in the working directory.

# 1. Inputs and parameters
input_files <- c(
  "IdenPlantGenu_Asm.csv",
  "IdenPlantSpec_Asm.csv",
  "IdenZooGenu_Asm.csv",
  "IdenZooSpec_Asm.csv",
  "SeqPlantOtu_Asm.csv",
  "SeqPlantSpec_Asm.csv",
  "SeqPlantGenu_Asm.csv",
  "SeqZooOtu_Asm.csv",
  "SeqZooSpec_Asm.csv",
  "SeqZooGenu_Asm.csv",
  "IdenPlantGenu_AsmCy.csv",
  "IdenPlantSpec_AsmCy.csv"
)
group_file <- "treatment.csv" # Sample names are in the first column; Treatment identifies the treatment group.
group_order <- c("C", "E", "P", "EP", "W", "WE", "WP", "WEP")
seed <- 123 
rand_number <- 1000 
n_worker <- 1

# 2. Load packages and sample groups
library(NST) 
library(vegan) 
treatment <- read.csv(group_file, row.names = 1, check.names = FALSE) 

# 3. Analyze each dataset
# Abundance tables: taxa in rows, samples in columns, unique taxon names in the first column; finite, nonnegative values without missing data.
for (input_file in input_files) {
  dataset <- sub("\\.csv$", "", input_file)
  abundance <- read.csv(input_file, row.names = 1, check.names = FALSE) 
  comm <- t(as.matrix(abundance)) 
  comm <- comm[rowSums(comm) > 0, , drop = FALSE]
  comm <- comm[, colSums(comm) > 0, drop = FALSE] 
  groups <- treatment[rownames(comm), "Treatment", drop = FALSE] 
  names(groups) <- "group" 
  groups$group <- droplevels(factor(groups$group, levels = group_order)) 
  comm_rel <- decostand(comm, method = "total") 

  set.seed(seed) 
  result <- tNST(
    comm = comm_rel, 
    group = groups, 
    dist.method = "jaccard", 
    abundance.weighted = TRUE, 
    rand = rand_number, 
    nworker = n_worker, 
    null.model = "PF", 
    dirichlet = TRUE 
  )
  pairwise <- result$index.pair.grp 
  write.csv(pairwise, paste0(dataset, "_tMST_pairwise.csv"), row.names = FALSE) }

# Fig. S4a: alpha diversity and paired comparisons of the MIC and SEQ assembly datasets.
# Inputs and outputs use the working directory; abundance tables have taxa in rows, samples in columns, and unique taxon names in the first column.

# 1. Inputs
input_files <- rbind(
  c("IdenPlantGenu_Asm.csv", "MIC", "Plant", "Genus"),
  c("IdenPlantSpec_Asm.csv", "MIC", "Plant", "Species"),
  c("IdenZooGenu_Asm.csv", "MIC", "Zoo", "Genus"),
  c("IdenZooSpec_Asm.csv", "MIC", "Zoo", "Species"),
  c("SeqPlantGenu_Asm.csv", "SEQ", "Plant", "Genus"),
  c("SeqPlantSpec_Asm.csv", "SEQ", "Plant", "Species"),
  c("SeqZooGenu_Asm.csv", "SEQ", "Zoo", "Genus"),
  c("SeqZooSpec_Asm.csv", "SEQ", "Zoo", "Species")
)
alpha <- data.frame()

# 2. Calculate three diversity indices from the filtered communities
for (i in seq_len(nrow(input_files))) { 
  abundance <- as.matrix(read.csv(input_files[i, 1], row.names = 1, check.names = FALSE)) 
  totals <- colSums(abundance) 
  relative <- sweep(abundance, 2, ifelse(totals > 0, totals, NA_real_), "/") 
  d <- data.frame(Sample = colnames(abundance), Method = input_files[i, 2], Organism = input_files[i, 3], Rank = input_files[i, 4],
                  Total = totals, Richness = colSums(abundance > 0), Shannon = -colSums(relative * log(pmax(relative, 1e-300))), Simpson = 1 - colSums(relative^2))
  alpha <- rbind(alpha, d) 
}

# 3. Pair MIC and SEQ observations from the same mesocosm
mic <- alpha[alpha$Method == "MIC", ]
seq <- alpha[alpha$Method == "SEQ", ]
pairs <- merge(mic, seq, by = c("Sample", "Organism", "Rank"), suffixes = c("_MIC", "_SEQ"))
pairs$Included <- pairs$Total_MIC > 0 & pairs$Total_SEQ > 0
write.csv(pairs, "Alpha_paired_samples.csv", row.names = FALSE)
pairs <- pairs[pairs$Included, ]
tests <- data.frame()

# 4. Two-sided paired Wilcoxon tests
for (organism in c("Plant", "Zoo")) {
  for (rank in c("Genus", "Species")) {
    d <- pairs[pairs$Organism == organism & pairs$Rank == rank, ]
    for (metric in c("Richness", "Shannon", "Simpson")) {
      mic_values <- d[[paste0(metric, "_MIC")]]
      seq_values <- d[[paste0(metric, "_SEQ")]]
      test <- wilcox.test(mic_values, seq_values, paired = TRUE, alternative = "two.sided", exact = FALSE, correct = TRUE) # Retain the original test settings without multiple-testing correction.
      result <- data.frame(Organism = organism, Rank = rank, Metric = metric, N = nrow(d), W = unname(test$statistic), P = test$p.value,
                           Median_SEQ_minus_MIC = median(seq_values - mic_values))
      tests <- rbind(tests, result)
    }
  }
}
tests$Significance <- ifelse(tests$P < 0.0001, "****", ifelse(tests$P < 0.001, "***", ifelse(tests$P < 0.01, "**", ifelse(tests$P < 0.05, "*", "ns"))))
write.csv(tests, "Alpha_paired_wilcoxon.csv", row.names = FALSE)

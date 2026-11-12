# Binomial logistic models of agreement in pMST classifications.
# Inputs and outputs use the working directory; abundance tables have taxa in rows, samples in columns, and unique taxon names in the first column.
library(mgcv)

# 1. Inputs
comparisons <- rbind(
  c("IdenPlantGenu", "SeqPlantGenu", "Phytoplankton", "Genus"),
  c("IdenPlantSpec", "SeqPlantSpec", "Phytoplankton", "Species"),
  c("IdenZooGenu", "SeqZooGenu", "Zooplankton", "Genus"),
  c("IdenZooSpec", "SeqZooSpec", "Zooplankton", "Species")
)
master <- data.frame()
sample_alpha <- list()

# 2. Calculate diversity from the abundance tables used for assembly analysis
for (dataset in as.vector(comparisons[, 1:2])) {
  abundance <- as.matrix(read.csv(paste0(dataset, "_Asm.csv"), row.names = 1, check.names = FALSE))
  totals <- colSums(abundance)
  relative <- sweep(abundance, 2, ifelse(totals > 0, totals, NA_real_), "/")
  sample_alpha[[dataset]] <- data.frame(Sample = colnames(abundance), Richness = colSums(abundance > 0),
    Shannon = -colSums(relative * log(pmax(relative, 1e-300))), Total = totals)
}

# 3. Match the same sample pairs across observation methods
for (i in seq_len(nrow(comparisons))) {
  mic <- read.csv(paste0(comparisons[i, 1], "_Asm_ProxyTree_pMST_pairwise.csv"), check.names = FALSE)
  seq <- read.csv(paste0(comparisons[i, 2], "_Asm_SeqTree_pMST_pairwise.csv"), check.names = FALSE)
  mic$sample_a <- pmin(mic$name1, mic$name2)
  mic$sample_b <- pmax(mic$name1, mic$name2)
  mic$Pair_ID <- paste(mic$sample_a, mic$sample_b, sep = "__")
  seq$Pair_ID <- paste(pmin(seq$name1, seq$name2), pmax(seq$name1, seq$name2), sep = "__")
  stopifnot(!anyDuplicated(mic$Pair_ID), !anyDuplicated(seq$Pair_ID))
  mic <- mic[, c("Pair_ID", "sample_a", "sample_b", "group", "MST.ij.bMNTD")]
  seq <- seq[, c("Pair_ID", "MST.ij.bMNTD")]
  names(mic)[5] <- "MIC_pMST"
  names(seq)[2] <- "SEQ_pMST"
  d <- merge(mic, seq, by = "Pair_ID", all = TRUE)
  d <- d[is.finite(d$MIC_pMST) & is.finite(d$SEQ_pMST), ]
  d$Organism <- comparisons[i, 3]
  d$Level <- comparisons[i, 4]
  alpha <- merge(sample_alpha[[comparisons[i, 1]]], sample_alpha[[comparisons[i, 2]]], by = "Sample", suffixes = c("_MIC", "_SEQ"))
  alpha$Richness_diff <- abs(alpha$Richness_MIC - alpha$Richness_SEQ)
  alpha$Shannon_diff <- abs(alpha$Shannon_MIC - alpha$Shannon_SEQ)
  a <- match(d$sample_a, alpha$Sample)
  b <- match(d$sample_b, alpha$Sample)
  stopifnot(all(alpha$Total_MIC[c(a, b)] > 0), all(alpha$Total_SEQ[c(a, b)] > 0))
  d$Richness_diff_abs_pair <- (alpha$Richness_diff[a] + alpha$Richness_diff[b]) / 2
  d$Shannon_diff_abs_pair <- (alpha$Shannon_diff[a] + alpha$Shannon_diff[b]) / 2
  d$same_binary_num <- as.integer((d$MIC_pMST > 0.5) == (d$SEQ_pMST > 0.5))
  stopifnot(all(complete.cases(d)), all(d$sample_a != d$sample_b))
  master <- rbind(master, d)
}

# 4. Use pooled scaling constants for all models
rich_mean <- mean(master$Richness_diff_abs_pair)
rich_sd <- sd(master$Richness_diff_abs_pair)
shan_mean <- mean(master$Shannon_diff_abs_pair)
shan_sd <- sd(master$Shannon_diff_abs_pair)
master$z_Richness_diff_abs_pair <- (master$Richness_diff_abs_pair - rich_mean) / rich_sd
master$z_Shannon_diff_abs_pair <- (master$Shannon_diff_abs_pair - shan_mean) / shan_sd
write.csv(data.frame(Variable = c("Richness difference", "Shannon difference"), Mean = c(rich_mean, shan_mean), SD = c(rich_sd, shan_sd)), "Agreement_scaling.csv", row.names = FALSE)
write.csv(master, "Agreement_pair_data.csv", row.names = FALSE)

# 5. Fit four models
master$Organism <- factor(master$Organism, levels = c("Zooplankton", "Phytoplankton")) # Reference group: zooplankton.
master$Level <- factor(master$Level, levels = c("Genus", "Species")) # Reference level: genus.
# Pooled model
d <- master
ids <- sort(unique(c(d$sample_a, d$sample_b))) # Each mesocosm shares one effect across taxonomic levels and organism groups.
incidence <- matrix(0, nrow(d), length(ids), dimnames = list(NULL, make.names(ids)))
incidence[cbind(seq_len(nrow(d)), match(d$sample_a, ids))] <- 1
incidence[cbind(seq_len(nrow(d)), match(d$sample_b, ids))] <- 1
stopifnot(all(rowSums(incidence) == 2)) # Each record must contain two distinct mesocosms.
d$sample_incidence <- incidence
penalty <- list(sample_incidence = list(diag(length(ids))))
model <- gam(same_binary_num ~ Organism + Level + z_Richness_diff_abs_pair + z_Shannon_diff_abs_pair + sample_incidence, data = d, family = binomial(link = "logit"), method = "REML", paraPen = penalty)
fixed_ids <- which(!grepl("^sample_incidence", names(coef(model))))
b <- coef(model)[fixed_ids]
V <- vcov(model, unconditional = TRUE)[fixed_ids, fixed_ids, drop = FALSE] # Retain the original unconditional covariance estimate.
se <- sqrt(diag(V))
coef_pooled <- data.frame(Model = "pooled", Term = names(b), Beta = as.numeric(b), SE = se, Beta_CI_low = b - qnorm(0.975) * se, Beta_CI_high = b + qnorm(0.975) * se, P = 2 * pnorm(abs(b / se), lower.tail = FALSE), OR = exp(b), OR_CI_low = exp(b - qnorm(0.975) * se), OR_CI_high = exp(b + qnorm(0.975) * se), row.names = NULL) # Report separate intervals for beta and odds ratios, without BH correction.

slopes <- setdiff(seq_along(b), which(names(b) == "(Intercept)"))
wald <- as.numeric(t(b[slopes]) %*% solve(V[slopes, slopes, drop = FALSE], b[slopes]))
vc <- gam.vcomp(model)
fit_pooled <- data.frame(Model = "pooled", N_records = nrow(d), N_mesocosms = length(ids), N_agree = sum(d$same_binary_num), N_disagree = sum(1 - d$same_binary_num), Converged = isTRUE(model$converged), Deviance_explained = summary(model)$dev.expl, AIC = AIC(model), Overall_Wald = wald, Df = length(slopes), Overall_P = pchisq(wald, df = length(slopes), lower.tail = FALSE), Shared_sample_SD = vc["sample_incidence", "std.dev"], Min_probability = min(fitted(model)), Max_probability = max(fitted(model)), Max_abs_beta = max(abs(b)))

residual_pooled <- data.frame(Model = "pooled", sample_a = d$sample_a, sample_b = d$sample_b, Organism = d$Organism, Level = d$Level, Treatment = d$group, Observed = d$same_binary_num, Fitted = fitted(model), Response_residual = residuals(model, type = "response"), Deviance_residual = residuals(model, type = "deviance"))
fixed_X <- model.matrix(~ Organism + Level + z_Richness_diff_abs_pair + z_Shannon_diff_abs_pair, data = d)[, -1, drop = FALSE]
fixed_X <- scale(fixed_X, center = TRUE, scale = TRUE)
vif <- diag(solve(cor(fixed_X)))
vif_pooled <- data.frame(Model = "pooled", Term = names(vif), VIF = unname(vif))
# Eukaryotic phytoplankton model
d <- master[master$Organism == "Phytoplankton", ]
ids <- sort(unique(c(d$sample_a, d$sample_b)))
incidence <- matrix(0, nrow(d), length(ids), dimnames = list(NULL, make.names(ids)))
incidence[cbind(seq_len(nrow(d)), match(d$sample_a, ids))] <- 1
incidence[cbind(seq_len(nrow(d)), match(d$sample_b, ids))] <- 1
stopifnot(all(rowSums(incidence) == 2))
d$sample_incidence <- incidence
penalty <- list(sample_incidence = list(diag(length(ids))))
model <- gam(same_binary_num ~ Level + z_Richness_diff_abs_pair + z_Shannon_diff_abs_pair + sample_incidence, data = d, family = binomial(link = "logit"), method = "REML", paraPen = penalty)
fixed_ids <- which(!grepl("^sample_incidence", names(coef(model))))
b <- coef(model)[fixed_ids]
V <- vcov(model, unconditional = TRUE)[fixed_ids, fixed_ids, drop = FALSE]
se <- sqrt(diag(V))
coef_phytoplankton <- data.frame(Model = "phytoplankton", Term = names(b), Beta = as.numeric(b), SE = se, Beta_CI_low = b - qnorm(0.975) * se, Beta_CI_high = b + qnorm(0.975) * se, P = 2 * pnorm(abs(b / se), lower.tail = FALSE), OR = exp(b), OR_CI_low = exp(b - qnorm(0.975) * se), OR_CI_high = exp(b + qnorm(0.975) * se), row.names = NULL)

slopes <- setdiff(seq_along(b), which(names(b) == "(Intercept)"))
wald <- as.numeric(t(b[slopes]) %*% solve(V[slopes, slopes, drop = FALSE], b[slopes]))
vc <- gam.vcomp(model)
fit_phytoplankton <- data.frame(Model = "phytoplankton", N_records = nrow(d), N_mesocosms = length(ids), N_agree = sum(d$same_binary_num), N_disagree = sum(1 - d$same_binary_num), Converged = isTRUE(model$converged), Deviance_explained = summary(model)$dev.expl, AIC = AIC(model), Overall_Wald = wald, Df = length(slopes), Overall_P = pchisq(wald, df = length(slopes), lower.tail = FALSE), Shared_sample_SD = vc["sample_incidence", "std.dev"], Min_probability = min(fitted(model)), Max_probability = max(fitted(model)), Max_abs_beta = max(abs(b)))

residual_phytoplankton <- data.frame(Model = "phytoplankton", sample_a = d$sample_a, sample_b = d$sample_b, Organism = d$Organism, Level = d$Level, Treatment = d$group, Observed = d$same_binary_num, Fitted = fitted(model), Response_residual = residuals(model, type = "response"), Deviance_residual = residuals(model, type = "deviance"))
fixed_X <- model.matrix(~ Level + z_Richness_diff_abs_pair + z_Shannon_diff_abs_pair, data = d)[, -1, drop = FALSE]
fixed_X <- scale(fixed_X, center = TRUE, scale = TRUE)
vif <- diag(solve(cor(fixed_X)))
vif_phytoplankton <- data.frame(Model = "phytoplankton", Term = names(vif), VIF = unname(vif))
# Zooplankton model
d <- master[master$Organism == "Zooplankton", ]
ids <- sort(unique(c(d$sample_a, d$sample_b)))
incidence <- matrix(0, nrow(d), length(ids), dimnames = list(NULL, make.names(ids)))
incidence[cbind(seq_len(nrow(d)), match(d$sample_a, ids))] <- 1
incidence[cbind(seq_len(nrow(d)), match(d$sample_b, ids))] <- 1
stopifnot(all(rowSums(incidence) == 2))
d$sample_incidence <- incidence
penalty <- list(sample_incidence = list(diag(length(ids))))
model <- gam(same_binary_num ~ Level + z_Richness_diff_abs_pair + z_Shannon_diff_abs_pair + sample_incidence, data = d, family = binomial(link = "logit"), method = "REML", paraPen = penalty)
fixed_ids <- which(!grepl("^sample_incidence", names(coef(model))))
b <- coef(model)[fixed_ids]
V <- vcov(model, unconditional = TRUE)[fixed_ids, fixed_ids, drop = FALSE]
se <- sqrt(diag(V))
coef_zooplankton <- data.frame(Model = "zooplankton", Term = names(b), Beta = as.numeric(b), SE = se, Beta_CI_low = b - qnorm(0.975) * se, Beta_CI_high = b + qnorm(0.975) * se, P = 2 * pnorm(abs(b / se), lower.tail = FALSE), OR = exp(b), OR_CI_low = exp(b - qnorm(0.975) * se), OR_CI_high = exp(b + qnorm(0.975) * se), row.names = NULL)

slopes <- setdiff(seq_along(b), which(names(b) == "(Intercept)"))
wald <- as.numeric(t(b[slopes]) %*% solve(V[slopes, slopes, drop = FALSE], b[slopes]))
vc <- gam.vcomp(model)
fit_zooplankton <- data.frame(Model = "zooplankton", N_records = nrow(d), N_mesocosms = length(ids), N_agree = sum(d$same_binary_num), N_disagree = sum(1 - d$same_binary_num), Converged = isTRUE(model$converged), Deviance_explained = summary(model)$dev.expl, AIC = AIC(model), Overall_Wald = wald, Df = length(slopes), Overall_P = pchisq(wald, df = length(slopes), lower.tail = FALSE), Shared_sample_SD = vc["sample_incidence", "std.dev"], Min_probability = min(fitted(model)), Max_probability = max(fitted(model)), Max_abs_beta = max(abs(b)))

residual_zooplankton <- data.frame(Model = "zooplankton", sample_a = d$sample_a, sample_b = d$sample_b, Organism = d$Organism, Level = d$Level, Treatment = d$group, Observed = d$same_binary_num, Fitted = fitted(model), Response_residual = residuals(model, type = "response"), Deviance_residual = residuals(model, type = "deviance"))
fixed_X <- model.matrix(~ Level + z_Richness_diff_abs_pair + z_Shannon_diff_abs_pair, data = d)[, -1, drop = FALSE]
fixed_X <- scale(fixed_X, center = TRUE, scale = TRUE)
vif <- diag(solve(cor(fixed_X)))
vif_zooplankton <- data.frame(Model = "zooplankton", Term = names(vif), VIF = unname(vif))
# Interaction model
d <- master
ids <- sort(unique(c(d$sample_a, d$sample_b)))
incidence <- matrix(0, nrow(d), length(ids), dimnames = list(NULL, make.names(ids)))
incidence[cbind(seq_len(nrow(d)), match(d$sample_a, ids))] <- 1
incidence[cbind(seq_len(nrow(d)), match(d$sample_b, ids))] <- 1
stopifnot(all(rowSums(incidence) == 2))
d$sample_incidence <- incidence
penalty <- list(sample_incidence = list(diag(length(ids))))
model <- gam(same_binary_num ~ Organism * z_Richness_diff_abs_pair + Level + z_Shannon_diff_abs_pair + sample_incidence, data = d, family = binomial(link = "logit"), method = "REML", paraPen = penalty)
fixed_ids <- which(!grepl("^sample_incidence", names(coef(model))))
b <- coef(model)[fixed_ids]
V <- vcov(model, unconditional = TRUE)[fixed_ids, fixed_ids, drop = FALSE]
se <- sqrt(diag(V))
coef_interaction <- data.frame(Model = "interaction", Term = names(b), Beta = as.numeric(b), SE = se, Beta_CI_low = b - qnorm(0.975) * se, Beta_CI_high = b + qnorm(0.975) * se, P = 2 * pnorm(abs(b / se), lower.tail = FALSE), OR = exp(b), OR_CI_low = exp(b - qnorm(0.975) * se), OR_CI_high = exp(b + qnorm(0.975) * se), row.names = NULL)

slopes <- setdiff(seq_along(b), which(names(b) == "(Intercept)"))
wald <- as.numeric(t(b[slopes]) %*% solve(V[slopes, slopes, drop = FALSE], b[slopes]))
vc <- gam.vcomp(model)
fit_interaction <- data.frame(Model = "interaction", N_records = nrow(d), N_mesocosms = length(ids), N_agree = sum(d$same_binary_num), N_disagree = sum(1 - d$same_binary_num), Converged = isTRUE(model$converged), Deviance_explained = summary(model)$dev.expl, AIC = AIC(model), Overall_Wald = wald, Df = length(slopes), Overall_P = pchisq(wald, df = length(slopes), lower.tail = FALSE), Shared_sample_SD = vc["sample_incidence", "std.dev"], Min_probability = min(fitted(model)), Max_probability = max(fitted(model)), Max_abs_beta = max(abs(b)))

residual_interaction <- data.frame(Model = "interaction", sample_a = d$sample_a, sample_b = d$sample_b, Organism = d$Organism, Level = d$Level, Treatment = d$group, Observed = d$same_binary_num, Fitted = fitted(model), Response_residual = residuals(model, type = "response"), Deviance_residual = residuals(model, type = "deviance"))
fixed_X <- model.matrix(~ Organism * z_Richness_diff_abs_pair + Level + z_Shannon_diff_abs_pair, data = d)[, -1, drop = FALSE]
fixed_X <- scale(fixed_X, center = TRUE, scale = TRUE)
vif <- diag(solve(cor(fixed_X)))
vif_interaction <- data.frame(Model = "interaction", Term = names(vif), VIF = unname(vif))
interaction_term <- grep(":", names(b))
write.csv(coef_interaction[grepl(":", coef_interaction$Term), ], "Agreement_interaction_test.csv", row.names = FALSE)
zoo_id <- which(names(b) == "z_Richness_diff_abs_pair")
plant_b <- b[zoo_id] + b[interaction_term]
plant_se <- sqrt(V[zoo_id, zoo_id] + V[interaction_term, interaction_term] + 2 * V[zoo_id, interaction_term]) # Include coefficient covariance in the standard error of the simple slope.
simple <- data.frame(Organism = c("Zooplankton", "Eukaryotic phytoplankton"), Beta = c(b[zoo_id], plant_b), SE = c(sqrt(V[zoo_id, zoo_id]), plant_se))
simple$Beta_CI_low <- simple$Beta - qnorm(0.975) * simple$SE
simple$Beta_CI_high <- simple$Beta + qnorm(0.975) * simple$SE
simple$P <- 2 * pnorm(abs(simple$Beta / simple$SE), lower.tail = FALSE)
write.csv(simple, "Agreement_interaction_richness_simple_slopes.csv", row.names = FALSE)
coefficients <- rbind(coef_pooled, coef_phytoplankton, coef_zooplankton, coef_interaction)
fits <- rbind(fit_pooled, fit_phytoplankton, fit_zooplankton, fit_interaction)
write.csv(coefficients, "Agreement_all_model_coefficients.csv", row.names = FALSE)
write.csv(fits, "Agreement_all_model_fit_summary.csv", row.names = FALSE)

vifs <- rbind(vif_pooled, vif_phytoplankton, vif_zooplankton, vif_interaction)
residual_data <- rbind(residual_pooled, residual_phytoplankton, residual_zooplankton, residual_interaction)
write.csv(vifs, "Agreement_VIF.csv", row.names = FALSE)
write.csv(residual_data, "Agreement_fitted_and_residuals.csv", row.names = FALSE)

# 6. Summarize residuals within fixed predicted-probability intervals
interval_labels <- c("[0.0, 0.2)", "[0.2, 0.4)", "[0.4, 0.6)", "[0.6, 0.8)", "[0.8, 1.0]")
residual_data$Interval <- cut(residual_data$Fitted, breaks = seq(0, 1, 0.2), right = FALSE, include.lowest = TRUE, labels = interval_labels)
residual_summary <- data.frame()
for (model_name in c("pooled", "phytoplankton", "zooplankton", "interaction")) {
  for (interval in interval_labels) {
    d <- residual_data[residual_data$Model == model_name & residual_data$Interval == interval, ]
    row <- data.frame(Model = model_name, Interval = interval, N = nrow(d), Observed_agreement = mean(d$Observed), Mean_prediction = mean(d$Fitted), Mean_response_residual = mean(d$Response_residual))
    row[is.na(row)] <- NA
    residual_summary <- rbind(residual_summary, row)
  }
}
write.csv(residual_summary, "Agreement_residual_intervals.csv", row.names = FALSE)

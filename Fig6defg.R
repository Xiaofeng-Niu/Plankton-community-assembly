# Fig. 6d-g: comparison of reference and randomized proxy-tree results.
# Read randomization outputs from the working directory containing this script; do not rerun null models.
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(gtable)

# 1. Inputs
input_files <- c("IdenPlantGenu_AsmCy_ProxyTree_randomization_pairwise.csv",
                 "IdenPlantSpec_AsmCy_ProxyTree_randomization_pairwise.csv")
dataset_labels <- c("Genus", "Species")
tree_number <- 999
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) stop(paste("缺少结果文件：", paste(missing_files, collapse = ", ")), call. = FALSE)
process_levels <- c("Homogeneous selection", "Heterogeneous selection", "Homogenizing dispersal", "Dispersal limitation", "Undominated")
process_colours <- c("Homogeneous selection" = "#6C94B1", "Heterogeneous selection" = "#94C5D1", "Homogenizing dispersal" = "#8FB794", "Dispersal limitation" = "#718457", "Undominated" = "#C0D09D")
pair_summary <- composition <- agreement <- sample_counts <- data.frame()

# 2. Continuous metrics, process proportions and classification agreement
for (i in seq_along(input_files)) {
  raw <- read.csv(input_files[i], check.names = FALSE)
  stopifnot(all(c("replicate", "group", "name1", "name2", "pMST", "betaNTI", "RCbray") %in% names(raw)))
  stopifnot(setequal(unique(raw$replicate), 0:tree_number), all(is.finite(as.matrix(raw[, c("pMST", "betaNTI", "RCbray")]))))
  raw$pair_id <- paste(pmin(raw$name1, raw$name2), pmax(raw$name1, raw$name2), sep = "|")
  stopifnot(!anyDuplicated(raw[, c("replicate", "pair_id")]), all(raw$name1 != raw$name2))
  reference <- raw[raw$replicate == 0, ]
  reference_index <- match(raw$pair_id, reference$pair_id)
  stopifnot(!anyNA(reference_index), all(table(raw$replicate) == nrow(reference)), all(raw$group == reference$group[reference_index]))
  stopifnot(all(raw$RCbray == reference$RCbray[reference_index]))
  raw$pMST_class <- ifelse(raw$pMST > 0.5, "Stochastic", "Deterministic")
  raw$beta_class <- ifelse(abs(raw$betaNTI) > 2, "Deterministic", "Stochastic")
  raw$process <- "Undominated"
  raw$process[raw$betaNTI > 2] <- "Heterogeneous selection"
  raw$process[raw$betaNTI < -2] <- "Homogeneous selection"
  raw$process[abs(raw$betaNTI) <= 2 & raw$RCbray > 0.95] <- "Dispersal limitation"
  raw$process[abs(raw$betaNTI) <= 2 & raw$RCbray < -0.95] <- "Homogenizing dispersal"
  reference <- raw[raw$replicate == 0, ]
  alt <- raw[raw$replicate > 0, ]
  ref_index <- match(alt$pair_id, reference$pair_id)
  alt$pMST_match <- alt$pMST_class == reference$pMST_class[ref_index]
  alt$beta_match <- alt$beta_class == reference$beta_class[ref_index]
  alt$five_match <- alt$process == reference$process[ref_index]

  continuous <- alt %>% group_by(pair_id) %>%
    summarise(pMST_median = median(pMST), pMST_lower = quantile(pMST, 0.025), pMST_upper = quantile(pMST, 0.975),
              beta_median = median(betaNTI), beta_lower = quantile(betaNTI, 0.025), beta_upper = quantile(betaNTI, 0.975), .groups = "drop")
  continuous <- left_join(reference[, c("pair_id", "pMST", "betaNTI")], continuous, by = "pair_id")
  continuous$Dataset <- dataset_labels[i]
  pair_summary <- rbind(pair_summary, continuous)

  per_tree <- alt %>% group_by(replicate) %>%
    summarise(pMST = mean(pMST_match), betaNTI = mean(beta_match), Five = mean(five_match), .groups = "drop")
  rates <- per_tree %>% pivot_longer(-replicate, names_to = "Method", values_to = "Rate") %>%
    group_by(Method) %>% summarise(Mean = mean(Rate), Lower = quantile(Rate, 0.025), Upper = quantile(Rate, 0.975), .groups = "drop")
  rates$Dataset <- dataset_labels[i]
  agreement <- rbind(agreement, rates)

  counts <- as.data.frame(table(Replicate = factor(raw$replicate, levels = 0:tree_number), Process = factor(raw$process, levels = process_levels)))
  counts$Proportion <- counts$Freq / nrow(reference)
  base_prop <- counts[counts$Replicate == "0", c("Process", "Proportion")]
  base_prop$Type <- "Reference"
  alt_prop <- counts[counts$Replicate != "0", ] %>% group_by(Process) %>% summarise(Proportion = mean(Proportion), .groups = "drop")
  alt_prop$Type <- "Randomized"
  comp <- rbind(base_prop, as.data.frame(alt_prop))
  comp$Dataset <- dataset_labels[i]
  composition <- rbind(composition, comp)
  sample_counts <- rbind(sample_counts, data.frame(Dataset = dataset_labels[i], Pairs = nrow(reference), Randomized_trees = tree_number))
}
pair_summary$Dataset <- factor(pair_summary$Dataset, levels = dataset_labels)
agreement$Dataset <- factor(agreement$Dataset, levels = rev(dataset_labels))
agreement$Method <- factor(agreement$Method, levels = c("pMST", "betaNTI", "Five"))
composition$Dataset <- factor(composition$Dataset, levels = rev(dataset_labels))
composition$Process <- factor(composition$Process, levels = process_levels)
composition$Y <- as.numeric(composition$Dataset) + ifelse(composition$Type == "Reference", 0.17, -0.17)

# 3. Plotting
figure_theme <- theme_classic(base_family = "Arial") +
  theme(axis.text = element_text(size = 18, colour = "black"), axis.title = element_text(size = 20),
        strip.text = element_text(size = 18, face = "bold"), strip.background = element_blank(),
        plot.tag = element_text(size = 20, face = "bold"), plot.tag.position = "topleft",
        legend.title = element_blank(), legend.text = element_text(size = 20), legend.position = "bottom",
        panel.spacing = unit(10, "mm"), plot.margin = margin(7, 7, 7, 7, "mm"))
p_mst <- ggplot(pair_summary, aes(pMST, pMST_median)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "#777777", linewidth = 0.35) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "#777777", linewidth = 0.35) +
  geom_errorbar(aes(ymin = pMST_lower, ymax = pMST_upper), width = 0.01, colour = "#777777", linewidth = 0.4) + # Quantile ranges across trees for each sample pair.
  geom_point(colour = "black", size = 1.1) +
  facet_wrap(~Dataset, ncol = 2, axes = "all", axis.labels = "all") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), expand = expansion(mult = 0.02)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), expand = expansion(mult = 0.02)) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "pMST from reference proxy tree", y = "Median pMST across\nrandomized proxy trees", tag = "(d)") + figure_theme
beta_limits <- range(c(-6, 6, pair_summary$betaNTI, pair_summary$beta_lower, pair_summary$beta_upper))
beta_breaks <- sort(unique(c(pretty(beta_limits, n = 5), -2, 0, 2)))
p_beta <- ggplot(pair_summary, aes(betaNTI, beta_median)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", colour = "#777777", linewidth = 0.35) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", colour = "#777777", linewidth = 0.35) +
  geom_errorbar(aes(ymin = beta_lower, ymax = beta_upper), width = 0.10, colour = "#777777", linewidth = 0.4) +
  geom_point(colour = "black", size = 1.1) +
  facet_wrap(~Dataset, ncol = 2, axes = "all", axis.labels = "all") +
  scale_x_continuous(breaks = beta_breaks, expand = expansion(mult = 0.02)) +
  scale_y_continuous(breaks = beta_breaks, expand = expansion(mult = 0.02)) +
  coord_fixed(xlim = beta_limits, ylim = beta_limits) +
  labs(x = "βNTI from reference proxy tree", y = "Median βNTI across\nrandomized proxy trees", tag = "(e)") + figure_theme
p_composition <- ggplot(composition, aes(Proportion, Y, fill = Process)) +
  geom_col(orientation = "y", width = 0.27, position = position_stack(reverse = TRUE)) + # For each dataset: reference tree above, mean across randomized trees below.
  scale_fill_manual(values = process_colours, drop = FALSE) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%"), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(dataset_labels), labels = rev(dataset_labels), expand = expansion(add = 0.5)) +
  labs(x = "Assembly proportion", y = NULL, tag = "(f)") + figure_theme +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))
agreement_min <- max(0, min(0.90, floor(min(agreement$Lower) * 20) / 20))
agreement_breaks <- sort(unique(c(agreement_min, 0.90, 0.95, 1)))
p_agreement <- ggplot(agreement, aes(Mean, Dataset, shape = Method, group = Method)) +
  geom_vline(xintercept = 0.95, linetype = "dashed", colour = "#777777", linewidth = 0.4) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper), orientation = "y", width = 0.12, position = position_dodge(width = 0.55), linewidth = 0.5) + # Quantile ranges of tree-level agreement.
  geom_point(position = position_dodge(width = 0.55), size = 3, colour = "black", fill = "white") +
  scale_shape_manual(values = c(pMST = 21, betaNTI = 22, Five = 23), labels = c(pMST = "pMST", betaNTI = "βNTI", Five = "Five-category")) +
  scale_x_continuous(breaks = agreement_breaks, labels = paste0(round(100 * agreement_breaks), "%"), expand = expansion(mult = 0.025)) +
  coord_cartesian(xlim = c(agreement_min, 1)) +
  labs(x = "Classification agreement", y = NULL, tag = "(g)") + figure_theme

# 4. Assemble and export
figure_name <- "Fig6defg"
panels <- list(p_mst, p_beta, p_composition, p_agreement)
positions <- rbind(c(1, 1), c(1, 2), c(2, 1), c(2, 2))
panel_width_mm <- c(100, 100, 180, 180)
panel_height_mm <- 100 # Each data region, including each facet, is 100 mm high.
ragg::agg_capture(width = 500, height = 500, units = "mm", res = 100)
grobs <- vector("list", length(panels))
width_mm <- height_mm <- numeric(length(panels))
for (j in seq_along(panels)) {
  g <- ggplotGrob(panels[[j]])
  cells <- g$layout[grepl("^panel", g$layout$name), ]
  g$widths[unique(cells$l)] <- unit(panel_width_mm[j], "mm")
  g$heights[unique(cells$t)] <- unit(panel_height_mm, "mm")
  grobs[[j]] <- g
  width_mm[j] <- convertWidth(sum(g$widths), "mm", valueOnly = TRUE)
  height_mm[j] <- convertHeight(sum(g$heights), "mm", valueOnly = TRUE)
}
column_widths <- as.numeric(tapply(width_mm, positions[, 2], max))
row_heights <- as.numeric(tapply(height_mm, positions[, 1], max))
layout_widths <- c(6, column_widths[1], 12, column_widths[2], 6)
layout_heights <- c(6, row_heights[1], 12, row_heights[2], 6)
layout <- gtable(widths = unit(layout_widths, "mm"), heights = unit(layout_heights, "mm"))
for (j in seq_along(grobs)) {
  cell <- grobTree(grobs[[j]], vp = viewport(width = unit(width_mm[j], "mm"), height = unit(height_mm[j], "mm")))
  layout <- gtable_add_grob(layout, cell, t = 2 * positions[j, 1], l = 2 * positions[j, 2], clip = "off")
}
canvas_width <- sum(layout_widths)
canvas_height <- sum(layout_heights)
dev.off()
svglite::svglite(paste0(figure_name, ".svg"), width = canvas_width / 25.4, height = canvas_height / 25.4, fix_text_size = FALSE)
grid.newpage()
grid.draw(layout)
dev.off()
ragg::agg_png(paste0(figure_name, ".png"), width = canvas_width, height = canvas_height, units = "mm", res = 300, background = "white")
grid.newpage()
grid.draw(layout)
dev.off()

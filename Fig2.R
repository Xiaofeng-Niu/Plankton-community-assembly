# Fig. 2: assembly-category proportions across methods and data sources
# Read pairwise results without rerunning null models.

library(ggplot2)
library(grid)
library(gtable)

# 1. Result files and panel assignments
inputs <- rbind(
  c("IdenPlantGenu_Asm_tMST_pairwise.csv", "a", "1", "Genus", "MIC", "tMST"),
  c("IdenPlantSpec_Asm_tMST_pairwise.csv", "a", "2", "Species", "MIC", "tMST"),
  c("SeqPlantGenu_Asm_tMST_pairwise.csv", "a", "3.5", "Genus", "SEQ", "tMST"),
  c("SeqPlantSpec_Asm_tMST_pairwise.csv", "a", "4.5", "Species", "SEQ", "tMST"),
  c("SeqPlantOtu_Asm_tMST_pairwise.csv", "a", "5.5", "OTU", "SEQ", "tMST"),
  c("IdenPlantGenu_Asm_ProxyTree_pMST_pairwise.csv", "b", "1", "Genus", "", "pMST"),
  c("IdenPlantSpec_Asm_ProxyTree_pMST_pairwise.csv", "b", "2", "Species", "", "pMST"),
  c("IdenPlantGenu_Asm_ProxyTree_betaNTI_RCbray_pairwise.csv", "c", "1", "Genus", "", "betaRC"),
  c("IdenPlantSpec_Asm_ProxyTree_betaNTI_RCbray_pairwise.csv", "c", "2", "Species", "", "betaRC"),
  c("IdenZooGenu_Asm_tMST_pairwise.csv", "d", "1", "Genus", "MIC", "tMST"),
  c("IdenZooSpec_Asm_tMST_pairwise.csv", "d", "2", "Species", "MIC", "tMST"),
  c("SeqZooGenu_Asm_tMST_pairwise.csv", "d", "3.5", "Genus", "SEQ", "tMST"),
  c("SeqZooSpec_Asm_tMST_pairwise.csv", "d", "4.5", "Species", "SEQ", "tMST"),
  c("SeqZooOtu_Asm_tMST_pairwise.csv", "d", "5.5", "OTU", "SEQ", "tMST"),
  c("IdenZooGenu_Asm_ProxyTree_pMST_pairwise.csv", "e", "1", "Genus", "", "pMST"),
  c("IdenZooSpec_Asm_ProxyTree_pMST_pairwise.csv", "e", "2", "Species", "", "pMST"),
  c("IdenZooGenu_Asm_ProxyTree_betaNTI_RCbray_pairwise.csv", "f", "1", "Genus", "", "betaRC"),
  c("IdenZooSpec_Asm_ProxyTree_betaNTI_RCbray_pairwise.csv", "f", "2", "Species", "", "betaRC")
)
colnames(inputs) <- c("File", "Panel", "X", "Label", "Block", "Metric")
missing_files <- unique(inputs[!file.exists(inputs[, "File"]), "File"])
if (length(missing_files)) stop(paste("缺少计算结果：", paste(missing_files, collapse = ", ")), call. = FALSE)

process_colours <- c(Stochastic = "#417654", Deterministic = "#376375", Undominated = "#C0D09D", "Dispersal limitation" = "#718457", "Homogenizing dispersal" = "#8FB794", "Heterogeneous selection" = "#94C5D1", "Homogeneous selection" = "#6C94B1")
process_levels <- names(process_colours)
process_labels <- c(Stochastic = "Stochastic process", Deterministic = "Deterministic process", Undominated = "Undominated", "Dispersal limitation" = "Dispersal limitation", "Homogenizing dispersal" = "Homogenizing dispersal", "Heterogeneous selection" = "Heterogeneous selection", "Homogeneous selection" = "Homogeneous selection")
box_colours <- c(OTU = "#5B83A6", Species = "#E99A41", Genus = "#69A363")
figure_theme <- theme_classic(base_family = "Arial") +
  theme(axis.text = element_text(size = 18, colour = "black"), axis.title = element_text(size = 20), plot.title = element_text(size = 22, face = "bold", hjust = 0.5), plot.subtitle = element_text(size = 18, hjust = 0.5), plot.tag = element_text(size = 20, face = "bold"), plot.tag.position = "topleft", strip.text = element_text(size = 18, face = "bold"), strip.background = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 20), legend.position = "none", plot.margin = margin(7, 7, 7, 7, "mm"), axis.ticks.length = unit(1.5, "mm"))

# 2. Summarize within-treatment sample pairs
plot_data <- data.frame()
pair_keys <- list()
for (k in seq_len(nrow(inputs))) {
  raw <- read.csv(inputs[k, "File"], check.names = FALSE)
  stopifnot(nrow(raw) > 0, all(c("group", "name1", "name2") %in% names(raw)))
  pair_key <- paste(pmin(raw$name1, raw$name2), pmax(raw$name1, raw$name2), sep = "|")
  stopifnot(!anyDuplicated(pair_key), all(raw$name1 != raw$name2))
  pair_keys[[inputs[k, "File"]]] <- sort(pair_key)
  metric <- inputs[k, "Metric"]
  if (metric == "tMST") {
    value <- raw$MST.ij.ruzicka
    process <- ifelse(value > 0.5, "Stochastic", "Deterministic")
    levels_used <- process_levels[1:2]
  } else if (metric == "pMST") {
    value <- raw$MST.ij.bMNTD
    process <- ifelse(value > 0.5, "Stochastic", "Deterministic")
    levels_used <- process_levels[1:2]
  } else {
    stopifnot(all(is.finite(raw$RCbray)))
    value <- raw$betaNTI
    process <- rep("Undominated", nrow(raw))
    process[value > 2] <- "Heterogeneous selection"
    process[value < -2] <- "Homogeneous selection"
    process[abs(value) <= 2 & raw$RCbray > 0.95] <- "Dispersal limitation"
    process[abs(value) <= 2 & raw$RCbray < -0.95] <- "Homogenizing dispersal"
    levels_used <- process_levels[3:7]
  }
  stopifnot(length(value) == nrow(raw), all(is.finite(value)))
  count <- table(factor(process, levels = levels_used))
  z <- data.frame(row.names = NULL, Panel = inputs[k, "Panel"], X = as.numeric(inputs[k, "X"]), Label = inputs[k, "Label"], Block = inputs[k, "Block"], Process = names(count), Count = as.numeric(count), Total = nrow(raw), Proportion = as.numeric(count) / nrow(raw))
  plot_data <- rbind(plot_data, z)
}
plot_data$Process <- factor(plot_data$Process, levels = process_levels)

legend_breaks <- c("Deterministic", "Stochastic", "Homogeneous selection", "Heterogeneous selection", "Homogenizing dispersal", "Dispersal limitation", "Undominated")
legend_rows <- 2

# 3. Draw six panels
d <- plot_data[plot_data$Panel == "a", ]
d <- d[order(d$X, d$Process), ]
d$Ymax <- ave(d$Proportion, d$X, FUN = cumsum)
d$Ymin <- d$Ymax - d$Proportion
bars <- d[!duplicated(d$X), c("X", "Label", "Block")]
blocks <- aggregate(X ~ Block, bars, range)
blocks$Left <- blocks$X[, 1] - 0.34
blocks$Right <- blocks$X[, 2] + 0.34
blocks <- blocks[nzchar(blocks$Block), , drop = FALSE]
p_a <- ggplot(d) +
  geom_hline(yintercept = 0.5, colour = "#D0D0D0", linewidth = 0.35, linetype = "dashed") +
  geom_rect(aes(xmin = X - 0.34, xmax = X + 0.34, ymin = Ymin, ymax = Ymax, fill = Process)) +
  geom_segment(data = blocks, aes(x = Left, xend = Right, y = -0.20, yend = -0.20), inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = blocks, aes(x = (Left + Right) / 2, y = -0.27, label = Block), inherit.aes = FALSE, family = "Arial", size = 18 / .pt) +
  scale_fill_manual(values = process_colours, limits = legend_breaks, breaks = legend_breaks, labels = process_labels[legend_breaks], drop = FALSE) +
  scale_x_continuous(breaks = bars$X, labels = bars$Label, expand = expansion(add = 0.55)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  labs(x = NULL, y = "Proportion of within-treatment sample pairs", title = "tMST · MIC + SEQ", subtitle = "Eukaryotic phytoplankton", tag = "(a)") +
  figure_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.margin = margin(8, 7, 34, 7, "mm"))

d <- plot_data[plot_data$Panel == "b", ]
d <- d[order(d$X, d$Process), ]
d$Ymax <- ave(d$Proportion, d$X, FUN = cumsum)
d$Ymin <- d$Ymax - d$Proportion
bars <- d[!duplicated(d$X), c("X", "Label", "Block")]
blocks <- aggregate(X ~ Block, bars, range)
blocks$Left <- blocks$X[, 1] - 0.34
blocks$Right <- blocks$X[, 2] + 0.34
blocks <- blocks[nzchar(blocks$Block), , drop = FALSE]
p_b <- ggplot(d) +
  geom_hline(yintercept = 0.5, colour = "#D0D0D0", linewidth = 0.35, linetype = "dashed") +
  geom_rect(aes(xmin = X - 0.34, xmax = X + 0.34, ymin = Ymin, ymax = Ymax, fill = Process)) +
  geom_segment(data = blocks, aes(x = Left, xend = Right, y = -0.14, yend = -0.14), inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = blocks, aes(x = (Left + Right) / 2, y = -0.21, label = Block), inherit.aes = FALSE, family = "Arial", size = 18 / .pt) +
  scale_fill_manual(values = process_colours, limits = legend_breaks, breaks = legend_breaks, labels = process_labels[legend_breaks], drop = FALSE) +
  scale_x_continuous(breaks = bars$X, labels = bars$Label, expand = expansion(add = 0.55)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  labs(x = NULL, y = NULL, title = "pMST · MIC", subtitle = "Proxy tree", tag = "(b)") +
  figure_theme + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), plot.margin = margin(8, 7, 34, 7, "mm"))

d <- plot_data[plot_data$Panel == "c", ]
d <- d[order(d$X, d$Process), ]
d$Ymax <- ave(d$Proportion, d$X, FUN = cumsum)
d$Ymin <- d$Ymax - d$Proportion
bars <- d[!duplicated(d$X), c("X", "Label", "Block")]
blocks <- aggregate(X ~ Block, bars, range)
blocks$Left <- blocks$X[, 1] - 0.34
blocks$Right <- blocks$X[, 2] + 0.34
blocks <- blocks[nzchar(blocks$Block), , drop = FALSE]
p_c <- ggplot(d) +
  geom_hline(yintercept = 0.5, colour = "#D0D0D0", linewidth = 0.35, linetype = "dashed") +
  geom_rect(aes(xmin = X - 0.34, xmax = X + 0.34, ymin = Ymin, ymax = Ymax, fill = Process)) +
  geom_segment(data = blocks, aes(x = Left, xend = Right, y = -0.14, yend = -0.14), inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = blocks, aes(x = (Left + Right) / 2, y = -0.21, label = Block), inherit.aes = FALSE, family = "Arial", size = 18 / .pt) +
  scale_fill_manual(values = process_colours, limits = legend_breaks, breaks = legend_breaks, labels = process_labels[legend_breaks], drop = FALSE) +
  scale_x_continuous(breaks = bars$X, labels = bars$Label, expand = expansion(add = 0.55)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  labs(x = NULL, y = NULL, title = "βNTI + RCbray · MIC", subtitle = "Proxy tree", tag = "(c)") +
  figure_theme + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), plot.margin = margin(8, 7, 34, 7, "mm"))

d <- plot_data[plot_data$Panel == "d", ]
d <- d[order(d$X, d$Process), ]
d$Ymax <- ave(d$Proportion, d$X, FUN = cumsum)
d$Ymin <- d$Ymax - d$Proportion
bars <- d[!duplicated(d$X), c("X", "Label", "Block")]
blocks <- aggregate(X ~ Block, bars, range)
blocks$Left <- blocks$X[, 1] - 0.34
blocks$Right <- blocks$X[, 2] + 0.34
blocks <- blocks[nzchar(blocks$Block), , drop = FALSE]
p_d <- ggplot(d) +
  geom_hline(yintercept = 0.5, colour = "#D0D0D0", linewidth = 0.35, linetype = "dashed") +
  geom_rect(aes(xmin = X - 0.34, xmax = X + 0.34, ymin = Ymin, ymax = Ymax, fill = Process)) +
  geom_segment(data = blocks, aes(x = Left, xend = Right, y = -0.20, yend = -0.20), inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = blocks, aes(x = (Left + Right) / 2, y = -0.27, label = Block), inherit.aes = FALSE, family = "Arial", size = 18 / .pt) +
  scale_fill_manual(values = process_colours, limits = legend_breaks, breaks = legend_breaks, labels = process_labels[legend_breaks], drop = FALSE) +
  scale_x_continuous(breaks = bars$X, labels = bars$Label, expand = expansion(add = 0.55)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  labs(x = NULL, y = "Proportion of within-treatment sample pairs", title = "tMST · MIC + SEQ", subtitle = "Zooplankton", tag = "(d)") +
  figure_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.margin = margin(8, 7, 34, 7, "mm"))

d <- plot_data[plot_data$Panel == "e", ]
d <- d[order(d$X, d$Process), ]
d$Ymax <- ave(d$Proportion, d$X, FUN = cumsum)
d$Ymin <- d$Ymax - d$Proportion
bars <- d[!duplicated(d$X), c("X", "Label", "Block")]
blocks <- aggregate(X ~ Block, bars, range)
blocks$Left <- blocks$X[, 1] - 0.34
blocks$Right <- blocks$X[, 2] + 0.34
blocks <- blocks[nzchar(blocks$Block), , drop = FALSE]
p_e <- ggplot(d) +
  geom_hline(yintercept = 0.5, colour = "#D0D0D0", linewidth = 0.35, linetype = "dashed") +
  geom_rect(aes(xmin = X - 0.34, xmax = X + 0.34, ymin = Ymin, ymax = Ymax, fill = Process)) +
  geom_segment(data = blocks, aes(x = Left, xend = Right, y = -0.14, yend = -0.14), inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = blocks, aes(x = (Left + Right) / 2, y = -0.21, label = Block), inherit.aes = FALSE, family = "Arial", size = 18 / .pt) +
  scale_fill_manual(values = process_colours, limits = legend_breaks, breaks = legend_breaks, labels = process_labels[legend_breaks], drop = FALSE) +
  scale_x_continuous(breaks = bars$X, labels = bars$Label, expand = expansion(add = 0.55)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  labs(x = NULL, y = NULL, title = "pMST · MIC", subtitle = "Proxy tree", tag = "(e)") +
  figure_theme + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), plot.margin = margin(8, 7, 34, 7, "mm"))

d <- plot_data[plot_data$Panel == "f", ]
d <- d[order(d$X, d$Process), ]
d$Ymax <- ave(d$Proportion, d$X, FUN = cumsum)
d$Ymin <- d$Ymax - d$Proportion
bars <- d[!duplicated(d$X), c("X", "Label", "Block")]
blocks <- aggregate(X ~ Block, bars, range)
blocks$Left <- blocks$X[, 1] - 0.34
blocks$Right <- blocks$X[, 2] + 0.34
blocks <- blocks[nzchar(blocks$Block), , drop = FALSE]
p_f <- ggplot(d) +
  geom_hline(yintercept = 0.5, colour = "#D0D0D0", linewidth = 0.35, linetype = "dashed") +
  geom_rect(aes(xmin = X - 0.34, xmax = X + 0.34, ymin = Ymin, ymax = Ymax, fill = Process)) +
  geom_segment(data = blocks, aes(x = Left, xend = Right, y = -0.14, yend = -0.14), inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = blocks, aes(x = (Left + Right) / 2, y = -0.21, label = Block), inherit.aes = FALSE, family = "Arial", size = 18 / .pt) +
  scale_fill_manual(values = process_colours, limits = legend_breaks, breaks = legend_breaks, labels = process_labels[legend_breaks], drop = FALSE) +
  scale_x_continuous(breaks = bars$X, labels = bars$Label, expand = expansion(add = 0.55)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  labs(x = NULL, y = NULL, title = "βNTI + RCbray · MIC", subtitle = "Proxy tree", tag = "(f)") +
  figure_theme + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), plot.margin = margin(8, 7, 34, 7, "mm"))

figure_name <- "Fig2"
panels <- list(p_a, p_b, p_c, p_d, p_e, p_f)
positions <- rbind(c(1, 1), c(1, 2), c(1, 3), c(2, 1), c(2, 2), c(2, 3))
panel_width_mm <- c(135, 75, 75, 135, 75, 75)
panel_height_mm <- c(100, 100, 100, 100, 100, 100) # Data-region height, excluding titles, axes and legends.
ragg::agg_capture(width = 300, height = 300, units = "mm", res = 100)
grobs <- vector("list", length(panels))
width_mm <- height_mm <- numeric(length(panels))
for (j in seq_along(panels)) {
  g <- ggplotGrob(panels[[j]])
  panel_cells <- g$layout[grepl("^panel", g$layout$name), ]
  g$widths[unique(panel_cells$l)] <- unit(panel_width_mm[j] / length(unique(panel_cells$l)), "mm")
  g$heights[unique(panel_cells$t)] <- unit(panel_height_mm[j] / length(unique(panel_cells$t)), "mm")
  grobs[[j]] <- g
  width_mm[j] <- convertWidth(sum(g$widths), "mm", valueOnly = TRUE)
  height_mm[j] <- convertHeight(sum(g$heights), "mm", valueOnly = TRUE)
}
legend_plot <- panels[[1]] + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = legend_rows, byrow = TRUE))
legend_table <- ggplotGrob(legend_plot)
legend_grob <- legend_table$grobs[[which(legend_table$layout$name == "guide-box-bottom")]]
legend_width <- convertWidth(sum(legend_grob$widths), "mm", valueOnly = TRUE)
legend_height <- convertHeight(sum(legend_grob$heights), "mm", valueOnly = TRUE)
column_widths <- as.numeric(tapply(width_mm, positions[, 2], max))
row_heights <- as.numeric(tapply(height_mm, positions[, 1], max))
column_count <- length(column_widths)
row_count <- length(row_heights)
layout_widths <- c(6, head(as.vector(rbind(column_widths, rep(8, column_count))), -1), 6)
layout_heights <- c(6, head(as.vector(rbind(row_heights, rep(10, row_count))), -1), 8, legend_height, 6)
extra_width <- max(0, legend_width + 12 - sum(layout_widths))
layout_widths[c(1, length(layout_widths))] <- layout_widths[c(1, length(layout_widths))] + extra_width / 2
layout <- gtable(widths = unit(layout_widths, "mm"), heights = unit(layout_heights, "mm"))
for (j in seq_along(grobs)) {
  cell <- grobTree(grobs[[j]], vp = viewport(width = unit(width_mm[j], "mm"), height = unit(height_mm[j], "mm")))
  layout <- gtable_add_grob(layout, cell, t = 2 * positions[j, 1], l = 2 * positions[j, 2], clip = "off")
}
layout <- gtable_add_grob(layout, legend_grob, t = 2 * row_count + 2, l = 2, r = 2 * column_count, clip = "off")
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

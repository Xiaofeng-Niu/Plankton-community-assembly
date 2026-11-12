# Fig. S2: distributions of pairwise tMST and pMST
# Use the data-release folder as the working directory; read pairwise results without rerunning null models.

library(ggplot2)
library(grid)
library(gtable)

# 1. Result files and panel assignments
inputs <- rbind(
  c("IdenPlantGenu_Asm_tMST_pairwise.csv", "a", "Eukaryotic phytoplankton", "Genus", "tMST"),
  c("IdenPlantSpec_Asm_tMST_pairwise.csv", "a", "Eukaryotic phytoplankton", "Species", "tMST"),
  c("IdenZooGenu_Asm_tMST_pairwise.csv", "a", "Zooplankton", "Genus", "tMST"),
  c("IdenZooSpec_Asm_tMST_pairwise.csv", "a", "Zooplankton", "Species", "tMST"),
  c("SeqPlantOtu_Asm_tMST_pairwise.csv", "b", "Eukaryotic phytoplankton", "OTU", "tMST"),
  c("SeqPlantSpec_Asm_tMST_pairwise.csv", "b", "Eukaryotic phytoplankton", "Species", "tMST"),
  c("SeqPlantGenu_Asm_tMST_pairwise.csv", "b", "Eukaryotic phytoplankton", "Genus", "tMST"),
  c("SeqZooOtu_Asm_tMST_pairwise.csv", "b", "Zooplankton", "OTU", "tMST"),
  c("SeqZooSpec_Asm_tMST_pairwise.csv", "b", "Zooplankton", "Species", "tMST"),
  c("SeqZooGenu_Asm_tMST_pairwise.csv", "b", "Zooplankton", "Genus", "tMST"),
  c("IdenPlantGenu_Asm_ProxyTree_pMST_pairwise.csv", "c", "Eukaryotic phytoplankton", "Genus", "pMST"),
  c("IdenPlantSpec_Asm_ProxyTree_pMST_pairwise.csv", "c", "Eukaryotic phytoplankton", "Species", "pMST"),
  c("IdenZooGenu_Asm_ProxyTree_pMST_pairwise.csv", "c", "Zooplankton", "Genus", "pMST"),
  c("IdenZooSpec_Asm_ProxyTree_pMST_pairwise.csv", "c", "Zooplankton", "Species", "pMST"),
  c("SeqPlantOtu_Asm_SeqTree_pMST_pairwise.csv", "d", "Eukaryotic phytoplankton", "OTU", "pMST"),
  c("SeqPlantSpec_Asm_SeqTree_pMST_pairwise.csv", "d", "Eukaryotic phytoplankton", "Species", "pMST"),
  c("SeqPlantGenu_Asm_SeqTree_pMST_pairwise.csv", "d", "Eukaryotic phytoplankton", "Genus", "pMST"),
  c("SeqZooOtu_Asm_SeqTree_pMST_pairwise.csv", "d", "Zooplankton", "OTU", "pMST"),
  c("SeqZooSpec_Asm_SeqTree_pMST_pairwise.csv", "d", "Zooplankton", "Species", "pMST"),
  c("SeqZooGenu_Asm_SeqTree_pMST_pairwise.csv", "d", "Zooplankton", "Genus", "pMST"),
  c("SeqPlantSpec_Asm_ProxyTree_pMST_pairwise.csv", "e", "Eukaryotic phytoplankton", "Species", "pMST"),
  c("SeqPlantGenu_Asm_ProxyTree_pMST_pairwise.csv", "e", "Eukaryotic phytoplankton", "Genus", "pMST"),
  c("SeqZooSpec_Asm_ProxyTree_pMST_pairwise.csv", "e", "Zooplankton", "Species", "pMST"),
  c("SeqZooGenu_Asm_ProxyTree_pMST_pairwise.csv", "e", "Zooplankton", "Genus", "pMST")
)
colnames(inputs) <- c("File", "Panel", "Community", "Resolution", "Metric")
missing_files <- unique(inputs[!file.exists(inputs[, "File"]), "File"])
if (length(missing_files)) stop(paste("缺少计算结果：", paste(missing_files, collapse = ", ")), call. = FALSE)

process_colours <- c(Stochastic = "#417654", Deterministic = "#376375", Undominated = "#C0D09D", "Dispersal limitation" = "#718457", "Homogenizing dispersal" = "#8FB794", "Heterogeneous selection" = "#94C5D1", "Homogeneous selection" = "#6C94B1")
process_levels <- names(process_colours)
process_labels <- c(Stochastic = "Stochastic process", Deterministic = "Deterministic process", Undominated = "Undominated", "Dispersal limitation" = "Dispersal limitation", "Homogenizing dispersal" = "Homogenizing dispersal", "Heterogeneous selection" = "Heterogeneous selection", "Homogeneous selection" = "Homogeneous selection")
box_colours <- c(OTU = "#5B83A6", Species = "#E99A41", Genus = "#69A363")
figure_theme <- theme_classic(base_family = "Arial") +
  theme(axis.text = element_text(size = 18, colour = "black"), axis.title = element_text(size = 20), plot.title = element_text(size = 22, face = "bold", hjust = 0.5), plot.subtitle = element_text(size = 18, hjust = 0.5), plot.tag = element_text(size = 20, face = "bold"), plot.tag.position = "topleft", strip.text = element_text(size = 18, face = "bold"), strip.background = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 20), legend.position = "none", plot.margin = margin(7, 7, 7, 7, "mm"), axis.ticks.length = unit(1.5, "mm"))

# 2. Extract pairwise metrics
plot_data <- data.frame()
for (k in seq_len(nrow(inputs))) {
  raw <- read.csv(inputs[k, "File"], check.names = FALSE)
  stopifnot(nrow(raw) > 0, all(c("group", "name1", "name2") %in% names(raw)))
  pair_key <- paste(pmin(raw$name1, raw$name2), pmax(raw$name1, raw$name2), sep = "|")
  stopifnot(!anyDuplicated(pair_key), all(raw$name1 != raw$name2))
  metric <- inputs[k, "Metric"]
  if (metric == "tMST") value <- raw$MST.ij.ruzicka
  if (metric == "pMST") value <- raw$MST.ij.bMNTD
  if (metric == "betaNTI") value <- raw$betaNTI
  stopifnot(length(value) == nrow(raw), all(is.finite(value)))
  z <- data.frame(row.names = NULL, Panel = inputs[k, "Panel"], Community = inputs[k, "Community"], Resolution = inputs[k, "Resolution"], Value = value, Pair = pair_key)
  plot_data <- rbind(plot_data, z)
}
plot_data$Community <- factor(plot_data$Community, levels = c("Eukaryotic phytoplankton", "Zooplankton"))
plot_data$Resolution <- factor(plot_data$Resolution, levels = c("OTU", "Species", "Genus"))
sample_counts <- aggregate(Pair ~ Panel + Community + Resolution, plot_data, length)
legend_breaks <- c("OTU", "Species", "Genus")
legend_rows <- 1

d <- plot_data[plot_data$Panel == "a", ]
p_a <- ggplot(d, aes(Resolution, Value)) +
  geom_rect(data = data.frame(Community = factor("Zooplankton", levels = levels(plot_data$Community))), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "#F0F1F2") +
  geom_boxplot(aes(fill = Resolution), width = 0.55, outlier.shape = NA, linewidth = 0.4, alpha = 0.85) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "#777777", linewidth = 0.4) +
  geom_point(colour = "black", shape = 16, size = 0.65, position = position_jitter(width = 0.16, height = 0, seed = 123)) +
  scale_fill_manual(values = box_colours, limits = legend_breaks, breaks = legend_breaks, drop = FALSE) +
  scale_x_discrete(limits = c("Species", "Genus"), drop = FALSE) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1.03)) +
  facet_grid(. ~ Community, labeller = as_labeller(c("Eukaryotic phytoplankton" = "Eukaryotic
phytoplankton", "Zooplankton" = "Zooplankton"))) +
  labs(x = NULL, y = "tMST", title = "Microscopy dataset", tag = "(a)") +
  figure_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0, "mm"), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4), axis.line = element_blank())

d <- plot_data[plot_data$Panel == "b", ]
p_b <- ggplot(d, aes(Resolution, Value)) +
  geom_rect(data = data.frame(Community = factor("Zooplankton", levels = levels(plot_data$Community))), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "#F0F1F2") +
  geom_boxplot(aes(fill = Resolution), width = 0.55, outlier.shape = NA, linewidth = 0.4, alpha = 0.85) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "#777777", linewidth = 0.4) +
  geom_point(colour = "black", shape = 16, size = 0.65, position = position_jitter(width = 0.16, height = 0, seed = 123)) +
  scale_fill_manual(values = box_colours, limits = legend_breaks, breaks = legend_breaks, drop = FALSE) +
  scale_x_discrete(limits = c("OTU", "Species", "Genus"), drop = FALSE) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1.03)) +
  facet_grid(. ~ Community, labeller = as_labeller(c("Eukaryotic phytoplankton" = "Eukaryotic
phytoplankton", "Zooplankton" = "Zooplankton"))) +
  labs(x = NULL, y = NULL, title = "Sequencing dataset", tag = "(b)") +
  figure_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0, "mm"), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4), axis.line = element_blank())

d <- plot_data[plot_data$Panel == "c", ]
p_c <- ggplot(d, aes(Resolution, Value)) +
  geom_rect(data = data.frame(Community = factor("Zooplankton", levels = levels(plot_data$Community))), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "#F0F1F2") +
  geom_boxplot(aes(fill = Resolution), width = 0.55, outlier.shape = NA, linewidth = 0.4, alpha = 0.85) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "#777777", linewidth = 0.4) +
  geom_point(colour = "black", shape = 16, size = 0.65, position = position_jitter(width = 0.16, height = 0, seed = 123)) +
  scale_fill_manual(values = box_colours, limits = legend_breaks, breaks = legend_breaks, drop = FALSE) +
  scale_x_discrete(limits = c("Species", "Genus"), drop = FALSE) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1.03)) +
  facet_grid(. ~ Community, labeller = as_labeller(c("Eukaryotic phytoplankton" = "Eukaryotic
phytoplankton", "Zooplankton" = "Zooplankton"))) +
  labs(x = NULL, y = "pMST", title = "Microscopy dataset\nProxy tree", tag = "(c)") +
  figure_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0, "mm"), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4), axis.line = element_blank())

d <- plot_data[plot_data$Panel == "d", ]
p_d <- ggplot(d, aes(Resolution, Value)) +
  geom_rect(data = data.frame(Community = factor("Zooplankton", levels = levels(plot_data$Community))), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "#F0F1F2") +
  geom_boxplot(aes(fill = Resolution), width = 0.55, outlier.shape = NA, linewidth = 0.4, alpha = 0.85) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "#777777", linewidth = 0.4) +
  geom_point(colour = "black", shape = 16, size = 0.65, position = position_jitter(width = 0.16, height = 0, seed = 123)) +
  scale_fill_manual(values = box_colours, limits = legend_breaks, breaks = legend_breaks, drop = FALSE) +
  scale_x_discrete(limits = c("OTU", "Species", "Genus"), drop = FALSE) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1.03)) +
  facet_grid(. ~ Community, labeller = as_labeller(c("Eukaryotic phytoplankton" = "Eukaryotic
phytoplankton", "Zooplankton" = "Zooplankton"))) +
  labs(x = NULL, y = NULL, title = "Sequencing dataset\nSequence-based tree", tag = "(d)") +
  figure_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0, "mm"), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4), axis.line = element_blank())

d <- plot_data[plot_data$Panel == "e", ]
p_e <- ggplot(d, aes(Resolution, Value)) +
  geom_rect(data = data.frame(Community = factor("Zooplankton", levels = levels(plot_data$Community))), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "#F0F1F2") +
  geom_boxplot(aes(fill = Resolution), width = 0.55, outlier.shape = NA, linewidth = 0.4, alpha = 0.85) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "#777777", linewidth = 0.4) +
  geom_point(colour = "black", shape = 16, size = 0.65, position = position_jitter(width = 0.16, height = 0, seed = 123)) +
  scale_fill_manual(values = box_colours, limits = legend_breaks, breaks = legend_breaks, drop = FALSE) +
  scale_x_discrete(limits = c("Species", "Genus"), drop = FALSE) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1.03)) +
  facet_grid(. ~ Community, labeller = as_labeller(c("Eukaryotic phytoplankton" = "Eukaryotic
phytoplankton", "Zooplankton" = "Zooplankton"))) +
  labs(x = NULL, y = NULL, title = "Sequencing dataset\nProxy tree", tag = "(e)") +
  figure_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0, "mm"), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4), axis.line = element_blank())

figure_name <- "FigS2"
panels <- list(p_a, p_b, p_c, p_d, p_e)
positions <- rbind(c(1, 1), c(1, 2), c(2, 1), c(2, 2), c(2, 3))
panel_width_mm <- c(125, 125, 125, 125, 125)
panel_height_mm <- c(100, 100, 100, 100, 100) # Data-region height, excluding titles, axes and legends.
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

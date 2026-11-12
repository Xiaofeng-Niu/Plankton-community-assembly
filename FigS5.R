# Fig. S5: assembly metrics for microscopy phytoplankton datasets including cyanobacteria
# Use the data-release folder as the working directory; read pairwise results without rerunning null models.

library(ggplot2)
library(grid)
library(gtable)

# 1. Result files and panel assignments
inputs <- rbind(
  c("IdenPlantGenu_AsmCy_tMST_pairwise.csv", "a", "Phytoplankton including cyanobacteria", "Genus", "tMST"),
  c("IdenPlantGenu_AsmCy_ProxyTree_pMST_pairwise.csv", "b", "Phytoplankton including cyanobacteria", "Genus", "pMST"),
  c("IdenPlantGenu_AsmCy_ProxyTree_betaNTI_RCbray_pairwise.csv", "c", "Phytoplankton including cyanobacteria", "Genus", "betaNTI"),
  c("IdenPlantSpec_AsmCy_tMST_pairwise.csv", "a", "Phytoplankton including cyanobacteria", "Species", "tMST"),
  c("IdenPlantSpec_AsmCy_ProxyTree_pMST_pairwise.csv", "b", "Phytoplankton including cyanobacteria", "Species", "pMST"),
  c("IdenPlantSpec_AsmCy_ProxyTree_betaNTI_RCbray_pairwise.csv", "c", "Phytoplankton including cyanobacteria", "Species", "betaNTI")
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
  if (metric == "betaNTI") {
    keep <- abs(raw$betaNTI) <= 2
    stopifnot(all(is.finite(raw$RCbray)))
    rc_panel <- c(c = "d")[inputs[k, "Panel"]]
    rc <- data.frame(row.names = NULL, Panel = rep(unname(rc_panel), sum(keep)), Community = rep(inputs[k, "Community"], sum(keep)), Resolution = rep(inputs[k, "Resolution"], sum(keep)), Value = raw$RCbray[keep], Pair = pair_key[keep])
    plot_data <- rbind(plot_data, rc)
  }
}
plot_data$Community <- factor(plot_data$Community, levels = "Phytoplankton including cyanobacteria")
plot_data$Resolution <- factor(plot_data$Resolution, levels = c("OTU", "Species", "Genus"))
sample_counts <- aggregate(Pair ~ Panel + Community + Resolution, plot_data, length)
legend_breaks <- c("OTU", "Species", "Genus")
legend_rows <- 1

legend_breaks <- c("Genus", "Species")
box_colours <- c(Genus = "#64A35C", Species = "#ED9838")
beta_limits <- range(c(-2.8, 2.8, plot_data$Value[plot_data$Panel == "c"]))
beta_limits <- beta_limits + c(-1, 1) * diff(beta_limits) * 0.03
beta_breaks <- c(-2, 0, 2)

d <- plot_data[plot_data$Panel == "a", ]
p_a <- ggplot(d, aes(Resolution, Value)) +
  geom_boxplot(aes(fill = Resolution), width = 0.55, outlier.shape = NA, linewidth = 0.4, alpha = 0.85) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "#777777", linewidth = 0.4) +
  geom_point(colour = "black", shape = 16, size = 0.65, position = position_jitter(width = 0.16, height = 0, seed = 123)) +
  scale_fill_manual(values = box_colours, limits = legend_breaks, breaks = legend_breaks, drop = FALSE) +
  scale_x_discrete(limits = c("Genus", "Species"), drop = FALSE) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1.03)) +
  labs(x = NULL, y = "tMST", title = "tMST", tag = "(a)") +
  figure_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0, "mm"), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4), axis.line = element_blank())

d <- plot_data[plot_data$Panel == "b", ]
p_b <- ggplot(d, aes(Resolution, Value)) +
  geom_boxplot(aes(fill = Resolution), width = 0.55, outlier.shape = NA, linewidth = 0.4, alpha = 0.85) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "#777777", linewidth = 0.4) +
  geom_point(colour = "black", shape = 16, size = 0.65, position = position_jitter(width = 0.16, height = 0, seed = 123)) +
  scale_fill_manual(values = box_colours, limits = legend_breaks, breaks = legend_breaks, drop = FALSE) +
  scale_x_discrete(limits = c("Genus", "Species"), drop = FALSE) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1.03)) +
  labs(x = NULL, y = "pMST", title = "pMST · Proxy tree", tag = "(b)") +
  figure_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0, "mm"), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4), axis.line = element_blank())

d <- plot_data[plot_data$Panel == "c", ]
p_c <- ggplot(d, aes(Resolution, Value)) +
  geom_boxplot(aes(fill = Resolution), width = 0.55, outlier.shape = NA, linewidth = 0.4, alpha = 0.85) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", colour = "#777777", linewidth = 0.4) +
  geom_point(colour = "black", shape = 16, size = 0.65, position = position_jitter(width = 0.16, height = 0, seed = 123)) +
  scale_fill_manual(values = box_colours, limits = legend_breaks, breaks = legend_breaks, drop = FALSE) +
  scale_x_discrete(limits = c("Genus", "Species"), drop = FALSE) +
  scale_y_continuous(breaks = beta_breaks, expand = c(0, 0)) +
  coord_cartesian(ylim = beta_limits) +
  labs(x = NULL, y = "βNTI", title = "βNTI · Proxy tree", tag = "(c)") +
  figure_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0, "mm"), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4), axis.line = element_blank())

d <- plot_data[plot_data$Panel == "d", ]
p_d <- ggplot(d, aes(Resolution, Value)) +
  geom_boxplot(aes(fill = Resolution), width = 0.55, outlier.shape = NA, linewidth = 0.4, alpha = 0.85) +
  geom_hline(yintercept = c(-0.95, 0.95), linetype = "dashed", colour = "#777777", linewidth = 0.4) +
  geom_point(colour = "black", shape = 16, size = 0.65, position = position_jitter(width = 0.16, height = 0, seed = 123)) +
  scale_fill_manual(values = box_colours, limits = legend_breaks, breaks = legend_breaks, drop = FALSE) +
  scale_x_discrete(limits = c("Genus", "Species"), drop = FALSE) +
  scale_y_continuous(breaks = c(-0.95, 0, 0.95), expand = c(0, 0)) +
  coord_cartesian(ylim = c(-1.03, 1.03)) +
  labs(x = NULL, y = "RCbray", title = "Conditional RCbray", tag = "(d)") +
  figure_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0, "mm"), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4), axis.line = element_blank())

figure_name <- "FigS5"
panels <- list(p_a, p_b, p_c, p_d)
positions <- rbind(c(1, 1), c(1, 2), c(2, 1), c(2, 2))
panel_width_mm <- c(100, 100, 100, 100)
panel_height_mm <- c(100, 100, 100, 100) # Data-region height, excluding titles, axes and legends.
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

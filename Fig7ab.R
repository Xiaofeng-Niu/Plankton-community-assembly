# Fig. 7a-b: five-process proportions within each treatment
# Use the data-release folder as the working directory; read pairwise results without rerunning null models.

library(ggplot2)
library(grid)
library(gtable)

# 1. Result files and panel assignments
inputs <- rbind(
  c("IdenPlantGenu_Asm_ProxyTree_betaNTI_RCbray_pairwise.csv", "Eukaryotic phytoplankton", "Genus", "MIC"),
  c("IdenPlantSpec_Asm_ProxyTree_betaNTI_RCbray_pairwise.csv", "Eukaryotic phytoplankton", "Species", "MIC"),
  c("IdenZooGenu_Asm_ProxyTree_betaNTI_RCbray_pairwise.csv", "Zooplankton", "Genus", "MIC"),
  c("IdenZooSpec_Asm_ProxyTree_betaNTI_RCbray_pairwise.csv", "Zooplankton", "Species", "MIC"),
  c("SeqPlantOtu_Asm_SeqTree_betaNTI_RCbray_pairwise.csv", "Eukaryotic phytoplankton", "OTU", "SEQ"),
  c("SeqPlantSpec_Asm_SeqTree_betaNTI_RCbray_pairwise.csv", "Eukaryotic phytoplankton", "Species", "SEQ"),
  c("SeqPlantGenu_Asm_SeqTree_betaNTI_RCbray_pairwise.csv", "Eukaryotic phytoplankton", "Genus", "SEQ"),
  c("SeqZooOtu_Asm_SeqTree_betaNTI_RCbray_pairwise.csv", "Zooplankton", "OTU", "SEQ"),
  c("SeqZooSpec_Asm_SeqTree_betaNTI_RCbray_pairwise.csv", "Zooplankton", "Species", "SEQ"),
  c("SeqZooGenu_Asm_SeqTree_betaNTI_RCbray_pairwise.csv", "Zooplankton", "Genus", "SEQ")
)
colnames(inputs) <- c("File", "Community", "Resolution", "Source")
missing_files <- unique(inputs[!file.exists(inputs[, "File"]), "File"])
if (length(missing_files)) stop(paste("缺少计算结果：", paste(missing_files, collapse = ", ")), call. = FALSE)

process_colours <- c(Stochastic = "#417654", Deterministic = "#376375", Undominated = "#C0D09D", "Dispersal limitation" = "#718457", "Homogenizing dispersal" = "#8FB794", "Heterogeneous selection" = "#94C5D1", "Homogeneous selection" = "#6C94B1")
process_levels <- names(process_colours)
process_labels <- c(Stochastic = "Stochastic process", Deterministic = "Deterministic process", Undominated = "Undominated", "Dispersal limitation" = "Dispersal limitation", "Homogenizing dispersal" = "Homogenizing dispersal", "Heterogeneous selection" = "Heterogeneous selection", "Homogeneous selection" = "Homogeneous selection")
box_colours <- c(OTU = "#5B83A6", Species = "#E99A41", Genus = "#69A363")
figure_theme <- theme_classic(base_family = "Arial") +
  theme(axis.text = element_text(size = 18, colour = "black"), axis.title = element_text(size = 20), plot.title = element_text(size = 22, face = "bold", hjust = 0.5), plot.subtitle = element_text(size = 18, hjust = 0.5), plot.tag = element_text(size = 20, face = "bold"), plot.tag.position = "topleft", strip.text = element_text(size = 18, face = "bold"), strip.background = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 20), legend.position = "none", plot.margin = margin(7, 7, 7, 7, "mm"), axis.ticks.length = unit(1.5, "mm"))

treatments <- c("C", "E", "P", "EP", "W", "WE", "WP", "WEP")
mirror_levels <- c("Homogeneous selection", "Heterogeneous selection", "Homogenizing dispersal", "Dispersal limitation", "Undominated")
legend_breaks <- mirror_levels
legend_rows <- 2
plot_data <- data.frame()

# 2. Summarize each dataset and treatment separately without intersecting MIC and SEQ pair sets
for (k in seq_len(nrow(inputs))) {
  raw <- read.csv(inputs[k, "File"], check.names = FALSE)
  stopifnot(nrow(raw) > 0, all(c("group", "name1", "name2") %in% names(raw)))
  pair_key <- paste(pmin(raw$name1, raw$name2), pmax(raw$name1, raw$name2), sep = "|")
  stopifnot(!anyDuplicated(pair_key), all(raw$name1 != raw$name2))
  stopifnot(all(is.finite(raw$betaNTI)), all(is.finite(raw$RCbray)), all(raw$group %in% treatments))
  process <- rep("Undominated", nrow(raw))
  process[raw$betaNTI > 2] <- "Heterogeneous selection"
  process[raw$betaNTI < -2] <- "Homogeneous selection"
  process[abs(raw$betaNTI) <= 2 & raw$RCbray > 0.95] <- "Dispersal limitation"
  process[abs(raw$betaNTI) <= 2 & raw$RCbray < -0.95] <- "Homogenizing dispersal"
  count <- as.data.frame(table(Treatment = factor(raw$group, levels = treatments), Process = factor(process, levels = mirror_levels)))
  names(count)[3] <- "Count"
  count$Total <- ave(count$Count, count$Treatment, FUN = sum)
  stopifnot(all(count$Total > 0))
  count$Percent <- 100 * count$Count / count$Total
  count$Community <- inputs[k, "Community"]
  count$Resolution <- inputs[k, "Resolution"]
  count$Source <- inputs[k, "Source"]
  count <- count[order(count$Treatment, count$Process), ]
  count$High <- ave(count$Percent, count$Treatment, FUN = cumsum)
  count$Low <- count$High - count$Percent
  count$Left <- ifelse(count$Source == "MIC", -count$High, count$Low)
  count$Right <- ifelse(count$Source == "MIC", -count$Low, count$High)
  count$Y <- ifelse(count$Resolution == "Genus", 27 - match(count$Treatment, treatments), ifelse(count$Resolution == "Species", 18 - match(count$Treatment, treatments), 9 - match(count$Treatment, treatments)))
  plot_data <- rbind(plot_data, count)
}
plot_data$Process <- factor(plot_data$Process, levels = mirror_levels)
sample_counts <- unique(plot_data[, c("Community", "Source", "Resolution", "Treatment", "Total")])

d <- plot_data[plot_data$Community == "Eukaryotic phytoplankton", ]
p_a <- ggplot(d, aes(fill = Process)) +
  geom_rect(aes(xmin = Left, xmax = Right, ymin = Y - 0.34, ymax = Y + 0.34)) +
  geom_vline(xintercept = 0, linewidth = 0.35) +
  geom_hline(yintercept = c(9, 18), colour = "#777777", linewidth = 0.25) +
  annotate("text", x = -125, y = c(22.5, 13.5, 4.5), label = c("Genus", "Species", "OTU"), angle = 90, family = "Arial", size = 18 / .pt) +
  annotate("text", x = -50, y = 4.5, label = "Not applicable", family = "Arial", size = 18 / .pt, colour = "#666666") +
  annotate("text", x = c(-50, 50), y = 27.3, label = c("MIC", "SEQ"), family = "Arial", size = 18 / .pt) +
  scale_fill_manual(values = process_colours, limits = legend_breaks, breaks = legend_breaks, drop = FALSE) +
  scale_x_continuous(breaks = seq(-100, 100, 50), labels = abs(seq(-100, 100, 50)), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(26:19, 17:10, 8:1), labels = rep(treatments, 3), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-100, 100), ylim = c(0.5, 26.5), clip = "off") +
  labs(x = "Proportion (%)", y = NULL, title = "Eukaryotic phytoplankton", tag = "(a)") +
  figure_theme + theme(plot.margin = margin(18, 8, 7, 22, "mm"))

d <- plot_data[plot_data$Community == "Zooplankton", ]
p_b <- ggplot(d, aes(fill = Process)) +
  geom_rect(aes(xmin = Left, xmax = Right, ymin = Y - 0.34, ymax = Y + 0.34)) +
  geom_vline(xintercept = 0, linewidth = 0.35) +
  geom_hline(yintercept = c(9, 18), colour = "#777777", linewidth = 0.25) +
  annotate("text", x = -125, y = c(22.5, 13.5, 4.5), label = c("Genus", "Species", "OTU"), angle = 90, family = "Arial", size = 18 / .pt) +
  annotate("text", x = -50, y = 4.5, label = "Not applicable", family = "Arial", size = 18 / .pt, colour = "#666666") +
  annotate("text", x = c(-50, 50), y = 27.3, label = c("MIC", "SEQ"), family = "Arial", size = 18 / .pt) +
  scale_fill_manual(values = process_colours, limits = legend_breaks, breaks = legend_breaks, drop = FALSE) +
  scale_x_continuous(breaks = seq(-100, 100, 50), labels = abs(seq(-100, 100, 50)), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(26:19, 17:10, 8:1), labels = rep(treatments, 3), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-100, 100), ylim = c(0.5, 26.5), clip = "off") +
  labs(x = "Proportion (%)", y = NULL, title = "Zooplankton", tag = "(b)") +
  figure_theme + theme(plot.margin = margin(18, 8, 7, 22, "mm"))

figure_name <- "Fig7ab"
panels <- list(p_a, p_b)
positions <- rbind(c(1, 1), c(1, 2))
panel_width_mm <- c(225, 225)
panel_height_mm <- c(270, 270) # Data-region height, excluding titles, axes and legends.
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

# Read calculated results from the working directory containing this script; do not refit models or rerun tests.
library(ggplot2) # Plotting.
library(grid)
library(gtable)
figure_theme <- theme_classic(base_family = "Arial") +
  theme(axis.text = element_text(size = 18, colour = "black"), axis.title = element_text(size = 20),
        plot.title = element_text(size = 18, hjust = 0.5), plot.tag = element_text(size = 20, face = "bold"),
        plot.tag.position = "topleft", legend.title = element_blank(), legend.text = element_text(size = 20),
        plot.margin = margin(7, 7, 7, 7, "mm"))

d <- read.csv("Alpha_paired_wilcoxon.csv", check.names = FALSE)
stopifnot(nrow(d) == 12, !anyDuplicated(d[, c("Organism", "Rank", "Metric")]))
d$X <- ifelse(d$Organism == "Plant", 0, 2) + match(d$Rank, c("Genus", "Species"))
d$Y <- match(d$Metric, c("Simpson", "Shannon", "Richness"))
d$Fill <- ifelse(d$Median_SEQ_minus_MIC < 0, "black", "white")
d$Ink <- ifelse(d$P < 0.05 & d$Median_SEQ_minus_MIC < 0, "white", "black")
p <- ggplot() +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = 0.5, ymax = 3.95, fill = "#EEEEEE") +
  geom_tile(data = d[d$P < 0.05, ], aes(X, Y, fill = Fill), width = 0.65, height = 0.53, colour = "black", linewidth = 0.5) +
  geom_text(data = d, aes(X, Y, label = Significance, colour = Ink), family = "Arial", fontface = "bold", size = 16 / 2.845) +
  annotate("text", x = c(1.5, 3.5), y = 3.65, label = c("Eukaryotic\nphytoplankton", "Zooplankton"), family = "Arial", size = 18 / 2.845) +
  scale_fill_identity() + scale_colour_identity() +
  scale_x_continuous(breaks = 1:4, labels = rep(c("Genus", "Species"), 2), expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:3, labels = c("Simpson", "Shannon", "Richness"), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.5, 4.5), ylim = c(0.5, 3.95)) +
  labs(x = NULL, y = NULL, tag = "(b)", caption = "Paired Wilcoxon; unadjusted P\nBlack: SEQ < MIC; white: SEQ > MIC") +
  figure_theme + theme(plot.caption = element_text(size = 16, hjust = 0))
figure_name <- "Fig4b"
panels <- list(p)
positions <- matrix(c(1, 1), nrow = 1)
panel_width_mm <- 160

ragg::agg_capture(width = 400, height = 500, units = "mm", res = 100)
grobs <- vector("list", length(panels))
width_mm <- height_mm <- numeric(length(panels))
for (j in seq_along(panels)) {
  g <- ggplotGrob(panels[[j]])
  cells <- g$layout[grepl("^panel", g$layout$name), ]
  g$widths[unique(cells$l)] <- unit(panel_width_mm, "mm")
  g$heights[unique(cells$t)] <- unit(100, "mm") # Data-region height: 100 mm.
  grobs[[j]] <- g
  width_mm[j] <- convertWidth(sum(g$widths), "mm", valueOnly = TRUE)
  height_mm[j] <- convertHeight(sum(g$heights), "mm", valueOnly = TRUE)
}

column_widths <- as.numeric(tapply(width_mm, positions[, 2], max))
row_heights <- as.numeric(tapply(height_mm, positions[, 1], max))
column_count <- length(column_widths)
row_count <- length(row_heights)
layout_widths <- c(6, head(as.vector(rbind(column_widths, rep(10, column_count))), -1), 6)
layout_heights <- c(6, head(as.vector(rbind(row_heights, rep(10, row_count))), -1), 6)

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

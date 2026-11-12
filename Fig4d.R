# Read calculated results from the working directory containing this script; do not refit models or rerun tests.
library(ggplot2) # Plotting.
library(grid)
library(gtable)
figure_theme <- theme_classic(base_family = "Arial") +
  theme(axis.text = element_text(size = 18, colour = "black"), axis.title = element_text(size = 20),
        plot.title = element_text(size = 18, hjust = 0.5), plot.tag = element_text(size = 20, face = "bold"),
        plot.tag.position = "topleft", legend.title = element_blank(), legend.text = element_text(size = 20),
        plot.margin = margin(7, 7, 7, 7, "mm"))

d <- read.csv("Agreement_all_model_coefficients.csv", check.names = FALSE)
d <- d[d$Model == "pooled" & d$Term != "(Intercept)", ]
term_order <- c("OrganismPhytoplankton", "LevelSpecies", "z_Richness_diff_abs_pair", "z_Shannon_diff_abs_pair")
stopifnot(nrow(d) == 4, setequal(d$Term, term_order))
d <- d[match(term_order, d$Term), ]
stopifnot(all(is.finite(as.matrix(d[, c("Beta", "Beta_CI_low", "Beta_CI_high", "P")]))), all(d$Beta_CI_low <= d$Beta), all(d$Beta_CI_high >= d$Beta))
d$Y <- 4:1
d$Ink <- ifelse(d$Term == "z_Richness_diff_abs_pair", "#287F8E", "#777777")
d$P_label <- ifelse(d$P < 0.001, "<0.001", sprintf("%.3f", d$P))
labels <- c("Organism\n(eukaryotic phytoplankton vs zooplankton)", "Taxonomic rank\n(species vs genus)", "Richness difference", "Shannon difference")
x_range <- range(c(0, d$Beta_CI_low, d$Beta_CI_high))
x_span <- max(diff(x_range), 1)
x_limits <- x_range + c(-0.08, 0.08) * x_span
p_x <- x_limits[2] + 0.12 * x_span
p <- ggplot(d) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#333333", linewidth = 0.4) +
  geom_errorbar(aes(xmin = Beta_CI_low, xmax = Beta_CI_high, y = Y, colour = Ink), orientation = "y", width = 0.13, linewidth = 0.7) +
  geom_point(aes(Beta, Y), colour = "black", size = 3) +
  geom_text(aes(x = p_x, y = Y, label = P_label), family = "Arial", size = 16 / 2.845) +
  annotate("text", x = p_x, y = 4.5, label = "italic(P)", parse = TRUE, family = "Arial", size = 18 / 2.845) +
  scale_colour_identity() +
  scale_x_continuous(breaks = pretty(x_limits, n = 5), expand = c(0, 0)) +
  scale_y_continuous(breaks = 4:1, labels = labels, expand = c(0, 0)) +
  coord_cartesian(xlim = x_limits, ylim = c(0.5, 4.8), clip = "off") +
  labs(x = "Regression coefficient (β)", y = NULL, tag = "(d)") + figure_theme +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(7, 28, 7, 7, "mm"))
figure_name <- "Fig4d"
panels <- list(p)
positions <- matrix(c(1, 1), nrow = 1)
panel_width_mm <- 130

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

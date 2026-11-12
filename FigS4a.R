# Read calculated results from the working directory containing this script; do not refit models or rerun tests.
library(ggplot2) # Plotting.
library(grid)
library(gtable)
figure_theme <- theme_classic(base_family = "Arial") +
  theme(axis.text = element_text(size = 18, colour = "black"), axis.title = element_text(size = 20),
        plot.title = element_text(size = 18, hjust = 0.5), plot.tag = element_text(size = 20, face = "bold"),
        plot.tag.position = "topleft", legend.title = element_blank(), legend.text = element_text(size = 20),
        plot.margin = margin(7, 7, 7, 7, "mm"))

pairs <- read.csv("Alpha_paired_samples.csv", check.names = FALSE)
tests <- read.csv("Alpha_paired_wilcoxon.csv", check.names = FALSE)
pairs <- pairs[pairs$Included, ]
stopifnot(!anyDuplicated(pairs[, c("Sample", "Organism", "Rank")]))
colours <- c(MIC = "#287F8E", SEQ = "#D58A45")
panels <- vector("list", 6)
panel_id <- 0
for (metric in c("Richness", "Shannon", "Simpson")) {
  top <- max(c(pairs[[paste0(metric, "_MIC")]], pairs[[paste0(metric, "_SEQ")]]))
  for (rank in c("Genus", "Species")) {
    panel_id <- panel_id + 1
    source <- pairs[pairs$Rank == rank, ]
    mic <- data.frame(Sample = source$Sample, Organism = source$Organism, Method = "MIC", Value = source[[paste0(metric, "_MIC")]])
    seq <- data.frame(Sample = source$Sample, Organism = source$Organism, Method = "SEQ", Value = source[[paste0(metric, "_SEQ")]])
    d <- rbind(mic, seq)
    d$Method <- factor(d$Method, levels = c("MIC", "SEQ"))
    d$X <- ifelse(d$Organism == "Plant", 0, 2) + as.numeric(d$Method)
    tt <- tests[tests$Metric == metric & tests$Rank == rank, ]
    stopifnot(nrow(tt) == 2, all(tt$N == as.integer(table(factor(source$Organism, levels = tt$Organism)))))
    tt$X <- ifelse(tt$Organism == "Plant", 1.5, 3.5)
    tt$Bracket <- top * 1.10
    tt$Star <- top * 1.20
    tt$Text <- top * 1.58
    tt$Label <- paste0(ifelse(tt$Organism == "Plant", "Eukaryotic\nphytoplankton", "Zooplankton"), "\nn = ", tt$N)
    p <- ggplot(d, aes(X, Value)) +
      annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf, fill = "#EEEEEE") +
      geom_boxplot(aes(fill = Method, group = X), width = 0.52, outlier.shape = NA, linewidth = 0.4, alpha = 0.8) +
      geom_point(colour = "black", size = 0.8, position = position_jitter(width = 0.13, height = 0, seed = 123)) +
      geom_segment(data = tt, aes(x = X - 0.5, xend = X + 0.5, y = Bracket, yend = Bracket), inherit.aes = FALSE, linewidth = 0.4) +
      geom_text(data = tt, aes(X, y = Star, label = Significance), inherit.aes = FALSE, family = "Arial", fontface = "bold", size = 16 / 2.845) +
      geom_text(data = tt, aes(X, y = Text, label = Label), inherit.aes = FALSE, family = "Arial", size = 18 / 2.845, lineheight = 1.05) +
      scale_fill_manual(values = colours, drop = FALSE) +
      scale_x_continuous(breaks = 1:4, labels = rep(c("MIC", "SEQ"), 2), limits = c(0.5, 4.5), expand = c(0, 0)) +
      scale_y_continuous(breaks = if (metric == "Simpson") seq(0, 1, 0.25) else waiver(), expand = c(0, 0)) +
      coord_cartesian(ylim = c(0, top * 1.85)) +
      labs(x = NULL, y = if (rank == "Genus") metric else NULL, title = if (metric == "Richness") rank else NULL,
           tag = if (panel_id == 1) "(a)" else NULL) + figure_theme + theme(legend.position = "none")
    panels[[panel_id]] <- p
  }
}
figure_name <- "FigS4a"
positions <- rbind(c(1, 1), c(1, 2), c(2, 1), c(2, 2), c(3, 1), c(3, 2))
panel_width_mm <- 140

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
common_widths <- grobs[[1]]$widths
common_heights <- grobs[[1]]$heights
for (j in seq_along(grobs)) {
  common_widths <- unit.pmax(common_widths, grobs[[j]]$widths)
  common_heights <- unit.pmax(common_heights, grobs[[j]]$heights)
}
for (j in seq_along(grobs)) {
  grobs[[j]]$widths <- common_widths
  grobs[[j]]$heights <- common_heights
  width_mm[j] <- convertWidth(sum(common_widths), "mm", valueOnly = TRUE)
  height_mm[j] <- convertHeight(sum(common_heights), "mm", valueOnly = TRUE)
}
legend_table <- ggplotGrob(panels[[1]] + theme(legend.position = "bottom"))
legend_grob <- legend_table$grobs[[which(legend_table$layout$name == "guide-box-bottom")]]
legend_height <- convertHeight(sum(legend_grob$heights), "mm", valueOnly = TRUE)
column_widths <- as.numeric(tapply(width_mm, positions[, 2], max))
row_heights <- as.numeric(tapply(height_mm, positions[, 1], max))
column_count <- length(column_widths)
row_count <- length(row_heights)
layout_widths <- c(6, head(as.vector(rbind(column_widths, rep(10, column_count))), -1), 6)
layout_heights <- c(6, head(as.vector(rbind(row_heights, rep(10, row_count))), -1), 6)
layout_heights <- c(head(layout_heights, -1), 8, legend_height, 6)
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

# Hex sticker for treeshapesim: four trees with increasing Aldous' beta (from a
# comb to a balanced tree), drawn with hexSticker. Run from the package root;
# writes man/figures/logo.png. Needs ggplot2, hexSticker and pkgload. The
# logo.png in the repository was made with data-raw/hex_preview_ggplot2.R,
# which draws the same sticker with ggplot2 only.
library(ggplot2)
library(hexSticker)
pkgload::load_all()

# segments of a rectangular tree drawn upwards, scaled to box (x0, x1, y0, y1)
tree_segments <- function(tr, box) {
  tr <- ape::reorder.phylo(tr, "cladewise")
  n <- ape::Ntip(tr); e <- tr$edge
  depth <- ape::node.depth.edgelength(tr)
  x <- numeric(n + tr$Nnode)
  tips <- e[e[, 2] <= n, 2]; x[tips] <- seq_along(tips)
  for (v in unique(ape::reorder.phylo(tr, "postorder")$edge[, 1]))
    x[v] <- mean(range(x[e[e[, 1] == v, 2]]))
  X <- box[1] + (x - 1) / (n - 1) * (box[2] - box[1])
  Y <- box[3] + depth / max(depth) * (box[4] - box[3])
  p <- e[, 1]; ch <- e[, 2]
  stems <- data.frame(x = X[ch], xend = X[ch], y = Y[p], yend = Y[ch])
  bars <- do.call(rbind, lapply(unique(p), function(v) {
    k <- ch[p == v]
    data.frame(x = min(X[k]), xend = max(X[k]), y = Y[v], yend = Y[v])
  }))
  rbind(stems, bars)
}

betas <- c(-1.95, -1, 0, 10)
pal <- grDevices::colorRampPalette(c("#E76F51", "#2A9D8F"))(4)
segs <- do.call(rbind, lapply(1:4, function(i) {
  set.seed(70 + i)
  x0 <- (i - 1) * 1.15
  cbind(tree_segments(sim_tree_shape(10, betas[i], 0), c(x0, x0 + 0.9, 0, 1.6)),
        tree = factor(i))
}))

p <- ggplot(segs) +
  geom_segment(aes(x, y, xend = xend, yend = yend, colour = tree),
               linewidth = 0.45, lineend = "square") +
  annotate("segment", x = 0.2, xend = 4.15, y = -0.35, yend = -0.35,
           colour = "#264653", linewidth = 0.4,
           arrow = grid::arrow(length = grid::unit(0.06, "inches"), type = "closed")) +
  annotate("text", x = 2.17, y = -0.65, label = "\u03b2",
           colour = "#264653", size = 3.2) +
  scale_colour_manual(values = pal, guide = "none") +
  coord_cartesian(xlim = c(-0.1, 4.4), ylim = c(-0.8, 1.7), expand = FALSE) +
  theme_void() + theme_transparent()

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
sticker(p, package = "treeshapesim",
        p_size = 17, p_color = "#264653", p_y = 0.47,
        s_x = 1, s_y = 1.08, s_width = 1.35, s_height = 0.85,
        h_fill = "#FFFFFF", h_color = "#264653", h_size = 1.6,
        url = "diogoprov.github.io/treeshapesim", u_size = 4.5, u_color = "#264653",
        filename = "man/figures/logo.png", dpi = 300)

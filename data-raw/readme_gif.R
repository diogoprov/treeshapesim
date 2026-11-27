# Animated figure for the README: trees simulated along a circuit in the
# (beta, c) parameter space. Run from the package root:
#   Rscript data-raw/readme_gif.R
# Requires pkgload and gifski.
pkgload::load_all(".", quiet = TRUE)
out <- "man/figures/treeshapesim.gif"
library(ape)

ink <- "#264653"; pal <- grDevices::colorRampPalette(c("#E76F51", "#E9C46A", "#2A9D8F"))
n <- 48
# circuit: imbalance (beta) and branching times (c) change one at a time
seg <- function(a, b, k) seq(a, b, length.out = k + 1)[-(k + 1)]
k <- 10
beta <- c(seg(-1.6, 4, k), rep(4, k), seg(4, -1.6, k), rep(-1.6, k))
cc   <- c(rep(-1.2, k), seg(-1.2, 1.2, k), rep(1.2, k), seg(1.2, -1.2, k))

frame <- function(i) {
  set.seed(2026)
  tr <- sim_tree_shape(n, beta = beta[i], c = cc[i])
  tr <- ladderize(tr)
  sh <- tree_shape(tr)
  depth <- node.depth.edgelength(tr)
  ecol <- pal(100)[pmax(1, ceiling(depth[tr$edge[, 2]] / max(depth) * 100))]
  layout(matrix(1:2, 1), widths = c(2.6, 1))
  par(mar = c(1.5, 1, 4.2, 0.5), bg = "white", family = "sans")
  plot(tr, show.tip.label = FALSE, edge.color = ecol, edge.width = 2.6,
       no.margin = FALSE, x.lim = c(0, 1.02))
  mtext("treeshapesim", side = 3, line = 2, adj = 0, cex = 2.1, font = 2, col = ink)
  mtext("Simulate trees with separately controlled imbalance and branching times",
        side = 3, line = 0.6, adj = 0, cex = 1, col = ink)
  # right panel: position in parameter space and shape statistics
  par(mar = c(5.5, 4.2, 6, 1.5))
  plot(beta, cc, type = "n", xlim = c(-2, 4.6), ylim = c(-1.6, 1.6), las = 1,
       xlab = "", ylab = expression("branching times  " * italic(c)),
       col.lab = ink, col.axis = ink, fg = ink, cex.lab = 1.1)
  mtext(expression("imbalance  " * beta), side = 1, line = 2.1, col = ink, cex = 0.95 * 1.1)
  polygon(c(-1.6, 4, 4, -1.6), c(-1.2, -1.2, 1.2, 1.2), border = "grey80", lty = 2)
  points(beta[seq_len(i)], cc[seq_len(i)], pch = 16, cex = 0.5, col = "grey70")
  points(beta[i], cc[i], pch = 21, bg = "#E76F51", col = ink, cex = 2.4)
  mtext(sprintf("β = %5.2f    c = %5.2f", beta[i], cc[i]), side = 3, line = 3.2, col = ink, cex = 1.05)
  mtext(sprintf("γ = %5.2f    Colless = %4.2f", sh[["gamma"]], sh[["colless"]]),
        side = 3, line = 1.8, col = ink, cex = 1.05)
  mtext("\u03b2 mainly sets the Colless index", side = 1, line = 3.4, col = "grey40", cex = 0.85)
  mtext("c mainly sets the \u03b3 statistic", side = 1, line = 4.3, col = "grey40", cex = 0.85)
}

gifski::save_gif(for (i in seq_along(beta)) frame(i), gif_file = out,
                 width = 1200, height = 675, res = 110, delay = 0.35, loop = TRUE, progress = FALSE)

# Preview of data-raw/hex_sticker.R without hexSticker (ggplot2 only): same subplot,
`%||%` <- function(a, b) if (is.na(a) || is.null(a)) b else a
# hexagon drawn by hand. Usage: Rscript data-raw/hex_preview_ggplot2.R out.png
library(ggplot2)
src <- readLines("data-raw/hex_sticker.R")
src <- src[!grepl("^library\\(hexSticker\\)|^pkgload::", src)]
src <- sub(" \\+ theme_transparent\\(\\)", "", src)
stop_at <- grep("^dir.create", src) - 1
source("R/simulate.R")
eval(parse(text = src[1:stop_at]))
a <- pi / 2 + seq(0, 2 * pi, length.out = 7)[-7]
hex <- data.frame(x = 1 + cos(a), y = 1 + sin(a))
g <- ggplot() +
  geom_polygon(aes(x, y), data = hex, fill = "#FFFFFF", colour = "#264653", linewidth = 1.3) +
  annotation_custom(ggplotGrob(p), xmin = 1 - 1.35 / 2, xmax = 1 + 1.35 / 2,
                    ymin = 1.08 - 0.85 / 2, ymax = 1.08 + 0.85 / 2) +
  annotate("text", 1, 0.47, label = "treeshapesim", colour = "#264653", size = 3.9, fontface = "bold") +
  annotate("text", 1.43, 0.33, label = "diogoprov.github.io/treeshapesim", colour = "#264653",
           size = 1.05, angle = 30) +
  coord_fixed(xlim = c(0.12, 1.88), ylim = c(-0.02, 2.02), expand = FALSE) + theme_void()
ggsave(commandArgs(TRUE)[1], g, width = 43.9, height = 50.8, units = "mm", dpi = as.numeric(commandArgs(TRUE)[2] %||% 600), bg = "transparent")

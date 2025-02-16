# Packages
library(readr)
library(sf)
library(RColorBrewer)

# Option
options(scipen = 10000)

# Working directory
# setwd("")

# Load data
comat_p <- read_delim("Data/1_Plant/CoMat_Plant.csv",
  delim = ";", col_name = TRUE,
  show_col_types = FALSE
)
comat_s <- read_delim("Data/2_Spectral/CoMat_Spectral.csv",
  delim = ";", col_name = TRUE,
  show_col_types = FALSE
)
shp <- st_read("Data/1_Plant/SHPs/Grid_Plant.shp")

# Figure 1 (Map of France)
box <- st_read("Data/0_Countries/Box.shp")
back <- st_read("Data/0_Countries/Background.shp")
cbn <- st_read("Data/0_Countries/CBN.shp")

pdf("Fig1.pdf", width = 7, height = 7, useDingbats = FALSE)

  par(mar = c(0, 0, 0, 0))
  
  plot(st_geometry(box), col = "lightgrey")
  plot(st_geometry(back), col = "darkgrey", add = TRUE)
  plot(st_geometry(cbn), col = "#5AA584", border = "black", add = TRUE)
  
  segments(st_bbox(box)[1] + 100000,
    st_bbox(box)[2] + 100000,
    st_bbox(box)[1] + 300000,
    st_bbox(box)[2] + 100000,
    lwd = 3, col = "grey2"
  )
  
  text(st_bbox(box)[1] + 200000,
    st_bbox(box)[2] + 60000,
    labels = "200 km", cex = 1.25, font = 2, col = "grey2"
  )

dev.off()

# Format data
comat_p <- data.frame(comat_p)
idp <- colnames(comat_p)[-1]
comat_p <- as.matrix(comat_p[, -1])

comat_s <- data.frame(comat_s)
ids <- colnames(comat_s)[-1]
comat_s <- as.matrix(comat_s[, -1])

shp <- shp[match(idp, shp$id), ]

# Stats
nbp <- apply(comat_p, 2, sum) # Number of observations (plant)
nbp0 <- apply(comat_p > 0, 2, sum) # Number of plant species

nbs <- apply(comat_s, 2, sum) # Number of pixels
nbs0 <- apply(comat_s > 0, 2, sum) # Number of spectral species

# Figure S1
pdf("FigS1.pdf", width = 8.84, height = 9.68, useDingbats = FALSE)
  
  layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8), 4, 2, byrow = FALSE), 
         width = c(1, 1, 1, 1), height = c(1, 0.2, 1, 0.2))
  
  # nbp
  q <- c(0, 500, 1000, 2000, 3000, 200000)
  qp <- cut(nbp, q, labels = FALSE)
  
  colo <- rev(brewer.pal(length(q) - 1, "Spectral"))
  
  par(mar = c(0, 0, 0, 0))
  plot(st_geometry(shp), col = colo[qp], border = NA)
  legend("topleft", inset = c(-0.11, -0.04), legend = "(a)", 
         bty = "n", cex = 2.8, xpd = TRUE, text.font = 1)
  
  par(mar = c(5, 6, 0, 6))
  image(seq(1, length(q) - 1, 1), 1, 
        matrix(data = seq(1, length(q) - 1, 1), 
               nrow = length(seq(1, length(q) - 1, 1)), ncol = 1), 
        col = colo, xlab = "", ylab = "", axes = FALSE)
  axis(1, at = (0.5 + seq(0, length(q) - 1, length.out = length(q))), 
       labels = c(q[-length(q)], "+"), las = 1, cex.axis = 1.5)
  mtext("Number of observations", 1, line = 3.8, cex = 1.7)
  box(lwd = 1.5)
  
  # nbs
  q <- c(0, 100, 110, 120, 130)
  qs <- cut(nbs, q, labels = FALSE)
  
  colo <- rev(brewer.pal(length(q) - 1, "Spectral"))
  
  par(mar = c(0, 0, 0, 0))
  plot(st_geometry(shp), col = colo[qs], border = NA)
  legend("topleft", inset = c(-0.11, -0.04), legend = "(c)", 
         bty = "n", cex = 2.8, xpd = TRUE, text.font = 1)
  
  par(mar = c(5, 6, 0, 6))
  image(seq(1, length(q) - 1, 1), 1, 
        matrix(data = seq(1, length(q) - 1, 1), 
               nrow = length(seq(1, length(q) - 1, 1)), ncol = 1), 
        col = colo, xlab = "", ylab = "", axes = FALSE)
  axis(1, at = (0.5 + seq(0, length(q) - 1, length.out = length(q))), 
       labels = q, las = 1, cex.axis = 1.5)
  mtext("Number of pixels", 1, line = 3.8, cex = 1.7)
  box(lwd = 1.5)
  
  # nbp0
  q <- seq(0, 1000, 200)
  qp0 <- cut(nbp0, q, labels = FALSE)
  
  colo <- rev(brewer.pal(length(q) - 1, "Spectral"))
  
  par(mar = c(0, 0, 0, 0))
  plot(st_geometry(shp), col = colo[qp0], border = NA)
  legend("topleft", inset = c(-0.11, -0.04), legend = "(b)", 
         bty = "n", cex = 2.8, xpd = TRUE, text.font = 1)
  
  par(mar = c(5, 6, 0, 6))
  image(seq(1, length(q) - 1, 1), 1, 
        matrix(data = seq(1, length(q) - 1, 1), 
               nrow = length(seq(1, length(q) - 1, 1)), ncol = 1), 
        col = colo, xlab = "", ylab = "", axes = FALSE)
  axis(1, at = (0.5 + seq(0, length(q) - 1, length.out = length(q))), 
       labels = q, las = 1, cex.axis = 1.5)
  mtext("Number of plant species", 1, line = 3.8, cex = 1.7)
  box(lwd = 1.5)
  
  # nbs0
  q <- c(0, 20, 40, 60, 80)
  qs0 <- cut(nbs0, q, labels = FALSE)
  
  colo <- rev(brewer.pal(length(q) - 1, "Spectral"))
  
  par(mar = c(0, 0, 0, 0))
  plot(st_geometry(shp), col = colo[qs0], border = NA)
  legend("topleft", inset = c(-0.11, -0.04), legend = "(d)", 
         bty = "n", cex = 2.8, xpd = TRUE, text.font = 1)
  
  par(mar = c(5, 6, 0, 6))
  image(seq(1, length(q) - 1, 1), 1, 
        matrix(data = seq(1, length(q) - 1, 1), 
               nrow = length(seq(1, length(q) - 1, 1)), ncol = 1), 
        col = colo, xlab = "", ylab = "", axes = FALSE)
  axis(1, at = (0.5 + seq(0, length(q) - 1, length.out = length(q))), 
       labels = q, las = 1, cex.axis = 1.5)
  mtext("Number of spectral species", 1, line = 3.8, cex = 1.7)
  box(lwd = 1.5)

dev.off()

# Filters 0, 1 & 2
f0 <- ((nbp > 0) & (nbs > 0))

f1p <- !(nbp < 200)
f1s <- !(nbs < 20)
f1 <- (f1p & f1s)

f2 <- (nbp0 > (1.5 * nbp^(1 / 1.5)))

sum(f0) / length(f0)
sum(f1) / length(f1)
sum(f1 & f2) / length(f2)

colo=c("#4285F4","#E94235","#A813DE")

# Figure S2
pdf("FigS2.pdf", width = 14.06250, height = 6.27632, useDingbats = FALSE)
  
  par(mfrow = c(1, 2))
  
  # a
  par(mar = c(5, 7, 1, 1))
  plot(nbp[f1p], nbp0[f1p], pch = 16, col = colo[1], 
       xlim = c(0, 1000), ylim = c(0, 500), axes = FALSE, xlab = "", ylab = "")
  par(new = TRUE)
  plot(nbp[!f1p], nbp0[!f1p], pch = 16, col = colo[2], 
       xlim = c(0, 1000), ylim = c(0, 500), axes = FALSE, xlab = "", ylab = "")
  abline(0, 1, col = "grey", lty = 2, lwd = 3)
  
  box(lwd = 1.5)
  axis(1, cex.axis = 1.5)
  axis(2, las = 1, cex.axis = 1.5)
  mtext("Number of observations", 1, line = 3.5, cex = 2.3)
  mtext("Number of plant species", 2, line = 4, cex = 2.3)
  
  legend("topleft", inset = c(-0.13, -0.05), legend = "(a)", 
         bty = "n", cex = 2.8, xpd = TRUE, text.font = 1)
  
  # b
  plot(nbs[f1s], nbs0[f1s], pch = 16, col = colo[1], 
       xlim = c(0, 130), ylim = c(0, 80), axes = FALSE, xlab = "", ylab = "")
  par(new = TRUE)
  plot(nbs[!f1s], nbs0[!f1s], pch = 16, col = colo[2], 
       xlim = c(0, 130), ylim = c(0, 80), axes = FALSE, xlab = "", ylab = "")
  abline(0, 1, col = "grey", lty = 2, lwd = 3)
  
  box(lwd = 1.5)
  axis(1, cex.axis = 1.5)
  axis(2, las = 1, cex.axis = 1.5)
  mtext("Number of pixels", 1, line = 3.5, cex = 2.3)
  mtext("Number of spectral species", 2, line = 4, cex = 2.3)
  
  legend("topleft", inset = c(-0.13, -0.05), legend = "(b)", 
         bty = "n", cex = 2.8, xpd = TRUE, text.font = 1)
  
dev.off()

# Figure S3
pdf("FigS3.pdf", width = 9.322917, height = 6.182488, useDingbats = FALSE)

  par(mar = c(5, 7, 1, 1))
  plot(nbp[f2], nbp0[f2], pch = 16, col = colo[1], 
       xlim = c(0, 5000), ylim = c(0, 800), axes = FALSE, xlab = "", ylab = "")
  par(new = TRUE)
  plot(nbp[!f2], nbp0[!f2], pch = 16, col = colo[3], 
       xlim = c(0, 5000), ylim = c(0, 800), axes = FALSE, xlab = "", ylab = "")
  
  box(lwd = 1.5)
  axis(1, cex.axis = 1.5)
  axis(2, las = 1, cex.axis = 1.5)
  mtext("Number of observations", 1, line = 3.5, cex = 2.3)
  mtext("Number of plant species", 2, line = 4, cex = 2.3)

dev.off()

# Figure S4
pdf("FigS4.pdf", width = 7, height = 7, useDingbats = FALSE)
  
  f <- rep(1, length(f1p))
  f[!f1] <- 2
  f[!f2] <- 3
  
  shp$f <- f
  
  par(mar = c(0, 0, 0, 0))
  plot(st_geometry(shp), col = colo[shp$f], border = NA)

dev.off()

# Export filters
if (!dir.exists("Analysis/Filters")) {
  dir.create("Analysis/Filters")
}
write_delim(data.frame(ID = shp$id, 
                       Filter0 = f0 * 1, 
                       Filter1 = f1 * 1, 
                       Filter2 = (f1 & f2) * 1),
  "Analysis/Filters/Filters.csv",
  delim = ";", col_name = TRUE
)

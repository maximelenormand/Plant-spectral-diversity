# Packages
library(readr)

# Option
options(scipen = 10000)

# Working directory
# setwd("")

# Load data
alpha_p <- read_delim("Analysis/Metrics/Metrics_Alpha_Plant.csv",
  delim = ";",
  col_name = TRUE, show_col_types = FALSE
)
alpha_s <- read_delim("Analysis/Metrics/Metrics_Alpha_Spectral.csv",
  delim = ";",
  col_name = TRUE, show_col_types = FALSE
)
beta_p <- read_delim("Analysis/Metrics/Metrics_Beta_Plant.csv",
  delim = ";",
  col_name = TRUE, show_col_types = FALSE
)
beta_s <- read_delim("Analysis/Metrics/Metrics_Beta_Spectral.csv",
  delim = ";",
  col_name = TRUE, show_col_types = FALSE
)
filters <- read_delim("Analysis/Filters/Filters.csv",
  delim = ";",
  col_name = TRUE, show_col_types = FALSE
)

# Match between sites in beta & filters 
idfiltersbeta1 <- match(beta_s$Site1, filters$ID)
idfiltersbeta2 <- match(beta_s$Site2, filters$ID)

# Extract and clean metrics
Rs <- alpha_s$Rs
Rs <- Rs[filters$Filter0 == 1] # Remove zero-zero cells (filter 0)
Rp <- alpha_p$Rp
Rp <- Rp[filters$Filter0 == 1]
Hs <- alpha_s$Hs
Hs <- Hs[filters$Filter0 == 1]
Hp <- alpha_p$Hp
Hp <- Hp[filters$Filter0 == 1]

SIMp <- beta_p$SIM
SIMp <- SIMp[filters$Filter0[idfiltersbeta1] == 1 &
             filters$Filter0[idfiltersbeta2] == 1] 
SIMs <- beta_s$SIM
SIMs <- SIMs[filters$Filter0[idfiltersbeta1] == 1 &
             filters$Filter0[idfiltersbeta2] == 1]
BCbalp <- beta_p$BCbal
BCbalp <- BCbalp[filters$Filter0[idfiltersbeta1] == 1 &
                 filters$Filter0[idfiltersbeta2] == 1]
BCbals <- beta_s$BCbal
BCbals <- BCbals[filters$Filter0[idfiltersbeta1] == 1 &
                 filters$Filter0[idfiltersbeta2] == 1]

gc()
rm(alpha_p, alpha_s, beta_p, beta_s, filters)
gc()

# Figure 2
pdf("Fig2.pdf", width = 13.32291, height = 8.642972, useDingbats = FALSE)
  
  par(mfrow = c(2, 2))
  
  colobox <- "steelblue3"
  colomu <- "#CC6666"
  
  # a
  minx <- 0
  maxx <- 70
  miny <- 0
  maxy <- 600
  bin <- cut(Rs, c(seq(0, 60, 10), 80), labels = FALSE, right = FALSE)
  posbin <- aggregate(Rs, list(bin), mean)
  
  mu <- aggregate(Rp, list(bin), mean)
  
  b <- boxplot(Rp ~ bin, plot = F)
  for (k in 1:length(table(bin))) {
    b$stats[1, k] <- quantile(Rp[bin == k], 0.1)
    b$stats[5, k] <- quantile(Rp[bin == k], 0.9)
  }
  
  par(mar = c(5.5, 7.5, 1, 1))
  plot(Rs, Rp,
    pch = 16, cex = 0.3, col = "grey",
    xlim = c(minx, maxx), ylim = c(miny, maxy), 
    axes = FALSE, xlab = "", ylab = ""
  )
  par(new = TRUE)
  bxp(b,
    at = posbin[, 2], outline = FALSE, boxcol = colobox, whiskcol = colobox,
    whisklty = "solid", whisklwd = 2, staplelwd = 2, boxwex = 4, 
    staplecol = colobox,
    medbg = colobox, boxfill = colobox, cex.axis = 1.5, las = 1,
    xlim = c(minx, maxx), ylim = c(miny, maxy), 
    axes = FALSE, xlab = "", ylab = ""
  )
  points(posbin[, 2], mu[, 2], col = colomu, cex = 3, pch = 16)
  
  box(lwd = 1.5)
  axis(1, cex.axis = 1.5)
  axis(2, las = 1, cex.axis = 1.5)
  mtext(expression(R^S), 1, line = 4, cex = 2.3)
  mtext(expression(R^P), 2, line = 4, cex = 2.3)
  
  legend("topleft",
    inset = c(-0.36, -0.12), legend = "(a)",
    bty = "n", cex = 2.6, xpd = TRUE, text.font = 1
  )
  
  # b
  minx <- -0.01
  maxx <- 0.7
  miny <- 0
  maxy <- 0.7
  
  bin <- cut(Hs, (0.7 / 0.4) * c(seq(0, 0.3, 0.05), 0.35, 0.45, 0.6),
    labels = FALSE, right = FALSE
  )
  posbin <- aggregate(Hs, list(bin), mean)
  
  mu <- aggregate(Hp, list(bin), mean)
  
  b <- boxplot(Hp ~ bin, plot = F)
  for (k in 1:length(table(bin))) {
    b$stats[1, k] <- quantile(Hp[bin == k], 0.1)
    b$stats[5, k] <- quantile(Hp[bin == k], 0.9)
  }
  
  par(mar = c(5.5, 7.5, 1, 1))
  plot(Hs, Hp,
    pch = 16, cex = 0.3, col = "grey",
    xlim = c(minx, maxx), ylim = c(miny, maxy), axes = FALSE, xlab = "", ylab = ""
  )
  par(new = TRUE)
  bxp(b,
    at = posbin[, 2], outline = FALSE, boxcol = colobox, whiskcol = colobox,
    whisklty = "solid", whisklwd = 2, staplelwd = 2, boxwex = 0.04, 
    staplecol = colobox,
    medbg = colobox, boxfill = colobox, cex.axis = 1.5, las = 1,
    xlim = c(minx, maxx), ylim = c(miny, maxy), 
    axes = FALSE, xlab = "", ylab = ""
  )
  points(posbin[, 2], mu[, 2], col = colomu, cex = 3, pch = 16)
  
  box(lwd = 1.5)
  axis(1, cex.axis = 1.5)
  axis(2, las = 1, cex.axis = 1.5)
  mtext(expression(H^S), 1, line = 4, cex = 2.3)
  mtext(expression(H^P), 2, line = 4, cex = 2.3)
  
  legend("topleft",
    inset = c(-0.36, -0.12), legend = "(b)",
    bty = "n", cex = 2.6, xpd = TRUE, text.font = 1
  )
  
  # c
  ind <- sample(length(SIMs), 20000) # Subsample the grey points
  
  minx <- 0
  maxx <- 1
  miny <- 0
  maxy <- 1
  bin <- cut(SIMs, c(0, 0.2, seq(0.4, 0.9, 0.1), 1.1), 
             labels = FALSE, right = FALSE)
  posbin <- aggregate(SIMs, list(bin), mean)
  gc()
  
  mu <- aggregate(SIMp, list(bin), mean)
  gc()
  
  b <- boxplot(SIMp ~ bin, plot = F)
  for (k in 1:length(table(bin))) {
    b$stats[1, k] <- quantile(SIMp[bin == k], 0.1)
    b$stats[5, k] <- quantile(SIMp[bin == k], 0.9)
  }
  gc()
  
  par(mar = c(5.5, 7.5, 1, 1))
  plot(SIMs[ind], SIMp[ind],
    pch = 16, cex = 0.3, col = "grey",
    xlim = c(minx, maxx), ylim = c(miny, maxy), 
    axes = FALSE, xlab = "", ylab = ""
  )
  par(new = TRUE)
  bxp(b,
    at = posbin[, 2], outline = FALSE, boxcol = colobox, whiskcol = colobox,
    whisklty = "solid", whisklwd = 2, staplelwd = 2, boxwex = 0.05, 
    staplecol = colobox,
    medbg = colobox, boxfill = colobox, cex.axis = 1.5, las = 1,
    xlim = c(minx, maxx), ylim = c(miny, maxy), 
    axes = FALSE, xlab = "", ylab = ""
  )
  points(posbin[, 2], mu[, 2], col = colomu, cex = 3, pch = 16)
  
  box(lwd = 1.5)
  axis(1, cex.axis = 1.5)
  axis(2, las = 1, cex.axis = 1.5)
  mtext(expression(beta[SIM]^S), 1, line = 4.5, cex = 2.3)
  mtext(expression(beta[SIM]^P), 2, line = 4, cex = 2.3)
  
  legend("topleft",
    inset = c(-0.36, -0.12), legend = "(c)",
    bty = "n", cex = 2.6, xpd = TRUE, text.font = 1
  )
  
  # d
  ind <- sample(length(BCbals), 20000) # Subsample the grey points
  
  minx <- 0
  maxx <- 1
  miny <- 0
  maxy <- 1
  bin <- cut(BCbals, c(0, seq(0.4, 0.9, 0.1), 1.1), 
             labels = FALSE, right = FALSE)
  posbin <- aggregate(BCbals, list(bin), mean)
  
  mu <- aggregate(BCbalp, list(bin), mean)
  
  b <- boxplot(BCbalp ~ bin, plot = F)
  for (k in 1:length(table(bin))) {
    b$stats[1, k] <- quantile(BCbalp[bin == k], 0.1)
    b$stats[5, k] <- quantile(BCbalp[bin == k], 0.9)
  }
  
  par(mar = c(5.5, 7.5, 1, 1))
  plot(BCbals[ind], BCbalp[ind],
    pch = 16, cex = 0.3, col = "grey",
    xlim = c(minx, maxx), ylim = c(miny, maxy), 
    axes = FALSE, xlab = "", ylab = ""
  )
  par(new = TRUE)
  bxp(b,
    at = posbin[, 2], outline = FALSE, boxcol = colobox, whiskcol = colobox,
    whisklty = "solid", whisklwd = 2, staplelwd = 2, boxwex = 0.05, 
    staplecol = colobox,
    medbg = colobox, boxfill = colobox, cex.axis = 1.5, las = 1,
    xlim = c(minx, maxx), ylim = c(miny, maxy), 
    axes = FALSE, xlab = "", ylab = ""
  )
  points(posbin[, 2], mu[, 2], col = colomu, cex = 3, pch = 16)
  
  box(lwd = 1.5)
  axis(1, cex.axis = 1.5)
  axis(2, las = 1, cex.axis = 1.5)
  mtext(expression(beta[BC - bal]^S), 1, line = 4.5, cex = 2.3)
  mtext(expression(beta[BC - bal]^P), 2, line = 4, cex = 2.3)
  
  legend("topleft",
    inset = c(-0.36, -0.12), legend = "(d)",
    bty = "n", cex = 2.6, xpd = TRUE, text.font = 1
  )

dev.off()

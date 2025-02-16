# Packages
library(readr)
library(gplots)
library(RColorBrewer)
library(TeachingDemos)

# Option
options(scipen = 10000)

# Working directory
#setwd("")

# Load data
alpha_p <- read_delim("Analysis/Metrics/Metrics_Alpha_Plant.csv",
  delim = ";", col_name = TRUE, show_col_types = FALSE
)
alpha_s <- read_delim("Analysis/Metrics/Metrics_Alpha_Spectral.csv",
  delim = ";", col_name = TRUE, show_col_types = FALSE
)
beta_p <- read_delim("Analysis/Metrics/Metrics_Beta_Plant.csv",
  delim = ";", col_name = TRUE, show_col_types = FALSE
)
beta_s <- read_delim("Analysis/Metrics/Metrics_Beta_Spectral.csv",
  delim = ";", col_name = TRUE, show_col_types = FALSE
)
filters <- read_delim("Analysis/Filters/Filters.csv",
  delim = ";", col_name = TRUE, show_col_types = FALSE
)
br <- read.csv2("Analysis/Bioregions/Bioregions.csv")
nbr <- max(br$BR)

print(c(length(filters$ID), length(br$ID), sum(filters$ID == br$ID)))

# Match between sites in beta & filters
idbeta1 <- match(beta_s$Site1, filters$ID)
idbeta2 <- match(beta_s$Site2, filters$ID)

# Extract and clean metrics
Hp <- list()
Hs <- list()
BCp <- list()
BCs <- list()
for (k in 1:nbr) {
  print(c(k, nbr))

  Hp[[k]] <- alpha_p$Hp[br$BR == k & filters$Filter0 == 1]
  Hs[[k]] <- alpha_s$Hs[br$BR == k & filters$Filter0 == 1]
  BCp[[k]] <- beta_p$BCbal[br$BR[idbeta1] == k &
    br$BR[idbeta2] == k &
    filters$Filter0[idbeta1] == 1 &
    filters$Filter0[idbeta2] == 1]
  BCs[[k]] <- beta_s$BCbal[br$BR[idbeta1] == k &
    br$BR[idbeta2] == k &
    filters$Filter0[idbeta1] == 1 &
    filters$Filter0[idbeta2] == 1]
}

int <- merge(1:nbr, 1:nbr)
int <- int[int[, 1] < int[, 2], ]
int <- int[order(int[, 1], int[, 2]), ]
nint <- dim(int)[1]
BCintp <- list()
BCints <- list()
for (k in 1:nint) {
  print(c(k, nint))

  BCintp[[k]] <- beta_p$BCbal[br$BR[idbeta1] == int[k, 1] &
    br$BR[idbeta2] == int[k, 2] &
    filters$Filter0[idbeta1] == 1 &
    filters$Filter0[idbeta2] == 1]
  BCints[[k]] <- beta_s$BCbal[br$BR[idbeta1] == int[k, 1] &
    br$BR[idbeta2] == int[k, 2] &
    filters$Filter0[idbeta1] == 1 &
    filters$Filter0[idbeta2] == 1]
}

Hptot <- alpha_p$Hp[filters$Filter0 == 1]
Hstot <- alpha_s$Hs[filters$Filter0 == 1]
BCptot <- beta_p$BCbal[filters$Filter0[idbeta1] == 1 & 
                       filters$Filter0[idbeta2] == 1]
BCstot <- beta_s$BCbal[filters$Filter0[idbeta1] == 1 & 
                       filters$Filter0[idbeta2] == 1]

gc()
rm(alpha_p, alpha_s, beta_p, beta_s)
gc()

# Figure 4
pdf("Fig4.pdf", width = 16.177083, height = 4.993949, useDingbats = FALSE)

  colobr <- brewer.pal(5, "Set2")
  
  par(mfrow = c(1, 3))
  
  # a
  mup <- NULL
  sdp <- NULL
  mus <- NULL
  sds <- NULL
  for (k in 1:nbr) {
    mup <- c(mup, mean(Hp[[k]]))
    sdp <- c(sdp, sd(Hp[[k]]))
    mus <- c(mus, mean(Hs[[k]]))
    sds <- c(sds, sd(Hs[[k]]))
  }
  
  par(mar = c(5.5, 7.5, 1, 1))
  plotCI(mus, mup, ui = mup + sdp, li = mup - sdp, err = "y", 
         axes = FALSE, xlab = "", ylab = "", type = "n", 
         pch = 16, cex = 2, col = colobr, lwd = 2, 
         xlim = c(0.3, 0.7), ylim = c(0.3, 0.7))
  par(new = TRUE)
  plotCI(mus, mup, ui = mus + sds, li = mus - sds, err = "x", 
         axes = FALSE, xlab = "", ylab = "", 
         pch = 16, cex = 2, col = colobr, lwd = 2, 
         xlim = c(0.3, 0.7), ylim = c(0.3, 0.7))
  
  box(lwd = 1.5)
  axis(1, cex.axis = 1.5)
  axis(2, las = 1, cex.axis = 1.5)
  mtext(expression(H^S), 1, line = 4.5, cex = 2.1)
  mtext(expression(H^P), 2, line = 4, cex = 2.1)
  
  legend("topleft", inset = c(-0.37, -0.08), legend = "(a)", bty = "n", cex = 3,
         xpd = TRUE, text.font = 1)
  
  legend("bottomleft", inset = c(0.01, 0), 
         legend = c("BR1", "BR2", "BR3", "BR4", "BR5"), col = colobr, pch = 16, 
         pt.cex = 2, lty = 1, lwd = 2, bty = "n", cex = 2, xpd = TRUE)
  
  # b
  mup <- NULL
  sdp <- NULL
  mus <- NULL
  sds <- NULL
  for (k in 1:nbr) {
    mup <- c(mup, mean(BCp[[k]]))
    sdp <- c(sdp, sd(BCp[[k]]))
    mus <- c(mus, mean(BCs[[k]]))
    sds <- c(sds, sd(BCs[[k]]))
  }
  
  par(mar = c(5.5, 7.5, 1, 1))
  plotCI(mus, mup, ui = mup + sdp, li = mup - sdp, err = "y", 
         axes = FALSE, xlab = "", ylab = "", type = "n", 
         pch = 16, cex = 2, col = colobr, lwd = 2, 
         xlim = c(0.6, 1), ylim = c(0.2, 1))
  par(new = TRUE)
  plotCI(mus, mup, ui = mus + sds, li = mus - sds, err = "x", 
         axes = FALSE, xlab = "", ylab = "", 
         pch = 16, cex = 2, col = colobr, lwd = 2, 
         xlim = c(0.6, 1), ylim = c(0.2, 1))
  
  box(lwd = 1.5)
  axis(1, cex.axis = 1.5)
  axis(2, las = 1, cex.axis = 1.5)
  mtext(expression(beta[BC - bal]^S), 1, line = 4.6, cex = 2.1)
  mtext(expression(beta[BC - bal]^P), 2, line = 3.8, cex = 2.1)
  
  legend("topleft", inset = c(-0.37, -0.08), legend = "(b)",
         bty = "n", cex = 3, xpd = TRUE, text.font = 1)
  
  # c
  mup <- NULL
  sdp <- NULL
  mus <- NULL
  sds <- NULL
  for (k in 1:nint) {
    mup <- c(mup, mean(BCintp[[k]]))
    sdp <- c(sdp, sd(BCintp[[k]]))
    mus <- c(mus, mean(BCints[[k]]))
    sds <- c(sds, sd(BCints[[k]]))
  }
  
  par(mar = c(5.5, 7.5, 1, 1))
  plot(mus, mup, axes = FALSE, xlab = "", ylab = "", type = "n", 
       xlim = c(0.88, 1))
  for (k in 1:nint) {
    mysymb <- function() {
      plot.new()
      polygon(c(0, 0, 1, 1, 0), c(0, 1, 1, 0, 0), col = colobr[int[k, 1]], 
              lty = 0)
      polygon(c(1, 1, 2, 2, 1), c(0, 1, 1, 0, 0), col = colobr[int[k, 2]], 
              lty = 0)
    }
    my.symbols(mus[k], mup[k], symb = mysymb, symb.plots = TRUE, inches = 0.2)
  }
  
  box(lwd = 1.5)
  axis(1, cex.axis = 1.5)
  axis(2, las = 1, cex.axis = 1.5)
  mtext(expression(beta[BC - bal]^S), 1, line = 4.6, cex = 2.1)
  mtext(expression(beta[BC - bal]^P), 2, line = 3.8, cex = 2.1)
  
  legend("topleft", inset = c(-0.37, -0.08), legend = "(c)", bty = "n", cex = 3, 
         xpd = TRUE, text.font = 1)
  
  legend("bottomright", inset = c(0.01, -0.01), fill = colobr, border = colobr, 
         legend = c("BR1", "BR2", "BR3", "BR4", "BR5"), 
         bty = "n", cex = 2, xpd = TRUE)
  
dev.off()

# Table S2
tex <- paste0("BR", int[, 1], " <-> BR", int[, 2], " & ", 
              round(mup, digits = 2), " (", round(sdp, digits = 2), ")", " & ", 
              round(mus, digits = 2), " (", round(sds, digits = 2), ")", "\\")
print(data.frame(tex), row.names = FALSE)

# Figure 5
pdf("Fig5.pdf", width = 16.177083, height = 4.993949, useDingbats = FALSE)
  
  layout(matrix(c(1, 2, 3), 1, 3), height = 1, width = c(1, 1, 1.2))
  # par(mfrow=c(1,3))
  
  # a
  cor <- NULL
  cormin <- NULL
  cormax <- NULL
  for (k in 1:nbr) {
    temp <- cor.test(Hp[[k]], Hs[[k]])
    cor <- c(cor, temp$estimate)
    cormin <- c(cormin, temp$conf.int[1])
    cormax <- c(cormax, temp$conf.int[2])
  }
  temp <- cor.test(Hptot, Hstot)
  cor <- c(cor, temp$estimate)
  cormin <- c(cormin, temp$conf.int[1])
  cormax <- c(cormax, temp$conf.int[2])
  
  par(mar = c(5.5, 7.5, 2, 1))
  plot(cor, 6:1, type = "n", axes = FALSE, xlab = "", ylab = "", 
       xlim = c(-0.15, 0.25), col = c(colobr, "steelblue3"), pch = 16, cex = 3)
  for (g in 1:6) {
    abline(g, 0, col = "grey", lty = 3)
  }
  par(new = TRUE)
  plotCI(cor, 6:1, ui = cormax, li = cormin, err = "x", xlim = c(-0.15, 0.25), 
         axes = FALSE, xlab = "", ylab = "", pch = 16, cex = 3,
         col = c(colobr, "steelblue3"), lwd = 3)
  
  axis(1, las = 1, cex.axis = 1.7)
  axis(2, at = 1:6, labels = rev(c("BR1", "BR2", "BR3", "BR4", "BR5", "All")), 
       las = 2, tick = FALSE, cex.axis = 2, font = 2, padj = 0.5)
  mtext("Correlation coefficient (H)", 1, line = 4, cex = 1.9)
  box(lwd = 1.5)
  
  legend("topleft", inset = c(-0.4, -0.13), legend = "(a)", bty = "n", cex = 3, 
         xpd = TRUE, text.font = 1)
  
  # b
  cor <- NULL
  cormin <- NULL
  cormax <- NULL
  for (k in 1:nbr) {
    temp <- cor.test(BCp[[k]], BCs[[k]])
    cor <- c(cor, temp$estimate)
    cormin <- c(cormin, temp$conf.int[1])
    cormax <- c(cormax, temp$conf.int[2])
  }
  temp <- cor.test(BCptot, BCstot)
  cor <- c(cor, temp$estimate)
  cormin <- c(cormin, temp$conf.int[1])
  cormax <- c(cormax, temp$conf.int[2])
  
  par(mar = c(5.5, 7.5, 2, 1))
  plot(cor, 6:1, type = "n", axes = FALSE, xlab = "", ylab = "", 
       xlim = c(0.1, 0.5), col = c(colobr, "steelblue3"), pch = 16, cex = 3)
  for (g in 1:6) {
    abline(g, 0, col = "grey", lty = 3)
  }
  par(new = TRUE)
  plotCI(cor, 6:1, ui = cormax, li = cormin, err = "x", xlim = c(0.1, 0.5), 
         axes = FALSE, xlab = "", ylab = "", pch = 16, cex = 3, 
         col = c(colobr, "steelblue3"), lwd = 3)
  
  axis(1, las = 1, cex.axis = 1.7)
  axis(2, at = 1:6, labels = rev(c("BR1", "BR2", "BR3", "BR4", "BR5", "All")), 
       las = 2, tick = FALSE, cex.axis = 2, font = 2, padj = 0.5)
  mtext(expression(paste("Correlation coefficient (", beta[BC - bal], ")")), 1, 
        line = 4.5, cex = 1.9)
  box(lwd = 1.5)
  
  legend("topleft", inset = c(-0.4, -0.13), legend = "(b)", bty = "n", cex = 3, 
         xpd = TRUE, text.font = 1)
  
  # c
  lab <- c(
    "BR1 <-> BR2", "BR1 <-> BR3", "BR1 <-> BR4", "BR1 <-> BR5",
    "BR2 <-> BR3", "BR2 <-> BR4", "BR2 <-> BR5",
    "BR3 <-> BR4", "BR3 <-> BR5",
    "BR4 <-> BR5"
  )
  cor <- NULL
  cormin <- NULL
  cormax <- NULL
  for (k in 1:nint) {
    temp <- cor.test(BCintp[[k]], BCints[[k]])
    cor <- c(cor, temp$estimate)
    cormin <- c(cormin, temp$conf.int[1])
    cormax <- c(cormax, temp$conf.int[2])
  }
  
  par(mar = c(5.5, 14, 2, 1))
  plot(cor, 10:1, type = "n", axes = FALSE, xlab = "", ylab = "", 
       xlim = c(0, 0.35), col = "grey", pch = 16, cex = 3)
  for (g in 1:10) {
    abline(g, 0, col = "grey", lty = 3)
  }
  par(new = TRUE)
  plotCI(cor, 10:1, ui = cormax, li = cormin, err = "x", 
         xlim = c(0, 0.35), axes = FALSE, xlab = "", ylab = "", 
         pch = 16, cex = 3, col = "grey", lwd = 3)
  
  axis(1, las = 1, cex.axis = 1.7)
  axis(2, at = 1:10, labels = rev(lab), las = 2, tick = FALSE, cex.axis = 2, 
       font = 2, padj = 0.5)
  mtext(expression(paste("Correlation coefficient (", beta[BC - bal], ")")), 1, 
        line = 4.5, cex = 1.9)
  box(lwd = 1.5)
  
  legend("topleft", inset = c(-0.6, -0.13), legend = "(c)", bty = "n", cex = 3, 
         xpd = TRUE, text.font = 1)

dev.off()

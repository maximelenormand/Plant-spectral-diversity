# Load packages
library(igraph)
library(sfsmisc)
library(readr)
library(sf)
library(RColorBrewer)
library(imager) # sudo apt-get install libx11-dev
library(scales)
# library(devtools)
# devtools::install_github('thomasp85/ggforce')
library(ggforce)

# Option
options(scipen = 10000)

# Working directory
#setwd("")

# Load shp
shp <- st_read("Data/1_Plant/SHPs/Grid_Plant.shp")

# Load, merge and clean CoMats
p <- read_delim("Data/1_Plant/CoMat_Plant.csv", delim = ";", 
                col_name = TRUE, show_col_types = FALSE)
s <- read_delim("Data/2_Spectral/CoMat_Spectral.csv", delim = ";", 
                col_name = TRUE, show_col_types = FALSE)
colnames(s) <- colnames(p)

mat <- rbind(p, s)
mat <- as.data.frame(mat)

idrow <- mat[, 1]
idcol <- colnames(mat)[-1]

mat <- as.matrix(mat[, -1])
mat[mat > 0] <- 1

print(c(length(idcol), dim(shp)[1], sum(sort(idcol) == sort(shp$id))))

# Remove empty cells
ind <- (apply(mat, 2, sum) > 0)
mat <- mat[, ind]
idcol <- idcol[ind]

# Load Bipartite
bip <- read.csv2("Data/3_Bipartite/Bipartite.csv", 
                 stringsAsFactors = FALSE)[, c(1, 2, 6:8)]
bip[, 3] <- as.numeric(bip[, 3])
bip[, 4] <- as.numeric(bip[, 4])
bip[, 5] <- as.numeric(bip[, 5])

# Build network
bip[, 4] <- pmin(bip[, 4], bip[, 5])
bip[, 5] <- 1 - (bip[, 4] / (bip[, 3] + bip[, 4])) # Bray-Curtis Turnover Similarity
bip <- bip[, c(1, 2, 5)]
net <- bip
colnames(net)[3] <- "weight"

idnodes <- c(net[, 1], net[, 2])
idnodes <- sort(idnodes[!duplicated(idnodes)])
com <- data.frame(ID = idnodes, Com = 0) # Empty community dataframe

# Choose threshold
w <- net[, 3]
w <- w[w > 0]

brea <- seq(0, 1, 0.01)
h <- hist(w, breaks = brea, plot = FALSE)
x <- h$mids
y <- h$density

pdf("FigS5.pdf", width = 9.6, height = 7.3, useDingbats = FALSE)
  par(mar = c(6, 6, 1, 1))
  plot(x, y, col = "steelblue3", pch = 16, cex = 2, log = "xy", 
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", axes = FALSE)
  abline(v = 0.2, lwd = 3, col = "#CC6666")
  mtext("Edge Weight", 1, line = 4, cex = 2)
  mtext("PDF", 2, line = 4, cex = 2)
  eaxis(1, at = c(10^-9, 10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 
                  10^-3, 10^-2, 10^-1, 10^0, 10, 10^2, 10^3), 
        cex.axis = 1.5, padj = 0.5)
  eaxis(2, at = c(10^-9, 10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 
                  10^-3, 10^-2, 10^-1, 10^0, 10, 10^2, 10^3),
        cex.axis = 1.5)
  box(lwd = 1.5)
dev.off()

net <- net[net[, 3] > 0.2, ] # Threshold = 0.2 based on FigS5

# Identify communities with Louvain 
net <- graph_from_data_frame(net, directed = FALSE)

set.seed(136390586) # MAY CHANGE ACCORDING TO VERSIONS/OS [SEE BACKUP BELOW]
comlou <- cluster_louvain(net)

table(comlou$membership) # Number of species per community
table(comlou$membership[substr(comlou$names, 1, 1) == "p"]) # Number of plant species per community
table(comlou$membership[substr(comlou$names, 1, 1) == "s"]) # Number of spectral species per community

comlou <- data.frame(comlou$names, as.numeric(comlou$membership))

com[, 2] <- comlou[match(com[, 1], comlou[, 1]), 2]
com[is.na(com)] <- 0

# Reclassify values
temp <- com
com[,2][temp[,2]==1]=5
com[,2][temp[,2]==2]=2
com[,2][temp[,2]==3]=3
com[,2][temp[,2]==4]=4
com[,2][temp[,2]==5]=1

# BACKUP IF PROBLEM WITH THE SEED... ###########################################
valobj = c(19,29,38,62,102)
val = as.numeric(table(comlou[,2][substr(comlou[,1], 1, 1) == "s"]))
val = sort(val)
if(length(val)!=5 | (length(val)==5 & sum(sort(val)==valobj)!=5)){
    com=read.csv2("Data/3_Bipartite/Communities_BACKUP.csv")
    print("BACKUP NEEDED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
}
################################################################################

# Export communities
if (!dir.exists("Analysis/Bioregions")) {
  dir.create("Analysis/Bioregions")
}
write.csv2(com, "Analysis/Bioregions/Communities.csv", row.names = FALSE)

# Reorder com
com <- com[match(idrow, com[, 1]), ]

print(c(length(idrow), length(com[, 1]), sum(idrow == com[, 1])))

# Weight species
w <- rep(1, length(idrow))

comp <- com[substr(idrow, 1, 1) == "p", 2] # Plant
wcomp <- aggregate(w[substr(idrow, 1, 1) == "p"], list(comp), sum)
wp <- wcomp[match(comp, wcomp[, 1]), 2]

coms <- com[substr(idrow, 1, 1) == "s", 2] # Spectral
wcoms <- aggregate(w[substr(idrow, 1, 1) == "s"], list(coms), sum)
ws <- wcoms[match(coms, wcoms[, 1]), 2]

w[substr(idrow, 1, 1) == "p"] <- wp
w[substr(idrow, 1, 1) == "s"] <- ws

w <- 1 / w

# Build Sandkey
colo <- brewer.pal(5, "Set2")

cat1 <- c("A", "B", "C", "D", "E")
cat2 <- c("F", "G", "H", "I", "J")
cat3 <- c("z1", "z2", "z3", "z4", "z5", "z6")

net <- bip[bip[, 3] > 0, ]

net[, 1] <- com[match(net[, 1], com[, 1]), 2]
net[, 2] <- com[match(net[, 2], com[, 1]), 2]
edges <- aggregate(net[, 3], list(net[, 1], net[, 2]), sum)
colnames(edges) <- c("N1", "N2", "Value")
edges <- edges[edges[, 1] > 0, ]
edges <- cbind(edges, id = 1:25, N3 = edges$N1)

edges <- rbind(cbind(edges, x = "N1", y = edges$N1), 
               cbind(edges, x = "N2", y = edges$N2))
edges$Value <- edges$Value / sum(edges$Value)
rownames(edges) <- 1:50
edges$N3[edges$N1 != edges$N2] <- 6

edges$N1 <- cat1[edges$N1]
edges$N2 <- cat2[edges$N2]
edges$y <- c(cat1[edges$y][1:25], cat2[edges$y][26:50])
edges$N3 <- cat3[edges$N3]

edges <- edges[order(edges$N3), ]

g <- ggplot(edges, aes(x, id = id, split = y, value = Value)) +
  scale_fill_manual(values = c(colo[c(1:5, 1:5, 1:5)], "grey")) +
  geom_parallel_sets(aes(fill = N3), alpha = 0.3, axis.width = 0.05, 
                     show.legend = FALSE) +
  geom_parallel_sets_axes(aes(fill = y), alpha = 0.3, axis.width = 0.05, 
                          show.legend = FALSE) +
  scale_x_discrete(expand = c(0.1, 0)) +
  theme_void()

# Export and load Sandkey
png("Sandkey.png", width = 15, height = 19, units = "cm", res = 600)
  g
dev.off()

g <- load.image("Sandkey.png")

# Figure S6 (plot Sandkey)
pdf("FigS6.pdf", width = 7, height = 7, useDingbats = FALSE)

  par(mar = c(0, 0, 0, 9))
    plot(g, xlab = "", ylab = "", axes = FALSE)
    legend("right", inset = c(-0.35, -0.011), fill = colo, border = colo, 
           legend = c("Com. 1", "Com. 2", "Com. 3", "Com. 4", "Com. 5"), 
           bty = "n", cex = 1.8, xpd = TRUE)
  
dev.off()

# Remove Sandkey file
file.remove("Sandkey.png")

# Color cells for bioregions
sp <- mat * w
sp <- aggregate(sp, list(com[, 2]), sum)
sp <- data.frame(idcol, t(sp[, -1]))

if (sum(com[, 2] == 0) > 0) {
  sp <- sp[, -2]
}

sp <- sp[apply(sp[, -1], 1, sum) > 0, ]

idsp <- sp[, 1]
sp <- sp[, -1]

sp <- 100 * sp / apply(sp, 1, sum)

idbio <- apply(sp, 1, function(x) {
  which(x == max(x))[1]
})
wbio <- apply(sp, 1, max)

colocom <- list()
colocom[[1]] <- alpha(colo[1], seq(0, 1, 0.01))
colocom[[2]] <- alpha(colo[2], seq(0, 1, 0.01))
colocom[[3]] <- alpha(colo[3], seq(0, 1, 0.01))
colocom[[4]] <- alpha(colo[4], seq(0, 1, 0.01))
colocom[[5]] <- alpha(colo[5], seq(0, 1, 0.01))

colocell <- NULL
for (i in 1:length(idbio)) {
  colocell <- c(colocell, colocom[[idbio[i]]][round(wbio[i])])
}

# Figure 3
pdf("Fig3.pdf", width = 10, height = 7, useDingbats = FALSE)
  
  par(mar = c(0, 0, 0, 9))
  plot(st_geometry(shp[match(as.character(idsp), shp$id), ]), 
       col = colocell, border = NA)
  legend("right", inset = c(-0.15, -0.011), fill = colo, border = colo, 
         legend = c("BR1", "BR2", "BR3", "BR4", "BR5"), 
         bty = "n", cex = 1.8, xpd = TRUE)
  
dev.off()

# Export Bioregions
temp <- data.frame(ID = idsp, BR = idbio)
br <- data.frame(ID = as.character(shp$id), BR = 0)
ind <- match(temp$ID, br$ID)
br[ind, 2] <- temp[, 2]
br <- br[order(br$ID), ]

write.csv2(br, "Analysis/Bioregions/Bioregions.csv", row.names = FALSE)

# Table S1
com <- read.csv2("Analysis/Bioregions/Communities.csv")
spcom <- read.csv2("Analysis/Bioregions/Bioregions.csv")
spcom <- spcom[spcom$BR > 0, ]
id <- substr(com[, 1], 1, 1)

tabcomp <- table(com[id == "p", 2])[-1]
tabcoms <- table(com[id == "s", 2])
tabcomsp <- table(spcom$BR)

tabcomp <- c(6650 - sum(tabcomp), tabcomp)
tabcoms <- c(250 - sum(tabcoms), tabcoms)
tabcomsp <- c(23060 - sum(tabcomsp), tabcomsp)

tex <- paste0(c("-", 1:5), " & ", tabcomp, " (", 
              round(100 * tabcomp / sum(tabcomp), digits = 2), "\\%)", 
              " & ", tabcoms, " (", 
              round(100 * tabcoms / sum(tabcoms), digits = 2), "\\%)", 
              " & ", tabcomsp, " (", 
              round(100 * tabcomsp / sum(tabcomsp), digits = 2), "\\%)", "\\")
print(data.frame(tex), row.names = FALSE)


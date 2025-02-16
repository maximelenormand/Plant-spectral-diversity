# Load packages
library(bioregion)
library(readr)

# Working directory
# setwd("")

# Import data
dat <- read_delim("Data/1_Plant/CoMat_Plant.csv",
  delim = ";",
  col_name = TRUE,
  show_col_types = FALSE
)

# Format comat
comat <- t(as.matrix(dat[, -1]))
rownames(comat) <- colnames(dat)[-1]
colnames(comat) <- as.character(as.matrix(dat[, 1]))

# Extract abc and ABC
abcABC <- similarity(comat, metric = c("abc", "ABC"))

# Export file
write_delim(abcABC, "Data/1_Plant/abcABC_Plant.csv",
  delim = ";",
  col_name = TRUE
)

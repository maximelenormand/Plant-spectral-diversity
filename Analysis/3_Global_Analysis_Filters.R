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
filters <- read_delim("Analysis/Filters/Filters.csv", delim = ";", 
                      col_name = TRUE, show_col_types = FALSE)

# Match between sites in beta & filters
idfiltersbeta1 <- match(beta_s$Site1, filters$ID)
idfiltersbeta2 <- match(beta_s$Site2, filters$ID)

gc()

# Correlations with Filter 0
Rs <- alpha_s$Rs
Rs <- Rs[filters$Filter0 == 1]
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

COR <- matrix(0, 4, 4)

temp <- cor.test(Rp, Rs, method = "pearson")
COR[1, 1] <- temp$estimate
COR[1, 2] <- temp$conf.int[1]
COR[1, 3] <- temp$conf.int[2]
COR[1, 4] <- temp$p.value

temp <- cor.test(Hp, Hs, method = "pearson")
COR[2, 1] <- temp$estimate
COR[2, 2] <- temp$conf.int[1]
COR[2, 3] <- temp$conf.int[2]
COR[2, 4] <- temp$p.value

temp <- cor.test(SIMp, SIMs, method = "pearson")
COR[3, 1] <- temp$estimate
COR[3, 2] <- temp$conf.int[1]
COR[3, 3] <- temp$conf.int[2]
COR[3, 4] <- temp$p.value

temp <- cor.test(BCbalp, BCbals, method = "pearson")
COR[4, 1] <- temp$estimate
COR[4, 2] <- temp$conf.int[1]
COR[4, 3] <- temp$conf.int[2]
COR[4, 4] <- temp$p.value

COR0 <- COR

gc()

# Correlations with Filter 1
Rs <- alpha_s$Rs
Rs <- Rs[filters$Filter1 == 1]
Rp <- alpha_p$Rp
Rp <- Rp[filters$Filter1 == 1]
Hs <- alpha_s$Hs
Hs <- Hs[filters$Filter1 == 1]
Hp <- alpha_p$Hp
Hp <- Hp[filters$Filter1 == 1]

SIMp <- beta_p$SIM
SIMp <- SIMp[filters$Filter1[idfiltersbeta1] == 1 & 
             filters$Filter1[idfiltersbeta2] == 1]
SIMs <- beta_s$SIM
SIMs <- SIMs[filters$Filter1[idfiltersbeta1] == 1 & 
             filters$Filter1[idfiltersbeta2] == 1]
BCbalp <- beta_p$BCbal
BCbalp <- BCbalp[filters$Filter1[idfiltersbeta1] == 1 & 
                 filters$Filter1[idfiltersbeta2] == 1]
BCbals <- beta_s$BCbal
BCbals <- BCbals[filters$Filter1[idfiltersbeta1] == 1 & 
                 filters$Filter1[idfiltersbeta2] == 1]

gc()

COR <- matrix(0, 4, 4)

temp <- cor.test(Rp, Rs, method = "pearson")
COR[1, 1] <- temp$estimate
COR[1, 2] <- temp$conf.int[1]
COR[1, 3] <- temp$conf.int[2]
COR[1, 4] <- temp$p.value

temp <- cor.test(Hp, Hs, method = "pearson")
COR[2, 1] <- temp$estimate
COR[2, 2] <- temp$conf.int[1]
COR[2, 3] <- temp$conf.int[2]
COR[2, 4] <- temp$p.value

temp <- cor.test(SIMp, SIMs, method = "pearson")
COR[3, 1] <- temp$estimate
COR[3, 2] <- temp$conf.int[1]
COR[3, 3] <- temp$conf.int[2]
COR[3, 4] <- temp$p.value

temp <- cor.test(BCbalp, BCbals, method = "pearson")
COR[4, 1] <- temp$estimate
COR[4, 2] <- temp$conf.int[1]
COR[4, 3] <- temp$conf.int[2]
COR[4, 4] <- temp$p.value

COR1 <- COR

gc()

# Correlations with Filter 2
Rs <- alpha_s$Rs
Rs <- Rs[filters$Filter2 == 1]
Rp <- alpha_p$Rp
Rp <- Rp[filters$Filter2 == 1]
Hs <- alpha_s$Hs
Hs <- Hs[filters$Filter2 == 1]
Hp <- alpha_p$Hp
Hp <- Hp[filters$Filter2 == 1]

SIMp <- beta_p$SIM
SIMp <- SIMp[filters$Filter2[idfiltersbeta1] == 1 & 
             filters$Filter2[idfiltersbeta2] == 1]
SIMs <- beta_s$SIM
SIMs <- SIMs[filters$Filter2[idfiltersbeta1] == 1 & 
             filters$Filter2[idfiltersbeta2] == 1]
BCbalp <- beta_p$BCbal
BCbalp <- BCbalp[filters$Filter2[idfiltersbeta1] == 1 & 
                 filters$Filter2[idfiltersbeta2] == 1]
BCbals <- beta_s$BCbal
BCbals <- BCbals[filters$Filter2[idfiltersbeta1] == 1 & 
                 filters$Filter2[idfiltersbeta2] == 1]

gc()

COR <- matrix(0, 4, 4)

temp <- cor.test(Rp, Rs, method = "pearson")
COR[1, 1] <- temp$estimate
COR[1, 2] <- temp$conf.int[1]
COR[1, 3] <- temp$conf.int[2]
COR[1, 4] <- temp$p.value

temp <- cor.test(Hp, Hs, method = "pearson")
COR[2, 1] <- temp$estimate
COR[2, 2] <- temp$conf.int[1]
COR[2, 3] <- temp$conf.int[2]
COR[2, 4] <- temp$p.value

temp <- cor.test(SIMp, SIMs, method = "pearson")
COR[3, 1] <- temp$estimate
COR[3, 2] <- temp$conf.int[1]
COR[3, 3] <- temp$conf.int[2]
COR[3, 4] <- temp$p.value

temp <- cor.test(BCbalp, BCbals, method = "pearson")
COR[4, 1] <- temp$estimate
COR[4, 2] <- temp$conf.int[1]
COR[4, 3] <- temp$conf.int[2]
COR[4, 4] <- temp$p.value

COR2 <- COR

gc()

# Table 1
temp=COR0[,4]
COR0[,4][temp<0.001]="(***)"
COR0[,4][temp>=0.001]="(ns)"
temp=COR2[,4]
COR2[,4][temp<0.001]="(***)"
COR2[,4][temp>=0.001]="(ns)"

tex=paste0(round(as.numeric(COR0[,1]),digits=3), 
           " [", round(as.numeric(COR0[,2]),digits=3),",", 
           round(as.numeric(COR0[,3]),digits=3),"] ",
           COR0[,4],
           " & ", 
           round(as.numeric(COR2[,1]),digits=3), 
           " [", round(as.numeric(COR2[,2]),digits=3),",", 
           round(as.numeric(COR2[,3]),digits=3),"] ", 
           COR2[,4],
           "\\")
print(data.frame(tex), row.names = FALSE)




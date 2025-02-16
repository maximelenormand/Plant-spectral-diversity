# Packages
library(readr)

# Option
options(scipen=10000)

# Working directory
#setwd("")

# Extract & export ALPHA diversity PLANT metrics ###############################
comat_p=read_delim("Data/1_Plant/CoMat_Plant.csv", delim=";", col_name=TRUE,
                   show_col_types = FALSE)
comat_p=t(as.matrix(comat_p[,-1]))
id=as.character(rownames(comat_p))

Rp=apply(comat_p>0,1,sum) # Rp

mat=comat_p               # Hp
mat=mat/apply(mat,1,sum)
mat=mat*log(mat)
mat[is.na(mat)]=0
Hp=apply(mat,1,sum)
Hp=-Hp/log(dim(mat)[2])

if(!dir.exists("Analysis/Metrics")){
  dir.create("Analysis/Metrics")
}
write_delim(data.frame(ID=id,Rp,Hp),
            "Analysis/Metrics/Metrics_Alpha_Plant.csv",
            delim=";", col_name=TRUE,)

# Extract & export ALPHA diversity SPECTRAL metrics ############################
comat_s=read_delim("Data/2_Spectral/CoMat_Spectral.csv", delim=";", 
                   col_name=TRUE,
                   show_col_types = FALSE)
comat_s=t(as.matrix(comat_s[,-1]))

Rs=apply(comat_s>0,1,sum)      # Rs

mat=comat_s                    # Hs
mat=mat/apply(mat,1,sum)
mat=mat*log(mat)
mat[is.na(mat)]=0
Hs=apply(mat,1,sum)
Hs=-Hs/log(dim(mat)[2])

write_delim(data.frame(ID=id,Rs,Hs),
            "Analysis/Metrics/Metrics_Alpha_Spectral.csv", 
            delim=";", col_name=TRUE)

gc()
rm(comat_p, comat_s, Rp, Rs, Hp, Hs, mat)
gc()

# Extract & export BETA diversity PLANT metrics ################################
abcABC_p=read_delim("Data/1_Plant/abcABC_Plant.csv", delim=";", col_name=TRUE,
                    show_col_types = FALSE)  

SOR=(abcABC_p$b+abcABC_p$c)/(2*abcABC_p$a+abcABC_p$b+abcABC_p$c)
SIM=pmin(abcABC_p$b,abcABC_p$c)/(abcABC_p$a+pmin(abcABC_p$b,abcABC_p$c))

BC=(abcABC_p$B+abcABC_p$C)/(2*abcABC_p$A+abcABC_p$B+abcABC_p$C)
BCbal=pmin(abcABC_p$B,abcABC_p$C)/(abcABC_p$A+pmin(abcABC_p$B,abcABC_p$C))

write_delim(data.frame(Site1=abcABC_p$Site1,Site2=abcABC_p$Site2,SOR,SIM,BC,BCbal),
            "Analysis/Metrics/Metrics_Beta_Plant.csv", 
            delim=";", col_name=TRUE)

gc()
rm(abcABC_p,SOR,SIM,BC,BCbal)
gc()

# Extract & export BETA diversity SPECTRAL metrics #############################
abcABC_s=read_delim("Data/2_Spectral/abcABC_Spectral.csv", delim=";", 
                    col_name=TRUE,
                    show_col_types = FALSE) 

SOR=(abcABC_s$b+abcABC_s$c)/(2*abcABC_s$a+abcABC_s$b+abcABC_s$c)
SIM=pmin(abcABC_s$b,abcABC_s$c)/(abcABC_s$a+pmin(abcABC_s$b,abcABC_s$c))

BC=(abcABC_s$B+abcABC_s$C)/(2*abcABC_s$A+abcABC_s$B+abcABC_s$C)
BCbal=pmin(abcABC_s$B,abcABC_s$C)/(abcABC_s$A+pmin(abcABC_s$B,abcABC_s$C))

write_delim(data.frame(Site1=abcABC_s$Site1,Site2=abcABC_s$Site2,SOR,SIM,BC,BCbal),
            "Analysis/Metrics/Metrics_Beta_Spectral.csv", 
            delim=";", col_name=TRUE)

gc()
rm(abcABC_s,SOR,SIM,BC,BCbal)
gc()






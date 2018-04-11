args = commandArgs(trailingOnly=TRUE)
# library(argparse)
library(data.table)
library(foreach)
library(doMC) # Unix only
registerDoMC(10)

beta = "1.0"
maf = '0.010.5'
snp_num = 500
sample_X = 50000
sample_Y = 20000
sim_num = 1000

rsqGX = 0.4
rsqXY = 0.05
rsqC = 0.2
rsqPG = as.numeric(args[1])

pleio_ratio=0.1
pleio_snp_num = round(snp_num*as.numeric(pleio_ratio), digits=0)

# message(rsqPG, "  Read files ...")
sim_para = paste0("beta",beta,"_maf",maf,"_snp",snp_num,"_sizeX",sample_X,"Y",sample_Y,"_sim",sim_num,"_pleio", pleio_ratio,"_rsqGX",rsqGX, "_rsqPG",rsqPG,"_rsqXY",rsqXY,"_rsqC", rsqC)
phenoX = fread(paste("gzip -dc", paste0(sim_para,"/phenoX.exp.xls.gz")), header=FALSE, skip=1, sep="\t")
phenoY = fread(paste("gzip -dc", paste0(sim_para,"/phenoY.out.xls.gz")), header=FALSE, skip=1, sep="\t")
phenoX2Y = fread(paste("gzip -dc", paste0(sim_para,"/phenoY.exp2out.xls.gz")), header=FALSE, skip=1, sep="\t")
genoX = fread(paste("gzip -dc", paste0(sim_para,"/genoG_standardised.exp.xls.gz")), header=FALSE, skip=1, sep="\t")
genoY = fread(paste("gzip -dc", paste0(sim_para,"/genoG_standardised.out.xls.gz")), header=FALSE, skip=1, sep="\t")
cx = fread(paste("gzip -dc", paste0(sim_para,"/c.exp.xls.gz")), header=FALSE, skip=1, sep="\t")
cy = fread(paste("gzip -dc", paste0(sim_para,"/c.out.xls.gz")), header=FALSE, skip=1, sep="\t")
rsqX_known = fread(paste("gzip -dc", paste0(sim_para,"/adj.RSQ_XG.exp.xls.gz")), header=FALSE, skip=1, sep="\t")
rsqY_known = fread(paste("gzip -dc", paste0(sim_para,"/adj.RSQ_YG.out.xls.gz")), header=FALSE, skip=1, sep="\t")


# Parallel test: 
# message("rsqGX = ")
# PctExp = c()
# PctExp <- foreach (j=1:sim_num) %dopar% {  
#   PctExpTMP = 0
#   for (i in 1:snp_num) {
#     fit = summary(lm(phenoX[,j] ~ genoX[,i]))
#     PctExpTMP = PctExpTMP+fit$adj.r.squared
#   }
#   PctExpTMP
# }
# rsqGX <- mean(as.numeric(PctExp))
rsqGX = 0 
for (i in 1:snp_num) { rsqGX = rsqGX + mean(as.numeric(rsqX_known[i,-1])) }
# message(rsqGX)


# message("rsqPX = ")
# PctExp = c()
# PctExp <- foreach (j=1:sim_num) %dopar% { 
#   PctExpTMP = 0
#   for (i in 1:pleio_snp_num) {
#     k = i + snp_num
#     fit = summary(lm(phenoX[,j] ~ genoX[,k]))
#     PctExpTMP = PctExpTMP+fit$adj.r.squared
#   }
#   PctExpTMP
# }
# rsqPX <- mean(as.numeric(PctExp))
rsqPX = 0 
for (i in 1:pleio_snp_num) { 
  k = i + snp_num
  rsqPX = rsqPX + mean(as.numeric(rsqX_known[k,-1])) 
}
# message(rsqPX)


# message("rsqPY = ")
# PctExp = c()
# PctExp <- foreach (j=1:sim_num) %dopar% { 
#   PctExpTMP = 0
#   for (i in 1:pleio_snp_num) {
#     k = i + snp_num
#     fit = summary(lm(phenoY[,j] ~ genoY[,k]))
#     PctExpTMP = PctExpTMP+fit$adj.r.squared
#   }
#   PctExpTMP
# }
# rsqPY <- mean(as.numeric(PctExp))
rsqPY = 0 
for (i in 1:pleio_snp_num) { 
  k = i + snp_num
  rsqPY = rsqPY + mean(as.numeric(rsqX_known[k,-1]))  
}
# message(rsqPY)


# message("rsqXY = ")
PctExp = c()
PctExp <- foreach (j=1:sim_num) %dopar% { 
    fit = summary(lm(phenoY[[j+1]] ~ phenoX2Y[[j+1]]))
    fit$adj.r.squared
}
rsqXY <- mean(as.numeric(PctExp))
# message(rsqXY)


# message("rsqCX = ")
PctExp = c()
PctExp <- foreach (j=1:sim_num) %dopar% { 
    fit = summary(lm(phenoX[[j+1]] ~ cx[[j+1]]))
    fit$adj.r.squared
}
rsqCX <- mean(as.numeric(PctExp))
# message(rsqCX)


# message("rsqCY = ")
PctExp = c()
PctExp <- foreach (j=1:sim_num) %dopar% { 
    fit = summary(lm(phenoY[[j+1]] ~ cy[[j+1]]))
    fit$adj.r.squared
}
rsqCY <- mean(as.numeric(PctExp))
# message(rsqCY)


message(paste(rsqPG,"-","RsqGX:",rsqGX, "RsqPX:",rsqPX, "RsqXY:",rsqXY, "RsqPY:",rsqPY, "RsqCX:", rsqCX, "RsqCY:",rsqCY))

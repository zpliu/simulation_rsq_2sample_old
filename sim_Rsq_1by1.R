# options(error=recover)
library(argparse)
library(foreach)
# library(doParallel)  # both Windows and Unix
library(doMC) # Unix only

#########################################################################
### MR Simulation: control Rsq 
###	Ref: YJ Nature Communication GSMR paper 

### Current function in this script:
###   + Simulate in Matrix format.
###   + Parallel computation;
###   + Include pleiotropic SNPs;
###     + Control Rsq(causal), Rsq(pleiotropy), Rsq(C) in exposure;
###     + Control Rsq(XY), Rsq(pleiotropy) in outcome
###   + Two different sample sizes for exposure and outcome

### X = g + pgX + c + eX, where
###   g = Wu, W is the causal SNP matrix, u is the effect (i.e. gamma in the script);
###   pgX = W(p)u(p), W(p) is the pleiotropic SNP matrix, u(p) is the effect;
###   SNP follows ~ Binomial(2,f), where f ~ Uniform(MAF_min, MAF_max);
###   each SNP effect u~ N(0, 1), i.e. both causal and pleiotropic SNP
###   c is the latent non-genetic confounding variable, c~N(0,var(c)), var(c)=var(g)*Rsq(c)/Rsq(g)
###   eX~N(0, var(eX)), var(eX) = var(g + pgX + c)*[1/(Rsq(g) + Rsq(pgX) + Rsq(c)) - 1]

### Y = bX + pgY + c + eY, where 
###   b is the causal effect (i.e. beta in the script)
###   pgY is the effect of pleiotropic SNP on outcome
###   c is the same latent non-genetic confounding variable for the same individual
###   eY~N(0, var(eY)), var(eY) = var(bX + pgY)*[1/(Rsq(XY) + Rsq(pgY) - 1] - var(c)*(1+2*beta), 
###           here, Rsq(pgY) is NOT equal to Rsq(pgX), Rsq(pgY) = var(pgY)*Rsq(XY)/var(bX)
#########################################################################

## Parse argument from command line
parser <- ArgumentParser()
parser$add_argument("--exposure_symbol", help="Symbol for exposure trait")
parser$add_argument("--outcome_symbol", help="Symbol for outcome trait")
parser$add_argument("--geno_symbol", help="Symbol for genotype")
parser$add_argument("--outdir", help="output directory")
parser$add_argument("--sample_X", type="integer", help="Sample size in exposure")
parser$add_argument("--sample_Y", type="integer", help="Sample size in outcome")
parser$add_argument("--sim_num", type="integer", help="Number of simulations")
parser$add_argument("--snp_num", type="integer", help="Number of causal SNPs in the simulation, excluding pleiotropic SNPs")
parser$add_argument("--min_maf", default=0.01, type="double", help="Minimum MAF:0.01 (default)")
parser$add_argument("--max_maf", default=0.5, type="double", help="Maximum MAF:0.5 (default)")
parser$add_argument("--rsq_GX", type="double", help="Variance explained by causal SNPs in exposure")
parser$add_argument("--rsq_PG", default=0, type="double", help="Variance explained by pleiotropic SNPs only in exposure: 0.0 (default)")
parser$add_argument("--rsq_XY", type="double", help="Variance explained by exposure in outcome")
parser$add_argument("--rsq_C", type="double", help="Variance explained by non-genetic confounding variable C in exposure")
parser$add_argument("--beta", type="double", help="Causal effect of exposure on outcome")
parser$add_argument("--thread", type="integer", help="multi-thread number")
parser$add_argument("--pleio_ratio", default=0, type="double", help="% Percentage of pleiotropic SNPs: 0.0 (default)")
args <- parser$parse_args()

registerDoMC(args$thread)  #change to your number of CPU cores

### Parsing parameters
exposure_symbol = args$exposure_symbol   # "Y2"
outcome_symbol = args$outcome_symbol  # "X2"
geno_symbol = args$geno_symbol  # "G2"
simuDIR = args$outdir
beta = args$beta

sample_X = args$sample_X 
sample_Y = args$sample_Y

sim_num = args$sim_num 
snp_num = args$snp_num 

MAF_min = args$min_maf # 0.1 # Allele frequency threshold
MAF_max = args$max_maf # 0.5

rsq_GX = args$rsq_GX
rsq_XY = args$rsq_XY
rsq_C = args$rsq_C
rsq_PG = args$rsq_PG
pleio_ratio = args$pleio_ratio

### Generate the causal and pleio SNPs effect
maf <- runif(snp_num, MAF_min, MAF_max)
var_gamma_X <- 1/(2*maf*(1-maf))
gamma_X <- rnorm(snp_num, 0, sqrt(var_gamma_X))

pleio_snp_num = 0
pleio_geno = c()
pleio_maf = c()
pleio_gamma_X = 0   # effect on exposure X
pleio_alpha_Y = 0   # effect on outcome Y
if (rsq_PG != 0) {
  pleio_snp_num <- round(snp_num*pleio_ratio, digits=0)
  pleio_maf <- runif(pleio_snp_num, MAF_min, MAF_max)

  var_pleio_effect <- 1/(2*pleio_maf*(1-pleio_maf))
  pleio_gamma_X <- rnorm(pleio_snp_num, 0, sqrt(var_pleio_effect))
  pleio_alpha_Y <- rnorm(pleio_snp_num, 0, sqrt(var_pleio_effect))
} 

sim_exposure <-function(sample_size) {
  res_out = list()
  # Generate causal SNPs
  geno_X <- sapply(1:snp_num, function(x) rbinom(sample_size, 2, p=maf[x]))
  geno_std_X <- sapply(1:snp_num, function(x) (geno_X[,x] - 2*maf[x])/sqrt(2*maf[x]*(1 - maf[x])))

  ### Generate the pleiotropic SNPs
  if (rsq_PG != 0) {
  	pleio_geno <- sapply(1:pleio_snp_num, function(x) rbinom(sample_size, 2, p=pleio_maf[x]))
  	pleio_geno_std <- sapply(1:pleio_snp_num, function(x) (pleio_geno[,x] - 2*pleio_maf[x])/sqrt(2*pleio_maf[x]*(1 - pleio_maf[x])))
  } 
    
  ### Generate exposure X for all individuals
  # g <- geno_X%*%gamma_X
  g <- geno_std_X%*%gamma_X
  var_g <- sd(g)^2     # var(zbzx) in GSMR paper

  pleio_g_X <- rep(0, times=sample_size)
  var_pleio_g_X <- 0
  if (rsq_PG != 0) {
  	# pleio_g_X <- pleio_geno%*%pleio_gamma_X
  	pleio_g_X <- pleio_geno_std%*%pleio_gamma_X
  	var_pleio_g_X <- sd(pleio_g_X)^2
  }
    
  var_c <- var_g*rsq_C/rsq_GX  ## This is important

  ### Obtain weight for pleiotropic effects
  pleio_weight_X = 0
  if (rsq_PG != 0) { pleio_weight_X = sqrt((rsq_PG*var_g)/(rsq_GX*var_pleio_g_X)) }

  var_epsilon_X <- (var_g + pleio_weight_X^2*var_pleio_g_X + var_c)*(1/(rsq_GX + rsq_PG + rsq_C) - 1) # from GSMR paper
  ## Equivalent to this: 
  ##  (1) var_epsilon_X <- (var_g+var_pleio_g_X)*(1/(rsq_GX+rsq_PG)-1) - var_c
  ##  (2) var_epsilon_X <- (var_g/rsq_GX)*(1-(rsq_GX+rsq_PG+rsq_C)) 
    
  epsilon_X <- foreach (j=1:sim_num, .combine=cbind) %dopar% { rnorm(sample_size, 0, sqrt(var_epsilon_X)) }
  cx <- foreach (j=1:sim_num, .combine=cbind) %dopar% { rnorm(sample_size, 0, sqrt(var_c)) }
  phenoX <- foreach (j=1:sim_num, .combine=cbind) %dopar% { g + pleio_weight_X*pleio_g_X + cx[,j] + epsilon_X[,j] }

  res_out[["var_c"]] <- var_c
  res_out[["cx"]] <- cx
  res_out[["phenoX"]] <- phenoX
  res_out[["geno_std"]] <- geno_std_X
  res_out[["pleio_geno_std"]] <- pleio_geno_std
  return(res_out)
}

### Generate exposure X for all individuals
message("Simulate phenotype X")
phenoX_s1 <- sim_exposure(sample_X)


### Generate outcome Y for all individuals
message("Simulate phenotype Y")
phenoX_s2 <- sim_exposure(sample_Y)

X_Y <- beta*phenoX_s2[["phenoX"]]
var_X_Y <- sd(X_Y)^2

pleio_g_Y <- rep(0, times=sample_Y)
var_pleio_g_Y <- 0
if (rsq_PG != 0) {
	# pleio_g_Y <- pleio_geno%*%pleio_alpha_Y
	pleio_g_Y <- phenoX_s2[["pleio_geno_std"]]%*%pleio_alpha_Y
	var_pleio_g_Y <- sd(pleio_g_Y)^2
}

### Obtain weight for pleiotropic effects
pleio_weight_Y = 0
if (rsq_PG != 0) { pleio_weight_Y = sqrt((rsq_PG*var_X_Y)/(rsq_XY*var_pleio_g_Y)) }

var_c <- phenoX_s2[["var_c"]]
var_epsilon_Y <- (var_X_Y + pleio_weight_Y^2*var_pleio_g_Y)*(1/(rsq_XY + rsq_PG) - 1) - var_c*(1+2*beta)

epsilon_Y <- foreach (j=1:sim_num, .combine=cbind) %dopar% { rnorm(sample_Y, 0, sqrt(var_epsilon_Y)) }
cy <- foreach (j=1:sim_num, .combine=cbind) %dopar% { rnorm(sample_Y, 0, sqrt(var_c)) }
phenoY <- foreach (j=1:sim_num, .combine=cbind) %dopar% { X_Y[,j] + pleio_weight_Y*pleio_g_Y + cy[,j] + epsilon_Y[,j] }



# -----------------
#  Output results
# -----------------
message("Writing results to files")
# Output log file
loginfo = paste(
                paste0("exposure_symbol=", exposure_symbol),
                paste0("outcome_symbol=", outcome_symbol),
                paste0("geno_symbol=", geno_symbol),
                paste0("sample_X=", sample_X),
                paste0("sample_Y=", sample_Y),
                paste0("sim_num=", sim_num),
                paste0("snp_num=", snp_num, "  # direct causal snp"),
                paste0("pleiotropic_snp_num=", pleio_snp_num),
                paste0("MAF_min=", MAF_min),
                paste0("MAF_max=", MAF_max),
                paste0("rsq_GX=", rsq_GX, " # variance explained by causal SNPs in exposure"),
                paste0("rsq_PG=", rsq_PG, " # variance explained by pleiotropic SNPs in exposure"),
                paste0("rsq_XY=", rsq_XY, " # variance explained by exposure in outcome"),
                paste0("rsq_C=", rsq_C, " # Variance explained by non-genetic confounding variable C"),
                paste0("beta=", beta, " # Causal effect"),
                paste0("Effec size of SNPs were generated from N(0,1/2f[1-f])"),
                sep="\n"
)
write.table(file=paste0(simuDIR, "/", geno_symbol,"_",exposure_symbol,"_",outcome_symbol,".log"), loginfo, col.names=FALSE, row.names=FALSE, quote=FALSE)

## Naming columns and rows
causal_snp_symbol = paste0(geno_symbol, "_", c(1:snp_num))
pleio_snp_symbol = c()
if (rsq_PG != 0) { pleio_snp_symbol = paste0(geno_symbol, "_pleio_", c(1:pleio_snp_num)) }

all_snps = c(causal_snp_symbol, pleio_snp_symbol)
all_simulations = paste0("sim", c(1:sim_num))

all_samples_X = c(1:sample_X)
all_samples_Y = c(1:sample_Y)

## Output genotype
# out_raw_geno_X = cbind(geno_X, pleio_geno)
# write.table(file=paste0(simuDIR, "/geno", geno_symbol, "_raw.xls"), out_raw_geno_X, col.names=all_snps, row.names=all_samples, quote=FALSE, sep="\t")
out_standardised_geno_X = cbind(phenoX_s1[["geno_std"]], phenoX_s1[["pleio_geno_std"]])
write.table(file=paste0(simuDIR, "/geno", geno_symbol, "_standardised.exp.xls"), out_standardised_geno_X, col.names=all_snps, row.names=all_samples_X, quote=FALSE, sep="\t")

out_standardised_geno_Y = cbind(phenoX_s2[["geno_std"]], phenoX_s2[["pleio_geno_std"]])
write.table(file=paste0(simuDIR, "/geno", geno_symbol, "_standardised.out.xls"), out_standardised_geno_Y, col.names=all_snps, row.names=all_samples_Y, quote=FALSE, sep="\t")

## Output phenotype
write.table(file=paste0(simuDIR, "/pheno", exposure_symbol, ".exp.xls"), phenoX_s1[["phenoX"]], col.names=all_simulations, row.names=all_samples_X, quote=FALSE, sep="\t")
write.table(file=paste0(simuDIR, "/pheno", outcome_symbol, ".out.xls"), phenoY, col.names=all_simulations, row.names=all_samples_Y, quote=FALSE, sep="\t")

## Output epsilon
# write.table(file=paste0(simuDIR, "/epsilon_", exposure_symbol, ".exp.xls"), epsilon_X, col.names=all_simulations, row.names=all_samples, quote=FALSE, sep="\t")
# write.table(file=paste0(simuDIR, "/epsilon_", outcome_symbol, ".out.xls"), epsilon_Y, col.names=all_simulations, row.names=all_samples, quote=FALSE, sep="\t")

## Output c
write.table(file=paste0(simuDIR, "/c.exp.xls"), phenoX_s1[["cx"]], col.names=all_simulations, row.names=all_samples_X, quote=FALSE, sep="\t")
write.table(file=paste0(simuDIR, "/c.out.xls"), cy, col.names=all_simulations, row.names=all_samples_Y, quote=FALSE, sep="\t")

## Output genotype ref and effects
causal_geno_info = cbind(maf, gamma_X, rep(0, times=snp_num))
pleio_geno_info = c()
if (rsq_PG != 0) { pleio_geno_info = cbind(pleio_maf, pleio_gamma_X, pleio_alpha_Y) }
geno_info = rbind(causal_geno_info, pleio_geno_info)
geno_info_col = c("MAF",paste0(exposure_symbol,".exp"), paste0(outcome_symbol,".out"))
write.table(file=paste0(simuDIR, "/geno_effect.xls"), geno_info, col.names=geno_info_col, row.names=all_snps, quote=FALSE, sep="\t")





## For test: Check is the Rsq is controlled
testbutton = 'NO'
if (testbutton == 'YES') {
sample_X = 1000
sample_Y = 500
sim_num = 1000
snp_num = 10
MAF_min = 0.01
MAF_max = 0.5
rsq_GX = 0.4
rsq_XY = 0.05
rsq_C = 0.2
rsq_PG=0.1
beta=1
pleio_ratio=0.1

exposure_symbol="X1"
outcome_symbol="Y1"
geno_symbol="G1"
# simuDIR=""

# phenoX1 = read.table(file="phenoX1.exp.xls", sep="\t", header=TRUE)
# phenoY1 = read.table(file="phenoY1.out.xls", sep="\t", header=TRUE)
# genoG1 = read.table(file="genoG1_raw.xls", sep="\t", header=TRUE)
# c = read.table(file="c_X1_Y1.xls", sep="\t", header=TRUE)


# Parallel test: 
PctExp = c()
PctExp <- foreach (j=1:sim_num) %dopar% {  
  PctExpTMP = 0
  for (i in 1:snp_num) {
    fit = summary(lm(phenoX[,j] ~ geno_std_X[,i]))
    PctExpTMP = PctExpTMP+fit$adj.r.squared
  }
  PctExpTMP
}
mean(as.numeric(PctExp))

PctExp = c()
PctExp <- foreach (j=1:sim_num) %dopar% { 
  PctExpTMP = 0
  for (i in 1:pleio_snp_num) {
    fit = summary(lm(phenoX[,j] ~ pleio_geno_std[,i]))
    PctExpTMP = PctExpTMP+fit$adj.r.squared
  }
  PctExpTMP
}
mean(as.numeric(PctExp))

PctExp = c()
PctExp <- foreach (j=1:sim_num) %dopar% { 
  PctExpTMP = 0
  for (i in 1:pleio_snp_num) {
    fit = summary(lm(phenoY[,j] ~ pleio_geno_std[,i]))
    PctExpTMP = PctExpTMP+fit$adj.r.squared
  }
  PctExpTMP
}
mean(as.numeric(PctExp))


PctExp = c()
PctExp <- foreach (j=1:sim_num) %dopar% { 
    fit = summary(lm(phenoY[,j] ~ phenoX[,j]))
    fit$adj.r.squared
}
mean(as.numeric(PctExp))

PctExp = c()
PctExp <- foreach (j=1:sim_num) %dopar% { 
    fit = summary(lm(phenoX[,j] ~ cx[,j]))
    fit$adj.r.squared
}
mean(as.numeric(PctExp))

## Test without parallel
# PctExp = c()
# for (j in 1:sim_num) { 
#   PctExpTMP = 0
#   for (i in 1:snp_num) {
#     fit = summary(lm(phenoX1[,j] ~ genoG1[,i]))
#     PctExpTMP = PctExpTMP+fit$adj.r.squared
#   }
#   PctExp = c(PctExp, PctExpTMP)
# }
# mean(PctExp)


# PctExp = c()
# for (j in 1:sim_num) { 
#   PctExpTMP = 0
#   for (i in 1:pleio_snp_num) {
#     k = i + snp_num
#     fit = summary(lm(phenoX1[,j] ~ genoG1[,k]))
#     PctExpTMP = PctExpTMP+fit$adj.r.squared
#   }
#   PctExp = c(PctExp, PctExpTMP)
# }
# mean(PctExp)
}
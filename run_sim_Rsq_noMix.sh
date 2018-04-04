#! /usr/bin/bash
################################################
## 1. Run the simulation;
## 2. Run MR analysis.

##  Condition:  Causal SNP + Pleiotropic SNP
#################################################

parallel="/home/zpliu/parallel/bin/parallel"
beta=$1         # Causal effect of exposure on outcome
min_maf=$2
max_maf=$3
snp_num=$4
sample_X=$5
sample_Y=$6
sim_num=$7
pleio_ratio=$8
rsq_GX=$9
rsq_PG=${10}
rsq_XY=${11}
rsq_C=${12}

thread=${13}

export sample_X; 
export sample_Y; 
export sim_num;
export snp_num;
export pleio_ratio;
export thread;
export rsq_GX;
export rsq_PG;
export rsq_XY;
export rsq_C;


function make_directory {
	if [ ! -d $1 ]; then 
		mkdir $1
	fi
}
export -f make_directory


## simulate two sub-phenotypes
function simulate {
	simDIR=$1
	Rscript sim_Rsq.R --outdir $simDIR \
						--sample_X $sample_X \
						--sample_Y $sample_Y \
						--sim_num $sim_num \
						--snp_num $snp_num \
						--pleio_ratio $pleio_ratio \
						--min_maf $min_maf \
						--max_maf $max_maf \
						--rsq_GX $rsq_GX \
						--rsq_PG $rsq_PG \
						--rsq_XY $rsq_XY \
						--rsq_C $rsq_C \
						--beta $beta \
						--thread $thread
}
export -f simulate

## MR for the mixed phenotypes
function TwoSampleMR {
	simDIR=$1
	Rscript TwoSampleMR.R  --beta_exp $simDIR/BETA_XG.exp.xls \
							--beta_out $simDIR/BETA_YG.out.xls \
							--se_exp $simDIR/SE_XG.exp.xls \
							--se_out $simDIR/SE_YG.out.xls \
							--pval_exp $simDIR/PVAL_XG.exp.xls \
							--pval_out $simDIR/PVAL_YG.out.xls \
							--sample_X $sample_X \
							--sample_Y $sample_Y \
							--sim_num $sim_num \
							--outprefix $simDIR/G_X_Y \
							--thread $thread
}
export -f TwoSampleMR


##--------------------------------
## Main function
##--------------------------------
simDIR="beta${beta}_maf${min_maf}${max_maf}_snp${snp_num}_sizeX${sample_X}Y${sample_Y}_sim${sim_num}_pleio${pleio_ratio}_rsqGX${rsq_GX}_rsqPG${rsq_PG}_rsqXY${rsq_XY}_rsqC${rsq_C}"
make_directory $simDIR
simulate $simDIR
TwoSampleMR $simDIR



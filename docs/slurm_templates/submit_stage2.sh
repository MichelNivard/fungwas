#!/bin/bash
#SBATCH --job-name=fungwas_stage2
#SBATCH --partition=short
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/stage2_%j.log

# FungWas Stage 2 SLURM Template
# Runs after Stage 1 completes

# Load modules
module load lang/R/4.5.1

# Configuration
STAGE1_DIR="results/stage1"
OUT_DIR="results/stage2"
WEIGHT_FILE="weights/my_weights.rds"

mkdir -p ${OUT_DIR}

# Merge Stage 1 outputs from all chromosomes
echo "Merging Stage 1 outputs..."

# Merge stage1 TSV files
zcat ${STAGE1_DIR}/chr*.stage1.tsv.gz | \
    awk 'NR==1 || !/^snp_id/' | \
    gzip > ${STAGE1_DIR}/all_chromosomes.stage1.tsv.gz

# Merge covariance files (simple concatenation)
cat ${STAGE1_DIR}/chr*.cov.gz > ${STAGE1_DIR}/all_chromosomes.cov.gz

# Run Stage 2
echo "Running Stage 2..."
Rscript -e "
library(fungwas)

# Load weights
weights <- readRDS('${WEIGHT_FILE}')
W <- weights\$W

# Run param_gwas with tau_cov SE mode
results <- param_gwas_from_file(
  '${STAGE1_DIR}/all_chromosomes.stage1.tsv.gz',
  W = W,
  se_mode = 'tau_cov',
  cov_file = '${STAGE1_DIR}/all_chromosomes.cov.gz'
)

# Save results
data.table::fwrite(results, '${OUT_DIR}/fungwas_results.tsv.gz', sep = '\t')
message('Stage 2 complete: ${OUT_DIR}/fungwas_results.tsv.gz')
"

echo "Stage 2 complete"

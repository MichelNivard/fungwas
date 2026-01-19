# Phenotype and Covariate Preparation

This document describes input requirements for FungWas analysis.

## Overview

FungWas Stage 1 takes a **raw phenotype** as input and computes RIF (Recentered Influence Function) values internally. You do NOT need to pre-compute RIF.

## Phenotype File Format

Standard REGENIE/PLINK-style format:

```
FID IID testosterone BMI height
1001 1001 15.3 24.5 175
1002 1002 12.1 26.2 168
1003 1003 NA 23.8 172
```

- **Space-delimited**
- First two columns: `FID`, `IID`
- Subsequent columns: phenotype(s)
- Missing values: `NA`

When running FungWas, specify which column to analyze:
```bash
fungwas-stage1 --pheno phenotypes.txt --pheno-col testosterone ...
```

## Covariate File Format

```
FID IID age sex PC1 PC2 PC3 PC4 PC5 ...
1001 1001 45 1 0.012 -0.034 0.056 -0.023 0.089 ...
1002 1002 52 0 -0.045 0.067 -0.012 0.045 -0.034 ...
```

Standard covariates:
- Age
- Sex
- Genetic principal components (PC1-PC20)
- Assessment center (if applicable)
- Genotyping batch

## Covariate Handling

FungWas uses the "residualize" approach:

1. Raw phenotype Y → RIF matrix computed
2. **RIF is residualized on covariates** (via QR decomposition)
3. **Genotypes are residualized on covariates** (separately)
4. RIF regressed on residualized genotypes

This differs from some GWAS software where Y is residualized before analysis.

## Quality Control Recommendations

Before running FungWas:

1. **Remove samples** with missing phenotype or covariates
2. **Apply any transformations** (e.g., log for skewed traits) before running
3. **Remove extreme outliers** if appropriate for your phenotype
4. **Ensure sample files match** between phenotype, covariate, and genetic data

```r
# Example QC in R
pheno <- read.table("phenotypes.txt", header=TRUE)
covar <- read.table("covariates.txt", header=TRUE)

# Keep complete cases
merged <- merge(pheno, covar, by=c("FID", "IID"))
complete <- merged[complete.cases(merged), ]

write.table(complete[, c("FID", "IID", "testosterone")], 
            "pheno_clean.txt", row.names=FALSE, quote=FALSE)
```

## Inverse Rank Normalization (Optional)

For heavily skewed phenotypes, consider INT before analysis:

```r
int_transform <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

pheno$testosterone_int <- int_transform(pheno$testosterone)
```

Then use `--pheno-col testosterone_int` in Stage 1.

## What FungWas Computes Internally

When you provide a raw phenotype, Stage 1:

1. Computes empirical quantiles at each tau level
2. Estimates density at each quantile (Silverman bandwidth)
3. Builds RIF matrix: `RIF = q + (tau - I(Y <= q)) / f(q)`
4. Residualizes RIF on covariates
5. Runs jackknife block regression for each SNP

This is equivalent to running `quantile_gwas(Y, G, taus)` in the R package.

# ============================================================
# Example: Load PLINK bed/bim/fam with bigsnpr
# Simulate phenotype
# Run Quantile GWAS with fungwas
# ============================================================

# Install if needed
# install.packages("bigsnpr")
# devtools::install_github("MichelNivard/fungwas")
devtools::document()
library(bigsnpr)
library(fungwasStage2)

set.seed(123)

# -------------------------
# 1. Load PLINK genotype data
# -------------------------
plink_prefix <- "/Users/michelnivard/Downloads/1KG.EAS.auto.snp.norm.nodup.split.rare002.common015.missing/1KG.EAS.auto.snp.norm.nodup.split.rare002.common015.missing"   # without extension
snp_readBed(paste0(plink_prefix, ".bed"))   # converts to .rds, RUNONCE
obj.bigSNP <- snp_attach(paste0(plink_prefix, ".rds"))



G   <- obj.bigSNP$genotypes     # FBM (genotype matrix)
G <- snp_fastImputeSimple(G)    # Quick impute missings...
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
SNP <- obj.bigSNP$map$marker.ID
N   <- nrow(G)
P   <- ncol(G)

cat("Loaded", N, "individuals and", P, "SNPs\n")

# -------------------------
# 2. Simulate phenotype
# -------------------------
set.seed(123)
causal <- sample(P, 5)
beta   <- rnorm(5, 0, 0.25)

genetic_score <- big_prodVec(G, beta, ind.col = causal)
y <- as.numeric(scale(genetic_score + rnorm(N)))

# -------------------------
# 3. Run Quantile GWAS (Stage 1)
# -------------------------
taus <- seq(0.1, 0.9, 0.1)

stage1 <- quantile_gwas(
  Y = y,
  G = as.matrix(G[,1:1000]),
  taus = taus,
  benchmark = FALSE,
  verbose = TRUE
)

# Inspect slopes
dim(stage1$Q_slope)   # T x P matrix of tau-slopes
head(stage1$Q_slope[,1:5])

# -------------------------
# 4. Parametric Mapping (Stage 2, example: vQTL model)
# -------------------------
W_vqtl <- make_weights_vqtl(taus, stage1$q_tau, mu = mean(y), sd = sd(y))

fit <- param_gwas(
  stage1,
  transform = "custom_W",
  transform_args = list(W = W_vqtl),
  se_mode = "dwls"
)

# Results
res <- data.frame(
  SNP = SNP[1:1000],
  CHR = CHR[1:1000],
  POS = POS[1:1000],
  beta_mu = fit$params["beta_mu", ],
  beta_sd = fit$params["beta_sigma2", ],
  se_mu = fit$SE_params["beta_mu", ],
  se_sd = fit$SE_params["beta_sigma2", ],
  Q       = fit$Q
)

head(res)

# -------------------------
# 5. Plot simple Manhattan for mean effect
# -------------------------
par(mar = c(4,4,2,1), bg = "white")
plot(res$POS, -log10(pchisq((res$beta_mu/res$SE_params["mu",])^2, df=1, lower.tail = FALSE)),
     col = res$CHR, pch = 20, cex = 0.6,
     xlab = "Genomic position", ylab = "-log10(p)",
     main = "Quantile GWAS (mu effect)")

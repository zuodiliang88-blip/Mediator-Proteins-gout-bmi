library(data.table)
library(dplyr)
library(coloc)

rm(list = ls())
gc()

dir.create("results", showWarnings = FALSE, recursive = TRUE)

pqtl_file <- "data/pqtl/example_protein_cis_1mb.txt"
gout_file <- "data/gwas/gout_summary_stats.txt.gz"
case_count <- 42034
window_bp <- 250000

# Read cis-pQTL data.
pqtl <- fread(pqtl_file)
pqtl <- pqtl[, c("CHR", "BP", "SNP", "A1", "A2", "EAF", "BETA", "SE", "P", "N")]
colnames(pqtl) <- c("CHR", "BP", "SNP", "A1", "A2", "freq", "beta", "se", "p", "n")
pqtl$varbeta <- pqtl$se^2
pqtl$z <- pqtl$beta / pqtl$se

lead_snp <- pqtl %>% arrange(p) %>% slice(1)
qtl_region <- pqtl %>%
  filter(CHR == lead_snp$CHR, between(BP, lead_snp$BP - window_bp, lead_snp$BP + window_bp)) %>%
  distinct(SNP, .keep_all = TRUE) %>%
  na.omit()

# Read gout GWAS data.
gwas <- fread(gout_file)
gwas$varbeta <- gwas$se^2
gwas$s <- case_count / gwas$n
gwas$z <- gwas$beta / gwas$se
gwas <- gwas %>% distinct(SNP, .keep_all = TRUE) %>% na.omit() %>% filter(p > 0)

shared_snps <- intersect(qtl_region$SNP, gwas$SNP)
qtl_region <- qtl_region[qtl_region$SNP %in% shared_snps, ]
gwas_region <- gwas[gwas$SNP %in% shared_snps, ]

dat <- merge(qtl_region, gwas_region, by = "SNP", suffixes = c("_pqtl", "_gwas"))

# Harmonise allele orientation.
swap <- dat$A1_pqtl == dat$A2_gwas & dat$A2_pqtl == dat$A1_gwas
dat$beta_gwas[swap] <- -dat$beta_gwas[swap]
dat$A1_gwas[swap] <- dat$A1_pqtl[swap]
dat$A2_gwas[swap] <- dat$A2_pqtl[swap]
dat <- dat[!duplicated(dat$SNP), ]

dat$MAF_pqtl <- ifelse(dat$freq_pqtl < 0.5, dat$freq_pqtl, 1 - dat$freq_pqtl)
dat$MAF_gwas <- ifelse(dat$freq_gwas < 0.5, dat$freq_gwas, 1 - dat$freq_gwas)

pqtl_form <- list(
  snp = dat$SNP,
  beta = dat$beta_pqtl,
  varbeta = dat$varbeta_pqtl,
  MAF = dat$MAF_pqtl,
  N = dat$n_pqtl,
  type = "quant"
)

gwas_form <- list(
  snp = dat$SNP,
  beta = dat$beta_gwas,
  varbeta = dat$varbeta_gwas,
  MAF = dat$MAF_gwas,
  N = dat$n_gwas,
  s = unique(dat$s)[1],
  type = "cc"
)

check_dataset(pqtl_form)
check_dataset(gwas_form)

coloc_res <- coloc.abf(dataset1 = pqtl_form, dataset2 = gwas_form)
write.table(as.data.frame(t(coloc_res$summary)), "results/06_protein_gout_coloc_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(coloc_res$results, "results/06_protein_gout_coloc_snp_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

library(TwoSampleMR)
library(ieugwasr)
library(plinkbinr)
library(dplyr)
library(MRPRESSO)

rm(list = ls())
gc()

dir.create("results", showWarnings = FALSE, recursive = TRUE)

bmi_file <- "data/gwas/bmi_top_snps.txt.gz"
gout_file <- "data/gwas/gout_summary_stats.txt.gz"
ld_reference <- "data/reference/1000G_EUR"
plink_bin <- get_plink_exe()

# Read BMI instruments.
bmi_exp <- read_exposure_data(
  filename = bmi_file,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq_Tested_Allele",
  effect_allele_col = "Tested_Allele",
  other_allele_col = "Other_Allele",
  pval_col = "P",
  samplesize_col = "N"
)

bmi_exp <- bmi_exp[bmi_exp$pval.exposure < 5e-8, ]

clump_input <- bmi_exp[, c("SNP", "pval.exposure")]
colnames(clump_input) <- c("rsid", "pval")
clump_input$pval <- as.numeric(clump_input$pval)

clumped <- ld_clump(
  dat = clump_input,
  plink_bin = plink_bin,
  bfile = ld_reference
)

bmi_exp <- bmi_exp[bmi_exp$SNP %in% clumped$rsid, ]
bmi_exp$F <- bmi_exp$beta.exposure^2 / bmi_exp$se.exposure^2
bmi_exp <- subset(bmi_exp, F >= 10)
bmi_exp$id.exposure <- "BMI"

# Read gout outcome data for selected instruments.
gout_out <- read_outcome_data(
  filename = gout_file,
  snps = bmi_exp$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "freq",
  pval_col = "p",
  samplesize_col = "n"
)

dat <- harmonise_data(exposure_dat = bmi_exp, outcome_dat = gout_out)
dat <- subset(dat, mr_keep == TRUE)
dat <- steiger_filtering(dat)
dat <- subset(dat, steiger_dir == TRUE)
rownames(dat) <- NULL

# Run MR-PRESSO and remove detected outliers when present.
dat_before_presso <- dat
presso_res <- run_mr_presso(dat_before_presso)
outlier_indices <- presso_res[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
presso_outlier_snps <- data.frame(SNP = character(0))

if (!is.null(outlier_indices)) {
  outlier_rows <- suppressWarnings(as.numeric(outlier_indices))
  outlier_rows <- outlier_rows[!is.na(outlier_rows)]
  outlier_rows <- outlier_rows[outlier_rows >= 1 & outlier_rows <= nrow(dat_before_presso)]
  if (length(outlier_rows) > 0) {
    presso_outlier_snps <- data.frame(SNP = dat_before_presso$SNP[outlier_rows])
    dat <- dat_before_presso[-outlier_rows, ]
  }
}

mr_res <- generate_odds_ratios(mr(dat))
het_res <- mr_heterogeneity(dat)
pleio_res <- mr_pleiotropy_test(dat)

write.table(mr_res, "results/01_bmi_gout_mr_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(het_res, "results/01_bmi_gout_heterogeneity.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pleio_res, "results/01_bmi_gout_pleiotropy.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(presso_outlier_snps, "results/01_bmi_gout_mrpresso_outliers.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)
library(plinkbinr)

rm(list = ls())
gc()

dir.create("results", showWarnings = FALSE, recursive = TRUE)

bmi_file <- "data/gwas/bmi_summary_stats.txt.gz"
outcome_dir <- "data/outcomes"
ld_reference <- "data/reference/1000G_EUR"
plink_bin <- get_plink_exe()

outcome_files <- list.files(outcome_dir, full.names = TRUE)

mr_results <- data.frame()
heterogeneity <- data.frame()
pleiotropy <- data.frame()

# Select independent BMI instruments.
bmi_exp <- read_exposure_data(
  filename = bmi_file,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "freq",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p",
  samplesize_col = "n"
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

# Run MR for each outcome trait.
for (i in seq_along(outcome_files)) {
  file_name <- outcome_files[i]
  message("Processing outcome ", i, " of ", length(outcome_files), ": ", basename(file_name))

  tryCatch({
    outcome_dat <- read_outcome_data(
      snps = bmi_exp$SNP,
      filename = file_name,
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

    outcome_id <- tools::file_path_sans_ext(basename(file_name))

    dat <- harmonise_data(exposure_dat = bmi_exp, outcome_dat = outcome_dat)
    dat$id.outcome <- outcome_id
    dat <- steiger_filtering(dat)
    dat <- subset(dat, mr_keep == TRUE)
    rownames(dat) <- NULL

    res <- generate_odds_ratios(mr(dat))
    res$id <- outcome_id

    het <- mr_heterogeneity(dat)
    het$id <- outcome_id

    pleio <- mr_pleiotropy_test(dat)
    pleio$id <- outcome_id

    mr_results <- rbind(mr_results, res)
    heterogeneity <- rbind(heterogeneity, het)
    pleiotropy <- rbind(pleiotropy, pleio)
  }, error = function(e) {
    message("Error in ", file_name, ": ", conditionMessage(e))
  })
}

write.table(mr_results, "results/02_bmi_to_multiple_traits_mr_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(heterogeneity, "results/02_bmi_to_multiple_traits_heterogeneity.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pleiotropy, "results/02_bmi_to_multiple_traits_pleiotropy.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

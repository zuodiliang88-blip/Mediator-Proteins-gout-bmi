library(plinkbinr)
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)
library(RadialMR)

rm(list = ls())
gc()

dir.create("results", showWarnings = FALSE, recursive = TRUE)

body_fat_file <- "data/gwas/body_fat_percentage.txt.gz"
protein_file_list <- "data/pqtl/protein_outcome_files.txt"
ld_reference <- "data/reference/1000G_EUR"
plink_bin <- get_plink_exe()

protein_files <- readLines(protein_file_list)

mr_results <- data.frame()
heterogeneity <- data.frame()
pleiotropy <- data.frame()

# Select independent instruments for body fat percentage.
body_fat_exp <- read_exposure_data(
  filename = body_fat_file,
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

body_fat_exp <- body_fat_exp[body_fat_exp$pval.exposure < 5e-8, ]
clump_input <- body_fat_exp[, c("SNP", "pval.exposure")]
colnames(clump_input) <- c("rsid", "pval")
clump_input$pval <- as.numeric(clump_input$pval)

clumped <- ld_clump(dat = clump_input, plink_bin = plink_bin, bfile = ld_reference)
body_fat_exp <- body_fat_exp[body_fat_exp$SNP %in% clumped$rsid, ]
body_fat_exp$F <- body_fat_exp$beta.exposure^2 / body_fat_exp$se.exposure^2
body_fat_exp <- subset(body_fat_exp, F >= 10)
body_fat_exp$id.exposure <- "body_fat_percentage"

# Test body fat percentage against each protein outcome.
for (file_name in protein_files) {
  tryCatch({
    protein_out <- read_outcome_data(
      snps = body_fat_exp$SNP,
      filename = file_name,
      sep = "\t",
      snp_col = "rsids",
      beta_col = "Beta",
      se_col = "SE",
      effect_allele_col = "effectAllele",
      other_allele_col = "otherAllele",
      eaf_col = "ImpMAF",
      pval_col = "Pval",
      samplesize_col = "N"
    )

    protein_id <- tools::file_path_sans_ext(basename(file_name))
    dat <- harmonise_data(exposure_dat = body_fat_exp, outcome_dat = protein_out)
    dat$id.outcome <- protein_id
    dat <- steiger_filtering(dat)
    dat <- subset(dat, mr_keep == TRUE)
    rownames(dat) <- NULL

    radial_input <- tsmr_to_rmr_format(dat)
    radial_ivw <- ivw_radial(radial_input, alpha = 0.05)
    outlier_snps <- radial_ivw$outliers$SNP[radial_ivw$outliers$p.value < 0.05]
    dat <- dat[!dat$SNP %in% outlier_snps, ]

    res <- generate_odds_ratios(mr(dat))
    res$id <- protein_id
    mr_results <- rbind(mr_results, res)

    het <- mr_heterogeneity(dat)
    het$I2 <- with(het, ifelse(Q == 0, 0, pmax(0, (Q - Q_df) / Q) * 100))
    het$id <- protein_id
    heterogeneity <- rbind(heterogeneity, het)

    pleio <- mr_pleiotropy_test(dat)
    pleio$id <- protein_id
    pleiotropy <- rbind(pleiotropy, pleio)
  }, error = function(e) {
    message("Error in ", file_name, ": ", conditionMessage(e))
  })
}

write.table(mr_results, "results/03_body_fat_protein_mr_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(heterogeneity, "results/03_body_fat_protein_heterogeneity.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pleiotropy, "results/03_body_fat_protein_pleiotropy.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

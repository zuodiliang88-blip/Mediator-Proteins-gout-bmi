library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)

rm(list = ls())
gc()

dir.create("results", showWarnings = FALSE, recursive = TRUE)

protein_list_file <- "data/pqtl/selected_protein_files.txt"
protein_dir <- "data/pqtl/clumped_instruments"
gout_file <- "data/gwas/gout_summary_stats.txt.gz"

protein_files <- file.path(protein_dir, readLines(protein_list_file))

mr_results <- data.frame()
heterogeneity <- data.frame()
pleiotropy <- data.frame()

gout_out <- read_outcome_data(
  filename = gout_file,
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

# Run MR for each selected protein instrument file.
for (file_name in protein_files) {
  tryCatch({
    protein_exp <- read_exposure_data(
      filename = file_name,
      sep = "\t",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      eaf_col = "EAF",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      pval_col = "P",
      samplesize_col = "N",
      gene_col = "Gene"
    )

    protein_exp$F <- protein_exp$beta.exposure^2 / protein_exp$se.exposure^2
    protein_exp <- subset(protein_exp, F >= 10)
    protein_exp$id.exposure <- protein_exp$gene.exposure

    dat <- harmonise_data(
      exposure_dat = protein_exp,
      outcome_dat = gout_out[gout_out$SNP %in% protein_exp$SNP, ]
    )

    dat$id.outcome <- "gout"
    dat <- steiger_filtering(dat)
    dat <- subset(dat, mr_keep == TRUE)

    protein_id <- tools::file_path_sans_ext(basename(file_name))

    if (nrow(dat) >= 1) {
      res <- generate_odds_ratios(mr(dat))
      res$id <- protein_id
      mr_results <- rbind(mr_results, res)
    }

    if (nrow(dat) >= 2) {
      het <- mr_heterogeneity(dat)
      het$I2 <- with(het, ifelse(Q == 0, 0, pmax(0, (Q - Q_df) / Q) * 100))
      het$id <- protein_id
      heterogeneity <- rbind(heterogeneity, het)
    }

    if (nrow(dat) >= 3) {
      pleio <- mr_pleiotropy_test(dat)
      pleio$id <- protein_id
      pleiotropy <- rbind(pleiotropy, pleio)
    }
  }, error = function(e) {
    message("Error in ", file_name, ": ", conditionMessage(e))
  })
}

write.table(mr_results, "results/02_protein_gout_mr_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(heterogeneity, "results/02_protein_gout_heterogeneity.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pleiotropy, "results/02_protein_gout_pleiotropy.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

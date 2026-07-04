library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)
library(plinkbinr)

rm(list = ls())
gc()

dir.create("results", showWarnings = FALSE, recursive = TRUE)

adiposity_file <- "data/gwas/epicardial_adipose_tissue.txt.gz"
exposure_label <- "epicardial_adipose_tissue"
protein_dir <- "data/pqtl/protein_outcomes"
ld_reference <- "data/reference/1000G_EUR"
plink_bin <- get_plink_exe()

protein_files <- list.files(protein_dir, full.names = TRUE)

mr_results <- data.frame()
heterogeneity <- data.frame()
pleiotropy <- data.frame()

# Select independent instruments for the adiposity trait.
# To analyse another adipose depot, replace adiposity_file and exposure_label.
adiposity_exp <- read_exposure_data(
  filename = adiposity_file,
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

adiposity_exp <- adiposity_exp[adiposity_exp$pval.exposure < 5e-8, ]
clump_input <- adiposity_exp[, c("SNP", "pval.exposure")]
colnames(clump_input) <- c("rsid", "pval")
clump_input$pval <- as.numeric(clump_input$pval)

clumped <- ld_clump(dat = clump_input, plink_bin = plink_bin, bfile = ld_reference)
adiposity_exp <- adiposity_exp[adiposity_exp$SNP %in% clumped$rsid, ]
adiposity_exp$F <- adiposity_exp$beta.exposure^2 / adiposity_exp$se.exposure^2
adiposity_exp <- subset(adiposity_exp, F >= 10)
adiposity_exp$id.exposure <- exposure_label

# Test the adiposity trait against each protein outcome.
for (file_name in protein_files) {
  tryCatch({
    protein_out <- read_outcome_data(
      snps = adiposity_exp$SNP,
      filename = file_name,
      sep = "\t",
      snp_col = "rsids",
      beta_col = "Beta",
      se_col = "SE",
      eaf_col = "effectAlleleFreq",
      effect_allele_col = "effectAllele",
      other_allele_col = "otherAllele",
      pval_col = "Pval",
      samplesize_col = "N"
    )

    protein_id <- tools::file_path_sans_ext(basename(file_name))
    dat <- harmonise_data(exposure_dat = adiposity_exp, outcome_dat = protein_out)
    dat$id.outcome <- protein_id
    dat <- steiger_filtering(dat)
    dat <- subset(dat, mr_keep == TRUE)

    res <- generate_odds_ratios(mr(dat))
    res$id <- protein_id

    het <- mr_heterogeneity(dat)
    het$id <- protein_id

    pleio <- mr_pleiotropy_test(dat)
    pleio$id <- protein_id

    mr_results <- rbind(mr_results, res)
    heterogeneity <- rbind(heterogeneity, het)
    pleiotropy <- rbind(pleiotropy, pleio)
  }, error = function(e) {
    message("Error in ", file_name, ": ", conditionMessage(e))
  })
}

write.table(mr_results, "results/04_adiposity_protein_mr_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(heterogeneity, "results/04_adiposity_protein_heterogeneity.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pleiotropy, "results/04_adiposity_protein_pleiotropy.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

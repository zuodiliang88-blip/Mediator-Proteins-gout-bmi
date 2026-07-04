library(TwoSampleMR)
library(ieugwasr)
library(plinkbinr)
library(dplyr)
library(MVMR)
library(MendelianRandomization)

rm(list = ls())
gc()

dir.create("results", showWarnings = FALSE, recursive = TRUE)

fat_file <- "data/gwas/fat_mass.txt.gz"
lean_file <- "data/gwas/lean_mass.txt.gz"
gout_file <- "data/gwas/gout_summary_stats.txt.gz"
ld_reference <- "data/reference/1000G_EUR"
plink_bin <- get_plink_exe()

p_thresh <- 5e-8
clump_kb <- 10000
clump_r2 <- 0.001

# Fat mass instruments.
fat_all <- read_exposure_data(fat_file, sep = "\t", snp_col = "SNP", beta_col = "beta", se_col = "se",
                              eaf_col = "freq", effect_allele_col = "A1", other_allele_col = "A2",
                              pval_col = "p", samplesize_col = "n")
fat_all$id.exposure <- "fat_mass"
fat_sig <- fat_all[fat_all$pval.exposure < p_thresh, ]
fat_snps <- fat_sig[, c("SNP", "pval.exposure")]
colnames(fat_snps) <- c("rsid", "pval")
fat_snps$pval <- as.numeric(fat_snps$pval)
fat_clumped <- ld_clump(fat_snps, clump_kb = clump_kb, clump_r2 = clump_r2, plink_bin = plink_bin, bfile = ld_reference)
fat_iv <- fat_all[fat_all$SNP %in% fat_clumped$rsid, ]

# Lean mass instruments.
lean_all <- read_exposure_data(lean_file, sep = "\t", snp_col = "SNP", beta_col = "beta", se_col = "se",
                               eaf_col = "freq", effect_allele_col = "A1", other_allele_col = "A2",
                               pval_col = "p", samplesize_col = "n")
lean_all$id.exposure <- "lean_mass"
lean_sig <- lean_all[lean_all$pval.exposure < p_thresh, ]
lean_snps <- lean_sig[, c("SNP", "pval.exposure")]
colnames(lean_snps) <- c("rsid", "pval")
lean_snps$pval <- as.numeric(lean_snps$pval)
lean_clumped <- ld_clump(lean_snps, clump_kb = clump_kb, clump_r2 = clump_r2, plink_bin = plink_bin, bfile = ld_reference)
lean_iv <- lean_all[lean_all$SNP %in% lean_clumped$rsid, ]

# Unified clumping across both exposure instrument sets.
fat_snps2 <- fat_iv[, c("SNP", "pval.exposure")]
colnames(fat_snps2) <- c("rsid", "pval")
lean_snps2 <- lean_iv[, c("SNP", "pval.exposure")]
colnames(lean_snps2) <- c("rsid", "pval")

combined_snps <- bind_rows(fat_snps2, lean_snps2) %>% distinct(rsid, .keep_all = TRUE)
final_snps <- ld_clump(combined_snps, clump_kb = clump_kb, clump_r2 = clump_r2, plink_bin = plink_bin, bfile = ld_reference)
snps_final <- final_snps$rsid

fat_final <- fat_all[fat_all$SNP %in% snps_final, ]
fat_final$id.exposure <- "fat_mass"
lean_final <- lean_all[lean_all$SNP %in% snps_final, ]
lean_final$id.exposure <- "lean_mass"
exposure_data <- rbind(fat_final, lean_final)

gout_out <- read_outcome_data(gout_file, snps = snps_final, sep = "\t", snp_col = "SNP", beta_col = "beta",
                              se_col = "se", effect_allele_col = "A1", other_allele_col = "A2",
                              eaf_col = "freq", pval_col = "p", samplesize_col = "n")

mv_data <- mv_harmonise_data(exposure_dat = exposure_data, outcome_dat = gout_out)
mv_input <- mr_mvinput(bx = mv_data$exposure_beta, bxse = mv_data$exposure_se,
                       by = mv_data$outcome_beta, byse = mv_data$outcome_se)

mvmr_ivw <- mr_mvivw(mv_input, model = "default", robust = TRUE, correl = FALSE, distribution = "normal")
mvmr_egger <- mr_mvegger(mv_input, correl = FALSE)
pleio <- pleiotropy_mvmr(mv_input, gencov = 0)
strength <- strength_mvmr(mv_input, gencov = 0)

capture.output(mvmr_ivw, file = "results/06_fat_lean_gout_mvmr_ivw.txt")
capture.output(mvmr_egger, file = "results/06_fat_lean_gout_mvmr_egger.txt")
capture.output(pleio, file = "results/06_fat_lean_gout_mvmr_pleiotropy.txt")
capture.output(strength, file = "results/06_fat_lean_gout_mvmr_strength.txt")

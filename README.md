# Reviewer Code

This directory contains reviewer-facing analysis scripts for the manuscript. The scripts use relative paths and placeholder input files so that no local server paths or user-specific directories are exposed.

## Directory Layout

- `R/01_bmi_gout_mr.R`: two-sample MR for BMI and gout.
- `R/02_bmi_to_multiple_traits_mr.R`: two-sample MR from BMI to multiple outcome traits. If the outcome folder contains protein GWAS files, this script is used for BMI-to-protein MR.
- `R/03_protein_gout_mr.R`: two-sample MR for selected proteins and gout.
- `R/04_adiposity_protein_mr.R`: MR from an adiposity trait to protein traits. The default exposure is epicardial adipose tissue; replace the exposure file and `exposure_label` for other adipose depots.
- `R/05_body_fat_protein_mr.R`: MR from body fat percentage to protein traits.
- `R/06_fat_lean_gout_mvmr.R`: multivariable MR for fat mass and lean mass on gout.
- `R/07_protein_gout_coloc.R`: coloc analysis for a protein locus and gout.
- `R/08_ukb_association_models.R`: UK Biobank Cox and linear association models.
- `scripts/09_smr_gtex.sh`: SMR analysis across GTEx tissues.

The observational association script uses `INHBB` and `INHBC` as the selected Olink proteins.

## Expected Input Folders

Place input files under the following relative directories before running:

- `data/gwas/`: GWAS summary statistics.
- `data/pqtl/`: protein QTL summary statistics.
- `data/eqtl/`: eQTL reference files for SMR.
- `data/ukb/`: individual-level phenotype, covariate, proteomics, and follow-up files.
- `data/reference/`: LD reference panel for clumping and SMR.
- `results/`: output directory.

## Notes

The scripts are intended to document the analysis workflow for review. Some source datasets are controlled-access and are not redistributed with this repository.

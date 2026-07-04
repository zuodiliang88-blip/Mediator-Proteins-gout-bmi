library(readr)
library(dplyr)
library(data.table)
library(survival)

rm(list = ls())
gc()

dir.create("results", showWarnings = FALSE, recursive = TRUE)

gout_file <- "data/ukb/gout_first_occurrence.csv"
bmi_file <- "data/ukb/bmi.csv"
protein_file <- "data/ukb/proteomics_olink_instance0.csv"
covariate_file <- "data/ukb/covariates.csv"
baseline_file <- "data/ukb/baseline_date.csv"
death_file <- "data/ukb/death_date.csv"
loss_file <- "data/ukb/loss_to_follow_up.csv"

end_of_followup <- as.Date("2022-10-30")
selected_proteins <- c("ADH1B", "AP3B1", "APOE", "ENTPD6", "GRP", "RBP1")

gout <- fread(gout_file)
bmi <- read_csv(bmi_file, show_col_types = FALSE)
colnames(bmi)[2] <- "BMI"
bmi <- na.omit(bmi)

protein <- fread(protein_file)
colnames(protein)[1] <- "eid"
for (col in 2:ncol(protein)) {
  names(protein)[col] <- sub(";.*", "", names(protein)[col])
}
protein <- data.frame(protein)
protein <- protein[, c("eid", selected_proteins)]

covariates <- read_csv(covariate_file, show_col_types = FALSE)
factor_vars <- c("sex", "Ethnicity", "house", "alcohol", "smoking", "education",
                 "hypertension", "heartFailure", "renalFailure", "asthma",
                 "dementia", "MI", "stroke", "copd", "diabetes", "healthy_diet")
covariates[factor_vars] <- lapply(covariates[factor_vars], factor)

baseline <- read_csv(baseline_file, show_col_types = FALSE)[, 1:2]
colnames(baseline)[2] <- "beginTime"
baseline$beginTime <- as.Date(baseline$beginTime)

death <- unique(read_csv(death_file, show_col_types = FALSE))
begin_death <- left_join(baseline, death, by = "eid")
begin_death[is.na(begin_death$date_of_death), "date_of_death"] <- as.Date("3000-01-01")

loss <- na.omit(read_csv(loss_file, show_col_types = FALSE))

followup <- left_join(begin_death, gout, by = "eid")
followup$minTime <- as.Date(followup$minTime)
followup[is.na(followup$minTime), "minTime"] <- as.Date("3000-01-01")
followup <- na.omit(followup)

controls <- filter(followup, minTime > end_of_followup)
controls$status <- 0
controls$minTime <- end_of_followup
cases <- anti_join(followup, controls, by = "eid")
cases$status <- 1

analysis_time <- full_join(controls, cases) %>%
  mutate(enddate = pmin(as.Date(date_of_death), as.Date(minTime)),
         time = as.numeric(difftime(enddate, beginTime, units = "days")) / 365.25) %>%
  filter(time > 0) %>%
  select(eid, status, time)

analysis_base <- inner_join(analysis_time, covariates, by = "eid") %>%
  anti_join(loss, by = "eid")

# BMI and incident gout.
bmi_gout <- inner_join(analysis_base, bmi, by = "eid")
bmi_cox <- coxph(Surv(time, status) ~ BMI + sex + age + Ethnicity + TDindex + MET +
                   alcohol + smoking + education + hypertension + renalFailure +
                   diabetes + healthy_diet + ua + hdl + ldl + CRP + glucose +
                   Cholesterol + HbA1c + triglyceride + creatinine,
                 data = bmi_gout)

bmi_coef <- data.frame(summary(bmi_cox)$coefficients)
bmi_ci <- data.frame(exp(confint(bmi_cox, level = 0.95)))
bmi_coef$exp_coef <- exp(bmi_coef[, 1])
bmi_coef$lower_bound <- bmi_ci[, 1]
bmi_coef$upper_bound <- bmi_ci[, 2]
bmi_coef$Case <- sum(bmi_gout$status == 1)
bmi_coef$total <- nrow(bmi_gout)
write.table(bmi_coef, "results/07_bmi_gout_cox.tsv", sep = "\t", quote = FALSE)

# BMI and protein linear models.
linear_results <- data.frame()
for (j in 2:ncol(protein)) {
  exposure <- protein[, c(1, j)]
  colnames(exposure)[2] <- "protein"
  exposure <- na.omit(exposure)

  dat <- inner_join(covariates, exposure, by = "eid")
  fit <- lm(BMI ~ protein + sex + age + Ethnicity + TDindex + MET + alcohol + smoking + education, data = dat)

  res <- data.frame(summary(fit)$coefficients)
  ci <- data.frame(confint(fit, level = 0.95))
  res$variable <- rownames(res)
  res$lower <- ci[, 1]
  res$upper <- ci[, 2]
  res$total <- nrow(dat)
  res$protein <- names(protein)[j]
  linear_results <- rbind(linear_results, res)
}
write.table(linear_results, "results/07_bmi_protein_linear.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Protein and incident gout Cox models.
protein_cox_results <- data.frame()
for (j in 2:ncol(protein)) {
  exposure <- protein[, c(1, j)]
  colnames(exposure)[2] <- "protein"
  exposure <- na.omit(exposure)

  dat <- inner_join(analysis_base, exposure, by = "eid")
  fit <- coxph(Surv(time, status) ~ protein + sex + age + Ethnicity + TDindex + MET +
                 alcohol + smoking + education + hypertension + renalFailure +
                 diabetes + healthy_diet + ua + hdl + ldl + CRP + glucose +
                 Cholesterol + HbA1c + triglyceride + creatinine,
               data = dat)

  res <- data.frame(summary(fit)$coefficients)
  ci <- data.frame(exp(confint(fit, level = 0.95)))
  res$exp_coef <- exp(res[, 1])
  res$lower_bound <- ci[, 1]
  res$upper_bound <- ci[, 2]
  res$Case <- sum(dat$status == 1)
  res$total <- nrow(dat)
  res$protein <- names(protein)[j]
  protein_cox_results <- rbind(protein_cox_results, res)
}
write.table(protein_cox_results, "results/07_protein_gout_cox.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

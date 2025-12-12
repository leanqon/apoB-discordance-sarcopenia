# ===================================================================
# ApoB/Non-HDL-C Discordance and Sarcopenia: NHANES 2011-2016 Analysis
# ===================================================================
#
# Description: Analysis of the association between apolipoprotein B and
#              non-HDL cholesterol discordance and sarcopenia risk in
#              middle-aged US adults
#
# Data Source: National Health and Nutrition Examination Survey (NHANES)
#              2011-2012, 2013-2014, 2015-2016
#              https://www.cdc.gov/nchs/nhanes/
#
# Authors: [Your Team]
# Date: 2025
# License: MIT
# ===================================================================

# Clear environment
rm(list = ls())
gc()

# ===================================================================
# 1. LOAD REQUIRED PACKAGES
# ===================================================================

required_packages <- c(
  "haven",        # Read XPT files
  "tidyverse",    # Data manipulation
  "survey",       # Complex survey design
  "tableone",     # Table 1 creation
  "ggplot2",      # Visualization
  "splines",      # Restricted cubic splines
  "broom",        # Tidy model outputs
  "patchwork",    # Combine plots
  "rms"           # Regression modeling strategies
)

# Install and load packages
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Handle namespace conflicts
if(requireNamespace("conflicted", quietly = TRUE)) {
  library(conflicted)
  conflict_prefer("select", "dplyr")
  conflict_prefer("filter", "dplyr")
}

# Set survey options
options(survey.lonely.psu = "adjust")
options(digits = 4)

cat("=== NHANES ApoB/Non-HDL-C Discordance and Sarcopenia Analysis ===\n")
cat("Analysis started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ===================================================================
# 2. SET WORKING DIRECTORY
# ===================================================================
# NOTE: Update this path to your NHANES data directory
# The directory should contain the downloaded XPT files from NHANES

# setwd("path/to/your/NHANES/data")  # UPDATE THIS PATH

# ===================================================================
# 3. DATA IMPORT FUNCTIONS
# ===================================================================

#' Read NHANES 2011-2012 cycle data
read_2011_2012 <- function() {
  cat("Reading 2011-2012 data...\n")

  # Demographics and examination data
  demo_g <- read_xpt("DEMO_G.xpt")
  bmx_g <- read_xpt("BMX_G.xpt")
  dxx_g <- read_xpt("DXX_G.xpt")

  # Laboratory data
  apob_g <- read_xpt("APOB_G.xpt")
  tchol_g <- read_xpt("TCHOL_G.xpt")
  hdl_g <- read_xpt("HDL_G.xpt")
  trigly_g <- read_xpt("TRIGLY_G.xpt")
  glu_g <- read_xpt("GLU_G.xpt")
  ghb_g <- read_xpt("GHB_G.xpt")
  biopro_g <- read_xpt("BIOPRO_G.xpt")

  # Questionnaire data
  diq_g <- read_xpt("DIQ_G.xpt")
  bpq_g <- read_xpt("BPQ_G.xpt")
  mcq_g <- read_xpt("MCQ_G.xpt")
  smq_g <- read_xpt("SMQ_G.xpt")
  alq_g <- read_xpt("ALQ_G.xpt")
  paq_g <- read_xpt("PAQ_G.xpt")

  # Merge all datasets
  data_g <- demo_g %>%
    left_join(bmx_g, by = "SEQN") %>%
    left_join(dxx_g, by = "SEQN") %>%
    left_join(apob_g, by = "SEQN") %>%
    left_join(tchol_g, by = "SEQN") %>%
    left_join(hdl_g, by = "SEQN") %>%
    left_join(trigly_g, by = "SEQN") %>%
    left_join(glu_g, by = "SEQN") %>%
    left_join(ghb_g, by = "SEQN") %>%
    left_join(biopro_g, by = "SEQN") %>%
    left_join(diq_g, by = "SEQN") %>%
    left_join(bpq_g, by = "SEQN") %>%
    left_join(mcq_g, by = "SEQN") %>%
    left_join(smq_g, by = "SEQN") %>%
    left_join(alq_g, by = "SEQN") %>%
    left_join(paq_g, by = "SEQN") %>%
    mutate(CYCLE = "2011-2012")

  return(data_g)
}

#' Read NHANES 2013-2014 cycle data
read_2013_2014 <- function() {
  cat("Reading 2013-2014 data...\n")

  demo_h <- read_xpt("DEMO_H.xpt")
  bmx_h <- read_xpt("BMX_H.xpt")
  dxx_h <- read_xpt("DXX_H.xpt")
  apob_h <- read_xpt("APOB_H.xpt")
  tchol_h <- read_xpt("TCHOL_H.xpt")
  hdl_h <- read_xpt("HDL_H.xpt")
  trigly_h <- read_xpt("TRIGLY_H.xpt")
  glu_h <- read_xpt("GLU_H.xpt")
  ins_h <- read_xpt("INS_H.xpt")
  ghb_h <- read_xpt("GHB_H.xpt")
  biopro_h <- read_xpt("BIOPRO_H.xpt")
  diq_h <- read_xpt("DIQ_H.xpt")
  bpq_h <- read_xpt("BPQ_H.xpt")
  mcq_h <- read_xpt("MCQ_H.xpt")
  smq_h <- read_xpt("SMQ_H.xpt")
  alq_h <- read_xpt("ALQ_H.xpt")
  paq_h <- read_xpt("PAQ_H.xpt")

  data_h <- demo_h %>%
    left_join(bmx_h, by = "SEQN") %>%
    left_join(dxx_h, by = "SEQN") %>%
    left_join(apob_h, by = "SEQN") %>%
    left_join(tchol_h, by = "SEQN") %>%
    left_join(hdl_h, by = "SEQN") %>%
    left_join(trigly_h, by = "SEQN") %>%
    left_join(glu_h, by = "SEQN") %>%
    left_join(ins_h, by = "SEQN") %>%
    left_join(ghb_h, by = "SEQN") %>%
    left_join(biopro_h, by = "SEQN") %>%
    left_join(diq_h, by = "SEQN") %>%
    left_join(bpq_h, by = "SEQN") %>%
    left_join(mcq_h, by = "SEQN") %>%
    left_join(smq_h, by = "SEQN") %>%
    left_join(alq_h, by = "SEQN") %>%
    left_join(paq_h, by = "SEQN") %>%
    mutate(CYCLE = "2013-2014")

  return(data_h)
}

#' Read NHANES 2015-2016 cycle data
read_2015_2016 <- function() {
  cat("Reading 2015-2016 data...\n")

  demo_i <- read_xpt("DEMO_I.xpt")
  bmx_i <- read_xpt("BMX_I.xpt")
  dxx_i <- read_xpt("DXX_I.xpt")
  apob_i <- read_xpt("APOB_I.xpt")
  tchol_i <- read_xpt("TCHOL_I.xpt")
  hdl_i <- read_xpt("HDL_I.xpt")
  trigly_i <- read_xpt("TRIGLY_I.xpt")
  glu_i <- read_xpt("GLU_I.xpt")
  ins_i <- read_xpt("INS_I.xpt")
  ghb_i <- read_xpt("GHB_I.xpt")
  biopro_i <- read_xpt("BIOPRO_I.xpt")
  diq_i <- read_xpt("DIQ_I.xpt")
  bpq_i <- read_xpt("BPQ_I.xpt")
  mcq_i <- read_xpt("MCQ_I.xpt")
  smq_i <- read_xpt("SMQ_I.xpt")
  alq_i <- read_xpt("ALQ_I.xpt")
  paq_i <- read_xpt("PAQ_I.xpt")

  data_i <- demo_i %>%
    left_join(bmx_i, by = "SEQN") %>%
    left_join(dxx_i, by = "SEQN") %>%
    left_join(apob_i, by = "SEQN") %>%
    left_join(tchol_i, by = "SEQN") %>%
    left_join(hdl_i, by = "SEQN") %>%
    left_join(trigly_i, by = "SEQN") %>%
    left_join(glu_i, by = "SEQN") %>%
    left_join(ins_i, by = "SEQN") %>%
    left_join(ghb_i, by = "SEQN") %>%
    left_join(biopro_i, by = "SEQN") %>%
    left_join(diq_i, by = "SEQN") %>%
    left_join(bpq_i, by = "SEQN") %>%
    left_join(mcq_i, by = "SEQN") %>%
    left_join(smq_i, by = "SEQN") %>%
    left_join(alq_i, by = "SEQN") %>%
    left_join(paq_i, by = "SEQN") %>%
    mutate(CYCLE = "2015-2016")

  return(data_i)
}

# ===================================================================
# 4. IMPORT AND COMBINE DATA
# ===================================================================

cat("Importing NHANES data...\n")
data_2011 <- read_2011_2012()
data_2013 <- read_2013_2014()
data_2015 <- read_2015_2016()

cat("Combining cycles...\n")
nhanes_combined <- bind_rows(data_2011, data_2013, data_2015)
cat("Total participants:", nrow(nhanes_combined), "\n\n")

# ===================================================================
# 5. VARIABLE CREATION AND DATA CLEANING
# ===================================================================

cat("Creating variables and cleaning data...\n")

nhanes_clean <- nhanes_combined %>%
  mutate(
    # Age and demographics
    age = RIDAGEYR,
    sex = case_when(RIAGENDR == 1 ~ "Male", RIAGENDR == 2 ~ "Female"),

    # Race/ethnicity
    race = case_when(
      RIDRETH3 == 1 ~ "Mexican American",
      RIDRETH3 == 2 ~ "Other Hispanic",
      RIDRETH3 == 3 ~ "Non-Hispanic White",
      RIDRETH3 == 4 ~ "Non-Hispanic Black",
      TRUE ~ "Other"
    ),

    # Education
    education = case_when(
      DMDEDUC2 %in% c(1, 2) ~ "<High school",
      DMDEDUC2 == 3 ~ "High school",
      DMDEDUC2 %in% c(4, 5) ~ ">High school",
      TRUE ~ NA_character_
    ),

    # Poverty Income Ratio
    PIR = INDFMPIR,

    # Anthropometrics
    height = BMXHT,
    weight = BMXWT,
    BMI = BMXBMI,
    waist_circumference = BMXWAIST,

    # Body composition (DXA)
    ASM_kg = (DXXLGBM + DXXRGBM) / 1000,  # kg
    ASM_BMI = ASM_kg / (BMI),

    # Sarcopenia definition (FNIH criteria)
    sarcopenia = case_when(
      sex == "Male" & ASM_BMI < 0.789 ~ 1,
      sex == "Female" & ASM_BMI < 0.512 ~ 1,
      TRUE ~ 0
    ),

    # Lipid panel
    apoB = LBXAPB,
    total_cholesterol = LBXTC,
    hdl_c = LBDHDD,
    triglycerides = LBXTR,
    non_hdl_c = total_cholesterol - hdl_c,

    # Glucose metabolism
    fasting_glucose = LBXGLU,
    HbA1c = LBXGH,

    # Kidney and liver function
    creatinine = LBXSCR,
    eGFR = case_when(
      sex == "Male" & creatinine <= 0.9 ~ 141 * (creatinine/0.9)^(-0.411) * 0.993^age,
      sex == "Male" & creatinine > 0.9 ~ 141 * (creatinine/0.9)^(-1.209) * 0.993^age,
      sex == "Female" & creatinine <= 0.7 ~ 144 * (creatinine/0.7)^(-0.329) * 0.993^age,
      sex == "Female" & creatinine > 0.7 ~ 144 * (creatinine/0.7)^(-1.209) * 0.993^age
    ),
    albumin = LBXSAL,

    # Comorbidities
    diabetes = case_when(
      DIQ010 == 1 | HbA1c >= 6.5 | fasting_glucose >= 126 ~ 1,
      TRUE ~ 0
    ),

    hypertension = case_when(
      BPQ020 == 1 | BPXSY1 >= 140 | BPXDI1 >= 90 ~ 1,
      TRUE ~ 0
    ),

    CVD = case_when(
      MCQ160B == 1 | MCQ160C == 1 | MCQ160D == 1 | MCQ160E == 1 | MCQ160F == 1 ~ 1,
      TRUE ~ 0
    ),

    # Lifestyle factors
    smoking = case_when(
      SMQ040 %in% c(1, 2) ~ "Current",
      SMQ040 == 3 ~ "Former",
      TRUE ~ "Never"
    ),

    alcohol_use = case_when(ALQ101 == 1 ~ 1, TRUE ~ 0),

    physical_activity = case_when(
      PAQ605 == 1 ~ 1,
      PAQ605 == 2 ~ 0,
      TRUE ~ NA_real_
    ),

    # Survey weights (adjusted for 3 cycles)
    weight_final = WTMEC2YR / 3,
    SDMVPSU = SDMVPSU,
    SDMVSTRA = SDMVSTRA
  ) %>%
  filter(
    # Age restriction: 40-59 years
    age >= 40 & age <= 59,

    # Complete data for key variables
    !is.na(sarcopenia),
    !is.na(apoB),
    !is.na(non_hdl_c),
    !is.na(BMI),

    # Exclude pregnancy
    RIDEXPRG != 1 | is.na(RIDEXPRG),

    # Biologically plausible values
    BMI >= 16 & BMI <= 50,
    apoB >= 30 & apoB <= 300,

    # Valid survey weights
    weight_final > 0
  )

cat("Final sample size:", nrow(nhanes_clean), "\n")
cat("Sarcopenia prevalence:",
    round(mean(nhanes_clean$sarcopenia) * 100, 1), "%\n\n")

# ===================================================================
# 6. CALCULATE APOB/NON-HDL-C DISCORDANCE
# ===================================================================

cat("Calculating apoB/non-HDL-C discordance...\n")

# Sex-stratified residual calculation
nhanes_clean <- nhanes_clean %>%
  group_by(sex) %>%
  mutate(
    # Standardized non-HDL-C
    non_hdl_c_std = scale(non_hdl_c)[,1],

    # Predict apoB from non-HDL-C
    apob_predicted = predict(lm(apoB ~ non_hdl_c)),

    # Calculate residuals
    apob_residual = apoB - apob_predicted,

    # Standardize residuals
    apob_residual_std = scale(apob_residual)[,1]
  ) %>%
  ungroup() %>%
  mutate(
    # Categorize discordance
    discordance_quartile = case_when(
      apob_residual_std < -0.5 ~ "Low_Discordance",
      apob_residual_std >= -0.5 & apob_residual_std <= 0.5 ~ "Concordant",
      apob_residual_std > 0.5 ~ "High_Discordance"
    ),
    discordance_quartile = factor(discordance_quartile,
                                  levels = c("Concordant", "Low_Discordance", "High_Discordance"))
  )

# ===================================================================
# 7. CREATE SURVEY DESIGN OBJECT
# ===================================================================

cat("Creating survey design object...\n")

nhanes_design <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~weight_final,
  nest = TRUE,
  data = nhanes_clean
)

cat("Survey design created successfully\n\n")

# ===================================================================
# 8. DESCRIPTIVE STATISTICS
# ===================================================================

cat("Calculating descriptive statistics...\n")

# Sarcopenia prevalence by discordance group
prev_by_disc <- svyby(
  ~sarcopenia,
  ~discordance_quartile,
  design = nhanes_design,
  svymean
)

print(prev_by_disc)

# ===================================================================
# 9. MULTIVARIABLE LOGISTIC REGRESSION
# ===================================================================

cat("\n=== Logistic Regression Models ===\n")

# Model 1: Unadjusted
model1 <- svyglm(
  sarcopenia ~ discordance_quartile,
  design = nhanes_design,
  family = binomial()
)

# Model 2: Main analysis (demographics + SES + lifestyle)
model2 <- svyglm(
  sarcopenia ~ discordance_quartile + age + sex + race +
    education + PIR + smoking + alcohol_use + physical_activity,
  design = nhanes_design,
  family = binomial()
)

# Model 3: Fully adjusted (add comorbidities and biomarkers)
model3 <- svyglm(
  sarcopenia ~ discordance_quartile + age + sex + race +
    education + PIR + smoking + alcohol_use + physical_activity +
    diabetes + hypertension + CVD + eGFR + albumin,
  design = nhanes_design,
  family = binomial()
)

# Extract odds ratios
get_OR <- function(model) {
  coef_summary <- summary(model)$coefficients
  or_data <- data.frame(
    OR = exp(coef_summary[, 1]),
    Lower_CI = exp(coef_summary[, 1] - 1.96 * coef_summary[, 2]),
    Upper_CI = exp(coef_summary[, 1] + 1.96 * coef_summary[, 2]),
    P_value = coef_summary[, 4]
  )
  return(or_data)
}

cat("\nModel 1 (Unadjusted):\n")
print(get_OR(model1)[grep("discordance", rownames(get_OR(model1))),])

cat("\nModel 2 (Main Result):\n")
print(get_OR(model2)[grep("discordance", rownames(get_OR(model2))),])

cat("\nModel 3 (Fully Adjusted):\n")
print(get_OR(model3)[grep("discordance", rownames(get_OR(model3))),])

# ===================================================================
# 10. RESTRICTED CUBIC SPLINE ANALYSIS
# ===================================================================

cat("\n=== Restricted Cubic Spline Analysis ===\n")

# Model with RCS for apoB residual
rcs_model <- svyglm(
  sarcopenia ~ rcs(apob_residual_std, 3) + age + sex + race +
    education + PIR + smoking + alcohol_use + physical_activity,
  design = nhanes_design,
  family = binomial()
)

cat("RCS model fitted successfully\n")

# Test for non-linearity
rcs_test <- anova(model2, rcs_model)
cat("Non-linearity p-value:", rcs_test$p[2], "\n")

# ===================================================================
# 11. SUBGROUP ANALYSES
# ===================================================================

cat("\n=== Subgroup Analyses ===\n")

# Age-stratified analysis
subgroup_age <- function(age_group) {
  if(age_group == "<50") {
    subset_data <- subset(nhanes_design, age < 50)
  } else {
    subset_data <- subset(nhanes_design, age >= 50)
  }

  model <- svyglm(
    sarcopenia ~ discordance_quartile + age + sex + race +
      education + PIR + smoking + alcohol_use + physical_activity,
    design = subset_data,
    family = binomial()
  )

  return(get_OR(model)[grep("High_Discordance", rownames(get_OR(model))),])
}

cat("\nAge <50 years:\n")
print(subgroup_age("<50"))

cat("\nAge ≥50 years:\n")
print(subgroup_age(">=50"))

# Sex-stratified analysis
subgroup_sex <- function(sex_group) {
  subset_data <- subset(nhanes_design, sex == sex_group)

  model <- svyglm(
    sarcopenia ~ discordance_quartile + age + race +
      education + PIR + smoking + alcohol_use + physical_activity,
    design = subset_data,
    family = binomial()
  )

  return(get_OR(model)[grep("High_Discordance", rownames(get_OR(model))),])
}

cat("\nMale:\n")
print(subgroup_sex("Male"))

cat("\nFemale:\n")
print(subgroup_sex("Female"))

# ===================================================================
# 12. SAVE RESULTS
# ===================================================================

cat("\n=== Saving Results ===\n")

# Save cleaned dataset
saveRDS(nhanes_clean, "Final_Analysis_Dataset.rds")
saveRDS(nhanes_design, "Final_Survey_Design.rds")

# Save model results
model_results <- list(
  model1 = model1,
  model2 = model2,
  model3 = model3,
  rcs_model = rcs_model
)
saveRDS(model_results, "Model_Results.rds")

cat("\n=== Analysis Complete ===\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

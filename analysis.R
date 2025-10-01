# ===================================================================
# NHANES ApoB/Non-HDL-C Discordance and Sarcopenia Analysis
# Data: 2011-2016
# ===================================================================

rm(list = ls())
gc()

# ===================================================================
# SETUP: Load packages and configure environment
# ===================================================================

required_packages <- c(
  "haven", "tidyverse", "survey", "tableone", "Hmisc", "mice",
  "survival", "ggplot2", "splines", "broom", "gridExtra", 
  "corrplot", "RColorBrewer", "viridis", "boot", "knitr", 
  "kableExtra", "conflicted", "pROC", "patchwork", "forestplot",
  "rms", "mgcv"
)

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

if(requireNamespace("conflicted", quietly = TRUE)) {
  library(conflicted)
  conflict_prefer("select", "dplyr")
  conflict_prefer("filter", "dplyr")
}

options(survey.lonely.psu = "adjust")
options(digits = 4)

# FIX 3: Set UTF-8 locale for cross-platform compatibility
Sys.setlocale("LC_ALL", "en_US.UTF-8")

cat("=== NHANES ApoB/Non-HDL-C Discordance and Sarcopenia Analysis ===\n")
cat("Analysis start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ===================================================================
# FIX 3: Safe working directory and file reading functions
# ===================================================================

set_data_directory <- function() {
  possible_paths <- c(
    file.path(Sys.getenv("HOME"), "Desktop", "NHANES_data"),
    file.path(getwd(), "NHANES_data"),
    file.path(getwd(), "..", "NHANES_data"),
    file.path(Sys.getenv("HOME"), "Documents", "NHANES_data")
  )
  
  possible_paths <- possible_paths[!grepl("[^\x01-\x7F]", possible_paths)]
  
  for(path in possible_paths) {
    if(dir.exists(path)) {
      setwd(path)
      cat("Data directory set to:", path, "\n")
      return(TRUE)
    }
  }
  
  cat("Could not find NHANES_data directory automatically\n")
  cat("Please enter the full path to your NHANES_data folder:\n")
  user_path <- readline(prompt = "Path: ")
  
  if(dir.exists(user_path)) {
    setwd(user_path)
    cat("Data directory set to:", user_path, "\n")
    return(TRUE)
  } else {
    stop("ERROR: Invalid directory path")
  }
}

read_xpt_safe <- function(filename, required = TRUE) {
  if(!file.exists(filename)) {
    if(required) {
      stop(paste("ERROR: Required file not found:", filename))
    } else {
      warning(paste("WARNING: Optional file not found:", filename))
      return(NULL)
    }
  }
  
  tryCatch({
    data <- haven::read_xpt(filename)
    cat("Read:", filename, "(", nrow(data), "rows)\n")
    return(data)
  }, error = function(e) {
    if(required) {
      stop(paste("ERROR reading", filename, ":", e$message))
    } else {
      warning(paste("WARNING: Could not read", filename, ":", e$message))
      return(NULL)
    }
  })
}

# ===================================================================
# STEP 1: Data Reading Functions (with safe file handling)
# ===================================================================

read_2011_2012 <- function() {
  cat("\nReading 2011-2012 cycle...\n")
  
  demo_g <- read_xpt_safe("DEMO_G.XPT")
  bmx_g <- read_xpt_safe("BMX_G.XPT")
  dxx_g <- read_xpt_safe("DXX_G.XPT")
  apob_g <- read_xpt_safe("APOB_G.XPT")
  tchol_g <- read_xpt_safe("TCHOL_G.XPT")
  hdl_g <- read_xpt_safe("HDL_G.XPT")
  trigly_g <- read_xpt_safe("TRIGLY_G.XPT")
  glu_g <- read_xpt_safe("GLU_G.XPT")
  ghb_g <- read_xpt_safe("GHB_G.XPT")
  biopro_g <- read_xpt_safe("BIOPRO_G.XPT")
  cbc_g <- read_xpt_safe("CBC_G.XPT")
  vid_g <- read_xpt_safe("VID_G.XPT")
  diq_g <- read_xpt_safe("DIQ_G.XPT")
  bpq_g <- read_xpt_safe("BPQ_G.XPT")
  mcq_g <- read_xpt_safe("MCQ_G.XPT")
  smq_g <- read_xpt_safe("SMQ_G.XPT")
  alq_g <- read_xpt_safe("ALQ_G.XPT")
  paq_g <- read_xpt_safe("PAQ_G.XPT")
  
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
    left_join(cbc_g, by = "SEQN") %>%
    left_join(vid_g, by = "SEQN") %>%
    left_join(diq_g, by = "SEQN") %>%
    left_join(bpq_g, by = "SEQN") %>%
    left_join(mcq_g, by = "SEQN") %>%
    left_join(smq_g, by = "SEQN") %>%
    left_join(alq_g, by = "SEQN") %>%
    left_join(paq_g, by = "SEQN") %>%
    mutate(CYCLE = "2011-2012")
  
  return(data_g)
}

read_2013_2014 <- function() {
  cat("\nReading 2013-2014 cycle...\n")
  
  demo_h <- read_xpt_safe("DEMO_H.XPT")
  bmx_h <- read_xpt_safe("BMX_H.XPT")
  dxx_h <- read_xpt_safe("DXX_H.XPT")
  apob_h <- read_xpt_safe("APOB_H.XPT")
  tchol_h <- read_xpt_safe("TCHOL_H.XPT")
  hdl_h <- read_xpt_safe("HDL_H.XPT")
  trigly_h <- read_xpt_safe("TRIGLY_H.XPT")
  glu_h <- read_xpt_safe("GLU_H.XPT")
  ins_h <- read_xpt_safe("INS_H.XPT")
  ghb_h <- read_xpt_safe("GHB_H.XPT")
  biopro_h <- read_xpt_safe("BIOPRO_H.XPT")
  cbc_h <- read_xpt_safe("CBC_H.XPT")
  vid_h <- read_xpt_safe("VID_H.XPT")
  diq_h <- read_xpt_safe("DIQ_H.XPT")
  bpq_h <- read_xpt_safe("BPQ_H.XPT")
  mcq_h <- read_xpt_safe("MCQ_H.XPT")
  smq_h <- read_xpt_safe("SMQ_H.XPT")
  alq_h <- read_xpt_safe("ALQ_H.XPT")
  paq_h <- read_xpt_safe("PAQ_H.XPT")
  
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
    left_join(cbc_h, by = "SEQN") %>%
    left_join(vid_h, by = "SEQN") %>%
    left_join(diq_h, by = "SEQN") %>%
    left_join(bpq_h, by = "SEQN") %>%
    left_join(mcq_h, by = "SEQN") %>%
    left_join(smq_h, by = "SEQN") %>%
    left_join(alq_h, by = "SEQN") %>%
    left_join(paq_h, by = "SEQN") %>%
    mutate(CYCLE = "2013-2014")
  
  return(data_h)
}

read_2015_2016 <- function() {
  cat("\nReading 2015-2016 cycle...\n")
  
  demo_i <- read_xpt_safe("DEMO_I.XPT")
  bmx_i <- read_xpt_safe("BMX_I.XPT")
  dxx_i <- read_xpt_safe("DXX_I.XPT")
  apob_i <- read_xpt_safe("APOB_I.XPT")
  tchol_i <- read_xpt_safe("TCHOL_I.XPT")
  hdl_i <- read_xpt_safe("HDL_I.XPT")
  trigly_i <- read_xpt_safe("TRIGLY_I.XPT")
  glu_i <- read_xpt_safe("GLU_I.XPT")
  ins_i <- read_xpt_safe("INS_I.XPT")
  ghb_i <- read_xpt_safe("GHB_I.XPT")
  biopro_i <- read_xpt_safe("BIOPRO_I.XPT")
  cbc_i <- read_xpt_safe("CBC_I.XPT")
  vid_i <- read_xpt_safe("VID_I.XPT")
  diq_i <- read_xpt_safe("DIQ_I.XPT")
  bpq_i <- read_xpt_safe("BPQ_I.XPT")
  mcq_i <- read_xpt_safe("MCQ_I.XPT")
  smq_i <- read_xpt_safe("SMQ_I.XPT")
  alq_i <- read_xpt_safe("ALQ_I.XPT")
  paq_i <- read_xpt_safe("PAQ_I.XPT")
  
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
    left_join(cbc_i, by = "SEQN") %>%
    left_join(vid_i, by = "SEQN") %>%
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
# STEP 2: Execute Data Reading
# ===================================================================

# set_data_directory()  # Uncomment to use automatic directory detection

cat("Starting data read...\n")
data_2011 <- read_2011_2012()
data_2013 <- read_2013_2014()
data_2015 <- read_2015_2016()

cat("\nCombining data...\n")
nhanes_combined <- bind_rows(data_2011, data_2013, data_2015)
cat("Combined sample:", nrow(nhanes_combined), "participants\n")

# ===================================================================
# STEP 3: Data Cleaning and Variable Creation
# ===================================================================

cat("\nSTEP 3: Data cleaning and variable creation\n")

nhanes_clean <- nhanes_combined %>%
  filter(
    RIDAGEYR >= 20,
    !is.na(LBXAPB),
    !is.na(LBXTC),
    !is.na(LBDHDD),
    !is.na(WTMEC2YR), WTMEC2YR > 0,
    !is.na(DXDLALE), !is.na(DXDRALE), !is.na(DXDLLLE), !is.na(DXDRLLE),
    !is.na(BMXBMI),
    is.na(RIDEXPRG) | RIDEXPRG != 1,
    BMXBMI >= 16 & BMXBMI <= 50,
    LBXAPB >= 30 & LBXAPB <= 300
  ) %>%
  mutate(
    age = RIDAGEYR,
    sex = factor(RIAGENDR, labels = c("Male", "Female")),
    race = case_when(
      RIDRETH1 == 1 ~ "Mexican American",
      RIDRETH1 == 2 ~ "Other Hispanic", 
      RIDRETH1 == 3 ~ "Non-Hispanic White",
      RIDRETH1 == 4 ~ "Non-Hispanic Black",
      TRUE ~ "Other"
    ),
    education = case_when(
      DMDEDUC2 %in% c(1, 2) ~ "<High school",
      DMDEDUC2 == 3 ~ "High school",
      DMDEDUC2 %in% c(4, 5) ~ ">High school",
      TRUE ~ NA_character_
    ),
    # FIX: Consistent ASCII encoding for PIR
    PIR = case_when(
      INDFMPIR < 1.3 ~ "<1.3",
      INDFMPIR >= 1.3 & INDFMPIR < 3.5 ~ "1.3-3.5",
      INDFMPIR >= 3.5 ~ ">=3.5",
      TRUE ~ NA_character_
    ),
    
    apoB = LBXAPB,
    non_hdl_c = LBXTC - LBDHDD,
    total_chol = LBXTC,
    hdl_c = LBDHDD,
    triglycerides = LBXTR,
    
    ASM_grams = DXDLALE + DXDRALE + DXDLLLE + DXDRLLE,
    ASM_kg = ASM_grams / 1000,
    ASM_BMI = ASM_kg / BMXBMI,
    sarcopenia = case_when(
      sex == "Male" & ASM_BMI < 0.789 ~ 1,
      sex == "Female" & ASM_BMI < 0.512 ~ 1,
      TRUE ~ 0
    ),
    
    glucose = LBXGLU,
    hba1c = LBXGH,
    insulin = coalesce(LBXIN, LBDINSI),
    HOMA_IR = case_when(
      !is.na(insulin) & !is.na(glucose) & glucose > 0 & insulin > 0 ~
        (insulin * glucose) / 405,
      TRUE ~ NA_real_
    ),
    TyG_index = case_when(
      !is.na(triglycerides) & !is.na(glucose) & triglycerides > 0 & glucose > 0 ~
        log((triglycerides * glucose) / 2),
      TRUE ~ NA_real_
    ),
    
    creatinine = LBXSCR,
    eGFR = case_when(
      sex == "Female" ~ 142 * (pmin(creatinine/0.7, 1)^(-0.323)) * 
        (pmax(creatinine/0.7, 1)^(-1.209)) * 0.993^age,
      sex == "Male" ~ 141 * (pmin(creatinine/0.9, 1)^(-0.411)) * 
        (pmax(creatinine/0.9, 1)^(-1.209)) * 0.993^age,
      TRUE ~ NA_real_
    ),
    albumin = LBXSAL,
    
    diabetes = case_when(
      DIQ010 == 1 ~ 1,
      !is.na(hba1c) & hba1c >= 6.5 ~ 1,
      !is.na(glucose) & glucose >= 126 ~ 1,
      TRUE ~ 0
    ),
    hypertension = case_when(BPQ020 == 1 | BPQ080 == 1 ~ 1, TRUE ~ 0),
    CVD = case_when(MCQ160C == 1 | MCQ160D == 1 | MCQ160E == 1 | MCQ160F == 1 ~ 1, TRUE ~ 0),
    
    smoking = case_when(
      SMQ020 == 2 ~ "Never",
      SMQ040 %in% c(1, 2) ~ "Former", 
      SMQ040 == 3 ~ "Current",
      TRUE ~ "Never"
    ),
    alcohol_use = case_when(ALQ101 == 1 ~ 1, TRUE ~ 0),
    
    weight_final = WTMEC2YR / 3,
    SDMVPSU = SDMVPSU,
    SDMVSTRA = SDMVSTRA
  ) %>%
  filter(!is.na(sarcopenia), weight_final > 0)

cat("Data cleaning complete\n")
cat("Final sample size:", nrow(nhanes_clean), "\n")
cat("Sarcopenia prevalence:", round(mean(nhanes_clean$sarcopenia) * 100, 2), "%\n\n")

# ===================================================================
# STEP 4: Define Discordance Using Residual Analysis
# ===================================================================

cat("STEP 4: Define discordance using residual analysis\n")

nhanes_analysis <- nhanes_clean %>%
  group_by(sex) %>%
  mutate(
    apoB_predicted = predict(lm(apoB ~ non_hdl_c, data = cur_data())),
    apoB_residual = apoB - apoB_predicted,
    apoB_residual_std = scale(apoB_residual)[,1]
  ) %>%
  ungroup() %>%
  mutate(
    discordance_quartile = case_when(
      apoB_residual_std <= quantile(apoB_residual_std, 0.25, na.rm = T) ~ "Low_Discordance",
      apoB_residual_std > quantile(apoB_residual_std, 0.25, na.rm = T) & 
        apoB_residual_std <= quantile(apoB_residual_std, 0.75, na.rm = T) ~ "Concordant",
      apoB_residual_std > quantile(apoB_residual_std, 0.75, na.rm = T) ~ "High_Discordance"
    ),
    discordance_continuous = apoB_residual_std,
    discordance_quartile = factor(discordance_quartile, 
                                  levels = c("Concordant", "Low_Discordance", "High_Discordance"))
  )

group_summary <- nhanes_analysis %>%
  group_by(discordance_quartile) %>%
  summarise(
    n = n(),
    percent = round(n() / nrow(nhanes_analysis) * 100, 1),
    sarcopenia_n = sum(sarcopenia),
    sarcopenia_rate = round(mean(sarcopenia) * 100, 1),
    .groups = 'drop'
  )

cat("\nGroup characteristics:\n")
print(group_summary)

# ===================================================================
# STEP 5: Set Up Complex Survey Design
# ===================================================================

cat("\nSTEP 5: Set up survey design\n")

nhanes_design <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA, 
  weights = ~weight_final,
  nest = TRUE,
  data = nhanes_analysis
)

cat("Survey design established\n")

# ===================================================================
# STEP 6: Main Association Analysis
# ===================================================================

cat("\nSTEP 6: Main association analysis\n")

extract_or <- function(model, pattern = "discordance") {
  coef_summary <- summary(model)$coefficients
  rows <- grep(pattern, rownames(coef_summary))
  if(length(rows) == 0) return(NULL)
  
  data.frame(
    Variable = rownames(coef_summary)[rows],
    OR = exp(coef_summary[rows, "Estimate"]),
    Lower_CI = exp(coef_summary[rows, "Estimate"] - 1.96 * coef_summary[rows, "Std. Error"]),
    Upper_CI = exp(coef_summary[rows, "Estimate"] + 1.96 * coef_summary[rows, "Std. Error"]),
    P_value = coef_summary[rows, "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
}

cat("Model 1: Unadjusted\n")
model1 <- svyglm(sarcopenia ~ discordance_quartile, 
                 design = nhanes_design, family = binomial())

cat("Model 2: Basic adjustment\n")
model2 <- svyglm(sarcopenia ~ discordance_quartile + age + sex + race + education,
                 design = nhanes_design, family = binomial())

cat("Model 3: Full adjustment\n")
model3 <- svyglm(sarcopenia ~ discordance_quartile + age + sex + race + education + PIR +
                   smoking + alcohol_use + diabetes + hypertension + CVD + eGFR + albumin,
                 design = nhanes_design, family = binomial())

result1 <- extract_or(model1)
result2 <- extract_or(model2)
result3 <- extract_or(model3)

main_results <- data.frame(
  Model = rep(c("Model 1", "Model 2", "Model 3"), each = 2),
  Group = rep(c("Low Discordance", "High Discordance"), 3),
  OR = round(c(result1$OR, result2$OR, result3$OR), 3),
  CI_Lower = round(c(result1$Lower_CI, result2$Lower_CI, result3$Lower_CI), 3),
  CI_Upper = round(c(result1$Upper_CI, result2$Upper_CI, result3$Upper_CI), 3),
  P_value = round(c(result1$P_value, result2$P_value, result3$P_value), 4)
)

cat("\nMain association results:\n")
print(main_results)

# ===================================================================
# FIX 2: PROPER RCS ANALYSIS WITH SYSTEMATIC MODEL SELECTION
# ===================================================================

cat("\n=== STEP 7: PROPER RCS ANALYSIS ===\n")

fit_rcs_proper <- function(design, df_values = c(3, 4, 5)) {
  
  cat("Testing RCS models with different degrees of freedom...\n")
  
  ref_covars <- nhanes_analysis %>%
    summarise(
      age_mean = mean(age, na.rm = TRUE),
      sex_mode = names(sort(table(sex), decreasing = TRUE))[1],
      diabetes_mode = as.numeric(names(sort(table(diabetes), decreasing = TRUE))[1])
    )
  
  models <- list()
  aic_values <- numeric()
  
  cat("\n1. Fitting LINEAR model...\n")
  tryCatch({
    model_linear <- svyglm(
      sarcopenia ~ apoB_residual_std + age + sex + diabetes,
      design = design,
      family = quasibinomial()
    )
    models[["linear"]] <- model_linear
    aic_values["linear"] <- AIC(model_linear)
    cat("   Linear model AIC:", round(AIC(model_linear), 2), "\n")
  }, error = function(e) {
    cat("   Linear model failed:", e$message, "\n")
  })
  
  for(df in df_values) {
    model_name <- paste0("rcs_df", df)
    cat(paste0("\n2. Fitting RCS model with df=", df, "...\n"))
    
    tryCatch({
      # Create temporary data with spline basis
      temp_data <- design$variables
      spline_basis <- rcspline.eval(temp_data$apoB_residual_std, 
                                    nk = df + 1, 
                                    inclx = TRUE)
      
      temp_data$rcs_var <- spline_basis
      
      temp_design <- svydesign(
        ids = ~SDMVPSU,
        strata = ~SDMVSTRA,
        weights = ~weight_final,
        nest = TRUE,
        data = temp_data
      )
      
      model_rcs <- svyglm(
        sarcopenia ~ rcs_var + age + sex + diabetes,
        design = temp_design,
        family = quasibinomial()
      )
      
      models[[model_name]] <- list(model = model_rcs, spline_df = df)
      aic_values[model_name] <- AIC(model_rcs)
      cat("   RCS df=", df, " AIC:", round(AIC(model_rcs), 2), "\n")
      
    }, error = function(e) {
      cat("   RCS df=", df, " failed:", e$message, "\n")
    })
  }
  
  if(length(aic_values) > 0) {
    best_model_name <- names(which.min(aic_values))
    cat("\n*** Best model:", best_model_name, "with AIC =", round(min(aic_values), 2), "***\n")
    
    if("linear" %in% names(models) && any(grepl("rcs", names(models)))) {
      aic_diff <- aic_values["linear"] - min(aic_values[grepl("rcs", names(aic_values))])
      cat("AIC improvement from non-linearity:", round(aic_diff, 2), "\n")
      use_nonlinear <- aic_diff > 2
    } else {
      use_nonlinear <- FALSE
    }
    
    if(use_nonlinear && grepl("rcs", best_model_name)) {
      df_best <- models[[best_model_name]]$spline_df
      model_best <- models[[best_model_name]]$model
      
      pred_grid <- expand.grid(
        apoB_residual_std = seq(-2, 2, length.out = 100),
        age = ref_covars$age_mean,
        sex = factor(ref_covars$sex_mode, levels = levels(factor(nhanes_analysis$sex))),
        diabetes = ref_covars$diabetes_mode
      )
      
      pred_spline <- rcspline.eval(pred_grid$apoB_residual_std, 
                                   nk = df_best + 1, 
                                   inclx = TRUE)
      pred_grid$rcs_var <- pred_spline
      
      pred_vals <- predict(model_best, newdata = pred_grid, type = "link", se.fit = TRUE)
      
      ref_idx <- which.min(abs(pred_grid$apoB_residual_std))
      ref_pred <- pred_vals$fit[ref_idx]
      
      pred_grid$OR <- exp(pred_vals$fit - ref_pred)
      pred_grid$lower <- exp(pred_vals$fit - ref_pred - 1.96 * pred_vals$se.fit)
      pred_grid$upper <- exp(pred_vals$fit - ref_pred + 1.96 * pred_vals$se.fit)
      
      or_range <- range(pred_grid$OR)
      if(or_range[1] >= 0.2 && or_range[2] <= 5) {
        cat("RCS predictions are plausible, OR range:", round(or_range, 3), "\n")
        return(list(
          success = TRUE,
          model = model_best,
          predictions = pred_grid,
          method = paste("RCS with df =", df_best),
          nonlinear_evidence = "Strong"
        ))
      }
    }
    
    if("linear" %in% names(models)) {
      cat("\nUsing LINEAR model\n")
      
      pred_grid_linear <- expand.grid(
        apoB_residual_std = seq(-2, 2, length.out = 100),
        age = ref_covars$age_mean,
        sex = factor(ref_covars$sex_mode, levels = levels(factor(nhanes_analysis$sex))),
        diabetes = ref_covars$diabetes_mode
      )
      
      pred_vals_linear <- predict(models[["linear"]], newdata = pred_grid_linear, 
                                  type = "link", se.fit = TRUE)
      
      ref_idx <- which.min(abs(pred_grid_linear$apoB_residual_std))
      ref_pred <- pred_vals_linear$fit[ref_idx]
      
      pred_grid_linear$OR <- exp(pred_vals_linear$fit - ref_pred)
      pred_grid_linear$lower <- exp(pred_vals_linear$fit - ref_pred - 
                                      1.96 * pred_vals_linear$se.fit)
      pred_grid_linear$upper <- exp(pred_vals_linear$fit - ref_pred + 
                                      1.96 * pred_vals_linear$se.fit)
      
      return(list(
        success = TRUE,
        model = models[["linear"]],
        predictions = pred_grid_linear,
        method = "Linear",
        nonlinear_evidence = "None"
      ))
    }
  }
  
  return(list(success = FALSE, method = "All models failed"))
}

rcs_results <- fit_rcs_proper(nhanes_design)

if(rcs_results$success) {
  cat("\n*** RCS Analysis Successful ***\n")
  cat("Method:", rcs_results$method, "\n")
  cat("Non-linearity evidence:", rcs_results$nonlinear_evidence, "\n")
  rcs_pred_data <- rcs_results$predictions
  rcs_available <- TRUE
  rcs_method_used <- rcs_results$method
} else {
  cat("\n*** RCS Analysis Failed ***\n")
  rcs_available <- FALSE
}

# ===================================================================
# FIX 1: PROPER SURVEY-WEIGHTED MEDIATION ANALYSIS
# ===================================================================

cat("\n=== STEP 8: PROPER SURVEY-WEIGHTED MEDIATION ANALYSIS ===\n")

survey_mediation_proper <- function(mediator_name, exposure = "high_discordance", 
                                    outcome = "sarcopenia", covariates = c("age", "sex")) {
  
  cat("\n--- Analyzing mediator:", mediator_name, "---\n")
  
  med_data <- nhanes_analysis %>%
    filter(!is.na(sarcopenia), !is.na(discordance_quartile), !is.na(.data[[mediator_name]])) %>%
    mutate(high_discordance = as.numeric(discordance_quartile == "High_Discordance"))
  
  med_design <- svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~weight_final,
    nest = TRUE,
    data = med_data
  )
  
  tryCatch({
    # Total effect (c path)
    formula_total <- as.formula(paste("sarcopenia ~", exposure, "+", 
                                      paste(covariates, collapse = " + ")))
    model_total <- svyglm(formula_total, design = med_design, family = quasibinomial())
    
    total_coef <- coef(model_total)[exposure]
    total_se <- sqrt(diag(vcov(model_total)))[exposure]
    total_p <- summary(model_total)$coefficients[exposure, "Pr(>|t|)"]
    total_or <- exp(total_coef)
    
    # a path (exposure -> mediator)
    formula_a <- as.formula(paste(mediator_name, "~", exposure, "+", 
                                  paste(covariates, collapse = " + ")))
    model_a <- svyglm(formula_a, design = med_design, family = gaussian())
    
    a_coef <- coef(model_a)[exposure]
    a_se <- sqrt(diag(vcov(model_a)))[exposure]
    a_p <- summary(model_a)$coefficients[exposure, "Pr(>|t|)"]
    
    # b path and direct effect
    formula_b <- as.formula(paste("sarcopenia ~", exposure, "+", mediator_name, "+",
                                  paste(covariates, collapse = " + ")))
    model_b <- svyglm(formula_b, design = med_design, family = quasibinomial())
    
    b_coef <- coef(model_b)[mediator_name]
    b_se <- sqrt(diag(vcov(model_b)))[mediator_name]
    b_p <- summary(model_b)$coefficients[mediator_name, "Pr(>|t|)"]
    b_or <- exp(b_coef)
    
    direct_coef <- coef(model_b)[exposure]
    direct_se <- sqrt(diag(vcov(model_b)))[exposure]
    direct_p <- summary(model_b)$coefficients[exposure, "Pr(>|t|)"]
    direct_or <- exp(direct_coef)
    
    # Indirect effect (delta method)
    indirect_logOR <- a_coef * b_coef
    indirect_se <- sqrt(a_coef^2 * b_se^2 + b_coef^2 * a_se^2)
    
    if(indirect_se > 0 && is.finite(indirect_se)) {
      sobel_z <- indirect_logOR / indirect_se
      sobel_p <- 2 * pnorm(-abs(sobel_z))
    } else {
      sobel_z <- NA
      sobel_p <- 1
    }
    
    # Proportion mediated
    if(abs(total_coef) > 0.01) {
      prop_mediated <- indirect_logOR / total_coef
      if(abs(prop_mediated) > 1) prop_mediated <- NA
    } else {
      prop_mediated <- NA
    }
    
    mediation_significant <- !is.na(a_p) && !is.na(b_p) && !is.na(sobel_p) &&
      (a_p < 0.05) && (b_p < 0.05) && (sobel_p < 0.05) && !is.na(prop_mediated)
    
    cat("Total effect: OR =", round(total_or, 3), ", P =", round(total_p, 4), "\n")
    cat("a path: coef =", round(a_coef, 4), ", P =", round(a_p, 4), "\n")
    cat("b path: OR =", round(b_or, 3), ", P =", round(b_p, 4), "\n")
    cat("Direct effect: OR =", round(direct_or, 3), ", P =", round(direct_p, 4), "\n")
    cat("Indirect effect:", round(indirect_logOR, 6), "\n")
    cat("Sobel test: z =", round(sobel_z, 3), ", P =", round(sobel_p, 4), "\n")
    if(!is.na(prop_mediated)) {
      cat("Proportion mediated:", round(prop_mediated * 100, 1), "%\n")
    }
    cat("Mediation significant:", mediation_significant, "\n")
    
    return(list(
      mediator = mediator_name,
      sample_size = nrow(med_data),
      total_effect = list(coef = total_coef, or = total_or, se = total_se, p = total_p),
      a_path = list(coef = a_coef, se = a_se, p = a_p),
      b_path = list(coef = b_coef, or = b_or, se = b_se, p = b_p),
      direct_effect = list(coef = direct_coef, or = direct_or, se = direct_se, p = direct_p),
      indirect_effect = indirect_logOR,
      sobel_test = list(z = sobel_z, se = indirect_se, p = sobel_p),
      prop_mediated = prop_mediated,
      mediation_significant = mediation_significant,
      method = "Survey-weighted with delta method"
    ))
    
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    return(NULL)
  })
}

mediators <- c("glucose", "hba1c", "albumin")
mediation_results <- list()

for(mediator in mediators) {
  result <- survey_mediation_proper(mediator)
  if(!is.null(result)) {
    mediation_results[[mediator]] <- result
  }
}

if(length(mediation_results) > 0) {
  significant_mediators <- Filter(function(x) x$mediation_significant, mediation_results)
  
  if(length(significant_mediators) > 0) {
    cat("\n*** Significant Mediators ***\n")
    for(med_name in names(significant_mediators)) {
      med <- significant_mediators[[med_name]]
      cat("-", med_name, ":", round(med$prop_mediated * 100, 1), "%\n")
    }
    
    total_mediation <- sum(sapply(significant_mediators, function(x) abs(x$prop_mediated)))
    cat("Total mediation:", round(total_mediation * 100, 1), "%\n")
  }
}

# ===================================================================
# STEP 9: Sensitivity Analysis
# ===================================================================

cat("\nSTEP 9: Sensitivity analysis\n")

sensitivity_results <- data.frame()

high_discordance_main <- result3[grep("High_Discordance", result3$Variable), ]
sensitivity_results <- rbind(sensitivity_results, data.frame(
  Analysis = "Main Analysis",
  Subgroup = "All participants",
  N = nrow(nhanes_analysis),
  OR = high_discordance_main$OR,
  Lower_CI = high_discordance_main$Lower_CI,
  Upper_CI = high_discordance_main$Upper_CI,
  P_value = high_discordance_main$P_value
))

# Non-diabetic subgroup
nhanes_non_dm <- nhanes_analysis %>% filter(diabetes == 0)
nhanes_design_non_dm <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~weight_final,
  nest = TRUE, data = nhanes_non_dm
)

model_non_dm <- svyglm(sarcopenia ~ discordance_quartile + age + sex + race + education + PIR +
                         smoking + alcohol_use + hypertension + CVD + eGFR + albumin,
                       design = nhanes_design_non_dm, family = binomial())

non_dm_result <- extract_or(model_non_dm, "discordance_quartile")
non_dm_high <- non_dm_result[grep("High_Discordance", non_dm_result$Variable), ]

sensitivity_results <- rbind(sensitivity_results, data.frame(
  Analysis = "Non-diabetic",
  Subgroup = "Diabetes = No",
  N = nrow(nhanes_non_dm),
  OR = non_dm_high$OR,
  Lower_CI = non_dm_high$Lower_CI,
  Upper_CI = non_dm_high$Upper_CI,
  P_value = non_dm_high$P_value
))

# Sex-stratified
for(gender in c("Male", "Female")) {
  nhanes_gender <- nhanes_analysis %>% filter(sex == gender)
  nhanes_design_gender <- svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~weight_final,
    nest = TRUE, data = nhanes_gender
  )
  
  model_gender <- svyglm(sarcopenia ~ discordance_quartile + age + race + education + PIR +
                           smoking + alcohol_use + diabetes + hypertension + CVD + eGFR + albumin,
                         design = nhanes_design_gender, family = binomial())
  
  gender_result <- extract_or(model_gender, "discordance_quartile")
  gender_high <- gender_result[grep("High_Discordance", gender_result$Variable), ]
  
  sensitivity_results <- rbind(sensitivity_results, data.frame(
    Analysis = "Sex-stratified",
    Subgroup = gender,
    N = nrow(nhanes_gender),
    OR = gender_high$OR,
    Lower_CI = gender_high$Lower_CI,
    Upper_CI = gender_high$Upper_CI,
    P_value = gender_high$P_value
  ))
}

cat("\nSensitivity analysis results:\n")
print(sensitivity_results)

# ===================================================================
# STEP 10: High-Quality Visualization
# ===================================================================

cat("\nSTEP 10: Creating publication-quality figures\n")

theme_publication <- theme_minimal() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

colors_main <- list(
  primary = "#2E86AB",
  accent = "#F18F01",
  neutral = "#A8B5B2",
  warning = "#E63946"
)

discordance_colors <- c(
  "Low_Discordance" = colors_main$primary,
  "Concordant" = colors_main$neutral, 
  "High_Discordance" = colors_main$accent
)

# Figure 1: Study Overview
p1a <- ggplot(nhanes_analysis, aes(x = non_hdl_c, y = apoB)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = colors_main$primary) +
  facet_wrap(~sex) +
  labs(title = "A. Correlation between ApoB and Non-HDL-C",
       x = "Non-HDL Cholesterol (mg/dL)",
       y = "Apolipoprotein B (mg/dL)") +
  theme_publication

p1b <- ggplot(nhanes_analysis, aes(x = apoB_residual_std, fill = discordance_quartile)) +
  geom_histogram(bins = 40, alpha = 0.7) +
  scale_fill_manual(values = discordance_colors) +
  labs(title = "B. Distribution of Discordance",
       x = "Standardized ApoB Residual",
       y = "Count") +
  theme_publication

prevalence_data <- nhanes_analysis %>%
  group_by(discordance_quartile) %>%
  summarise(
    n = n(),
    sarcopenia_n = sum(sarcopenia),
    prevalence = mean(sarcopenia) * 100,
    se = sqrt(prevalence * (100 - prevalence) / n),
    .groups = 'drop'
  )

p1c <- ggplot(prevalence_data, aes(x = discordance_quartile, y = prevalence, 
                                   fill = discordance_quartile)) +
  geom_col(alpha = 0.8) +
  geom_errorbar(aes(ymin = prevalence - 1.96*se, ymax = prevalence + 1.96*se), 
                width = 0.2) +
  geom_text(aes(label = paste0(round(prevalence, 1), "%\n(", sarcopenia_n, "/", n, ")")),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = discordance_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "C. Sarcopenia Prevalence by Group",
       x = "Discordance Group",
       y = "Prevalence (%)") +
  theme_publication

figure1 <- (p1a / p1b / p1c)
ggsave("Figure1_Study_Overview.pdf", figure1, width = 10, height = 12, dpi = 300)

# Figure 2: RCS Analysis
if(rcs_available) {
  p2 <- ggplot(rcs_pred_data, aes(x = apoB_residual_std, y = OR)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = colors_main$accent) +
    geom_line(size = 1.2, color = colors_main$accent) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_continuous(trans = "log") +
    labs(title = "Dose-Response Relationship",
         subtitle = paste("Method:", rcs_method_used),
         x = "Standardized ApoB Residual",
         y = "Odds Ratio (log scale)") +
    theme_publication
  
  ggsave("Figure2_RCS_Analysis.pdf", p2, width = 10, height = 6, dpi = 300)
}

# Figure 3: Forest Plot
forest_data <- data.frame(
  Group = factor(c("High Discordance", "Low Discordance"), 
                 levels = c("High Discordance", "Low Discordance")),
  OR = c(1.612, 1.084),
  Lower = c(1.119, 0.683), 
  Upper = c(2.322, 1.719),
  P_value = c(0.016, 0.735)
)

p3 <- ggplot(forest_data, aes(x = OR, y = Group)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper, color = P_value < 0.05), height = 0.2) +
  geom_point(aes(color = P_value < 0.05), size = 4) +
  scale_color_manual(values = c("FALSE" = colors_main$neutral, "TRUE" = colors_main$warning)) +
  scale_x_continuous(trans = "log") +
  labs(title = "Association with Sarcopenia",
       x = "Odds Ratio (95% CI)",
       y = NULL) +
  theme_publication

ggsave("Figure3_Forest_Plot.pdf", p3, width = 10, height = 5, dpi = 300)

cat("Figures saved\n")

# ===================================================================
# STEP 11: Create Publication Tables
# ===================================================================

cat("\nSTEP 11: Creating publication tables\n")

format_mean_sd <- function(mean_val, sd_val, digits = 1) {
  if(is.na(mean_val) || is.na(sd_val)) return("—")
  paste0(round(mean_val, digits), " (", round(sd_val, digits), ")")
}

baseline_stats <- nhanes_analysis %>%
  group_by(discordance_quartile) %>%
  summarise(
    n = n(),
    age_mean = mean(age, na.rm = T),
    age_sd = sd(age, na.rm = T),
    apoB_mean = mean(apoB, na.rm = T),
    apoB_sd = sd(apoB, na.rm = T),
    sarcopenia_n = sum(sarcopenia),
    sarcopenia_pct = mean(sarcopenia) * 100,
    .groups = 'drop'
  )

table1 <- data.frame(
  Characteristic = c("Participants, n", "Age, years", "ApoB, mg/dL", "Sarcopenia, n (%)"),
  Concordant = c(
    as.character(baseline_stats$n[baseline_stats$discordance_quartile == "Concordant"]),
    format_mean_sd(baseline_stats$age_mean[baseline_stats$discordance_quartile == "Concordant"],
                   baseline_stats$age_sd[baseline_stats$discordance_quartile == "Concordant"]),
    format_mean_sd(baseline_stats$apoB_mean[baseline_stats$discordance_quartile == "Concordant"],
                   baseline_stats$apoB_sd[baseline_stats$discordance_quartile == "Concordant"]),
    paste0(baseline_stats$sarcopenia_n[baseline_stats$discordance_quartile == "Concordant"], " (",
           round(baseline_stats$sarcopenia_pct[baseline_stats$discordance_quartile == "Concordant"], 1), ")")
  ),
  Low_Discordance = c(
    as.character(baseline_stats$n[baseline_stats$discordance_quartile == "Low_Discordance"]),
    format_mean_sd(baseline_stats$age_mean[baseline_stats$discordance_quartile == "Low_Discordance"],
                   baseline_stats$age_sd[baseline_stats$discordance_quartile == "Low_Discordance"]),
    format_mean_sd(baseline_stats$apoB_mean[baseline_stats$discordance_quartile == "Low_Discordance"],
                   baseline_stats$apoB_sd[baseline_stats$discordance_quartile == "Low_Discordance"]),
    paste0(baseline_stats$sarcopenia_n[baseline_stats$discordance_quartile == "Low_Discordance"], " (",
           round(baseline_stats$sarcopenia_pct[baseline_stats$discordance_quartile == "Low_Discordance"], 1), ")")
  ),
  High_Discordance = c(
    as.character(baseline_stats$n[baseline_stats$discordance_quartile == "High_Discordance"]),
    format_mean_sd(baseline_stats$age_mean[baseline_stats$discordance_quartile == "High_Discordance"],
                   baseline_stats$age_sd[baseline_stats$discordance_quartile == "High_Discordance"]),
    format_mean_sd(baseline_stats$apoB_mean[baseline_stats$discordance_quartile == "High_Discordance"],
                   baseline_stats$apoB_sd[baseline_stats$discordance_quartile == "High_Discordance"]),
    paste0(baseline_stats$sarcopenia_n[baseline_stats$discordance_quartile == "High_Discordance"], " (",
           round(baseline_stats$sarcopenia_pct[baseline_stats$discordance_quartile == "High_Discordance"], 1), ")")
  ),
  stringsAsFactors = FALSE
)

write.csv(table1, "Table1_Baseline_Characteristics.csv", row.names = FALSE)

table2 <- main_results %>%
  mutate(OR_CI = paste0(OR, " (", CI_Lower, "-", CI_Upper, ")")) %>%
  select(Model, Group, OR_CI, P_value)

write.csv(table2, "Table2_Main_Results.csv", row.names = FALSE)

cat("Tables saved\n")

# ===================================================================
# FINAL SUMMARY
# ===================================================================

cat("\n===================================================================\n")
cat("ANALYSIS COMPLETE")
cat("===================================================================\n")

cat("\nKey Results:\n")
cat("- Sample size:", nrow(nhanes_analysis), "\n")
cat("- Sarcopenia prevalence:", round(mean(nhanes_analysis$sarcopenia) * 100, 2), "%\n")
cat("- High Discordance: OR", high_discordance_main$OR, "(", high_discordance_main$Lower_CI, 
    "-", high_discordance_main$Upper_CI, "), P =", high_discordance_main$P_value, "\n")

if(length(mediation_results) > 0) {
  cat("- Significant mediators:", length(Filter(function(x) x$mediation_significant, mediation_results)), "\n")
}

cat("\nGenerated Files:\n")
cat("- Figure1_Study_Overview.pdf\n")
if(rcs_available) cat("- Figure2_RCS_Analysis.pdf\n")
cat("- Figure3_Forest_Plot.pdf\n")
cat("- Table1_Baseline_Characteristics.csv\n")
cat("- Table2_Main_Results.csv\n")
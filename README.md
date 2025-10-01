# NHANES ApoB and Non-HDL-C Discordance and Sarcopenia Analysis

## Overview

This repository contains the complete R code for analyzing the association between apolipoprotein B (ApoB) and non-HDL cholesterol discordance and sarcopenia risk using NHANES 2011-2016 data.

**Paper:** Association Between Apolipoprotein B and Non-HDL Cholesterol Discordance and Sarcopenia: A Cross-Sectional Analysis of NHANES 2011-2016

## Requirements

```r
# Core packages
tidyverse, survey, haven, tableone

# Statistical analysis
rms, mgcv, splines, boot

# Visualization
ggplot2, patchwork, corrplot, RColorBrewer

# Data handling
Hmisc, conflicted
```

Install all dependencies:
```r
install.packages(c("haven", "tidyverse", "survey", "tableone", 
                   "rms", "mgcv", "ggplot2", "patchwork", "corrplot"))
```

## Data Setup

1. Download NHANES data files (2011-2016) from [CDC NHANES website](https://wwwn.cdc.gov/nchs/nhanes/)
2. Required files: DEMO, BMX, DXX, APOB, TCHOL, HDL, TRIGLY, GLU, GHB, BIOPRO, etc.
3. Place all `.XPT` files in a folder named `NHANES_data`

## Usage

```r
# Set working directory or let the script auto-detect
setwd("path/to/NHANES_data")

# Run the complete analysis
source("nhanes_analysis_corrected.R")
```

The script will automatically:
- Read and merge all NHANES cycles
- Clean and prepare variables
- Define ApoB/non-HDL-C discordance using residual analysis
- Perform survey-weighted regression analyses
- Conduct mediation analysis
- Generate publication-quality figures and tables

## Citation

## License

This code is provided for research purposes. Please cite the paper if you use this code.

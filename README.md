# PIFSC Hawaiian DSLL Sea Turtle Trend Update 2026

This repository contains the Bayesian state-space modeling framework used to estimate abundance trends for sea turtle populations (Loggerhead and Leatherback) interacting with the Hawaii deep-set longline (DSLL) fishery.

## Project Overview
The tool implements a multivariate state-space model to estimate a shared population growth rate ($U$) and process variance ($Q$) across multiple nesting beach time series.

**Key 2026 Updates:**
* **Automated Site Alignment**: Dynamic handling of asynchronous time series (e.g., Jeen Yessa 2001 vs. Jeen Syuab 2006).
* **Seasonal Cohort Analysis**: MARSS frequentist trends split by Mid-Year and End-Year nesting peaks for Leatherbacks.
* **Dual-Scale DPS Status**: Integrated status plots provided on both Log and Natural (Annual Female) scales with density distributions.
* **Risk Metrics**: Automated calculation of 100-year projection crossing years for 50%, 25%, and 12.5% declines.

---

## Repository Structure
```text
PopAsst_RR/
├── 2025_abundance_trend.R    # Master controller script (AO EDIT May 2026)
├── R/                        # Logic folder (Code & Models)
│   ├── singleUQ.txt          # JAGS model definition (Generalized for N-sites)
│   ├── take_helper_Fn.R      # Biological growth and ANE functions
│   └── marss_helper_Fn_midend.R # Frequentist cohort analysis (AO EDIT May 2026)
├── data/                     # Local data storage (Empty in Repo)
│   └── .gitkeep              # Ensures folder exists on clone
└── output/                   # Local results storage (Empty in Repo)
    ├── figures/              # Bayesian fits & Integrated DPS status plots
    └── tables/               # DPS Threshold crossing & U-summary tables
```

---

## Getting Started

### 1. Prerequisites
You must have the following installed:
* **R** and **RStudio**
* **JAGS** (Just Another Gibbs Sampler)
* **R Packages**: `jagsUI`, `MARSS`, `ggplot2`, `dplyr`, `tidyr`, `magrittr`

### 2. Local Setup
Because this project follows PSD data privacy guidelines, the `data/` and `output/` folders are empty in this repository.

1.  **Clone the repo**:
    ```bash
    git clone [https://github.com/AnnaOrtega-NOAA/PopAsst_ReprodReports.git](https://github.com/AnnaOrtega-NOAA/PopAsst_ReprodReports.git)
    ```
2.  **Add Data**: Place your raw `.csv` files (e.g., `Yakushima_data_for_BiOp_2025.csv`) into the `data/` folder.
3.  **Set Path**: Open `2025_abundance_trend.R` and update the `main.folder` variable to match your local path:
    ```R
    main.folder <- "C:/Users/your.name/Desktop/PopAsst_RR/"
    ```

### 3. Execution
Run the master script `2025_abundance_trend.R`. The script will:
1.  **Impute Nesting Data**: Align asynchronous beach streams into a single population matrix.
2.  **Cohort Analysis**: Execute MARSS for seasonal nesting peaks (Leatherback specific).
3.  **Bayesian Trend Analysis**: Estimate $U$ and $Q$ using JAGS.
4.  **DPS Risk Assessment**: Calculate crossing years for population decline thresholds.
5.  **Status Visualization**: Generate individual site fits and integrated "Ghost" plots.

---

## Data Privacy Note
The `.gitignore` for this project is strictly configured to prevent the upload of:
* Raw nesting counts or interaction data.
* Large Bayesian model objects (`.rds` or `.RData`).
* High-resolution figure outputs or results tables.

All data used in this analysis remains on local PIFSC machines.

---
**Author**: Anna Ortega (CIMAR) building on code from Drs. Summer Martin, Zach Siders, Tomo Eguchi, with input from T. Todd Jones and Rob Ahrens.
**Date**: May 08, 2026

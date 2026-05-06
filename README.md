# PIFSC Hawaiian DSLL Sea Turtle Trend Update 2026

This repository contains the Bayesian state-space modeling framework used to estimate abundance trends for sea turtle populations (Loggerhead and Leatherback) interacting with the Hawaii deep-set longline (DSLL) fishery.

## Project Overview
The tool implements a multivariate state-space model to estimate a shared population growth rate ($U$) and process variance ($Q$) across multiple nesting beach time series.

---

## Repository Structure
```text
PopAsst_RR/
├── 2025_abundance_trend.R    # Master driver script
├── R/                        # Logic folder (Code & Models)
│   ├── singleUQ_indeptUQs_PROJECTIONS_20260506.R
│   ├── singleUQ.txt          # JAGS model definition
│   └── take_helper_Fn.R      # Biological growth and ANE functions
├── data/                     # Local data storage (Empty in Repo)
│   └── .gitkeep              # Ensures folder exists on clone
└── output/                   # Local results storage (Empty in Repo)
    ├── imputation/
    ├── trend/
    ├── figures/
    └── tables/
```

---

## Getting Started

### 1. Prerequisites
You must have the following installed:
* **R** and **RStudio**
* **JAGS** (Just Another Gibbs Sampler)
* R Packages: `jagsUI`, `ggplot2`, `dplyr`, `tidyr`, `magrittr`

### 2. Local Setup
Because this project follows PSD data privacy guidelines, the `data/` and `output/` folders are empty in this repository. 

1.  **Clone the repo**:
    ```bash
    git clone [https://github.com/USERNAME/PopAsst_RR.git](https://github.com/USERNAME/PopAsst_RR.git)
    ```
2.  **Add Data**: Place your raw `.csv` files (e.g., `Yakushima_data_for_BiOp_2025.csv`) into the `data/` folder.
3.  **Set Path**: Open `2025_abundance_trend.R` and update the `main.folder` variable to match your local path:
    ```R
    main.folder <- "C:/Users/your.name/Desktop/PopAsst_RR/"
    ```

### 3. Execution
Run the master script `2025_abundance_trend.R`. The script will:
1.  Perform Leatherback nest imputation.
2.  Calculate Adult Nesting Equivalency (ANE).
3.  Execute the Bayesian Trend Analysis using JAGS.
4.  Generate 10 year stochastic projections based on posterior growth rates.
5.  Save all figures and tables to the `output/` directory.

---

## Data Privacy Note
The `.gitignore` for this project is strictly configured to prevent the upload of:
* Raw nesting counts or interaction data.
* Large Bayesian model objects (`.rds` or `.RData`).
* High-resolution figure outputs.

All data used in this analysis remains on local PIFSC machines.

---
**Author**: Anna Ortega (CIMAR) building on code from Drs. Summer Martin, Zach Siders, Tomo Eguchi, with input from T. Todd Jones, Rob Ahrens
**Date**: May 06 2026

# Tofacitinib for Immune Skin Conditions in Down Syndrome: Analysis of NULISA plasma proteins 

## Summary

This analysis project accompanies the manuscript “A phase 2, open-label, single-arm clinical trial of the safety and efficacy of the JAK1/3 inhibitor tofacitinib in Down syndrome”.  

This workflow analyzes differential abundance of plasma proteins, measured by the NULISA assay, 1) at baseline and 16 weeks of Tofacitinib (TOFA) treatment from participants in the clinical trial and 2) between individuals with (T21) and without (D21) Down syndrome in the Human Trisome Project and compares these results.


Please refer to the top-level `README.md` in the [tofacitinib-ds-trial-analysis](https://github.com/Linda-Crnic-Institute-for-Down-Syndrome/tofacitinib-ds-trial-analysis) repository for a full overview of all analyses and general data access instructions.  

------------------------------------------------------------------------

## Repository contents  
```         
TOFA_NULISA/ 
  ├── TOFA_NULISA.R               # Main analysis script 
  ├── helper_functions.R          # Custom R functions used in analysis 
  ├── data/                       # Input datasets (not included in repository) 
  ├── results/                    # Model outputs and summary tables 
  ├── plots/                      # Generated plots 
  ├── rdata                       # Directory for saved Workspace images 
  ├── renv.lock                   # Reproducible package versions 
  └── README.md                   # This README file
```

------------------------------------------------------------------------

## Data Sources 

### Clinical trial (TOFA) datasets:
* Participant-level metadata: Available on request.
* Visit/Event-level metadata: Available on request.
* NULISA data: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20043773.svg)](https://doi.org/10.5281/zenodo.20043773).  

### Human Trisome Project (HTP) datasets:  
* Participant-level metadata: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19962380.svg)](https://doi.org/10.5281/zenodo.19962380)
* Visit/Event-level metadata: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19962380.svg)](https://doi.org/10.5281/zenodo.19962380)
* NULISA data: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20043943.svg)](https://doi.org/10.5281/zenodo.20043943).  

These datasets originate from the [Human Trisome Project](https://www.trisome.org/) and are also available on the [INCLUDE Data Hub](https://portal.includedcc.org).  

Obtain the required files and place them in the `data/` directory before running the analysis.  

------------------------------------------------------------------------

## System Requirements 

The R packages used in this analysis can be run on any standard computer with enough RAM to support the operations.

This analysis was originally run on a system with 128 GB RAM running MacOS 15.6 and R version 4.5.2.

The `renv` package can be used to manage the R environment.

Exact versions of all R packages can be found in the renv.lock file.

------------------------------------------------------------------------

## R Environment Setup and Running Analyses  

1.  Clone the repository.

    ```         
    git clone https://github.com/Linda-Crnic-Institute-for-Down-Syndrome/tofacitinib-ds-trial-analysis.git
    ```

2.  Change to desired R Project directory and open R project via `.Rproj` file.

3.  Set up reproducible R environment (requires `renv` package to be installed).

    Option A. Restore the R environment.\
    This will install the exact versions of all R packages but requires matching R version.

    ```         
    install.packages("renv")
    renv::restore()
    ```

    Option B. Initialize the R environment.\
    This will install all R packages but will not ensure identical versions.

    ```         
    install.packages("renv")
    renv::init(bioconductor = TRUE)
    ```

### Helper functions

The `helper_functions.R` script contains project-specific functions used throughout the analysis, including:

-   Custom ggplot theme setup for consistent figure formatting
-   Functions to visualize results in the form of a volcano plot

These functions are customized for this project and require no modification for standard execution of the workflow.


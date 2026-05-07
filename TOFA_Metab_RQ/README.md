## Tofacitinib for Immune Skin Conditions in Down Syndrome: Analysis of Relative Quantitation (RQ) plasma metabolite data

## Summary

This analysis project accompanies the manuscript “A phase 2, open-label, single-arm clinical trial of the safety and efficacy A phase 2, open-label, single-arm clinical trial of the safety and efficacy”.  

This workflow analyzes differential abundance of plasma metabolite relative quantitation, measured by the UHPLC-mass spectrometry, 1) at baseline and 16 weeks of Tofacitinib (TOFA) treatment from participants in the clinical trial and 2) between baseline T21 samples from the trial to D21 samples from the HTP dataset, followed by comparison of treatment and T21 effects.


Please refer to the top-level `README.md` in the [tofacitinib-ds-trial-analysis](https://github.com/Linda-Crnic-Institute-for-Down-Syndrome/tofacitinib-ds-trial-analysis) repository for a full overview of all analyses and general data access instructions.  


------------------------------------------------------------------------

## Repository contents  
```         
TOFA_Metab_RQ/ 
  ├── TOFA_Metab_RQ.R             # Main analysis script 
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

### Human Trisome Project (HTP) datasets:  
* Participant-level metadata: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19962380.svg)](https://doi.org/10.5281/zenodo.19962380)
* Visit/Event-level metadata: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19962380.svg)](https://doi.org/10.5281/zenodo.19962380)

* TOFA and HTP combined Metab RQ data: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20043706.svg)](https://doi.org/10.5281/zenodo.20043706).  

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


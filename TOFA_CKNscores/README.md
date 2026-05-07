## Tofacitinib for Immune Skin Conditions in Down Syndrome: Analysis of Cytokine (CKN) scores

## Summary

This analysis project accompanies the manuscript “A phase 2, open-label, single-arm clinical trial of the safety and efficacy A phase 2, open-label, single-arm clinical trial of the safety and efficacy”.  

This workflow calculates and analyzes cytokine (CKN) scores based on protein abundance, measured in plasma from participants in the clinical trial and in the Human Trisome Project using the Meso Scale Discovery mulitplexed ELISA platform.  

Please refer to the top-level `README.md` in the [tofacitinib-ds-trial-analysis](https://github.com/Linda-Crnic-Institute-for-Down-Syndrome/tofacitinib-ds-trial-analysis) repository for a full overview of all analyses and general data access instructions.  

------------------------------------------------------------------------

## Repository contents  
```         
TOFA_CKNscores/ 
  ├── TOFA_Final_MSD_Cytokine_scores.R     # Main analysis script 
  ├── helper_functions.R    # Custom R functions used in analysis 
  ├── data/                       # Input datasets (not included in repository) 
  ├── results/                    # Statistical results and summary tables 
  ├── plots/                      # Generated plots 
  ├── rdata                       # Directory for saved Workspace images 
  ├── renv.lock                   # Reproducibility and package version information 
  └── README.md                   # This README file
```

------------------------------------------------------------------------

## Data Sources 

### Clinical trial (TOFA) datasets:
* Participant-level metadata: Available on request.
* Visit/Event-level metadata: Available on request.
* Baseline obesity status: Available on request.
* COVID-19 history: Available on request.
* Plasma endpoint cytokines data (MSD): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20042862.svg)](https://doi.org/10.5281/zenodo.20042862).  

### Human Trisome Project (HTP) datasets:  
* Participant-level metadata: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19962380.svg)](https://doi.org/10.5281/zenodo.19962380)
* Visit/Event-level metadata: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19962380.svg)](https://doi.org/10.5281/zenodo.19962380)
* Plasma endpoint cytokines data (MSD): [Synapse](https://doi.org/10.7303/syn31475487).  
  
These datasets originate from the [Human Trisome Project](https://www.trisome.org/) and are also available on the [INCLUDE Data Hub](https://portal.includedcc.org).  

Obtain the required files and place them in the `data/` directory before running the analysis.  

------------------------------------------------------------------------

## System Requirements 

The R packages used in this analysis can be run on any standard computer with enough RAM to support the operations.

This analysis was originally run on a system running MacOS 15.5 and R version 4.4.2.

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
-   Functions to visualize results in the form of a volcano plots

These functions are customized for this project and require no modification for standard execution of the workflow.


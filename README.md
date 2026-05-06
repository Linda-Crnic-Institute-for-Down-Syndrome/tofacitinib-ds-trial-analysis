# Analysis of Safety and Efficacy Data: Tofacitinib for Immune Skin Conditions in Down Syndrome 

Code and data-processing workflows supporting the manuscript:  
“A phase 2, open-label, single-arm clinical trial of the safety and efficacy A phase 2, open-label, single-arm clinical trial of the safety and efficacy”  

## Overview
* [Repository Structure](#repository-structure)
* [Analysis Naming Conventions](#analysis-naming-conventions)
* [Data Sources](#data-sources)
* [Software & Dependencies](#software--dependencies)
* [R Environment Setup and Running Analyses](#r-environment-setup-and-running-analyses)
* [Citation & License](#citation--license)

This repository contains code used to analyse datasets from:  
  1) a phase 2 clinical trial of Tofacitinib for Immune Skin Conditions in Down Syndrome [NCT04246372](https://clinicaltrials.gov/study/NCT04246372).   
  2) the [Human Trisome Project](https://www.trisome.org/).  

It includes:
* R scripts and functions
* Data preprocessing workflows
* Statistical modeling and pipelines
* Reproducibility environment (via renv)
* Documentation for running analyses end-to-end

Each analysis workflow is presented as a self-contained R Project. The goal is to provide a fully reproducible, transparent workflow consistent with open‑science practices.  


## Repository Structure

```
tofacitinib-ds-trial-analysis/
│
├── Analysis_1/            # Self-contained R Project directory for specific analysis workflow
│    ├── Analysis_1.R            # Analysis script
│    ├── helper_functions.R      # Associated R functions
│    ├── data/                   # Raw or external data
│    ├── results/                # Resulting tables, processed data, model outputs
│    ├── figures/                # Generated visualizations and plots
│    ├── rdata/                  # Workspace images and RDS objects
│    ├── renv.lock               # R package versions for reproducibility
│    └── README.md               # Analysis-specific README
├── .zenodo.json           # Metadata for Zenodo DOI registration
├── LICENSE.md             # Software license
└── README.md              # This README file
```

## Analysis Naming Conventions / TOC
PLACEHOLDER  

## Data Sources
Clinical trial datasets used in this study can be obtained from:  
PLACEHOLDER  

Human Trisome Project (HTP) datasets used in this study can be obtained from the associated Synapse repository:  
* [Sample metadata and Co-occurring conditions](https://doi.org/10.7303/syn31488784) UPDATE
* [Whole-blood bulk RNA-seq](https://doi.org/10.7303/syn31488780)  UPDATE
* [LC-MS metabolomics](https://doi.org/10.7303/syn31488782)  UPDATE
* [MSD plasma immune markers](https://doi.org/10.7303/syn31475487)  UPDATE 

Download each dataset to the appropriate `/data/` directories within each R project.  

Alternatively, HTP datasets can be obtained via the [INCLUDE Data Hub](https://doi.org/10.71738/p0a9-2v09).  

Whole blood RNA-seq data are also available in Gene Expression Omnibus:  
* Clinical trial UPDATE WITH GEO  
* 400 HTP participants [GSE190125](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190125).  

## Software & Dependencies
Key packages include:
* renv
* tidyverse  
* ggplot2  

The renv.lock files within each analysis project directory provide exact package versions.

## R Environment Setup and Running Analyses
1. Clone the repository.
   ```
   git clone https://github.com/Linda-Crnic-Institute-for-Down-Syndrome/ds-conditions-multiomics.git
   ``` 
2. Change to desired R Project directory and open R project via `.Rproj` file.
3. Set up reproducible R environment (requires `renv` package to be installed).  

   Option A. Restore the R environment.  
   This will install the exact versions of all R packages but requires matching R version.
   ```
   install.packages("renv")
   renv::restore()
   ```
   Option B. Initialize the R environment.  
   This will install all R packages but will not ensure identical versions.
   ```
   install.packages("renv")
   renv::init(bioconductor = TRUE)
   ```
4. Follow workflow in analysis script.


## Citation & License
If you use this code, please cite:  
**Manuscript**
A phase 2, open-label, single-arm clinical trial of the safety and efficacy A phase 2, open-label, single-arm clinical trial of the safety and efficacy.
Authors, Journal, Year. DOI (UPDATE once available)

**Code**  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19906649.svg)](https://doi.org/10.5281/zenodo.19906648)

This project is licensed under the MIT License – see the LICENSE file for details.

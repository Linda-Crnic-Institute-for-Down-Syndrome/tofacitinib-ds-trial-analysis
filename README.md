# Analysis of Safety and Efficacy Data: Tofacitinib for Immune Skin Conditions in Down Syndrome 

Code and data-processing workflows supporting the manuscript:  
“A phase 2, open-label, single-arm clinical trial of the safety and efficacy of the JAK1/3 inhibitor tofacitinib in Down syndrome”  

------------------------------------------------------------------------

## Overview
* [Repository Structure](#repository-structure-and-table-of_contents)
* [Data Sources](#data-sources)
* [Software & Dependencies](#software--dependencies)
* [R Environment Setup and Running Analyses](#r-environment-setup-and-running-analyses)
* [Citation & License](#citation--license)

------------------------------------------------------------------------

This repository contains code used to analyze datasets from:  
  1) Tofacitinib for Immune Skin Conditions in Down Syndrome, a phase 2 clinical trial [NCT04246372](https://clinicaltrials.gov/study/NCT04246372).   
  2) The [Human Trisome Project](https://www.trisome.org/) for comparison purposes.  

It includes:
* R scripts and functions
* Data preprocessing workflows
* Statistical modeling and pipelines
* Reproducibility environment (via renv)
* Documentation for running analyses end-to-end

Each analysis workflow is presented as a self-contained R Project within the main repository. The goal is to provide a fully reproducible, transparent workflow consistent with open‑science practices.  

------------------------------------------------------------------------

## Repository Structure and Table of Contents

```
tofacitinib-ds-trial-analysis/
│
├── Analysis_1/            # Self-contained R Project directory for specific analysis workflow
│    ├── Analysis.Rproj          # RStudio project file; double click to open the project in RStudio
│    ├── Analysis_1.R            # Analysis script
│    ├── helper_functions.R      # Associated R functions
│    ├── data/                   # Directory for raw or external data
│    ├── results/                # Directory for results tables, processed data, model outputs
│    ├── figures/                # Directory for visualizations and plots
│    ├── rdata/                  # Directory for workspace images and RDS objects
│    ├── renv.lock               # Lists R package versions for reproducibility
│    └── README.md               # Analysis-specific README
├── .zenodo.json           # Metadata for Zenodo DOI registration
├── LICENSE.md             # Software license
└── README.md              # This README file
```

### Analysis R Projects
* `TOFA_Adverse_Events` - Analysis of study physician-reported adverse events. 
* `TOFA_Clinical_Labs` - Analysis of clinical laboratory values.
* `TOFA_IFNscores` - Analysis of RNA-seq-based Interferon (IFN) scores.
* `TOFA_CKNscores` - Analysis of plasma protein-based Cytokine (CKN) scores (MSD platform).
* `TOFA_Metab_AQ` - Analysis of plasma endpoint metabolites Absolute Quantitation (AQ) (UHPLC-MS).
* `TOFA_Skin_scores` - Analysis of Dermatological (skin) scores.
* `TOFA_Neurocognitive` - Analysis of Neurocognitive assessment scores.
* `TOFA_RNAseq_DESeq2` - Analysis of whole blood transcriptomics (RNA-seq).
* `TOFA_Olink` - Analysis of plasma proteomics (Olink platform).
* `TOFA_NULISA` - Analysis of plasma proteomics (NULISA platform).

------------------------------------------------------------------------

## Data Sources
Download each dataset to the appropriate `/data/` directories within each R project.  

### Clinical trial (TOFA) datasets:
* Participant-level metadata: Available on request.
* Visit/Event-level metadata: Available on request.
* Baseline obesity status: Available on request.
* COVID-19 history: Available on request.
* Adverse events reporting data: Available on request.
* PAXgene whole blood RNAseq data (RPKMs): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19954573.svg)](https://doi.org/10.5281/zenodo.19954573) and GEO: [GSE33018](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330188).  
* Plasma endpoint cytokines data (MSD): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20042862.svg)](https://doi.org/10.5281/zenodo.20042862).  
* Plasma metabolite AQ data (UHPLC-MS): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20046361.svg)](https://doi.org/10.5281/zenodo.20046361).  
* Plasma metabolite RQ data (UHPLC-MS): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20043706.svg)](https://doi.org/10.5281/zenodo.20043706).  
* Dermatological (skin) scores: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20077742.svg)](https://doi.org/10.5281/zenodo.20077742)
* Neurocognitive Assessment scores: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20080323.svg)](https://doi.org/10.5281/zenodo.20080323).
* Plasma proteomics data (Olink): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19962923.svg)](https://doi.org/10.5281/zenodo.19962923).  
* Plasma proteomics data (NULISA): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20043773.svg)](https://doi.org/10.5281/zenodo.20043773).  


### Human Trisome Project (HTP) datasets:  
* Participant-level metadata: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19962380.svg)](https://doi.org/10.5281/zenodo.19962380)
* Visit/Event-level metadata: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19962380.svg)](https://doi.org/10.5281/zenodo.19962380)
* PAXgene whole blood RNAseq data (RPKMs): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20044079.svg)](https://doi.org/10.5281/zenodo.20044079) and GEO: [GSE190125](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190125).  
* Plasma cytokines data (MSD): [Synapse 31475487](https://doi.org/10.7303/syn31475487). 
* Plasma metabolite AQ data (UHPLC-MS): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20074289.svg)](https://doi.org/10.5281/zenodo.20074289).   
* Plasma proteomics data (Olink): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20046326.svg)](https://doi.org/10.5281/zenodo.20046326).  
* Plasma proteomics data (NULISA: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20043943.svg)](https://doi.org/10.5281/zenodo.20043943).  


These datasets originate from the Linda Crnic Institute for Down Syndrome's [Human Trisome Project](https://www.trisome.org/) and are also available on the [INCLUDE Data Hub](https://portal.includedcc.org) [DOI: 10.71738/p0a9-2v09](https://doi.org/10.71738/p0a9-2v09).  

------------------------------------------------------------------------

## Software & Dependencies
* [R](https://cran.r-project.org/)  
* [RStudio](https://posit.co/download/rstudio-desktop)  
Key packages include:
* renv
* tidyverse  
* ggplot2  

The renv.lock files within each analysis project directory contains a full list of packages and versions.

------------------------------------------------------------------------

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

------------------------------------------------------------------------

## Citation & License
If you use this code, please cite:  
**Manuscript**
A phase 2, open-label, single-arm clinical trial of the safety and efficacy A phase 2, open-label, single-arm clinical trial of the safety and efficacy.
Authors, Journal, Year. DOI (UPDATE once available)

**Code**  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19906649.svg)](https://doi.org/10.5281/zenodo.19906648)

This project is licensed under the MIT License – see the LICENSE file for details.

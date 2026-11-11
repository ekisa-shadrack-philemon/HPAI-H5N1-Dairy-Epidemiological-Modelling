
# HPAI H5N1 in Dairy Cattle: Longitudinal Analysis of Herd-Level Milk Production, Cell Count and Cow-Level Rumination 

[![R Version](https://img.shields.io/badge/R-v4.2%2B-blue.svg)](https://www.r-project.org/)
[![Python Version](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Project Overview

This repository contains the analytical code, figures, and workflow documentation for a longitudinal epidemiological investigation of a commercial dairy herd in Minnesota, USA infected by **Highly Pathogenic Avian Influenza (HPAI H5N1)**.

The study evaluates how cow-level behavioral abnormalities (rumination), milk quality indicators (cell count), and herd-level milk production disruptions relate temporally to the H5N1 outbreak. Clinical signs were first observed on **19 June 2024**, with laboratory confirmation on **9 July 2024**.

The analysis combines longitudinal statistical modelling, time-series methods, and Change Point analysis to investigate the potential of milk production, cell count and rumination patterns as indicators of HPAI H5N1 infection in dairy cattle.

---

---

## Related Manuscript

A manuscript arising from this research analysis, titled:

**"The Temporal Relationship Between Highly Pathogenic Avian Influenza (HPAI H5N1) Outbreak and Changes in Milk Production"**

is currently **under peer review**.

The manuscript builds on the longitudinal analysis presented in this repository, examining the temporal association between the HPAI H5N1 outbreak and changes in dairy herd milk production.

Further details and publication information will be provided following the completion of the peer-review process.

## Key Methodological Features

- **Pre-Clinical Baseline and Trend Fitting:** Modelled baseline lactation dynamics using non-linear thin-plate regression splines to account for seasonal variation in the normal milk production cycle.

- **Interrupted Time-Series (ITS) GAMM:** Quantified  milk-yield losses, maximum single-day declines, and sustained recovery trajectories following clinical onset.

- **Autoregressive Residual Structure:** Accounted for temporal autocorrelation in daily milk production data using an explicit correlation structure, `nlme::corARMA(p = 2, q = 0)`.

- **Cow-Level Behavioral Tracking:** Aggregated high-frequency sensor records into daily means and identified abnormal rumination relative to individual cow pre-clinical baselines using a standardized threshold (Z ≤ -2).

- **Change-Point and Statistical Diagnostics:** Assessed temporal stability and changes using AMOC change-point procedures, Welch's t-tests, bootstrap confidence intervals, and residual diagnostics.

---

## Repository Architecture

```text
HPAI-H5N1-Dairy-Epidemiological-Modelling/
├── code/                                    # Analytical scripts (Python preprocessing & R statistical modeling)
│   ├── H5N1_Final_ITS-GAMM_AR(2)_Analysis.R # Main R script for herd-level GAMM and AR(2) and change point
│   └── exploratory_analysis.py              # Python scripts for data cleaning, aggregation, and rumination analysis
├── Tables/                                  # Summary tables and statistical outputs
├── results/                                 # Processed model outputs, CSV summaries.
├── Figures/
│   └── PNG/                                 # Plots
└── README.md                                # Research documentation and setup guide
```
---

## Software and tools

The project uses the following analytical tools:

- **R:** Statistical modelling, generalized additive mixed models (GAMMs), interrupted time-series analysis, and diagnostic procedures.
- **Python:** Data preprocessing, exploratory data analysis, cow-level analysis and visualization.
- **Key R Packages:** `mgcv`, `nlme`, `dplyr`, `ggplot2`, and other statistical analysis packages.
- **Python Libraries:** `pandas`, `numpy`, `matplotlib`, and `seaborn`.

---

## Data Availability and Access

The datasets used in this research are **secondary data obtained from external sources** and are not publicly distributed through this repository.

Access to the underlying datasets is restricted and subject to the relevant data ownership, permissions, and access conditions.

Researchers interested in accessing the data for verification, replication, or further academic investigation should contact the **Academic Supervisor** to discuss the appropriate access procedures and obtain the necessary authorization.

---

## Project Team & Authors

This research report was submitted to the **Department of Mathematics, Gulu University**, in partial fulfillment of the requirements for the award of a Bachelor of Science Education (Physical).

### Research Team

| Name | Registration Number | Role |
|---|---|---|
| **Ekisa Shadrack Philemon** | 23/U/0027/SEG | Lead Analyst |
| **Batwala George** | 23/U/0017/SEG | Collaborator |
| **Tumwesigye Osward** | 23/U/0110/SEG | Collaborator |

### Academic Supervisor

**Dr. Ssematimba Amos**  
Senior Lecturer, Department of Mathematics  
Gulu University, Uganda.

---

## Research Report

The complete research report is titled:

**HPAI H5N1 in Dairy Cattle: Longitudinal Analysis of Herd-Level Milk Production, Cell Count and Cow-Level Rumination**

It documents the research background, methodology, statistical analyses, Results, Discussion, and conclusions.

---

## License

This project is licensed under the MIT License. See the [MIT License](https://opensource.org/licenses/MIT) for details.

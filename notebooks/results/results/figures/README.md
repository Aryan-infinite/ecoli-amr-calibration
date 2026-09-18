# AMR Prediction Under Geographic and Population Shift

## Overview

Antimicrobial resistance (AMR) prediction models can perform well on the data used to develop them, but their performance may change when applied to genetically or geographically different bacterial populations.

This project investigates the generalizability and calibration of machine learning models for antimicrobial resistance prediction in *Escherichia coli*, with a focus on whether model performance remains consistent across genetically distinct bacterial populations.

The analysis uses publicly available genomic and antimicrobial susceptibility data and evaluates model performance across different population groups.

> **Project status:** Analysis in progress.
## Research Question

**Do machine learning models for antimicrobial resistance prediction maintain their performance and calibration when evaluated across genetically distinct populations of *Escherichia coli*?**
## Study Focus

* **Organism:** *Escherichia coli*
* **Antibiotic:** Ciprofloxacin
* **Data:** Publicly available bacterial genomic and antimicrobial susceptibility data
* **Approach:** Machine learning-based antimicrobial resistance prediction
* **Main focus:** Model generalizability and calibration across genetically distinct bacterial populations
## Motivation

Machine learning has increasingly been explored as a tool for predicting antimicrobial resistance from bacterial genomic data. However, bacterial populations can differ substantially in their genetic backgrounds and population structures.

A model trained on one set of bacterial populations may therefore perform differently when applied to another population. Evaluating this variation is important for understanding how reliably genomic AMR prediction models can generalize beyond their training data.

This project examines that problem using *E. coli* and ciprofloxacin resistance, focusing on performance and calibration across genetically distinct populations.
## Methods

The project uses publicly available *E. coli* genomic and antimicrobial susceptibility data to construct and evaluate machine learning models for ciprofloxacin resistance prediction.

The analysis includes:

1. **Data preparation**
   Genomic and phenotype data are processed and filtered to construct the study cohort.

2. **Population structure analysis**
   Bacterial isolates are examined according to their genetic population structure to assess differences between groups.

3. **Feature construction**
   Genomic features are prepared for use in machine learning models.

4. **Model development**
   Machine learning models are trained to predict ciprofloxacin resistance.

5. **Generalizability assessment**
   Model performance is evaluated across genetically distinct populations to examine whether predictive performance changes under population shift.

6. **Calibration assessment**
   Predicted probabilities are evaluated to determine how well they correspond to observed resistance outcomes across populations.

The specific preprocessing procedures, models, evaluation metrics, and final analytical decisions will be documented after the analysis is completed.
## Project Status

The project is currently in the analysis stage.

The data preparation and computational workflow are being developed in Google Colab, with the project repository maintained on GitHub for version control and reproducibility.

Final model results, calibration analyses, figures, and conclusions will be added after the analytical workflow has been completed and verified.
## Data

The project uses publicly available genomic and antimicrobial susceptibility data for *Escherichia coli*.

The raw dataset is not fully stored in this repository because some files exceed GitHub's file-size limits. The repository therefore contains documentation describing the data and its use, while large raw files are kept outside the repository.

The exact data sources, filtering criteria, cohort construction, and preprocessing steps will be documented as the analysis is finalized.
## Reproducibility

The computational analysis is being conducted in Python using Google Colab. The analysis notebooks and supporting materials will be maintained in this repository through Git version control.

The final repository will document the data sources, preprocessing steps, model development, evaluation procedures, and analytical workflow required to reproduce the study.
## Repository Structure

```text
amr-calibration-geographic-shift/
├── data/
│   └── raw/
├── notebooks/
├── results/
│   └── figures/
├── .gitignore
├── LICENSE
└── README.md
```

The repository structure will be expanded as the analysis progresses and additional scripts, tables, figures, and documentation are finalized.
## Citation

This repository contains an ongoing research project. A formal citation will be added when the study and its accompanying manuscript are finalized.


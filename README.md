# Major Research Project: Decoding Neurodegeneration: Leveraging Machine Learning Approaches to Classify Single Cells and Identify Transcriptomic Features in Alzheimer’s Disease and ALS

Welcome to the **Major Research Project** repository. This repository contains all the scripts and data necessary to reproduce the analyses presented in the thesis. The README file serves as a bridge between the methods section of the thesis and the scripts available in this GitHub repository, guiding you through the workflow and linking each methodological step to the corresponding scripts.

## Table of Contents

1. [Introduction](#introduction)
2. [Data](#data)
   - [Data Acquisition](#data-acquisition)
   - [Data Preprocessing](#data-preprocessing)
3. [Feature Extraction Methods](#feature-extraction-methods)
   - [Random Feature Selection](#random-feature-selection)
   - [Most Variable Features Selection](#most-variable-features-selection)
   - [Principal Component Analysis (PCA)](#principal-component-analysis-pca)
   - [Differential Expression (DE) and Weighted Gene Co-expression Network Analysis (WGCNA)](#differential-expression-de-and-weighted-gene-co-expression-network-analysis-wgcna)
4. [Model Training and Evaluation](#model-training-and-evaluation)
5. [Data Combination before Analysis](#data-combination-before-analysis)
6. [Differential Expression Analysis](#differential-expression-analysis)
7. [Enrichment Analysis](#enrichment-analysis)
   - [Gene Ontology (GO) Biological Process Enrichment](#gene-ontology-go-biological-process-enrichment)
   - [KEGG Pathway Enrichment](#kegg-pathway-enrichment)
   - [Semantic Similarity Analysis](#semantic-similarity-analysis)
8. [Visualisation](#visualisation)
9. [Supplemental Notes](#supplemental-notes)

---

## Introduction

This repository encompasses the complete analytical workflow for studying the classification of single cells and assessing its possibilites as a candidate biomarker discovery tool in neurodegenerative diseases, specifically Alzheimer's Disease (AD), C9ALS, and SALS. The analysis leverages single-cell RNA sequencing (scRNA-seq) data to identify key genes and pathways associated with these conditions.

---

## Data

### Data Acquisition

For this study, two single-cell RNA sequencing (scRNA-seq) datasets were utilised:

1. **Alzheimer's Disease (AD) Dataset**:
   - **Description**: Single-nucleus RNA sequencing (snRNA-seq) data from the dorsolateral prefrontal cortex (DLPFC) of 15 individuals (7 AD patients and 8 controls).
   - **Source**: Originally described by Anderson et al. (2023)[12], available on Gene Expression Omnibus (GEO) under accession number [GSE214911](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214911).
   - **Data File**: `SeuratObject.rds`

2. **Amyotrophic Lateral Sclerosis (ALS) Dataset**:
   - **Description**: scRNA-seq data from approximately 380,000 nuclei isolated from the primary motor cortex of 64 individuals (17 controls and 47 ALS patients).
   - **Subtypes**: C9ALS and SALS.
   - **Source**: Originally described by Pineda et al. (2021)[25][8], accessible via Synapse ID [syn53421586](https://www.synapse.org/#!Synapse:syn53421586).
   - **Data Files**: `C9ALS_Seurat_Downsampled.rds`, `SALS_Seurat_Downsampled.rds`

### Data Preprocessing

All data preprocessing and analyses were conducted using the [Seurat](https://satijalab.org/seurat/) package (version 5.1.0). The preprocessing steps include:

- Loading datasets and standardising metadata.
- Removing non-RNA assays (e.g., ATAC-seq data) to focus on RNA data.
- Dividing the ALS dataset into C9ALS and SALS subsets.
- Identifying and removing outlier samples to balance group sizes.

**Relevant Scripts:**
- `{SeparateSeuratObject.Rmd}`
- `{RenameCelltypes.R}`

---

## Feature Extraction Methods

To manage the high dimensionality of scRNA-seq data, various feature extraction methods were employed:

### Random Feature Selection

**Purpose**: To evaluate classifier performance using randomly selected gene subsets, serving as a baseline.

**Scripts:**
- `{Classify500randomfeatures.R}`
- `{Classify1000randomfeatures.R}`
- `{Classify2000randomfeatures.R}`
- `{Classify3000randomfeatures.R}`

*Note: Separate scripts were created to facilitate simultaneous execution.*

### Most Variable Features Selection

**Purpose**: To select genes with the highest variability across the dataset, assuming they are more informative for classification.

**Scripts:**
- `{ClassifyMostVariableFeaturesParallel.R}`

### Principal Component Analysis (PCA)

**Purpose**: To reduce data dimensionality by transforming original genes into principal components capturing major variance sources.

**Scripts:**
- `{ClassifyTopPCsParallel.R}`

### Differential Expression (DE) and Weighted Gene Co-expression Network Analysis (WGCNA)

**Purpose**: To identify genes differentially expressed between conditions and co-expressed gene modules associated with disease status.

**Scripts:**
- `{ClassifyDEWGCNAParallel.R}`
- `{DE.R}`
- `{WGCNA.R}`

---

## Model Training and Evaluation


**Purpose**: To assess model performance and prevent data leakage using a five-fold cross-validation strategy.
**Purpose**: To train logistic regression models with L1 regularisation for feature selection.
**Purpose**: To train random forest classifiers for robust feature importance estimation.
**Purpose**: To evaluate models using metrics such as Accuracy, AUC, Precision, Recall, and F1 Score.

**Scripts:**
- `{TrainClassifier.R}`
- `{PlotROC.R}`

---

## Data Combination before Analysis

**Purpose**: To extract and combine all classification results and feature importance scores from trained models in two large dataframes for further analysis

**Scripts:**
- `{CreateCompleteResultsFile.Rmd}`
- `{CreateCompleteCoefsImportances.Rmd}`

---

## Differential Expression Analysis

**Purpose**: To identify differentially expressed genes using DESeq2 on pseudobulk RNA-seq data and compare these with classifier-identified important genes.

**Scripts:**
- `{DifferentialExpressionAnalysis.Rmd}`

---

## Enrichment Analysis

### Gene Ontology (GO) Biological Process Enrichment

**Purpose**: To identify biological processes enriched among important genes using the clusterProfiler package.

**Scripts:**
- `{EnrichmentAnalysis.Rmd}`

### KEGG Pathway Enrichment

**Purpose**: To identify KEGG pathways enriched among important genes.

**Scripts:**
- `{EnrichmentAnalysis.Rmd}`

### Semantic Similarity Analysis

**Purpose**: To reduce redundancy in enriched GO terms using semantic similarity measures and visualize results with treemaps.

**Scripts:**
- `{EnrichmentAnalysis.Rmd}`

---

## Visualisation

**Purpose**: To create visual representations of the data and analysis results, including UMAP plots, volcano plots, bar plots, and ROC curves.

**Scripts:**
- `{CreateUMAPplots.Rmd}`
- `{CreateFigures.Rmd}`
- `{PlotROC.R}`
- `{GetTrendlineData.R}`

---

## Supplemental Notes

### R Version, Data Visualization, and Use of Tidyverse Packages

All analyses were conducted using R (version 4.3.2) with various [tidyverse](https://www.tidyverse.org/) packages for data manipulation and visualization

---

## Code Availability

All scripts used in this project are available in this GitHub repository: [https://github.com/user485584/MajorResearchProject](https://github.com/user485584/MajorResearchProject)

**Scripts Overview:**

- **AnalyseResults.Rmd**: Sandbox file for exploratory analysis (can be ignored).
- **Classify500randomfeatures.R**, **Classify1000randomfeatures.R**, **Classify2000randomfeatures.R**: Scripts for random feature selection methods.
- **ClassifyDEWGCNAParallel.R**: DEWGCNA feature extraction method.
- **ClassifyMostVariableFeaturesParallel.R**: Variable features extraction method.
- **ClassifyTopPCsParallel.R**: Top principal components feature extraction method.
- **CreateCompleteResultsFile.Rmd**: Concatenates all classification predictions/results.
- **CreateCompleteCoefsImportances.Rmd**: Concatenates all feature importances across models.
- **CreateFigures.Rmd**: Generates multiple figures for the thesis.
- **CreatePseudobulk.R**: Generates pseudobulk counts necessary for DE and WGCNA analyses.
- **CreateUMAPplots.Rmd**: Creates UMAP plots (Figure 1 in thesis).
- **DE.R**: Performs Differential Expression Analysis (part of DEWGCNA feature extraction method).
- **DifferentialExpressionAnalysis.Rmd**: Conducts DE analysis for all datasets and compares with feature importances.
- **EnrichmentAnalysis.Rmd**: Performs enrichment analysis on feature importances.
- **GetTrendlineData.R**: Generates data for regression plots and statistics.
- **PlotROC.R**: Creates ROC plots for trained models.
- **PrepareDfClassifier.R**: Prepares data frames for classification models.
- **RenameCelltypes.R**: Standardises cell type names across datasets.
- **SaveObjects.R**: Saves objects throughout the analysis pipeline.
- **SeparateSeuratObject.Rmd**: Separates and downsamples ALS datasets.
- **TrainClassifier.R**: Trains and tests classifiers.
- **WGCNA.R**: Performs Weighted Gene Co-expression Network Analysis (part of DEWGCNA feature extraction method).

---
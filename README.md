# 🧬 Introduction
## Overview

The purpose of this project is to explore the similarities, differences, and overlapping features—including differentially expressed genes (DEGs), Gene Ontology (GO) terms, biological pathways, and protein-protein interaction (PPI) networks—between three mutant mouse models of Alzheimer’s disease (AD). This analysis is part of a larger project aimed at understanding overlapping molecular features between Alzheimer’s Disease and spaceflight-related neurodegeneration.

The immediate focus of this project is on establishing which rodent AD models best recapitulate the neurodegenerative conditions associated with space exposure. Achieving this objective is crucial before extending the analysis to include additional AD mouse models, different tissue types, and mice subjected to actual spaceflight. If successful, this project will serve as a springboard for broader investigations into the molecular and functional connections between Alzheimer’s pathology and space-associated neurodegeneration.

This repository is structured into five main sections:
- **Introduction:** Provides an overview of the project and the data sources used, as well as scientific background for the project.  
- **Project Staging:** Details project dependencies, data loading and preperation, exploratory data analysis (EDA), quality control, filtering, and normalization steps.
- **Analyses:** Includes differential expression analysis, functional enrichment analysis, and PPI network analysis. Summarizations of the results for each analysis are also provided here.
- **Results and Discussion:** Synthesizes the results from all analyses and discusses their broader biological context and implications.
- **Reproducibility:** Explains how to reproduce all the analyses using the provided Docker images and scripts.

## Scientific Background

Alzheimer’s Disease (AD) is a progressive neurodegenerative disorder characterized by hallmark pathologies such as amyloid-β plaques, tau protein tangles, synaptic dysfunction, and neuroinflammation. While mice do not naturally develop Alzheimer’s Disease, transgenic mouse models have been developed to recapitulate key pathological features of the disease, enabling researchers to study its underlying mechanisms and evaluate potential therapeutic interventions.

In this project, we focus on three widely used AD mouse models:
- 5xFAD: This model carries five familial AD mutations across APP (amyloid precursor protein) and PSEN1 (presenilin-1) genes. It exhibits aggressive amyloid-β plaque deposition, neuronal loss, and gliosis at an early age, making it an ideal model for studying amyloid pathology.
- 3xTG-AD: This triple transgenic model carries mutations in APP, PSEN1, and MAPT (tau). It develops both amyloid plaques and tau tangles, as well as cognitive impairments, providing a more comprehensive representation of AD pathology.
- PS3O1S: This model harbors mutations in presenilin-1 and other associated genes. It exhibits a slower progression of amyloid-β accumulation and a more signficant accumulation of abnormal tau protein aggregates characteristic of tauopathies like Alzheimer's disease. 

This project is significant not only for Alzheimer’s research but also for understanding the phenomenon of space-associated neurodegeneration, often referred to as "space brain." Extended exposure to microgravity and cosmic radiation has been linked to changes in brain structure, neuroinflammation, oxidative stress, and cognitive impairments in astronauts. These effects share striking parallels with neurodegenerative diseases like Alzheimer’s. To identify mutual features between AD and space-associated neurodegeneration, it is essential to determine which rodent models best reflect the "space brain" condition. By comparing the molecular and functional features of these AD models, we aim to pinpoint those that most closely align with the neurodegenerative changes observed in space-flown mice.

## Data Sources and Rationale for Selection
The data for this project was obtained from three separate studies, each providing valuable insights into different transgenic mouse models of Alzheimer’s disease (AD). These studies are as follows:
- Systematic Phenotyping and Characterization of the 5xFAD mouse model of Alzheimer’s Disease ([GSE168137](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168137))
- Transcriptional, behavioural and biochemical profiling in the 3xTg-AD mouse model reveals a specific signature of amyloid deposition and functional decline in Alzheimer’s disease ([GSE161904](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi))
- CNS cell-type specific gene profiling of aging P301S tau transgenic mice, which are accessible via the Gene Expression Omnibus (GEO) via the following accension numbers ([GSE118523](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi))

In this project we'll analyze bulk RNA-sequencing data from the cortex of 5xFAD, 3xTG-AD, and PS3O1S transgenic mice and their corresponding wild type controls. To ensure robust and meaningful comparisons, several criteria were applied when selecting the above datasets. For example, all three AD mouse models were vbred on hte same background strain (C57BL/6J) to minimize genetic variability, RNA-sequencing data was exclusively derived from the cortex to ensure tissue consistency, and only mice aged 7–10 months were included, as this age range corresponds to the onset or progression of AD-like pathology in these models. Additionally, I've ensured each data sets has enough replicates for each model to give statistically significant DEG results (>=3/condition) and that each data set has high sequencing depth, minimal technical bias, and proper normalization in downstream analysis to ensure accurate comparisons.

By meeting these criteria, we aimed to optimize the comparability between datasets and reduce confounding factors that could skew downstream analyses, such as pathway enrichment or protein-protein interaction network construction. While additional data could theoretically enhance the robustness of the analysis, its inclusion often introduces tradeoffs. For example, some datasets included RNA-seq data from AD mice that were either much younger (~2 months) or older (~14 months) than the selected age range. This variability could obscure the biological signals specific to the age range of interest, leading to less reliable DEG comparisons. Additionally, other datasets provided RNA-seq data from different brain regions, such as the hippocampus, rather than the cortex. Including these datasets would add unwanted variability, as gene expression patterns can vary significantly between brain regions, even within the same model and experimental condition.

# 🧬 Project Staging

This section is broken into three sub-sections:
- **Project Dependencies:** Import the necessary libraries for downstream analyses. 
- **Load, Inspect, and Prepare Data:** Steps to load raw datasets, inspect data structure, and prepare it for downstream processing.
- **Quality Control, Filtering, and Normalization:** Applying quality control measures and data transformations to generate high-quality, reliable data for down stream analysis.

## Dependencies
To start, we'll import the following Python libraries, which will be used in various steps of the analysis. Ensure that these libraries are installed in your Python environment before proceeding.

```python
import os 
import subprocess
import gzip
import pandas as pd #v2.2.2
import mygene #v3.2.2
import numpy as np #v1.26.4
import matplotlib.pyplot as plt
import seaborn as sns #v0.13.2
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import gseapy as gp #v1.1.4
from matplotlib_venn import venn3
import requests #v2.32.3
import networkx as nx #v3.4.2
from io import StringIO

```

## Load, Inspect, and Prepare Data

## Quality Control, Filtering, and Normalization

# 🧬 Analyses
## Differential Expression Analysis (DEA)
## Functional Enrichment Analysis
## PPI Network Analysis

# 🧬 Results and Discussion
## Results
 - analysis resuts
## Discussion 
-  broader context
## Next Steps

# 🧬 Reproducibility

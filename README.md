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

In this project we'll analyze bulk RNA-sequencing data from the cortex of 5xFAD, 3xTG-AD, and PS3O1S transgenic mice and their corresponding wild type controls. To ensure robust and meaningful comparisons, several criteria were applied when selecting the above datasets. All three AD mouse models were bred on the same background strain (C57BL/6J) to minimize genetic variability, RNA-sequencing data was exclusively derived from the cortex to ensure tissue consistency, and only mice aged 7–10 months were included, as this age range corresponds to the onset or progression of AD-like pathology in these models. Additionally, each dataset has enough replicates for to give statistically significant DEG results, has high sequencing depth, and minimal technical bias. By meeting these criteria, we aimed to optimize the comparability between datasets and reduce confounding factors that could skew downstream analyses.

While additional data could theoretically enhance the robustness of the analysis, its inclusion often introduces tradeoffs. For example, some datasets included RNA-seq data from AD mice that were either much younger (~2 months) or older (~14 months) than the selected age range. This variability could obscure the biological signals specific to the age range of interest, leading to less reliable DEG comparisons, and as a result these datasets were excluded. Additionally, other datasets provided RNA-seq data from different brain regions, such as the hippocampus, rather than the cortex. Including these datasets would add unwanted variability, as gene expression patterns can vary significantly between brain regions, even within the same model and experimental condition.

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
Next, we'll  retrieve our RNA-sequencing data for each of the three AD mouse models and their corresponding wild type controls. 

### Load, Inspect, and Prepare Data for 5xFAD Mouse Model
First, we'll start by definng a function to download a gzipped file from a URL, unzip it, and theb load it into pandas data frame, as demonstrated int the code block below. 

```python
def download_and_load_data(url, output_filename, sep="\t", column_filter=None):
    # Step 1: Download the file using wget
    print(f"Downloading {output_filename} from {url}...")
    subprocess.run(["wget", "-O", output_filename + ".gz", url], check=True)

    # Step 2: Gunzip the file
    print(f"Unzipping {output_filename}.gz...")
    with gzip.open(output_filename + ".gz", "rb") as gz_file:
        with open(output_filename, "wb") as out_file:
            out_file.write(gz_file.read())

    # Step 3: Load the data into a pandas DataFrame
    print(f"Loading {output_filename} into a pandas DataFrame...")
    df = pd.read_csv(output_filename, sep=sep, index_col=0)

    # Optional: Filter columns based on the keyword
    if column_filter:
        print(f"Filtering columns with keyword '{column_filter}'...")
        filtered_columns = [col for col in df.columns if column_filter in col]
        df = df[filtered_columns]

    return df

# Load data for 5xFAD mouse model (8mo only)
url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE168137&format=file&file=GSE168137%5FcountList%2Etxt%2Egz"
output_filename = "GSE168137_countList.txt"
column_keyword = "cortex_8mon"
countlist_5xFAD = download_and_load_data(url, output_filename, column_filter=column_keyword)

# Rename columns
new_column_names = ["5xFAD_cortex_8mon_Female_295", "5xFAD_cortex_8mon_Female_312","5xFAD_cortex_8mon_Female_339", "5xFAD_cortex_8mon_Female_341", "5xFAD_cortex_8mon_Female_342", "5xFAD_cortex_8mon_Male_299", "5xFAD_cortex_8mon_Male_300", "5xFAD_cortex_8mon_Male_307", "5xFAD_cortex_8mon_Male_387", "5xFAD_cortex_8mon_Male_390",  "BL6_cortex_8mon_Female_322", "BL6_cortex_8mon_Female_338", "BL6_cortex_8mon_Female_340",  "BL6_cortex_8mon_Female_348", "BL6_cortex_8mon_Female_351", "BL6_cortex_8mon_Male_389", "BL6_cortex_8mon_Male_396", "BL6_cortex_8mon_Male_399", "BL6_cortex_8mon_Male_410", "BL6_cortex_8mon_Male_412"]
countlist_5xFAD.columns = new_column_names

# Drop ensemble version ID from gene_id's
countlist_5xFAD.index = countlist_5xFAD.index.str.split('.').str[0]

# View first 5 rows of data
countlist_5xFAD.head()
```
<img width="1422" alt="Screenshot 2025-01-07 at 10 22 18 AM" src="https://github.com/user-attachments/assets/df7bd76d-90f2-45fb-8a38-481c806ec8a4" />

The data above is indexed by Ensemble gene ID (ENSMUSG) with 20 columns of RNA-sequencing expression data (i.e., counts). Before performing downstream analyses we'll want to convert the Ensemble gene IDs to gene names, as demonstrated in next code block.

```python
# create MyGeneInfo object
mg = mygene.MyGeneInfo()

# get the ensembl id from index
ensembl_ids = countlist_5xFAD.index.tolist()

# query the gene symbols for the ensemble ids and onvert result to dataframe
gene_info = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='mouse')
gene_df = pd.DataFrame(gene_info)

# remove duplicate ensemble ids and rows where symbol is missing or duplicated
gene_df = gene_df.dropna(subset=['symbol']).drop_duplicates(subset='query')

# map gene symbols back to original dataframe and move gene_name column to front column
countlist_5xFAD['Gene_Name'] = countlist_5xFAD.index.map(gene_df.set_index('query')['symbol'])
cols = ['Gene_Name'] + [col for col in countlist_5xFAD.columns if col != 'Gene_Name']
countlist_5xFAD = countlist_5xFAD[cols]
```
Now, we'll display the results to ensure our data frame includes gene names along with ensemble IDs. 

<img width="1315" alt="Screenshot 2025-01-07 at 11 23 49 AM" src="https://github.com/user-attachments/assets/5b3ae70e-421f-46b6-943c-7294eb930a14" />

Next, we'll check for missing data and perform basic data exploration to understand the distribution and variability of RNA sequencing counts across the samples before performing any downstream analysis. First, let's check for missing value:

```python
# check for missing values
print(countlist_5xFAD.isnull().sum())
```

<img width="247" alt="Screenshot 2025-01-07 at 11 40 02 AM" src="https://github.com/user-attachments/assets/4cdcd638-b2ec-48f1-a0cf-5fab0b4378ac" />

Notably, the dataset has no null (missing) values. Next, we'll explore the distribution and variability in our dataset, as demonstrated in the code block below:

```python
# Drop the Gene Name column from countlist_5xFAD for counting
countlist_no_name = countlist_5xFAD.iloc[:, 1:]

# Calculate total counts per sample and log transform counts
total_counts = countlist_no_name.sum(axis=0)
log_counts = countlist_no_name.apply(lambda x: np.log2(x + 1))

# Create subplots
fig, axes = plt.subplots(1, 2, figsize=(18, 6))

# Plot total counts per sample
axes[0].bar(countlist_no_name.columns, total_counts, color='skyblue')  # Use countlist_no_name columns here
axes[0].set_ylabel('Total Counts')
axes[0].set_title('Total Counts per Sample')
axes[0].tick_params(axis='x', rotation=85)

# Plot log transformed counts per sample
log_counts.boxplot(ax=axes[1])
axes[1].set_ylabel('Log2(Counts + 1)')
axes[1].set_title('Log Transformed Counts per Sample')
axes[1].tick_params(axis='x', rotation=85)

plt.tight_layout()
plt.show()
```

<img width="1336" alt="Screenshot 2025-01-07 at 11 24 56 AM" src="https://github.com/user-attachments/assets/cf0e96a2-e09c-4080-bfc7-d43f2606cc0a" />

The chart on the left titled, "Total Counts Per Sample," helps visualize the overall sequencing depth across the samples. Ideally, the bars, representing the total counts, should be of similar height, indicating that sequencing depth is consistent across samples. However, our data suggests uneven sequencing depth across samples, which can impat downstream analyses since samples with higher depth may show higher expression counts, not due to biological differences but due to technical variability. This is a common issue in RNA-seq analysis, and will be addressed via data normalization later on.

The rightmost chart, titled "Log Transformed Total Per Sample," helps assess the variability and distribution of gene expression counts across samples. Log-transforming counts reduces skewness by dampening large values, allowing for more comparable distributions across samples. Here, the boxes (indicating the interquartile range) and whiskers (indicating variability outside the quartiles) appear similar across samples. This suggests that while total counts differ, log transformation has minimized variability among samples, partially masking depth differences and indicating a more consistent distribution across samples.

> Note: You can find the processing count list for the 5xFAD mice (+ WT controls) [here](https://github.com/evanpeikon/mouse_AD_models/tree/main/data/count_lists). 

### Load, Inspect, and Prepare Data for 3xTG-AD Mouse Model

In the code blocks below, I'm going to load, clean, and inspect the data for the 3xTG-AD mouse model. In the last sub-section I explained the rationale behind each of the preceding steps in detail, so in this section I'll provide minimal additional commentary as only minor modifications are made to the code from the previous sub-section. 

```
# Load data for 3xTG-AS mouse model (8mo only)
url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE161904&format=file&file=GSE161904%5FRaw%5Fgene%5Fcounts%5Fcortex%2Etxt%2Egz"
output_filename = "GSE161904_Raw_gene_counts_cortex.txt"
column_keyword = "G3" #keyword for 8mo old data
countlist_3xTG_AD = download_and_load_data(url, output_filename, column_filter=column_keyword)

# Rename columns
new_column_names = ['3xTG_AD_Cortex_R1', '3xTG_AD_Cortex_R3', '3xTG_AD_Cortex_R4', 'WT_Cortex_R10', 'WT_Cortex_R17', 'WT_Cortex_R9']
countlist_3xTG_AD.columns = new_column_names

# Add column name to index
countlist_3xTG_AD.index.set_names('gene_id', inplace=True)

# create MyGeneInfo object
mg = mygene.MyGeneInfo()

# get the ensembl id from index
ensembl_ids = countlist_3xTG_AD.index.tolist()

# query the gene symbols for the ensemble ids and onvert result to dataframe
gene_info = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='mouse')
gene_df = pd.DataFrame(gene_info)

# remove duplicate ensemble ids and rows where symbol is missing or duplicated
gene_df = gene_df.dropna(subset=['symbol']).drop_duplicates(subset='query')

# map gene symbols back to original dataframe and move gene_name column to front column
countlist_3xTG_AD['Gene_Name'] = countlist_3xTG_AD.index.map(gene_df.set_index('query')['symbol'])
cols = ['Gene_Name'] + [col for col in countlist_3xTG_AD.columns if col != 'Gene_Name']
countlist_3xTG_AD = countlist_3xTG_AD[cols]

# view first 5 rows of data to make sure gene_names are added
countlist_3xTG_AD.head()
```

<img width="836" alt="Screenshot 2025-01-07 at 11 30 16 AM" src="https://github.com/user-attachments/assets/5ea86ce9-f1bc-42e3-ae8a-9890b18cd6fb" />

Now that we've added gene names to our data frame, we'll check for missing values. 

```python
# check for missing values
print(countlist_3xTG_AD.isnull().sum())
```

<img width="172" alt="Screenshot 2025-01-07 at 11 40 23 AM" src="https://github.com/user-attachments/assets/d95588a0-632e-458f-8b8c-5ba91a305ca1" />

Again, the dataset has no null (missing) values, so, we'll move on to exploring the distribution and variability in our dataset. 

```python
# Drop the Gene Name column from countlist_3xTG_AD for counting
countlist_no_name = countlist_3xTG_AD.iloc[:, 1:]

# Calculate total counts per sample and log transform counts
total_counts = countlist_no_name.sum(axis=0)
log_counts = countlist_no_name.apply(lambda x: np.log2(x + 1))

# Create subplots
fig, axes = plt.subplots(1, 2, figsize=(18, 6))

# Plot total counts per sample
axes[0].bar(countlist_no_name.columns, total_counts, color='skyblue')  # Use countlist_no_name columns here
axes[0].set_ylabel('Total Counts')
axes[0].set_title('Total Counts per Sample')
axes[0].tick_params(axis='x', rotation=85)

# Plot log transformed counts per sample
log_counts.boxplot(ax=axes[1])
axes[1].set_ylabel('Log2(Counts + 1)')
axes[1].set_title('Log Transformed Counts per Sample')
axes[1].tick_params(axis='x', rotation=85)

plt.tight_layout()
plt.show()
```

> Note: You can find the processing count list for the 3xTG-AD mice (+ WT controls) [here](https://github.com/evanpeikon/mouse_AD_models/tree/main/data/count_lists). 

### Load, Inspect, and Prepare Data for PS3O1S Mouse Model

In this section we'll retrieve the data for the 10mo old PS301S mice and their corresponding wild type controls. 

```python
# Load data for 3xTG-AS mouse model (8mo only)
url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE118523&format=file&file=GSE118523%5F20161109%5Fold%5Fwt%5Ftg%2Ecsv%2Egz"
output_filename = "GSE118523_20161109_old_wt_tg.csv"
column_keyword = ""
countlist_PS3O1S = download_and_load_data(url, output_filename, sep=',', column_filter=column_keyword)

# View first 5 rows of data
countlist_PS3O1S.head()
```

<img width="1136" alt="Screenshot 2025-01-07 at 11 43 00 AM" src="https://github.com/user-attachments/assets/769c05ef-54ae-48a8-867a-a2522f420c70" />

Notably, upon downloading the data for the PS3O1S mice we can see that the researchers did not provide a countlist in GEO and instead this CSV contains the output from differential expression analysis. As a result, we'll save this file as DEA_PS3O1S for reference at a later stage of this project.

```
# rename countlist_PS3O1S to DEA_PS3O1S
DEA_PS3O1S = countlist_PS3O1S
```

> Note: You can find the file containing the differential expression analysis results from the PS3O1S mice [here](https://github.com/evanpeikon/mouse_AD_models/tree/main/data/DEA_results). 

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

# 🧬 Introduction
## Overview

The purpose of this project is to explore the similarities, differences, and overlapping features—including differentially expressed genes (DEGs), Gene Ontology (GO) terms, biological pathways, and protein-protein interaction (PPI) networks—between three mutant mouse models of Alzheimer’s disease (AD). This analysis is part of a larger project aimed at understanding overlapping molecular features between Alzheimer’s Disease and spaceflight-related neurodegeneration.

The immediate focus of this project is on establishing which rodent AD models best recapitulate the neurodegenerative conditions associated with space exposure. Achieving this objective is crucial before extending the analysis to include additional AD mouse models, different tissue types, and mice subjected to actual spaceflight. If successful, this project will serve as a springboard for broader investigations into the molecular and functional connections between Alzheimer’s pathology and space-associated neurodegeneration.

This repository is structured into four main sections:
- **Introduction:** Provides an overview of the project and the data sources used, as well as scientific background for the project.  
- **Project Staging:** Details project dependencies, data loading and preperation, exploratory data analysis (EDA), quality control, filtering, and normalization steps.
- **Analyses:** Includes differential expression analysis, functional enrichment analysis, and PPI network analysis. Summarizations of the results for each analysis are also provided here.
- **Results and Discussion:** Synthesizes the results from all analyses and discusses their broader biological context and implications.

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
- Transcriptional, behavioural and biochemical profiling in the 3xTg-AD mouse model reveals a specific signature of amyloid deposition and functional decline in Alzheimer’s disease ([GSE161904](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161904))
- CNS cell-type specific gene profiling of aging P301S tau transgenic mice, which are accessible via the Gene Expression Omnibus (GEO) via the following accension numbers ([GSE118523](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118523))
  
In this project we'll analyze bulk RNA-sequencing data from the cortex of 5xFAD, 3xTG-AD, and PS3O1S transgenic mice and their corresponding wild type controls. To ensure robust and meaningful comparisons, several criteria were applied when selecting the above datasets. All three AD mouse models were bred on the same background strain (C57BL/6J) to minimize genetic variability, RNA-sequencing data was exclusively derived from the cortex to ensure tissue consistency, and only mice aged 7–10 months were included, as this age range corresponds to the onset or progression of AD-like pathology in these models. Additionally, each dataset has enough replicates for to give statistically significant DEG results, has high sequencing depth, and minimal technical bias. By meeting these criteria, we aimed to optimize the comparability between datasets and reduce confounding factors that could skew downstream analyses.

While additional data could theoretically enhance the robustness of the analysis, its inclusion often introduces tradeoffs. For example, some datasets included RNA-seq data from AD mice that were either much younger (~2 months) or older (~14 months) than the selected age range. This variability could obscure the biological signals specific to the age range of interest, leading to less reliable DEG comparisons, and as a result these datasets were excluded. Additionally, other datasets provided RNA-seq data from different brain regions, such as the hippocampus, rather than the cortex. Including these datasets would add unwanted variability, as gene expression patterns can vary significantly between brain regions, even within the same model and experimental condition.

> Note: You can find a non-technical write up of this project, and the results of my analysis, [here](https://decodingbiology.substack.com/p/unraveling-molecular-pathways-in).

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

```python
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

```python
# rename countlist_PS3O1S to DEA_PS3O1S
DEA_PS3O1S = countlist_PS3O1S
```

> Note: You can find the file containing the differential expression analysis results from the PS3O1S mice [here](https://github.com/evanpeikon/mouse_AD_models/tree/main/data/DEA_results). 

## Quality Control, Filtering, and Normalization

The next step in our analysis is to filter out genes with low expression levels across all samples, which can introduce noise in the data. By filtering these out, we can make our results more reliable and improve statistical power, making detecting real biological differences between conditions easier. Additionally, filtering out genes with low expression counts decreases computational load by reducing the number of genes in the dataset, making future downstream analyses faster. To determine the optimal filtering criteria, we'll plot the number of genes retained with different filtering criteria, then use a function named filter_normalize to filter and normalize our data.

### Quality Control, Filtering, and Normalization for 5xFAD Mouse Model

```python
# plot the number of genes retained as a function of differnet CPM thresholds
def plot_genes_retained_by_cpm(data, min_samples=2):
    # convert raw counts to CPM to normalize the data
    cpm = data.apply(lambda x: (x / x.sum()) * 1e6) #convert raw counts to CPM to normalize
    # define a range of CPM thresholds to test, from 0 to 5 with increments of 0.1
    thresholds = np.arange(0, 5, 0.1)
    # initialize list to store the # of genes retained for ea/ threshold
    genes_retained = []

    # loop through ea/ threshold value to determine the # of genes retained
    for min_cpm in thresholds:
        # create mask where CPM > min_cpm in at least min_samples samples
        mask = (cpm > min_cpm).sum(axis=1) >= min_samples
        # count # of genes that meet the criteria and append to the list
        genes_retained.append(mask.sum())

    # plot # of genes retained as a function of CPM threshold
    plt.figure(figsize=(10, 6))
    plt.plot(thresholds, genes_retained, marker='o', color='green')
    plt.axvline(x=1.0, color='red', linestyle='--', label='CPM = 1')
    plt.xlabel('Threshold (CPM)')
    plt.ylabel('Num Genes Retained')
    plt.legend()
    plt.show()
```
```python
# Drop the Gene Name column from countlist_5xFAD for counting
countlist_no_name = countlist_5xFAD.iloc[:, 1:]

# call plot_genes_retained_by_cpm function
plot_genes_retained_by_cpm(countlist_no_name)
```

<img width="651" alt="Screenshot 2025-01-07 at 11 58 07 AM" src="https://github.com/user-attachments/assets/998c0f5e-e549-4510-9770-7ad8dde486c9" />

Based on the data in the chart above, we'll filter genes with an expression threshold of <0.70 CPM. For many bulk RNA-seq datasets, a CPM threshold of 1 is a common filtering point, but 0.70 is slightly more lenient is justifiable given the distribution of our data. Now, In the code block below, well perform filtering and normalization of our data, then visualize the results:

```python
def filter_normalize(data, min_cpm=1.0, min_samples=2):
    # Separate the gene_name column
    gene_names = data.iloc[:, 0]  # First column is gene_name
    raw_counts = data.iloc[:, 1:]  # Remaining columns are raw counts

    # Convert raw counts to CPM
    cpm = raw_counts.apply(lambda x: (x / x.sum()) * 1e6, axis=0)

    # Filter genes based on CPM thresholds
    mask = (cpm > min_cpm).sum(axis=1) >= min_samples  # Keep genes with CPM > min_cpm in at least min_samples
    filtered_counts = raw_counts[mask]
    filtered_gene_names = gene_names[mask]

    # Compute geometric mean of non-zero values for each gene
    geometric_means = filtered_counts.apply(lambda row: np.exp(np.log(row[row > 0]).mean()), axis=1)

    # Calculate size factors by dividing each gene expression by its geometric mean
    size_factors = filtered_counts.div(geometric_means, axis=0).median(axis=0)

    # Normalize data by dividing each gene expression by the size factors
    normalized_counts = filtered_counts.div(size_factors, axis=1)

    # Add back the gene_name column
    normalized_data = pd.concat([filtered_gene_names, normalized_counts], axis=1)

    # Return normalized data as a DataFrame
    return normalized_data
```
```python
# Apply the function to filter and normalize the data
filtered_normalized_5xFAD = filter_normalize(countlist_5xFAD, min_cpm=0.70)

# Plot the distribution of data after normalization
fig, axes = plt.subplots(1, 2, figsize=(18, 6))

# Total normalized counts per sample
total_counts_normalized = filtered_normalized_5xFAD.iloc[:, 1:].sum(axis=0)  # Exclude gene_name column
axes[0].bar(filtered_normalized_5xFAD.columns[1:], total_counts_normalized, color='lightcoral')
axes[0].set_ylabel('Total Normalized Counts')
axes[0].set_title('Total Counts per Sample (Normalized)')
axes[0].tick_params(axis='x', rotation=85)

# Log-transformed normalized counts per sample
log_normalized_data = filtered_normalized_5xFAD.iloc[:, 1:].apply(lambda x: np.log2(x + 1), axis=0)  # Exclude gene_name column
log_normalized_data.boxplot(ax=axes[1])
axes[1].set_ylabel('Log2(Normalized Counts + 1)')
axes[1].set_title('Log Transformed Counts per Sample (Normalized)')
axes[1].tick_params(axis='x', rotation=85)

plt.tight_layout()
plt.show()
```

<img width="1223" alt="Screenshot 2025-01-07 at 12 34 42 PM" src="https://github.com/user-attachments/assets/d5a79758-50e2-49c1-9991-5732f09e43ac" />

> Note: You can find the filtered and normalized data for the 5xFAD mouse model (+ WT controls) [here](https://github.com/evanpeikon/mouse_AD_models/tree/main/data/filtered_normalized). 

### Quality Control, Filtering, and Normalization for 3xTG-AD Mouse Model

Now, we'll repeat the steps fron the last sub-section for the 3xTG-AD model. 

```python
# Apply the function to filter and normalize the data
filtered_normalized_3xTG_AD = filter_normalize(countlist_3xTG_AD, min_cpm=0.70, min_samples=1)

# Plot the distribution of data after normalization
fig, axes = plt.subplots(1, 2, figsize=(18, 6))

# Total normalized counts per sample
total_counts_normalized = filtered_normalized_3xTG_AD.iloc[:, 1:].sum(axis=0)  # Exclude gene_name column
axes[0].bar(filtered_normalized_3xTG_AD.columns[1:], total_counts_normalized, color='lightcoral')
axes[0].set_ylabel('Total Normalized Counts')
axes[0].set_title('Total Counts per Sample (Normalized)')
axes[0].tick_params(axis='x', rotation=85)

# Log-transformed normalized counts per sample
log_normalized_data = filtered_normalized_3xTG_AD.iloc[:, 1:].apply(lambda x: np.log2(x + 1), axis=0)  # Exclude gene_name column
log_normalized_data.boxplot(ax=axes[1])
axes[1].set_ylabel('Log2(Normalized Counts + 1)')
axes[1].set_title('Log Transformed Counts per Sample (Normalized)')
axes[1].tick_params(axis='x', rotation=85)

plt.tight_layout()
plt.show()
```

<img width="1221" alt="Screenshot 2025-01-07 at 12 39 37 PM" src="https://github.com/user-attachments/assets/51f457c7-417a-4f36-a709-3760ff80dbb6" />

> Note: You can find the filtered and normalized data for the 3xTG-AD mouse model (+ WT controls) [here](https://github.com/evanpeikon/mouse_AD_models/tree/main/data/filtered_normalized). 

# 🧬 Analyses

This section is broken into three sub-sections:
- **Differential Expression Analysis:** Identifies genes with significant changes in expression between Alzheimer’s disease models and wild-type controls.
- **Functional Enrichment Analysis:** Explores the biological processes, molecular functions, and cellular components associated with differentially expressed genes.
- **PPI Network Analysis:** Constructs protein-protein interaction networks to highlight key hub proteins and assess connectivity among differentially expressed genes.

## Differential Expression Analysis (DEA)
Now that we've loaded our data and performed quality control and normalization, we can perform differential expression analysis to identify differentially expressed genes. In the case of the 5xFAD and 3xTG-AD models, we'll perform pairwise comparisons between the groups of wild type and transgenic mice, apply multiple testing corrections, and identify then differentially expressed gene. For the PS3O1S model, we already have a csv with the outputs of differential expression analysis so we'll simply identify genes that are up or down regualted.

### Differential Expression Analysis for 5xFAD Mouse Model

In the code block below, we'll perform pairwise analyses, apply multiple testing corrections, and identify differentially expressed genes in the 5xFAD mice.

```python
# Separate the groups
treated_columns = ["5xFAD_cortex_8mon_Female_295", "5xFAD_cortex_8mon_Female_312", "5xFAD_cortex_8mon_Female_339", "5xFAD_cortex_8mon_Female_341", "5xFAD_cortex_8mon_Female_342", "5xFAD_cortex_8mon_Male_299", "5xFAD_cortex_8mon_Male_300", "5xFAD_cortex_8mon_Male_307", "5xFAD_cortex_8mon_Male_387", "5xFAD_cortex_8mon_Male_390"]

control_columns = ["BL6_cortex_8mon_Female_322", "BL6_cortex_8mon_Female_338", "BL6_cortex_8mon_Female_340", "BL6_cortex_8mon_Female_348", "BL6_cortex_8mon_Female_351", "BL6_cortex_8mon_Male_389", "BL6_cortex_8mon_Male_396", "BL6_cortex_8mon_Male_399", "BL6_cortex_8mon_Male_410", "BL6_cortex_8mon_Male_412"]

# Initialize a list to store results
results = []

# Perform differential expression analysis for each gene
for gene in filtered_normalized_5xFAD.index:
    # Extract treated and control group data for the current gene
    treated = pd.to_numeric(filtered_normalized_5xFAD.loc[gene, treated_columns], errors='coerce')
    control = pd.to_numeric(filtered_normalized_5xFAD.loc[gene, control_columns], errors='coerce')

    # Drop NaN values if present
    treated = treated.dropna()
    control = control.dropna()

    # Skip genes where either group is empty after coercion
    if treated.empty or control.empty:
        continue

    # Calculate mean expression levels for control and treated groups
    mean_control = np.mean(control)
    mean_treated = np.mean(treated)

    # Compute log2 fold change, adding 1 to avoid log(0)
    log2fc = np.log2((mean_treated + 1) / (mean_control + 1))

    # Perform t-test to compare control and treated groups
    t_stat, p_val = ttest_ind(treated, control)

    # Append results for the current gene to the results list
    results.append({"gene": gene,"Gene_Name": filtered_normalized_5xFAD.loc[gene, "Gene_Name"],"log2fc": log2fc,"t_stat": t_stat,"p_val": p_val})

# Convert results list to DataFrame for easier manipulation
results_df = pd.DataFrame(results)

# Convert p-values to numeric values, coercing errors to NaN if conversion fails
results_df['p_val'] = pd.to_numeric(results_df['p_val'], errors='coerce')
results_df = results_df.dropna(subset=['p_val'])

# Extract p-values as a numpy array for multiple testing correction
pvals = results_df['p_val'].values

# Apply multiple testing correction using Benjamini-Hochberg method
results_df['p_adj'] = multipletests(results_df['p_val'], method='fdr_bh')[1]

# Calculate the absolute value of log2 fold change
results_df['abs_log2fc'] = results_df['log2fc'].abs()

# Filter results to get differentially expressed genes (DEGs) with p_adj < 0.05 and |log2fc| > 1
deg_5xFAD = results_df[(results_df['p_adj'] < 0.075) & (results_df['abs_log2fc'] > 0.75)]
```

The code above performs a differential expression analysis on gene expression data, and the final output, deg_5xFAD, is a DataFrame containing the genes that are significantly differentially expressed between the WT and AD samples.

Notably, the filtering criteria for DEGS is a log2fc > 0.75 and an adjusted p-value of <0.075. These are a little less stringent than the standard criteria of log2fc > 1 and a adjusted p-value of 0.05. A log2FC of 0.75 corresponds to a fold change of about 1.68 (or ~68% increase/decrease in expression), which is biologically significant in many contexts but less strict than a log2FC of 1 (two-fold change). It includes genes with moderate changes in expression, capturing more potential candidates than stricter criteria. Additionally a adjusted pvalue are used to control the false discovery rate (FDR). A threshold of 0.075 is less strict than 0.05, meaning it allows for a slightly higher rate of false positives. This threshold is useful when working with smaller sample sizes, as is the case in this scenario.

Now that we have a dataframe of differentially expressed genes, we can view the distribution of Absolute log fold changes, as demonstrated in the code block below.

```python
# view distribution of scores
sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))
sns.histplot(deg_5xFAD['abs_log2fc'], bins=50, kde=True)
plt.title('Distribution of Absolute log2FC from Differential Expression Analysis')
plt.xlabel('abs_log2fc')
plt.ylabel('Frequency (# of genes)')
plt.show()
```

<img width="663" alt="Screenshot 2025-01-07 at 12 59 37 PM" src="https://github.com/user-attachments/assets/01beaadf-0985-42f5-bcaf-20362bf881ec" />

Notably, the image above displays the total number of genes at each absolute log-2 fold change. We can see that the distribution is "right skewed", which means that most of the data is concentrated on the left side of the distribution, with a longer tail extending to the right. In essence, this means that the majority of genes have a absolute log2 fold change closer to 1. 

We'll also use a volcano plot to better understand our data. A volcano plot is a type of scatter plot commonly used in genomics and other areas of bioinformatics to visualize the results of differential expression analysis and help identify statistically significant changes in gene expression between different conditions. In the code block below, I'll demonstrate how to create a volcano plot using our data frame of filtered differentially expressed genes. The first plot will be from our results_df DataFrame before filtering, then the second plot will be from our deg_5xFAD DataFrame, a after filtering out genes that do not meet our selection criteria:

```python
# Create a figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)

# Subplot 1: Scatterplot for the full results_df
sns.scatterplot(data=results_df, x='log2fc', y='p_adj', hue='log2fc', palette='viridis', alpha=0.9, edgecolor=None, ax=axes[0])
axes[0].axhline(y=0.075, color='red', linestyle='-', linewidth=1)
axes[0].axvline(x=0.75, color='blue', linestyle='-', linewidth=1)
axes[0].axvline(x=-0.75, color='blue', linestyle='-', linewidth=1)
axes[0].set_xlabel('log2 Fold Change')
axes[0].set_ylabel('Adjusted P-value')
axes[0].legend(title='log2 Fold Change', loc='upper left')
axes[0].set_title('All Results')

# Subplot 2: Scatterplot for the filtered DEGs (deg_5xFAD)
sns.scatterplot(data=deg_5xFAD, x='log2fc', y='p_adj', hue='log2fc', palette='viridis', alpha=0.9, edgecolor=None, ax=axes[1])
axes[1].axhline(y=0.075, color='red', linestyle='-', linewidth=1)
axes[1].axvline(x=0.75, color='blue', linestyle='-', linewidth=1)
axes[1].axvline(x=-0.75, color='blue', linestyle='-', linewidth=1)
axes[1].set_xlabel('log2 Fold Change')
axes[1].set_ylabel('Adjusted P-value')
axes[1].legend(title='log2 Fold Change', loc='upper left')
axes[1].set_title('Filtered DEGs')

# Adjust layout for better spacing
plt.tight_layout()

# Show the plot
plt.show()
```

<img width="1195" alt="Screenshot 2025-01-07 at 1 04 04 PM" src="https://github.com/user-attachments/assets/37957ed9-de96-4dc0-bbef-f0e5b278fd53" />

As you can see in the image above, our volcano plot combines two critical pieces of information for each gene: the magnitude of change (fold change) and the statistical significance (p-value) of that change. Specifically, the x-axis on this graph shows the log2 fold change between the control and experimental samples in our pairwise analysis. A positive value indicates an upregulation of a gene in the experimental group compared to the control, and a negative value represents downregulation of a gene in the experimental group compared to the control. Additionally, the y-axis shows the significance of said change in gene expression.

Thus, when viewing this graph, we are most interested in the two boxes formed in the lower left and lower right corners, which represent down-regulated and up-regulated genes with high statistical significance. This type of chart also allows us to see how changing our criteria for defining differentially expressed genes can impact the total number of genes in our dataset.

> Note: You can find the file with DEGS for the 5xFAD model [here](https://github.com/evanpeikon/mouse_AD_models/tree/main/data/DEG). 

### Differential Expression Analysis for 3xTG-AD Mouse Model

Now, we'll perform pairwise analyses, apply multiple testing corrections, and identify differentially expressed genes for the 3xTg-AD mice.

```python
# Separate the groups
treated_columns = ['3xTG_AD_Cortex_R1','3xTG_AD_Cortex_R3','3xTG_AD_Cortex_R4']

control_columns = ['WT_Cortex_R10' ,'WT_Cortex_R17','WT_Cortex_R9' ]

# Initialize a list to store results
results = []

# Perform differential expression analysis for each gene
for gene in filtered_normalized_3xTG_AD.index:
    # Extract treated and control group data for the current gene
    treated = pd.to_numeric(filtered_normalized_3xTG_AD.loc[gene, treated_columns], errors='coerce')
    control = pd.to_numeric(filtered_normalized_3xTG_AD.loc[gene, control_columns], errors='coerce')

    # Drop NaN values if present
    treated = treated.dropna()
    control = control.dropna()

    # Skip genes where either group is empty after coercion
    if treated.empty or control.empty:
        continue

    # Calculate mean expression levels for control and treated groups
    mean_control = np.mean(control)
    mean_treated = np.mean(treated)

    # Compute log2 fold change, adding 1 to avoid log(0)
    log2fc = np.log2((mean_treated + 1) / (mean_control + 1))

    # Perform t-test to compare control and treated groups
    t_stat, p_val = ttest_ind(treated, control)

    # Append results for the current gene to the results list
    results.append({"gene": gene,"Gene_Name": filtered_normalized_3xTG_AD.loc[gene, "Gene_Name"],"log2fc": log2fc,"t_stat": t_stat,"p_val": p_val})

# Convert results list to DataFrame for easier manipulation
results_df = pd.DataFrame(results)

# Convert p-values to numeric values, coercing errors to NaN if conversion fails
results_df['p_val'] = pd.to_numeric(results_df['p_val'], errors='coerce')
results_df = results_df.dropna(subset=['p_val'])

# Extract p-values as a numpy array for multiple testing correction
pvals = results_df['p_val'].values

# Apply multiple testing correction using Benjamini-Hochberg method
results_df['p_adj'] = multipletests(results_df['p_val'], method='fdr_bh')[1]

# Calculate the absolute value of log2 fold change
results_df['abs_log2fc'] = results_df['log2fc'].abs()

# Filter results to get differentially expressed genes (DEGs) with p_adj < 0.05 and |log2fc| > 1
deg_3xTG_AD = results_df[(results_df['p_adj'] < 0.075) & (results_df['abs_log2fc'] > 0.75)]
```
```python
# view distribution of scores
sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))
sns.histplot(deg_3xTG_AD['abs_log2fc'], bins=50, kde=True)
plt.title('Distribution of Absolute log2FC from Differential Expression Analysis')
plt.xlabel('abs_log2fc')
plt.ylabel('Frequency (# of genes)')
plt.show()
```

<img width="669" alt="Screenshot 2025-01-07 at 1 07 54 PM" src="https://github.com/user-attachments/assets/a58925d5-d065-4f3d-a413-5b7755cb658e" />

As you can see in the chart above, there are less total DEGs in the 3xTG-AD mice as compared to the 5xFAD mice, and althought the DEGs for the 3xTG-AD mice also follow a "right skewed" pattern, the peak absolute log2 fold changes are significantly higher than for the 5xFAD mice. Now, we'll create a scatter plot to visualize genes that are significantly up and down regualted.  

```python
# Create a figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)

# Subplot 1: Scatterplot for the full results_df
sns.scatterplot(data=results_df, x='log2fc', y='p_adj', hue='log2fc', palette='viridis', alpha=0.9, edgecolor=None, ax=axes[0])
axes[0].axhline(y=0.075, color='red', linestyle='-', linewidth=1)
axes[0].axvline(x=0.75, color='blue', linestyle='-', linewidth=1)
axes[0].axvline(x=-0.75, color='blue', linestyle='-', linewidth=1)
axes[0].set_xlabel('log2 Fold Change')
axes[0].set_ylabel('Adjusted P-value')
axes[0].legend(title='log2 Fold Change', loc='upper left')
axes[0].set_title('All Results')

# Subplot 2: Scatterplot for the filtered DEGs (deg_3xTG-AD)
sns.scatterplot(data=deg_3xTG_AD, x='log2fc', y='p_adj', hue='log2fc', palette='viridis', alpha=0.9, edgecolor=None, ax=axes[1])
axes[1].axhline(y=0.075, color='red', linestyle='-', linewidth=1)
axes[1].axvline(x=0.75, color='blue', linestyle='-', linewidth=1)
axes[1].axvline(x=-0.75, color='blue', linestyle='-', linewidth=1)
axes[1].set_xlabel('log2 Fold Change')
axes[1].set_ylabel('Adjusted P-value')
axes[1].legend(title='log2 Fold Change', loc='upper left')
axes[1].set_title('Filtered DEGs')

# Adjust layout for better spacing
plt.tight_layout()

# Show the plot
plt.show()
```

<img width="1200" alt="Screenshot 2025-01-07 at 1 10 57 PM" src="https://github.com/user-attachments/assets/567270ef-dcfc-42c3-b0f4-f6c6cd35711e" />

> Note: You can find the file with DEGS for the 3xTG-AD model [here](https://github.com/evanpeikon/mouse_AD_models/tree/main/data/DEG). 

### Differential Expression Analysis for PS3O1S Mouse Model

Finally, we'll filter out DEGS from the dataset DEA_PS3O1S, downloaded from GEO, containing the results from a differential expression analysis and visualize the results. 

```python
# Add a column for absolute log2 fold change
DEA_PS3O1S['abs_log2fc'] = DEA_PS3O1S['log2fc'].abs()

# Filter for DEGs based on log2fc and pval thresholds
deg_PS3O1S = DEA_PS3O1S[(DEA_PS3O1S['abs_log2fc'] > 0.75) & (DEA_PS3O1S['pval'] < 0.075)]
```
```python
# Create a figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)

# Subplot 1: Scatterplot for the full results_df
sns.scatterplot(data=DEA_PS3O1S, x='log2fc', y='pval', hue='log2fc', palette='viridis', alpha=0.9, edgecolor=None, ax=axes[0])
axes[0].axhline(y=0.05, color='red', linestyle='-', linewidth=1)
axes[0].axvline(x=1, color='blue', linestyle='-', linewidth=1)
axes[0].axvline(x=-1, color='blue', linestyle='-', linewidth=1)
axes[0].set_xlabel('log2 Fold Change')
axes[0].set_ylabel('Adjusted P-value')
axes[0].legend(title='log2 Fold Change', loc='upper left')
axes[0].set_title('All Results')

# Subplot 2: Scatterplot for the filtered DEGs (deg_5xFAD)
sns.scatterplot(data=deg_PS3O1S, x='log2fc', y='pval', hue='log2fc', palette='viridis', alpha=0.9, edgecolor=None, ax=axes[1])
axes[1].axhline(y=0.05, color='red', linestyle='-', linewidth=1)
axes[1].axvline(x=1, color='blue', linestyle='-', linewidth=1)
axes[1].axvline(x=-1, color='blue', linestyle='-', linewidth=1)
axes[1].set_xlabel('log2 Fold Change')
axes[1].set_ylabel('Adjusted P-value')
axes[1].legend(title='log2 Fold Change', loc='upper left')
axes[1].set_title('Filtered DEGs')

# Adjust layout for better spacing
plt.tight_layout()

# Show the plot
plt.show()
```

<img width="1194" alt="Screenshot 2025-01-07 at 1 12 12 PM" src="https://github.com/user-attachments/assets/a8ecf9b2-2e2d-426d-9c5a-cd9af2d06f62" />

> Note: You can find the file with DEGS for the PS3O1S model [here](https://github.com/evanpeikon/mouse_AD_models/tree/main/data/DEG). 

### All Group Comparison of Differentially Expressed Genes (DEGs)

Now that we have identified differentially expressed genes for all of our AD mouse models, we can perform the following comparisons:
- Quantify absolute number of differentially expressed genes (DEGs), up-regulated DEGS, and down-regulated DEGs for each group.
- Quantify the number of shared versus unique DEGS across groups to identify overlapping DEGS.

```python
# Define groups and datasets
groups = ['5xFAD', '3xTG-AD', 'PS3O1S']
datasets = [deg_5xFAD, deg_3xTG_AD, deg_PS3O1S]

# Initialize counts for total, up-regulated, and down-regulated DEGs
total_counts = []
up_counts = []
down_counts = []

# Calculate counts for each group
for dataset in datasets:
    total_counts.append(len(dataset))
    up_counts.append(len(dataset[dataset['log2fc'] > 0]))
    down_counts.append(len(dataset[dataset['log2fc'] < 0]))

# Create 2x2 subplots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot total DEGs
axes[0, 0].bar(groups, total_counts, color='steelblue')
axes[0, 0].set_title('Total DEGs')
axes[0, 0].set_ylabel('Number of DEGs')

# Plot up-regulated DEGs
axes[0, 1].bar(groups, up_counts, color='green')
axes[0, 1].set_title('Up-regulated DEGs')

# Plot down-regulated DEGs
axes[1, 0].bar(groups, down_counts, color='red')
axes[1, 0].set_title('Down-regulated DEGs')

# Use gene_name column for comparison
genes_5xFAD = set(deg_5xFAD['Gene_Name'])  # Adjust column name as per your dataset
genes_3xTG_AD = set(deg_3xTG_AD['Gene_Name'])
genes_PS3O1S = set(deg_PS3O1S['gene_name'])

# Create the Venn diagram in the last subplot
axes[1, 1].axis('off')  # Turn off the axis for the Venn diagram subplot
venn = venn3(genes_5xFAD, genes_3xTG_AD, genes_PS3O1S],set_labels=('5xFAD', '3xTG-AD', 'PS3O1S'))

# Customize the diagram (optional)
for text in venn.set_labels:
    if text:  # Check if the label exists
        text.set_fontsize(14)
for text in venn.subset_labels:
    if text:  # Check if the subset label exists
        text.set_fontsize(12)

# Add a title for the figure
fig.suptitle('DEGs Analysis and Venn Diagram of DEGs', fontsize=16)

# Adjust layout and show the plots
plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust for suptitle
plt.show()
```

<img width="903" alt="Screenshot 2025-01-07 at 1 14 38 PM" src="https://github.com/user-attachments/assets/db729963-58bd-4f8b-b6c2-39e1dcb6f3dc" />

The differential expression analysis of the three Alzheimer's disease mouse models—5xFAD, 3xTG-AD, and PS3O1S—reveals interesting insights into the gene expression patterns of each model.

- 5xFAD Model: A total of 444 differentially expressed genes (DEGs) were identified, with the overwhelming majority being up-regulated. This suggests that the 5xFAD model, which typically reflects a more aggressive form of Alzheimer's, shows widespread activation of genes in response to disease pathology, possibly in an attempt to compensate for neurodegenerative changes.
- 3xTG-AD Model: In the 3xTG-AD model, 164 DEGs were detected, with approximately 75% of them being up-regulated. This model is characterized by the presence of amyloid and tau pathologies, and the up-regulation of genes likely reflects cellular responses to these stresses, particularly in the context of amyloid plaques and tau tangles. The fewer number of DEGs compared to 5xFAD suggests a milder or earlier disease stage.
- PS3O1S Model: The PS3O1S model, with 318 DEGs, shows a high percentage (85%) of up-regulated genes, indicating a strong genetic response to disease-related insults, possibly driven by tau pathology. This model's higher proportion of up-regulated genes may reflect compensatory mechanisms during the neurodegenerative process.

Additionally, as you can see in the chart above, there is minimal overlap in DEGs between the models. For example, there are only 7 overlapping genes between 5xFAD and 3xTG-AD suggest that while these models share some common disease mechanisms, their gene expression profiles differ significantly, possibly due to differences in the timing or nature of neurodegenerative pathology. Additionally, there are 8 overlapping genes between 5xFAD and PS3O1S, further suggesting some shared molecular pathways related to the disease processes, but still distinct enough to produce largely unique gene expression patterns in each model. Finally, there is only one overlapping gene between 3xTG-AD and PS3O1S, which indicates very little commonality in the genetic responses between these two models. This may reflect divergent disease processes in the two models or differences in the genes and pathways most affected by tau and amyloid pathologies.

Overall, the low overlap of DEGs between models reflects the complexity and heterogeneity of Alzheimer's disease and suggests that different mouse models may recapitulate different aspects of the disease. This could be due to variations in the timing of pathology onset, the specific proteins expressed (amyloid vs. tau), and the resulting downstream effects on gene expression.

## Functional Enrichment Analysis

Now that we've identified a number of differentially expressed genes, we'll perform a series of functional enrichment analyses, which allow us to identify and interpret the biological process, molecular functions, cellular components, and pathways that are overrepresented or significant in our list of differenrentially expressed genes.

### Gene Ontology (GO) Analysis

The first analysis we'll explore is Gene Ontology (GO) analysis, which categorizes differentially expressed genes according to their associated biological processes, cellular components, and molecular functions. This categorization is based on a structured, hierarchical vocabulary known as the Gene Ontology, which systematically describes gene functions.

While differential expression analysis identifies genes that are up- or down-regulated in response to an intervention, treatment, or drug regimen, GO analysis takes this a step further by linking these genes to broader biological contexts. By grouping genes into functional categories, GO analysis can reveal which biological processes, molecular functions, or cellular components are impacted, offering a more detailed understanding of the mechanisms through which an intervention, treatment, or drug exerts its effects.

First, we'll look at the top biological process for each of our mouse models.

```python
# Define the gene lists for each model (DEGs) here
gene_list_5xFAD = deg_5xFAD['Gene_Name'].dropna().astype(str).tolist()
gene_list_3xTG_AD = deg_3xTG_AD['Gene_Name'].dropna().astype(str).tolist()
gene_list_PS3O1S = deg_PS3O1S['gene_name'].dropna().astype(str).tolist()

# Perform GO enrichment analysis for Biological Process (BP), Molecular Function (MF), and Cellular Component (CC)
enrichment_5xFAD_BP = gp.enrichr(gene_list_5xFAD, gene_sets=['GO_Biological_Process_2018'], organism='mouse')
enrichment_3xTG_AD_BP = gp.enrichr(gene_list_3xTG_AD, gene_sets=['GO_Biological_Process_2018'], organism='mouse')
enrichment_PS3O1S_BP = gp.enrichr(gene_list_PS3O1S, gene_sets=['GO_Biological_Process_2018'], organism='mouse')

# Extract the Biological Process results for each
enrichment_5xFAD_BP_df = enrichment_5xFAD_BP.results
enrichment_3xTG_AD_BP_df = enrichment_3xTG_AD_BP.results
enrichment_PS3O1S_BP_df = enrichment_PS3O1S_BP.results
```
```python
enrichment_5xFAD_BP_df.head()
```

<img width="1224" alt="Screenshot 2025-01-07 at 1 34 49 PM" src="https://github.com/user-attachments/assets/08f02d6c-85a6-4c35-8ad8-8216db62ee80" />

```python
enrichment_3xTG_AD_BP_df.head()
```

<img width="1182" alt="Screenshot 2025-01-07 at 1 35 14 PM" src="https://github.com/user-attachments/assets/5f178c2c-b498-4bc2-9a02-58b7adf28615" />

```python
enrichment_PS3O1S_BP_df.head()
```

<img width="1225" alt="Screenshot 2025-01-07 at 1 35 31 PM" src="https://github.com/user-attachments/assets/22eb252a-4178-4311-bc75-95aace97db98" />

This same analysis was then repeated for molecular functions (MF) and cellular components (CC), which I have not included here, but can be viewed in the .ipynb file for this analysis. After running these analyses, I identified overlapping biological processes, molecular functions, and cellular components for our three AD mouse models, as demonstrated below:

```python
# Find overlapping Biological Process terms across the three datasets
bp_5xFAD_terms = set(enrichment_5xFAD_BP_df['Term'])
bp_3xTG_AD_terms = set(enrichment_3xTG_AD_BP_df['Term'])
bp_PS3O1S_terms = set(enrichment_PS3O1S_BP_df['Term'])

# Find overlapping Molecular Function terms across the three datasets
mf_5xFAD_terms = set(enrichment_5xFAD_MF_df['Term'])
mf_3xTG_AD_terms = set(enrichment_3xTG_AD_MF_df['Term'])
mf_PS3O1S_terms = set(enrichment_PS3O1S_MF_df['Term'])

# Find overlapping Cellular Component terms across the three datasets
cc_5xFAD_terms = set(enrichment_5xFAD_CC_df['Term'])
cc_3xTG_AD_terms = set(enrichment_3xTG_AD_CC_df['Term'])
cc_PS3O1S_terms = set(enrichment_PS3O1S_CC_df['Term'])

# Create subplots plots 
plt.figure(figsize=(18, 6))
plt.subplot(1, 3, 1)  # (number of rows, number of columns, plot index)
venn3([bp_5xFAD_terms, bp_3xTG_AD_terms, bp_PS3O1S_terms], set_labels=('5xFAD', '3xTG-AD', 'PS3O1S'))
plt.title('Overlapping Biological Process Terms')
plt.subplot(1, 3, 2)  # (number of rows, number of columns, plot index)
venn3([mf_5xFAD_terms, mf_3xTG_AD_terms, mf_PS3O1S_terms], set_labels=('5xFAD', '3xTG-AD', 'PS3O1S'))
plt.title('Overlapping Molecular Function Terms')
plt.subplot(1, 3, 3)  # (number of rows, number of columns, plot index)
venn3([cc_5xFAD_terms, cc_3xTG_AD_terms, cc_PS3O1S_terms], set_labels=('5xFAD', '3xTG-AD', 'PS3O1S'))
plt.title('Overlapping Cellular Component Terms')
plt.tight_layout()
plt.show()
```

<img width="1203" alt="Screenshot 2025-01-07 at 1 38 03 PM" src="https://github.com/user-attachments/assets/d7acc92f-5b85-4070-b349-e418dcf0fadd" />

The analysis of the differentially expressed genes (DEGs) in the Alzheimer's Disease mouse models revealed very little overlap between the models. However, when performing Gene Ontology (GO) enrichment analysis for Biological Process (BP), Molecular Function (MF), and Cellular Component (CC), there was substantial overlap in the top enriched terms across the three models. The biologial signiicance of this finding is as explained below:

- Lack of DEG overlap: The minimal overlap in DEGs suggests that the three AD mouse models (5xFAD, 3xTG-AD, and PS3O1S) may exhibit distinct molecular mechanisms or biological pathways in response to Alzheimer's-related pathology. Each model might be influenced by different genetic backgrounds or variations in how they manifest Alzheimer's disease, leading to unique sets of differentially expressed genes.
- Overlap in GO terms: Despite the lack of overlap in DEGs, the substantial overlap in the enriched GO terms for Biological Process (BP), Molecular Function (MF), and Cellular Component (CC) suggests that these models share common pathways or cellular processes that are activated in response to Alzheimer's disease. This could reflect conserved biological responses to neurodegeneration, such as neuroinflammatory pathways, synaptic dysfunction, and cellular stress responses, which are common features of Alzheimer's pathology.
- Implication for therapeutic development: The overlap in GO terms could highlight key biological processes and pathways that might be targeted in therapeutic strategies for Alzheimer's disease. While the specific genes driving these processes differ across models, the fact that similar functional categories (e.g., immune response, cell signaling, and synaptic activity) are involved in all models points to potential universal therapeutic targets or biomarkers for Alzheimer's.

In summary, while the differential gene expression profiles vary across models, the shared functional annotations in BP, MF, and CC indicate that certain fundamental biological processes are universally affected by Alzheimer's pathology. 

### Pathway Analysis

Gene Ontology (GO) analysis is useful for identifying broad biological changes associated with gene expression, but it may not always pinpoint specific pathways affected by drug treatment. To address this limitation, we will also perform pathway analysis. Pathway analysis focuses on identifying signaling and metabolic pathways that are enriched among differentially expressed genes, providing insights into how these genes interact within specific biological pathways.

While GO analysis offers a general overview of biological processes and cellular components, pathway analysis provides a more detailed perspective. It maps differentially expressed genes to established biological pathways, such as those documented in KEGG or Reactome databases. This approach clarifies how genes collaborate within biological systems and highlights key pathways altered by the drug treatment. This detailed understanding is crucial for unraveling complex biological mechanisms and identifying potential therapeutic targets.

```python
# Perform pathway analysis for each dataset
enrichment_5xFAD = gp.enrichr(gene_list_5xFAD, gene_sets=['KEGG_2016'], organism='mouse')
enrichment_3xTG_AD = gp.enrichr(gene_list_3xTG_AD, gene_sets=['KEGG_2016'], organism='mouse')
enrichment_PS3O1S = gp.enrichr(gene_list_PS3O1S, gene_sets=['KEGG_2016'], organism='mouse')

# Retrieve and display results for each analysis
enrichment_df_5xFAD = enrichment_5xFAD.results
enrichment_df_3xTG_AD = enrichment_3xTG_AD.results
enrichment_df_PS3O1S = enrichment_PS3O1S.results
```

As with the previous sub-sections, we can look at the top overexpressed pathways in each individual AD mouse model, which you can view in the .ipynb for this analysis. In the code block below we'll identify overlapping pathways between our mouse models.

```python
# Find overlapping Biological Process terms across the three datasets
bp_5xFAD_terms = set(enrichment_df_5xFAD['Term'])
bp_3xTG_AD_terms = set(enrichment_df_3xTG_AD['Term'])
bp_PS3O1S_terms = set(enrichment_df_PS3O1S['Term'])

# Venn Diagram for Biological Process
plt.figure(figsize=(8, 8))
venn3([bp_5xFAD_terms, bp_3xTG_AD_terms, bp_PS3O1S_terms], set_labels=('5xFAD', '3xTG-AD', 'PS3O1S'))
plt.title('Overlapping Pathway Terms')
plt.show()
```

<img width="463" alt="Screenshot 2025-01-07 at 1 45 10 PM" src="https://github.com/user-attachments/assets/89885329-e17e-4835-94f4-5b04dc132562" />

The high number of unique pathways for each model indicates distinct biological responses or molecular mechanisms of Alzheimer’s disease in each model. This suggests that each model may reflect different aspects of the disease process or variations in genetic background that drive divergent pathway activation. However, the overlap in pathways, especially the 44 shared across all models, points to core biological processes that are commonly disrupted in Alzheimer's pathology, such as inflammation, neurodegeneration, and cellular stress. These shared pathways represent potential therapeutic targets that could have broader efficacy across different Alzheimer's disease models. The shared pathways between 5xFAD and PS3O1S (94) and the smaller overlap with 3xTG-AD (12) may also suggest a closer relationship between the 5xFAD and PS3O1S models in terms of pathway activation, which may indicate similar disease mechanisms or stages.

## PPI Network Analysis

In this last analysis section, we'll map the differentially expressed genes (DEGs) from each model to a common protein-protein interaction (PPI) network to understand how DEGs interact with one another and to identify shared key proteins, hubs, or biological processes across different models. 

```
# Function to fetch PPI data for a list of genes from STRING database
def fetch_string_ppi(gene_list, species='9606'):
    base_url = "https://string-db.org/api/tsv/network"
    genes = "\n".join(gene_list)
    params = {'identifiers': genes,'species': species, 'limit': 1 }
    response = requests.post(base_url, data=params)
    if response.status_code == 200:
        return response.text
    else:
        return None

# Fetch PPI network for your gene lists
gene_list_all = gene_list_5xFAD + gene_list_3xTG_AD + gene_list_PS3O1S
ppi_data = fetch_string_ppi(gene_list_all, species='10090')

# Parse the PPI data into a pandas DataFrame
ppi_df = pd.read_csv(StringIO(ppi_data), sep="\t")

# Filter interactions with a score above 0.7
ppi_df_filtered = ppi_df[ppi_df['score'] > 0.7]

# Create an empty graph
G = nx.Graph()

# Add edges to the graph from the filtered PPI data
for index, row in ppi_df_filtered.iterrows():
    protein1 = row[0]  # ENSEMBL ID of the first protein
    protein2 = row[1]  # ENSEMBL ID of the second protein
    protein_name = row[2]  # This could be a gene/protein name (not used here, but you could use it for node labels)
    G.add_edge(protein1, protein2, weight=row['score'])

# Visualize the network
plt.figure(figsize=(12, 12))
nx.draw_networkx(G, node_size=50, with_labels=False, font_size=10, width=1, alpha=0.7)
plt.title('PPI Network for DEGs')
plt.show()
```

<img width="483" alt="Screenshot 2025-01-07 at 1 50 57 PM" src="https://github.com/user-attachments/assets/af4baa23-80dd-48bd-85d2-ecc765421997" />

Now that we've created a protein-protein interaction network for the DEGs we can identify the most important proteins by degree and betweeness centrality, as demonstrated below.

```python
# Calculate degree centrality for each node
degree_centrality = nx.degree_centrality(G)

# Sort nodes by degree centrality (most important hubs)
sorted_degree_centrality = sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)

# Display the top 10 most important proteins
print("Top 10 most important proteins based on degree centrality:")
for protein, centrality in sorted_degree_centrality[:10]:
    print(f"{protein}: {centrality}")

# Calculate betweenness centrality
betweenness_centrality = nx.betweenness_centrality(G)

# Sort nodes by betweenness centrality
sorted_betweenness = sorted(betweenness_centrality.items(), key=lambda x: x[1], reverse=True)

# Display the top 10 proteins by betweenness centrality
print("\nTop 10 proteins based on betweenness centrality:")
for protein, centrality in sorted_betweenness[:10]:
    print(f"{protein}: {centrality}")
```

<img width="377" alt="Screenshot 2025-01-07 at 1 52 58 PM" src="https://github.com/user-attachments/assets/72a7184a-9ee9-4076-aabf-a9adae905737" />

Interpreting a Protein-Protein Interaction (PPI) network for differentially expressed genes (DEGs) from your Alzheimer’s disease mouse models involves examining both the structural properties of the network and the biological roles of key proteins within it. This analysis helps uncover potential regulators, signaling hubs, and critical pathways altered in your models.

Degree centrality highlights the most connected nodes (proteins) within the PPI network. These nodes, often referred to as "hubs," represent proteins involved in numerous interactions and are crucial for maintaining the structural and functional integrity of the network. High-degree proteins are typically key regulators or scaffolding proteins that influence a wide range of cellular processes. In our analysis, the top proteins based on degree centrality include TNF (tumor necrosis factor), PTPRC, STAT1, CXCL10, CD86, IRF7, FCGR3, CCL5, TYROBP, and IFIT3. These proteins are heavily interconnected and may play central roles in the biological processes disrupted in Alzheimer’s disease models. For example:
- TNF is a well-known inflammatory cytokine involved in neuroinflammation, a hallmark of Alzheimer’s disease. Its high connectivity suggests a central role in mediating the immune response and inflammation in your models.
- PTPRC (CD45) is a key regulator of immune cell signaling and could indicate alterations in microglial or immune cell activity in the brain.
- STAT1 is a transcription factor involved in interferon signaling and immune response, suggesting its role in driving inflammatory cascades.
- CXCL10 and CCL5 are chemokines that regulate immune cell recruitment to the brain, potentially implicating changes in neuroimmune interactions.
- TYROBP (DAP12) is a microglial signaling adapter protein linked to Alzheimer’s pathology, emphasizing its relevance in disease progression.

The presence of these highly connected proteins suggests that inflammation, immune regulation, and microglial activation are central themes in the disrupted biological processes across our models.

Betweenness centrality, on the other hand, measures how often a protein acts as a bridge connecting different parts of the network. Proteins with high betweenness centrality are critical for integrating or coordinating communication between distinct biological modules, often representing key signaling nodes. In our analysis, the top proteins based on betweenness centrality include TNF, P42225, PTPRC, CD74, TYROBP, PLCG1, BRCA1, VAV1, P20029, and ICAM1. These proteins are likely to play pivotal roles in connecting various pathways or processes. For example:
- TNF appears again, reinforcing its importance not just as a hub but also as a signaling bridge coordinating inflammatory pathways.
- CD74 is involved in antigen presentation and may highlight interactions between neuroinflammation and adaptive immunity.
- TYROBP further supports its role as a key microglial signaling protein linking neuroinflammation and phagocytosis pathways.
- PLCG1 (phospholipase C gamma 1) is a key player in intracellular signaling pathways, including those downstream of immune and growth factor receptors, which could suggest disruptions in microglial or neuronal signaling.
- ICAM1 is a cell adhesion molecule involved in immune cell trafficking and neurovascular unit integrity, highlighting potential vascular contributions to Alzheimer’s pathology.

The overlap in high-degree and high-betweenness proteins, such as TNF and TYROBP, suggests that these proteins not only serve as central hubs but also act as key connectors between distinct biological processes. This dual role highlights their critical importance in the dysregulated networks of Alzheimer’s disease. Functionally, these proteins likely mediate crosstalk between neuroinflammatory pathways, microglial activation, and immune signaling, which are well-documented contributors to Alzheimer’s pathology. Ultimately, the PPI network analysis reveals key proteins that could serve as biomarkers or therapeutic targets. For example, targeting TNF or modulating TYROBP activity may offer avenues for interventions aimed at reducing neuroinflammation and restoring network balance in Alzheimer’s disease.

# 🧬 Results and Discussion

## Results

The results of this analysis illuminate the molecular heterogeneity and shared biological themes across three Alzheimer’s disease (AD) mouse models: 5xFAD, 3xTG-AD, and PS3O1S. Each model presents distinct gene expression patterns while converging on common processes central to AD pathology. These findings provide a nuanced understanding of the molecular underpinnings of AD and highlight potential therapeutic targets.

The 5xFAD model, characterized by aggressive amyloid-β plaque deposition, exhibited a robust transcriptional response with 444 differentially expressed genes (DEGs), most of which were up-regulated. This response aligns with the model’s pronounced neuroinflammatory phenotype and rapid disease progression. In contrast, the 3xTG-AD model, which features dual amyloid and tau pathologies, revealed 164 DEGs, approximately 75% up-regulated. The smaller number of DEGs suggests a milder or earlier-stage transcriptional profile consistent with its slower pathology development and cognitive impairments. Meanwhile, the PS3O1S model, emphasizing tau pathology with a slower amyloid-β accumulation, showed 318 DEGs, with 85% up-regulated. This pattern reflects the model’s focus on tau-driven neurodegeneration and compensatory genetic responses.

Despite the distinct DEG profiles, there was minimal overlap across models, with only seven genes shared between 5xFAD and 3xTG-AD, eight between 5xFAD and PS3O1S, and one between 3xTG-AD and PS3O1S. This limited overlap underscores the unique molecular mechanisms represented by each model. However, gene ontology (GO) enrichment analysis revealed substantial convergence at the level of biological processes, molecular functions, and cellular components. All three models exhibited enrichment in pathways related to neuroinflammation, cellular stress responses, and synaptic dysfunction. Immune signaling and receptor binding were common molecular functions, while terms associated with the extracellular matrix, synaptic structures, and immune complexes appeared frequently in the cellular components category. These shared pathways underscore the fundamental processes disrupted in Alzheimer’s pathology.

Protein-protein interaction (PPI) network analysis further highlighted the roles of key proteins in mediating disease processes. Central hub proteins such as TNF, PTPRC, STAT1, CXCL10, and TYROBP emerged as pivotal regulators of neuroinflammation and immune responses. Proteins like TNF and TYROBP not only exhibited high-degree centrality, indicating their extensive connections within the network, but also high-betweenness centrality, reflecting their roles as bridges between distinct biological modules. These findings emphasize the importance of inflammatory and immune pathways in Alzheimer’s pathology and suggest potential targets for therapeutic intervention.

## Discussion and Next Steps

The results highlight the complexity and heterogeneity of Alzheimer’s disease, as modeled by 5xFAD, 3xTG-AD, and PS3O1S. Each model recapitulates distinct aspects of the disease. The 5xFAD model captures the aggressive amyloid pathology and associated immune responses characteristic of advanced AD. The 3xTG-AD model offers a more balanced view of amyloid and tau pathologies but with a milder transcriptional signature, reflecting its intermediate disease phenotype. The PS3O1S model emphasizes tau-driven mechanisms and compensatory genetic responses, presenting a slower disease progression. The minimal overlap in DEGs highlights the molecular diversity across these models, while the shared GO terms and PPI network hubs indicate conserved biological processes central to AD pathology.

These findings also align with prior research identifying neuroinflammation and microglial activation as critical contributors to Alzheimer’s disease. Key regulators like TNF and TYROBP represent promising therapeutic targets, particularly for interventions aimed at modulating immune responses or restoring synaptic function. By focusing on these conserved pathways, researchers can develop strategies to address the shared aspects of AD pathology across different contexts.

Looking ahead, several avenues for future research emerge. Expanding the analysis to include additional AD mouse models, such as APP/PS1 and TauP301S, would provide a broader understanding of molecular diversity and shared features. Exploring gene expression patterns in non-brain tissues, such as the spinal cord, liver, or blood, could uncover systemic effects and peripheral biomarkers of AD. Integration with transcriptomic data from space-flown mice offers a unique opportunity to identify shared features between Alzheimer’s models and space-associated neurodegeneration, enhancing our understanding of cross-species and cross-environmental influences.

Temporal analysis of gene expression changes over the disease course would capture dynamic regulatory events, offering insights into early versus late-stage molecular mechanisms. Experimental validation of key findings, such as perturbing TNF or TYROBP, could elucidate their roles in neuroinflammation and synaptic dysfunction. Additionally, incorporating proteomics and metabolomics data would complement transcriptomic analyses, providing a comprehensive view of molecular changes in AD.


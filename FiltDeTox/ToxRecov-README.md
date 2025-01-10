# **FiltDeTox Module**

## **Description**

## **Description**

The **FiltDeTox** module is the final step in the FiltDeTox pipeline. Using the R script `ToxRecov.R`, it integrates outputs from the **HHMER module** (new in this version, identifying ORFs mapped to toxin families) and the **ToxinKeyMatch** module. This enhanced version introduces refined filtering rules to classify sequences more accurately, incorporating outputs from the HHMER-based analysis and additional criteria such as **`Rating`**, **`toxin_keywords`**, **`pfam_ToxinKeywords`**, and **`ToxinFamily`**. 

The module categorizes sequences into toxin candidates, unlikely toxins, non-toxins, and Secreted Cysteine-Rich sequences Without Annotation (SCRs-WA). It also generates comprehensive visualizations (dendrograms, dot plots, pie charts, and domain-wise statistics) and FASTA files for both mature and precursor sequences, supporting downstream bioinformatics and functional analyses.

---

## **Dependencies**

Before running the **FiltDeTox** R script, ensure the following R packages are installed. You can install them using the commands below.

### **List of R Packages**:
1. **dplyr**: For data manipulation tasks such as filtering, grouping, and summarizing the data.
2. **tidyr**: For reshaping and tidying data, such as splitting or unnesting columns.
3. **ggplot2**: For creating visualizations, including the dot plot and pie charts.
4. **stringdist**: To calculate distances between ORFs or Gene_IDs based on binary or Levenshtein distances.
5. **ggtree**: To visualize hierarchical trees (dendrograms) generated from clustering.
6. **ape**: For working with phylogenetic trees and handling hierarchical clustering results.
7. **cowplot**: For combining multiple plots (the dendrogram and dot plot) into one final layout.
8. **RColorBrewer**: For defining custom color palettes used in the plots.

### **Installing Required Packages**

To install the required packages, open an R session and run the following commands:

```
# r
# Install packages from CRAN
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("stringdist")
install.packages("ape")
install.packages("cowplot")
install.packages("RColorBrewer")

# Install Bioconductor manager to install ggtree:
install.packages("BiocManager")
BiocManager::install("ggtree")

# Alternatively, install ggtree using devtools (optional):
install.packages("devtools")
devtools::install_github("YuLab-SMU/ggtree")
```
# How to Run the FiltDeTox Module

## Pre-requisites:

### Input Files:
Ensure that the following input files are available in the FiltDeTox directory:

-	`../ToxinKeyMatch/combined_output_keywords.tsv`: **Output from the ToxinKeyMatch module.**
-	`../ToxinKeyMatch/ToxProt_domain_Keywords.tsv`: **Pfam domain keywords associated with toxins.**
-	`../hhmer_tx_VenomZone/hhmer_Tx_orf_mapping.csv`: **Output from the HHMER-based module containing the summary of ORFs mapped to toxin families.**
-	`../DeTox_output_Ss_candidate_toxins.fasta`:  **DeTox output (example data).**

### R Environment:
Make sure R or RStudio is installed and the required R packages have been installed (as described above).

---

## Steps to Run the Module:

1. **Open RStudio (or any R environment)** and set the working directory to the FiltDeTox folder:

```
# r
setwd("/path/to/FiltDeTox/")
```

2. Load the Required Libraries: Ensure all necessary libraries are loaded before running the script:

```
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringdist)
library(ggtree)
library(ape)
library(cowplot)
library(RColorBrewer)
```
## Run the FiltDeTox Classification and Filtering:

Open the ´ToxRecov.R´ script and run the code to classify and filter toxin candidates. The script generates several output files, including plots, statistics and the corresponding FASTA files resulting from the filtering and classification process:

```
FiltDeTox/            # Final classification and filtering of toxin candidates
│   ├── Dendrogram_and_DotPlot_ORFs_TPM_by_ORF.pdf     # Plot of toxin candidates
│   ├── Dendrogram_and_DotPlot_ORFs_TPM_by_ORF.png     # High-resolution plot image
│   ├── FiltDeTox_Stats.tsv                 # Summary statistics of filtered sequences
│   ├── Full_Classified_Data.tsv            # Complete classification of sequences
│   ├── NestedPie_FiltDeTox_Flag_Sorted.png # Sorted pie chart visualization
│   ├── NestedPie_FiltDeTox.png             # Summary pie chart
│   ├── Non_Toxins_mature.fasta             # FASTA file of non-toxin mature peptides
│   ├── Non_Toxins_precursor.fasta          # FASTA file of non-toxin precursors
│   ├── Non_Toxins.tsv                      # Classification of non-toxin sequences
│   ├── Pfam_Domain_Summary_with_ORFs_Genes.tsv # Pfam domain and ORF summary
│   ├── SCRs_WA.fasta                     # Secreted cysteine-rich “mature” sequences (SCRs-WA)
│   ├── SCRs_WA.tsv                         # Detailed information on SCRs-WA sequences
│   ├── SCRs-WA_precursor.fasta             # Precursor sequences of SCRs-WA
│   ├── Toxins_Candidate_Rating_PieChart.png # Pie chart for toxin candidate ratings
│   ├── Toxins_Candidates_mature.fasta    # FASTA file of toxin candidates (mature peptides)
│   ├── Toxins_Candidates_precursor.fasta   # FASTA file of toxin candidates (precursors)
│   ├── Toxins_Candidates.tsv               # High-confidence toxin candidates
│   ├── ToxRecov.R                          # R script for classification and filtering
│   ├── ToxRecov-README.md                  # Instructions and details for running ToxRecov.R
│   ├── Unlikely_Toxins_mature.fasta        # FASTA of unlikely toxins (mature peptides)
│   ├── Unlikely_Toxins_precursor.fasta     # FASTA of unlikely toxins (precursors)
│   ├── Unlikely_Toxins.tsv                 # Classification of unlikely toxins
```


     # Step-by-Step Filtering Logic

The filtering process systematically classifies ORFs (Open Reading Frames) into four categories:

1. **Toxins-Candidates**
2. **Unlikely-Toxins**
3. **SCRs-WA**
4. **Non-Toxins**

The classification is based on several columns, including:
- **`wolfpsort_prediction`**: Indicates secretion predictions (e.g., "extr" for extracellular).
- **`ToxinFamily`**: Matches identified by HMMER.
- **`Rating`**: Flags specific sequence properties (e.g., `*`, `BD`, `SBD`).
- **`ORF_precursor_length`**: Length of the precursor sequence.
- **`toxin_keywords`** and **`pfam_ToxinKeywords`**: Boolean indicators for toxin-related keywords or domains.
- **`hit_descr`**: Descriptions of hits to proteins.

---

## **Rules**

### 1. Non-Toxins Based on Secretion Prediction
- **Condition**: If `wolfpsort_prediction` does not contain "extr" or "E.R.", classify as **Non-Toxins**.

---

### 2. ToxinFamily and Precursor Length
- **Condition**: If `ToxinFamily` is identified:
  - If `ORF_precursor_length > 500`, classify as **Unlikely-Toxins**.
  - If `Rating` contains `*`, classify as **Unlikely-Toxins**.
  - Otherwise, classify as **Toxins-Candidates**.

---

### 3. Strong Toxin Candidates Based on Flags
- **Condition**: If `Rating` contains `SBCD`, `SBCDT`, `SBC`, or `SBCT`:
  - If `ORF_precursor_length > 500`, classify as **Unlikely-Toxins**.
  - Otherwise, classify as **Toxins-Candidates**.

---

### 4. Sequences Matching Both Toxin Keywords and Pfam Domains
- **Condition**: If both `toxin_keywords` and `pfam_ToxinKeywords` are `TRUE`:
  - If `ORF_precursor_length > 500`, classify as **Non-Toxins**.
  - Otherwise, classify as **Unlikely-Toxins**.

---

### 5. Specific Flags (BD, B, SBD)
- **Condition**: If `Rating` contains `*BD`, `*B`, `BD`, `B`, or `SBD`:
  - If `ORF_precursor_length > 500`, classify as **Non-Toxins**.
  - Otherwise, classify as **Unlikely-Toxins**.

---

### 6. SCRs-WA Sequences
- **Condition**: If `Rating` is `SC` and `hit_descr` contains "uncharacterized" or is empty:
  - Classify as **SCRs-WA**.
- Otherwise, classify as **Non-Toxins**.

---

### 7. Default Classification
- **Condition**: If no conditions are met, classify as **Non-Toxins**.

---

## **Summary of Categories**

1. **Toxins-Candidates**:
   - Valid toxin family matches with appropriate lengths.

2. **Unlikely-Toxins**:
   - Sequences with `*` flags or long precursors.

3. **SCRs-WA**:
   - Secreted cysteine-rich, uncharacterized sequences.

4. **Non-Toxins**:
   - Default classification for sequences without sufficient evidence.

---

### **Happy Toxin Identification!**

# **FiltDeTox Module**

## **Description**

The **FiltDeTox** module is the third step in the FiltDeTox pipeline. It takes the output from the **ToxinKeyMatch** module and further refines the classification of sequences based on their **`Rating`**, **`toxin_keywords`**, **`pfam_ToxinKeywords`**, and other criteria. The module applies a sophisticated classification process using R, which generates various output files for high-confidence toxin candidates, unlikely toxins, non-toxins, and secreted cysteine-rich proteins.

Additionally, this module generates several visualizations, including dendrograms, dot plots, pie charts, and detailed domain-wise statistics.

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

- `../ToxinKeyMatch/combined_output_keywords.tsv`: This is the output from the ToxinKeyMatch module.
- `../ToxinKeyMatch/ToxProt_domain_Keywords.tsv`: Contains Pfam domain keywords associated with toxins.
- `SCRs_WA.tsv`: This file contains secreted cysteine-rich proteins for further classification.

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

      Open the ToxRecov.R script and run the code to classify and filter the toxin candidates.
      The script will generate several output files, including:

          - `Toxins_Candidates.tsv`: Contains high-confidence toxin candidates.
          - `Unlikely_Toxins.tsv`: Contains sequences classified as unlikely toxins.
          - `SCRs_WA.tsv`: Contains secreted cysteine-rich sequences.
          - `Non_Toxins.tsv`: Contains sequences classified as non-toxins.
          - `SCRs-WA.fasta`: Secreted cysteine-rich “mature” sequences (SCRs-WA).
          - `SCRs-WA_precursor.fasta`: Precursor sequences of SCRs-WA in FASTA format.

      ### Generate Plots and Statistics: The script also generates various visualizations:

      - `NestedPie_FiltDeTox.png`: A nested pie chart showing the classification distribution.
      - `NestedPie_FiltDeTox_Flag_Sorted.png`: A sorted flag distribution chart.
      - `Toxins_Candidate_Rating_PieChart.png`: A pie chart visualizing toxin candidate ratings.
      - `Dendrogram_and_DotPlot_ORFs_TPM_by_ORF.png`: A dendrogram and dot plot of toxin candidates based on Pfam domains and TPM.

      5. **Run the Visualizations:** You can manually run the code sections for generating the visualizations after the classification is complete.



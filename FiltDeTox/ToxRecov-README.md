# **FiltDeTox Module**

## **Description**

The **FiltDeTox** module is the final step in the FiltDeTox pipeline. Using the R script `ToxRecov.R`, it integrates outputs from the **HHMER module** (identifying ORFs mapped to toxin families) and the **ToxinKeyMatch** module. This enhanced version introduces refined filtering rules to classify sequences more accurately, incorporating outputs from the HHMER-based analysis and additional criteria such as **`Rating`**, **`toxin_keywords`**, **`pfam_ToxinKeywords`**, and **`ToxinFamily`**.

The module categorizes sequences into toxin candidates, unlikely toxins, non-toxins, and Secreted Cysteine-Rich sequences Without Annotation (SCRs-WA). It also generates:

- summary tables and statistics,  
- stacked barplots for FiltDeTox classes and DeTox flags,  
- nested donut pie charts (inner ring = FiltDeTox classification; outer ring = DeTox flags),  
- a dendrogram and dot plot summarising Pfam domains,  
- FASTA files for both mature and precursor sequences for all four categories, and  
- additional outputs for **full-length (type:complete)** precursor sequences per category,

supporting downstream bioinformatics and functional analyses.

---

## **Dependencies**

Before running the **FiltDeTox** R script, ensure the following R packages are installed.  

### **List of R Packages**

- **dplyr**: data manipulation (filtering, grouping, joins).  
- **tidyr**: reshaping and tidying tables (e.g. `separate_rows`).  
- **ggplot2**: visualisation (barplots, pie charts, dot plots).  
- **ape**: clustering and phylogenetic tree handling.  
- **cowplot**: combining multiple ggplot figures into a single layout.  
- **RColorBrewer**: qualitative colour palettes for plots.  
- **grid**: layout utilities (e.g. `unit()` for legend sizing).  
- **ggnewscale**: multiple independent colour/fill scales in a single figure (used for nested pies).  
- **ggtree**: visualisation of dendrograms as tree objects.  
- **Biostrings**: reading/writing FASTA files and handling amino acid sequences (full-length detection and export).

### **Installing Required Packages**

To install the required packages, open an R session and run the following commands:

```r
# Install packages from CRAN
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("ape")
install.packages("cowplot")
install.packages("RColorBrewer")
install.packages("grid")
install.packages("ggnewscale")

# Install Bioconductor manager and required Bioconductor packages
install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("Biostrings")

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
library(ggtree)
library(ape)
library(cowplot)
library(RColorBrewer)
library(grid)
library(ggnewscale)
library(Biostrings)
```
## Run the FiltDeTox Classification and Filtering:

Open the ´ToxRecov.R´ script and run the code to classify and filter toxin candidates. The script generates several output files, including plots, statistics and the corresponding FASTA files resulting from the filtering and classification process:

```
FiltDeTox/                     # Final classification and filtering of toxin candidates
│   ├── StackedBar_FiltDeTox_Classification_Sorted.png        # Stacked bar: FiltDeTox classes (PNG)
│   ├── StackedBar_FiltDeTox_Classification_Sorted.pdf        # Stacked bar: FiltDeTox classes (PDF)
│   ├── StackedBar_DeTox_Flags_Sorted.png                     # Stacked bar: DeTox flags (PNG)
│   ├── StackedBar_DeTox_Flags_Sorted.pdf                     # Stacked bar: DeTox flags (PDF)
│   ├── NestedPie_FiltDeTox.png                               # Nested pie (unsorted flags, PNG)
│   ├── NestedPie_FiltDeTox.pdf                               # Nested pie (unsorted flags, PDF)
│   ├── NestedPie_FiltDeTox_Flag_Sorted.png                   # Nested pie (flags sorted by abundance, PNG)
│   ├── NestedPie_FiltDeTox_Flag_Sorted.pdf                   # Nested pie (flags sorted by abundance, PDF)
│   ├── Toxins_Candidate_Rating_PieChart.png                  # Rating composition within Toxins-Candidates (PNG)
│   ├── Toxins_Candidate_Rating_PieChart.pdf                  # Rating composition within Toxins-Candidates (PDF)
│   ├── Dendrogram_and_DotPlot_ORFs_TPM_by_ORF.png            # Pfam domain clustering and TPM dot plot (PNG)
│   ├── Dendrogram_and_DotPlot_ORFs_TPM_by_ORF.pdf            # Pfam domain clustering and TPM dot plot (PDF)
│   ├── FiltDeTox_Stats.tsv                                   # Summary statistics of counts per FiltDeTox class
│   ├── Full_Classified_Data.tsv                              # Complete table with FiltDeTox_Classification and all fields
│   ├── Pfam_Domain_Summary_with_ORFs_Genes.tsv               # Pfam domain summary (ORFs, genes, TPM)
│   ├── Non_Toxins.tsv                                        # Classification table: Non-Toxins
│   ├── Non_Toxins_mature.fasta                               # FASTA: Non-Toxins, mature peptides
│   ├── Non_Toxins_precursor.fasta                            # FASTA: Non-Toxins, precursor sequences
│   ├── Toxins_Candidates.tsv                                 # Classification table: Toxins-Candidates
│   ├── Toxins_Candidates_mature.fasta                        # FASTA: Toxins-Candidates, mature peptides
│   ├── Toxins_Candidates_precursor.fasta                     # FASTA: Toxins-Candidates, precursor sequences
│   ├── Unlikely_Toxins.tsv                                   # Classification table: Unlikely-Toxins
│   ├── Unlikely_Toxins_mature.fasta                          # FASTA: Unlikely-Toxins, mature peptides
│   ├── Unlikely_Toxins_precursor.fasta                       # FASTA: Unlikely-Toxins, precursor sequences
│   ├── SCRs_WA.tsv                                           # Classification table: SCRs-WA sequences
│   ├── SCRs_WA_mature.fasta                                  # FASTA: SCRs-WA, mature peptides
│   ├── SCRs_WA_precursor.fasta                               # FASTA: SCRs-WA, precursor sequences
│   ├── Toxins_Candidates_full-length.tsv                     # Toxins-Candidates with full.length tag (complete/partial)
│   ├── Toxins_Candidates_precursor_full-length-seqs.fasta    # Full-length precursor sequences: Toxins-Candidates
│   ├── Toxins_Candidates_full-length_stats.txt               # Summary of full-length counts: Toxins-Candidates
│   ├── Unlikely_Toxins_full-length.tsv                       # Unlikely-Toxins with full.length tag
│   ├── Unlikely_Toxins_precursor_full-length-seqs.fasta      # Full-length precursor sequences: Unlikely-Toxins
│   ├── Unlikely_Toxins_full-length_stats.txt                 # Summary of full-length counts: Unlikely-Toxins
│   ├── SCRs_WA_full-length.tsv                               # SCRs-WA with full.length tag
│   ├── SCRs_WA_precursor_full-length-seqs.fasta              # Full-length precursor sequences: SCRs-WA
│   ├── SCRs_WA_full-length_stats.txt                         # Summary of full-length counts: SCRs-WA
│   ├── Non_Toxins_full-length.tsv                            # Non-Toxins with full.length tag
│   ├── Non_Toxins_precursor_full-length-seqs.fasta           # Full-length precursor sequences: Non-Toxins
│   ├── Non_Toxins_full-length_stats.txt                      # Summary of full-length counts: Non-Toxins
│   ├── ToxRecov.R                                            # R script for classification, summaries, plots and exports
│   ├── ToxRecov-README.md                                    # Instructions and details for running ToxRecov.R
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

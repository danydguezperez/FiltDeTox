# FiltDeTox: A Tool to Enhance and Filter Animal Toxins Identification from RNAseq Analyses

## Tool Description

**FiltDeTox** is an advanced bioinformatics pipeline designed to enhance the identification and classification of toxin-related sequences from RNAseq data. The tool integrates outputs from **TransDecoder**, **DeTox**, **BLASTp** (e.g., against UniRef90), and **Pfam domain annotations** to provide a streamlined and precise method for identifying animal toxins.

The pipeline uses a toxin keyword classification system, applying filters to minimize false positives and refine the dataset into high-confidence toxin candidates. This helps researchers focus on the most relevant sequences while discarding likely non-toxins. The pipeline is modular, meaning each component can be run independently or as part of the complete workflow.

## Modules Overview

### 1. TransDeTox
This module merges the outputs from TransDecoder BLASTp and DeTox, identifying candidate toxin sequences by matching sequence identifiers. By combining BLASTp annotations and DeTox predictions, it provides a refined dataset of sequences with potential toxin-related properties. It uses multiprocessing for efficient handling of large datasets and generates a combined output of BLASTp and DeTox data.

### 2. ToxinKeyMatch
This module classifies sequences based on toxin-related Pfam domains and curated toxin-related keywords. It applies keyword matching and domain filtering based on curated lists of metazoan toxins and venom proteins (filtered from UniProt), allowing for more precise identification of potential toxin sequences.

### 3. FiltDeTox
The final module further refines the dataset by applying an additional filtering step based on binary classifications generated in previous modules. It retains only high-confidence toxin candidates, eliminating sequences classified as non-toxins or unlikely toxins based on ToxinKeyMatch and DeTox flags.

## Execution and Customization

- **Independent Module Execution:**
  Each module can be executed independently, allowing users to optimize specific stages of the workflow or customize individual parts to fit their needs.
  
- **Pipeline Automation:**
  A shell script (`FiltDeTox.sh`) can automate the entire process, executing all modules in sequence and generating the final output.
  
- **Seamless Updates:**
  FiltDeTox allows continuous improvement with easy updates to toxin-related keywords, Pfam domains, and filtering criteria.

## Future Improvements

FiltDeTox is designed for continuous improvement, with Modules 2 (ToxinKeyMatch) and Module 3 (FiltDeTox) being adaptable. Users can update the keyword lists, Pfam domains, and filtering logic as new toxin-related data becomes available, making the tool more precise and versatile over time.

## Pipeline Structure and Enclosed Files


```plaintext
FiltDeTox/                # Main directory containing the entire pipeline
├── TransDeTox/           # First module: Processes and merges BLASTp and DeTox results.
│   ├── blastp.outfmt6.w_pct_hit_length     # BLASTp output file with hit length percentage.
│   ├── combined_output.tsv                 # Merged output of TransDecoder BLASTp and DeTox results.
│   ├── DeTox_output_Ss_SE_toxins.tsv       # DeTox output from SE assembly of S. savaglia (example data).
│   ├── matched_content.tsv                 # Intermediate file matching BLASTp and DeTox outputs.
│   ├── TransDeTox.py                       # Python script for processing and merging results from BLASTp and DeTox.
│   └── TransDeTox-README.md                # README file with instructions on the TransDeTox module.
│
├── ToxinKeyMatch/        # This module enhances toxin identification through keywords matching.
│   ├── combined_output_keywords.tsv        # Output file after keyword and domain matching.
│   ├── ToxinKeyMatch.py                    # Python script for keyword and Pfam domain matching.
│   ├── ToxinKeyMatch-README.md             # README file with instructions for the ToxinKeyMatch module.
│   ├── toxins_keywords.csv                 # List of toxin-related keywords for matching.
│   ├── ToxProt_domain_Keywords.tsv         # Pfam domain keywords for identifying toxins.
│   └── Extract_ToxProtDomain/              # Subfolder for extracting Pfam domains.
│       ├── Extract_ToxProtDomain.md        # Documentation for the Pfam domain extraction process.
│       ├── Extract_ToxProtDomain.py        # Python script for extracting Pfam domains.
│       └── uniprotkb_taxonomy_id_33208_AND_cc_tiss_2024_04_18_pfam.xlsx # data from ToxProt.
│
├── FiltDeTox/            # Final classification and filtering of toxin candidates.
│   ├── Dendrogram_and_DotPlot_ORFs_TPM_by_ORF.pdf     # Plot of toxin candidates based on ORFs and TPM.
│   ├── FiltDeTox_Stats.tsv                             # Summary statistics of filtered sequences.
│   ├── Full_Classified_Data.tsv                        # Complete classification of sequences after filtering.
│   ├── NestedPie_FiltDeTox.pdf                         # Nested pie chart showing classifications.
│   ├── Non_Toxins.tsv                                  # Sequences classified as non-toxins.
│   ├── Pfam_Domain_Summary_with_ORFs_Genes.tsv         # Pfam domain summary with ORF and gene details.
│   ├── SCRs-WA.fasta                                   # Secreted cysteine-rich “mature” sequences (SCRs-WA).
│   ├── SCRs_WA.tsv                                     # Detailed information on SCRs-WA sequences.
│   ├── SCRs-WA_precursor.fasta                         # Precursor sequences of SCRs-WA in FASTA format.
│   ├── Toxins_Candidate_Rating_PieChart.pdf            # Pie chart of ratings for toxin candidates.
│   ├── Toxins_Candidates.tsv                           # High-confidence toxin candidates.
│   ├── ToxRecov.R                                      # R script for final classification and filtering steps.
│   ├── ToxRecov-README.md                              # Instructions and details for running ToxRecov.R.
│   ├── Unlikely_Toxins.tsv                             # Sequences classified as unlikely toxins.
│   ├── ToxProt_domain_Keywords_Stats.tsv               # Summary statistics of ToxProt domain matches.
│
├── FiltDeTox.sh         # Shell script to run the entire FiltDeTox pipeline.
```

## Automatic Execution Using FiltDeTox.sh Pipeline
### Conditions (Mandatory):

- **Input Files**:
        DeTox_output_toxins.tsv: DeTox output file (e.g., DeTox_output_Ss_SE_toxins.tsv)
        blastp.outfmt6.w_pct_hit_length: BLASTp output file from full-length transcript analysis
- **Files for Keyword Matching:**
        ToxProt_domain_Keywords.tsv: This file is generated during the process. To update it, download a new version and place it in the appropriate folder.
        toxins_keywords.csv: A curated list of toxin-related keywords. Users can customize and update the content as necessary.

> Note: Since generating or updating keyword lists is a specialized task, it is recommended to run
> this step separately. The toxins keyword list can be reviewed and edited by experts according to
> specific research needs.

## Dependencies
### Python-based Scripts:

- **Python 3:** Required for running the TransDeTox.py and ToxinKeyMatch.py scripts.
- **pandas:** For efficient data manipulation and processing.
- **multiprocessing:** To handle large datasets efficiently.

Install the required dependencies with:
`bash`
`pip install pandas`

## R-based Scripts (ToxRecov.R):
To run the FiltDeTox module using R, ensure **R** or **RStudio** is installed. Install the necessary R packages with the following commands:

```
# Install packages from CRAN
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("stringdist")
install.packages("ape")
install.packages("cowplot")
install.packages("RcolorBrewer")

# Install Bioconductor manager to install ggtree:
install.packages("BiocManager")
BiocManager::install("ggtree")

```
# Alternatively, install ggtree using devtools (optional):
```
install.packages("devtools")
devtools::install_github("YuLab-SMU/ggtree")
```

Ensure the working directory is set correctly before running any of the modules, particularly for FiltDeTox:

`R`
`setwd("/path/to/FiltDeTox/")`

This ensures that the R script saves the output files in the correct directory.

## Initial Shell Script (FiltDeTox.sh)

The FiltDeTox.sh script is designed to link all three modules of the pipeline (TransDeTox, ToxinKeyMatch, and FiltDeTox), allowing users to run the entire process seamlessly. Once all dependencies are installed, users can execute this script to run the pipeline from start to finish.
Note: The integrated pipeline has only been tested on **Linux** systems.

## How to Run the Shell Script:

  1. **Clone the Repository:**
Download the FiltDeTox repository or copy all necessary files into
a directory:

`git clone https://github.com/danydguezperez/FiltDeTox.git`

  2. **Navigate to the Repository Folder**
Once the cloning is complete, navigate into the project folder using the cd command:

`cd FiltDeTox`

  3. **Set Permissions**:
Ensure the shell script is executable:

     ```
     #bash
       chmod +x FiltDeTox.sh
     ``` 
     
  5. **Run the Script**:
Run the full pipeline using the shell script:

     ```
     #bash
       ./FiltDeTox.sh
     ```

This script will automatically call each module in sequence:

  1. **TransDeTox.py**
  2. **ToxinKeyMatch.py**
  3. **ToxRecov.R**
     
Upon completion, the final output files, figures, and statistics will be saved in their respective folders. Users can check the output for toxin candidates, summary statistics, and generated plots.

>Note:
>All the scripts within the FiltDeTox pipeline can be run standalone. If you experience issues, particularly with the ToxRecov.R script, detailed >explanations and instructions are provided in the respective README.md files in each module directory.

Happy Toxin Identification!
#done!!!


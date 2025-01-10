# **hhmer_tx_VenomZone Module**

## **Description**
The `hhmer_tx_VenomZone` module is designed to identify and classify toxin-related sequences from pre-aligned toxin families using Hidden Markov Models (HMMs). It automates the HMMER workflow, building HMM profiles from pre-aligned sequences and searching them against a specified FASTA file. This module streamlines toxin identification, providing ORF-to-family mappings, sequence extractions, and statistical summaries for downstream analyses.

##  Module Structure

```
FiltDeTox/                # Main directory containing the entire pipeline
├── hhmer_tx_VenomZone/   # module for HMMER-based toxin family identification
│   ├── tx_VenomZone_aln/     # Pre-aligned toxin family sequences and script
│   │   ├── *_aln.fasta       # Pre-aligned sequences for each toxin family
│   │   ├── *_aln.hmm         # HMM profiles built using hmmbuild
│   │   ├── *_aln_hmmsearch.out  # HMMER search results
│   │   ├── *_aln_grep.txt    # ORFs matched for each toxin family
│   │   ├── *_aln_grep.fasta  # Extracted sequences for matched ORFs
│   │   ├── summary_hhmer.txt # Log file summarizing HMMER processes
│   │   └── run_hmmbuild_and_hmmsearch.sh  # Script to run HMMER steps
│   ├── hhmer_Tx_fam_hits.csv              # Summary of ORFs per toxin family
│   ├── hhmer_Tx_orf_mapping.csv           # Mapping of ORFs to toxin families
│   ├── hhmer_Tx_fam_hits_stats.txt        # Statistics on ORFs per toxin family
│   ├── hhmer_tx_VenomZone-README.md # README file for the HMMER-based module
│   └── hhmerTxMatch.py     # Python script for processing and summarizing HMMER results
…
├── DeTox_output_Ss_candidate_toxins.fasta     # DeTox output (example data)
├── FiltDeTox_v2.0.sh     # Shell script to run the entire FiltDeTox pipeline
```
---

## **Conditions and Requirements**

### **Input Files**
1. **Pre-aligned sequences**: Files located in the `tx_VenomZone_aln` folder with the pattern `*_aln.fasta`. These files represent pre-aligned sequences for each toxin family.
2. **FASTA file**: A candidate toxin file located two levels up from the script directory (e.g., `../DeTox_output_Ss_candidate_toxins.fasta`).

### **Pre-requisites**
- **HMMER tools**: `hmmbuild` and `hmmsearch` must be installed and available in the PATH.
- **`seqkit`**: Required for extracting matching sequences.
- **Python 3**: Required for processing and summarizing results using `hhmerTxMatch.py`.
- **pandas**: Install using:

`pip install pandas`

**Script Summary**:
This Bash script automates the **HMMER** workflow to analyze pre-aligned sequences of **cnidarian toxins** retrieved from the **VenomZone** database. These pre-aligned sequences, grouped by toxin family, are stored in the folder **tx_VenomZone_aln**. The script uses these alignments to build HMM profiles and searches them against a **`.fasta file`** located **two levels up** from the working directory, specifically in the root of the **FiltDeTox** directory. This script can be used to profile other animal toxins families from other taxon, by naming the file with containing the aligned FASTA sequences with this format `*_aln.fasta`, placed in the proper directory. 

**Main Steps the Script Performs**:
1.	**Setup Directories**:
-	Identifies and configures key directories:
-	**tx_VenomZone_aln** (working directory with pre-aligned toxin sequences).
-	The FiltDeTox root directory two levels up, containing the target **`.fasta file`** for hmmsearch.
2.	**Clean Previous Results**:
-	Deletes all files from the working directory **except** `.sh` scripts and pre-aligned `*_aln.fasta` files, ensuring a clean start.
3.	**Build HMM Profiles**:
-	For each pre-aligned **toxin family** in the **tx_VenomZone_aln** folder (files ending with **`.fasta`**), the script:
-	Builds an HMM profile using hmmbuild.
-	Logs progress in a summary file **summary_hhmer.txt**.
4.	**Run hmmsearch**:
-	Searches each HMM profile against the provided **`.fasta file`** in the root of the **`FiltDeTox`** directory.
-	Extracts matched ORFs and writes them to corresponding **`_grep.txt files`**.
-	Logs families with **"No matches found"** where no ORFs are detected.
5.	**Generate Results in Tabular Format**:
-	Concatenates all _grep.txt files into a single summary table named **hhmer_Tx_hits.csv**.
-	Columns represent toxin families, and rows contain matched ORFs (or **"No matches found"**).
6.	**Copy and Rename Results**:
-	Copies the **hhmer_Tx_hits.csv** file **one level up** from the script's directory.
-	Renames it to **hhmer_Tx_fam_hits.csv**.
-	Cleans up the column headers by removing the `_aln` suffix for a cleaner presentation.
7.	**Extract Sequences**:
-	Extracts matching sequences for each toxin family listed in **`_grep.txt` files** from the target **`.fasta file**`**.
-	Saves these sequences into individual **`_grep.fasta` files**.
8.	**Generate ORF Statistics (New Step)**:
-	Counts the number of matched ORFs per toxin family.
-	Saves the statistics in a report file named **hhmer_Tx_fam_hits_stats.txt**, formatted as:
  
```
Acrorhagin_I_fam: 5 ORFs
Actinoporin_fam_HALT_subfam: 3 ORFs
```

**Main Outputs**:
1.	**summary_hhmer.txt**: Log file summarizing hmmbuild and hmmsearch progress.
2.	**hhmer_Tx_hits.csv**: Tabular summary of matched ORFs per toxin family.
3.	**hhmer_Tx_fam_hits.csv**: A cleaner version of the results, placed one level up with simplified column headers.
4.	**hhmer_Tx_fam_hits_stats.txt**: A report of ORF counts per toxin family.
5.	**_grep.fasta** files: Extracted FASTA sequences of matched ORFs for each toxin family.

**Summary**:
The script efficiently automates the workflow to build HMM profiles, search toxin family sequences, summarize ORFs, and extract matched sequences for downstream analysis. It ensures clean outputs, concise logs, and comprehensive reporting, making it suitable for toxin family classification and analysis in cnidarian sequences.

### Script `hhmerTxMatch.py`

This Python script processes HMMER search results and produces a summarized mapping of ORFs (Open Reading Frames) to their corresponding toxin families. Additionally, it generates a final report with the total count of unique ORFs.

**Key Steps of the Script**:

1.	**Input Data**:
-	**hhmer_Tx_fam_hits.csv**: The input file containing ORFs grouped under different toxin family columns.
-	This file must be in the same directory as the script.

2.	**Processing ORF-to-Toxin Mapping**:
-	Reads the input file (**hhmer_Tx_fam_hits.csv**).
-	For each toxin family (column in the file), it extracts the ORFs.
-	Skips any empty values or invalid rows.
-	Combines families for duplicate ORFs:
-	If an ORF belongs to multiple families, the family names are joined by semicolons (;).

3.	**Output 1**:
-	Creates a new file **hhmer_Tx_orf_mapping.csv**.
-	This file contains two columns:
-	ORF: Unique ORFs from the input.
-	ToxinFamily: Corresponding toxin family (or families) for each ORF.

4.	**Appending ORF Statistics**:
-	Calculates the **total unique ORFs**:
-	Excludes rows with **"No matches found"**.
-	Ensures the header or invalid rows are skipped.
-	Appends this count as a final line to the **hhmer_Tx_fam_hits_stats.txt** file.
-	Example of the appended line:

```
Total Unique ORFs hhmer Tx hits: 15
```

**Outputs**:
1.	**hhmer_Tx_orf_mapping.csv**:
-	A clean file mapping ORFs to their toxin families.
2.	**hhmer_Tx_fam_hits_stats.txt**:
-	Contains counts of ORFs per toxin family (previously generated).
-	Adds a final line summarizing the **total unique ORFs** found.

**Execution Flow**:
1.	The script processes the input file **hhmer_Tx_fam_hits.csv**.
2.	It generates **hhmer_Tx_orf_mapping.csv** with ORFs mapped to toxin families.
3.	Appends the total unique ORF count to **hhmer_Tx_fam_hits_stats.txt**.
---
## Integration with FiltDeTox Pipeline
The outputs of this module feed into the next steps of the FiltDeTox pipeline:
-	**hhmer_Tx_orf_mapping.csv**: Used by the **FiltDeTox** module to merge with the outputs from other modules for classification and filtering.
This module ensures accurate and clean identification of toxin-related ORFs, paving the way for further analysis in the FiltDeTox pipeline.


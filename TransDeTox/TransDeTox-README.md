# **TransDeTox Module**

## **Description**

The **TransDeTox** module is the first step in the FiltDeTox pipeline. This module processes and merges the results from two sources:
1. **TransDecoder BLASTp output**: Results from BLASTp searches of translated RNAseq data.
2. **DeTox output**: Predictions from the DeTox tool, flagging potential toxin-related sequences.

The **TransDeTox.py** script matches sequence identifiers (contigs) between these two datasets and produces a combined output. This merged file serves as input for subsequent filtering and classification in the FiltDeTox pipeline.

---

## **How the Script Works**

### 1. **Input Files**

- **`blastp.outfmt6.w_pct_hit_length`**: The BLASTp output from TransDecoder. It contains alignment results of RNAseq data, detailing which sequences matched known proteins in a database like UniProt. The file should have tab-separated fields, including the query sequence (`qseqid`).
  
- **`DeTox_output_Ss_SE_toxins.tsv`**: The output from the DeTox tool, which identifies potential toxin-related sequences from RNAseq data. The file should include fields such as `contig` and `ID`, where `contig` refers to the sequence identifier.

### 2. **Workflow**

The script performs the following steps:

#### **Step 1: Clean the BLASTp Output**
- The script reads the BLASTp file (`blastp.outfmt6.w_pct_hit_length`). 
- If the first line in the file begins with a comment (`#`), it removes that comment. This prepares the file for further processing.
- The cleaned file is saved as `matched_content.tsv`.

#### **Step 2: Load Data into DataFrames**
- The script reads the cleaned BLASTp file (`matched_content.tsv`) and the DeTox file (`DeTox_output_Ss_SE_toxins.tsv`) into two `pandas` DataFrames (`blastp_df` and `detox_df`, respectively).
  
#### **Step 3: Parallelized Matching**
- Using the `multiprocessing` module, the script identifies matching rows between the BLASTp and DeTox outputs based on the `contig` identifier.
- For each row in the DeTox file, it searches for matching sequences in the BLASTp file by looking for sequences where `qseqid` in the BLASTp data contains the `contig` value from the DeTox data.
- This matching is done in parallel using all available CPU cores to handle large datasets efficiently.

#### **Step 4: Save Matched Data**
- The matched content is saved as a new DataFrame (`matched_df`), which is written to the file `matched_content.tsv`.
  
#### **Step 5: Combine DeTox and BLASTp Data**
- The script merges the DeTox file (`DeTox_output_Ss_SE_toxins.tsv`) with the matched BLASTp content using the `ID` column as the key.
- The final merged data is saved as `combined_output.tsv`.

---

## **Output Files**

The script generates two output files:

1. **`matched_content.tsv`**: An intermediate file containing the merged data between BLASTp and DeTox based on matching sequence identifiers (contigs).
2. **`combined_output.tsv`**: The final combined output containing all relevant fields from the DeTox and BLASTp outputs. This file is used as input for further filtering and classification in the next steps of the FiltDeTox pipeline.

---

## **Dependencies**

### **Python 3**:
This script requires Python 3 or a higher version to run.

### **Python Libraries**:
The following Python libraries are required to run the script:
- **pandas**: For data manipulation and reading/writing CSV/TSV files.
- **multiprocessing**: To parallelize the matching of BLASTp and DeTox data, speeding up the process.

You can install the necessary libraries by running the following commands:

```
# bash
pip3 install pandas
```
# How to Run the Script
## Pre-requisites:
1. **Input Files:** Ensure that the following input files are available in the same directory as the script:

- blastp.outfmt6.w_pct_hit_length
- DeTox_output_Ss_SE_toxins.tsv
2. **Python Environment:** Ensure that Python 3 and the required libraries are installed.

## Running the Script:
Open a terminal and navigate to the directory containing the script and the input files.
Run the script using the following command:


```
# bash
cd path/to/TransDeTox.py
python3 TransDeTox.py
```
### Output:
The script will generate the following output files in the same directory:

- `matched_content.tsv`
- `combined_output.tsv`

These files can be used as input for the next module in the FiltDeTox pipeline.


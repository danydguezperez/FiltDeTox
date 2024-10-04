# **Extract_ToxProtDomain.py**

### **General Description**
The `Extract_ToxProtDomain.py` script processes an input file from **UniProt**, focusing specifically on the **"Pfam"** domain annotations. It extracts toxin-related **Pfam domains**, which are then used to generate a keyword file for identifying toxins in RNAseq data analysis. Additionally, the script generates a statistical report on the extracted domains.

This script is an integral part of the **FiltDeTox** pipeline, providing the necessary Pfam domain keywords for further filtering in **Module 2: ToxinKeyMatch**.

---

### **Dependencies**
To run the script, the following dependencies are required:

- **Python 3**
- **pandas**: For reading and manipulating Excel files.
- **re**: For handling regular expressions, which are used to extract Pfam domains.

You can install the necessary Python libraries by running:

```
bash
pip install pandas
```

---

### **Input File**

1. **uniprotkb_taxonomy_id_33208_AND_cc_tiss_2024_04_18_pfam.xlsx** 
   - This Excel file contains UniProt data for metazoan toxins and venom-related proteins.
   - The file includes several columns, but the script primarily focuses on the **"Pfam"** column to extract domain annotations.
   - Columns included in the file: **Entry**, **Reviewed**, **Entry Name**, **Protein names**, **Gene Names**, **Organism**, **Length**, **Domain [FT]**, **Protein families**, and **Pfam**.
   - This file is located in the `Extract_ToxProt_domain` folder within the **ToxinKeyMatch** module.

---

### **Output Files**

1. **ToxProt_domain_Keywords.tsv** 
   - A TSV file containing **unique Pfam domain terms** extracted from the "Pfam" column.
   - Each row in the file represents a single unique term related to toxin or venom protein domains.
   - This file will be used in **Module 2: ToxinKeyMatch** for further toxin identification and filtering based on domain matches.

2. **ToxProt_domain_Keywords_Stats.tsv** 
   - A TSV file that provides **statistical information** about the extraction process, including:
     - Total number of lines analyzed from the input file.
     - Total number of Pfam domain terms extracted.
     - Number of unique terms.
     - Number of duplicate terms.

---

### **What the Script Does**

1. **Read the Excel File**:
   - The script reads the input Excel file using **pandas** and checks for the presence of the **"Pfam"** column to ensure the correct data is processed.

2. **Extract Middle Terms**:
   - Using a regular expression pattern, the script extracts the **middle terms** from the entries in the "Pfam" column. Each entry is expected to follow the pattern: `PFXXXXX; term; 1..`. 
   - The extracted terms represent domain annotations related to toxins or venom proteins.

3. **Calculate Statistics**:
   - The script calculates:
     - The total number of lines analyzed in the DataFrame.
     - The total number of terms extracted from the "Pfam" column.
     - The number of unique terms (terms without duplicates).
     - The number of duplicate terms (terms that appeared more than once).

4. **Save Extracted Terms and Statistics**:
   - The **unique Pfam domain terms** are saved in `ToxProt_domain_Keywords.tsv`.
   - The **statistical report** is saved in `ToxProt_domain_Keywords_Stats.tsv`.

---

### **Example Usage**

1. **Prepare the Input File**: 
   Ensure that the input file `uniprotkb_taxonomy_id_33208_AND_cc_tiss_2024_04_18_pfam.xlsx` is located in the same directory as the script or provide the correct file path.

2. **Run the Script**: 
   To execute the script, use the following command in the terminal:

   ```
   # bash
   cd path/to/Extract_ToxProtDomain.py
   python3 Extract_ToxProtDomain.py
   ```

3. **Generated Output**: 
   After running the script, the following output files will be created:
   - `ToxProt_domain_Keywords.tsv`: Contains the unique toxin-related Pfam domains.
   - `ToxProt_domain_Keywords_Stats.tsv`: Provides statistics on the extracted Pfam domains.

---

### **Conclusion**

The `Extract_ToxProtDomain.py` script efficiently extracts toxin-related Pfam domains and generates key statistics. The extracted domains are used in downstream steps, particularly in **Module 2: ToxinKeyMatch**, to enhance the identification and classification of potential toxins based on domain annotations.

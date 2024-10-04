import pandas as pd
import re

# Function to read Excel file and handle potential issues
def read_excel_file(file_path):
    try:
        df = pd.read_excel(file_path)
        print("File read successfully.")
        return df
    except Exception as e:
        print(f"Error reading the file: {e}")
        return None

# Path to the Excel file
file_path = "uniprotkb_taxonomy_id_33208_AND_cc_tiss_2024_04_18_pfam.xlsx"

# Read the Excel file into a DataFrame
df = read_excel_file(file_path)

if df is not None:
    # Print the column headers to verify the column name
    print("Column Headers:", df.columns)

    # Try to find the correct column name for "Pfam"
    column_name = None
    for col in df.columns:
        if col.strip().lower() == 'pfam':
            column_name = col
            break

    if column_name:
        print(f"Using column name: {column_name}")

        # Replace NaN values in the "Pfam" column with an empty string
        df[column_name] = df[column_name].fillna('')

        # Initialize a set to store unique terms and a list to store all extracted terms
        unique_terms = set()
        all_terms = []

        # Define a regular expression pattern to extract the middle term
        pattern = r'PF\d+; ([^;]+);'

        # Extract and collect the unique terms from the "Pfam" column
        for entry in df[column_name]:
            matches = re.findall(pattern, entry)
            all_terms.extend(matches)
            unique_terms.update(matches)

        # Write the unique terms to the output file
        with open("ToxProt_domain_Keywords.tsv", "w") as f:
            f.write("Unique Terms\n")
            for term in unique_terms:
                f.write(f"{term}\n")

        # Calculate statistics
        total_lines = len(df)
        total_terms_extracted = len(all_terms)
        number_of_unique_values = len(unique_terms)
        number_of_duplicates = total_terms_extracted - number_of_unique_values

        # Write the statistics to a new TSV file
        with open("ToxProt_domain_Keywords_Stats.tsv", "w") as f:
            f.write("Statistic\tValue\n")
            f.write(f"Total Lines Analyzed\t{total_lines}\n")
            f.write(f"Number of Terms Extracted\t{total_terms_extracted}\n")
            f.write(f"Number of Duplicates\t{number_of_duplicates}\n")
            f.write(f"Number of Unique Values\t{number_of_unique_values}\n")

        print("Statistics extracted and saved to 'ToxProt_domain_Keywords_Stats.tsv'.")
    else:
        print(f"Column 'Pfam' not found in the DataFrame.")
else:
    print("Failed to read the DataFrame.")


import pandas as pd
import re

def load_keywords(file_path):
    # Load keywords directly into a Pandas DataFrame and select the first column
    df = pd.read_csv(file_path, header=None)
    keywords_series = df.iloc[:, 0]  # Selecting the first column as a Series
    keywords = set(keywords_series.str.lower().str.strip())  # Convert to lowercase for case-insensitive matching
    return keywords

def check_keywords(row, keywords):
    # Combine all cell values into one large string and convert to lower case
    combined_text = ' '.join(row.astype(str)).lower()
    # Check if any keyword is in the combined text as a whole word
    return any(re.search(r'\b' + re.escape(keyword) + r'\b', combined_text) for keyword in keywords)

def check_pfam_domains(row, keywords):
    # Check if any keyword is in the 'pfam domains' column as a whole word
    pfam_text = str(row['pfam domains']).lower()
    return any(re.search(r'\b' + re.escape(keyword) + r'\b', pfam_text) for keyword in keywords)

def add_keywords_column(data_path, keywords_path_1, keywords_path_2, output_path):
    # Load keywords from both files
    keywords_1 = load_keywords(keywords_path_1)
    keywords_2 = load_keywords(keywords_path_2)

    data = pd.read_csv(data_path, sep='\t', header=0)

    # Apply the check_keywords function to each row for the first set of keywords
    data['toxin_keywords'] = data.apply(lambda row: check_keywords(row, keywords_1), axis=1)
    data['toxin_keywords'] = data['toxin_keywords'].map({True: 'TRUE', False: 'FALSE'})

    # Apply the check_pfam_domains function to each row for the second set of keywords
    data['pfam_ToxinKeywords'] = data.apply(lambda row: check_pfam_domains(row, keywords_2), axis=1)
    data['pfam_ToxinKeywords'] = data['pfam_ToxinKeywords'].map({True: 'TRUE', False: 'FALSE'})

    # Save the updated DataFrame
    data.to_csv(output_path, sep='\t', index=False)
    print(f"Updated data saved to '{output_path}'.")

# Usage example, adjust the file paths as needed
if __name__ == "__main__":
    combined_output_path = 'TransDeTox/combined_output.tsv'
    keywords_path_1 = 'ToxinKeyMatch/toxins_keywords.csv'  # First set of keywords
    keywords_path_2 = 'ToxinKeyMatch/ToxProt_domain_Keywords.tsv'  # Second set of keywords
    output_path = 'ToxinKeyMatch/combined_output_keywords.tsv'

    add_keywords_column(combined_output_path, keywords_path_1, keywords_path_2, output_path)


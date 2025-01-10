import pandas as pd
import os

def process_hhmer_file(hhmer_file, output_file):
    """Process HMMER results to create ORF-to-toxin-family mapping."""
    # Load the HMMER data
    hhmer_data = pd.read_csv(hhmer_file, sep='\t', header=0)

    # Prepare lists for ORFs and their corresponding families
    orf_list = []
    family_list = []

    for family in hhmer_data.columns:
        for orf in hhmer_data[family].dropna():  # Skip NaN values
            if orf.strip():  # Ensure it's not an empty string
                orf_list.append(orf.strip())
                family_list.append(family)  # Keep the original family name without modifying it

    # Create a DataFrame from the lists
    orf_family_df = pd.DataFrame({'ORF': orf_list, 'ToxinFamily': family_list})

    # Combine families for duplicate ORFs
    combined_orf_family_df = (
        orf_family_df.groupby('ORF')['ToxinFamily']
        .apply(lambda families: '; '.join(sorted(set(families))))
        .reset_index()
    )

    # Save to a CSV file
    combined_orf_family_df.to_csv(output_file, index=False)
    print(f"Processed ORF-to-toxin-family mapping saved to '{output_file}'.")

def append_unique_orf_count(orf_mapping_file, stats_file):
    """Append the total unique ORF count to the stats file."""
    # Load ORF mapping file
    orf_mapping_df = pd.read_csv(orf_mapping_file)

    # Exclude rows with 'No matches found' and count unique ORFs
    filtered_orfs = orf_mapping_df[
        (orf_mapping_df['ORF'] != 'No matches found') & 
        (orf_mapping_df['ORF'].str.strip() != 'ORF')  # Exclude header or invalid rows
    ]['ORF'].unique()

    total_unique_orfs = len(filtered_orfs)

    # Append the total count to the stats file
    with open(stats_file, 'a') as f:
        f.write(f"Total Unique ORFs hhmer Tx hits: {total_unique_orfs}\n")
    print(f"Total Unique ORFs count appended to '{stats_file}'.")

# Usage example
if __name__ == "__main__":
    # Define paths relative to the script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # File paths
    hhmer_file = os.path.join(script_dir, 'hhmer_Tx_fam_hits.csv')
    output_file = os.path.join(script_dir, 'hhmer_Tx_orf_mapping.csv')
    stats_file = os.path.join(script_dir, 'hhmer_Tx_fam_hits_stats.txt')

    # Run the processing script
    process_hhmer_file(hhmer_file, output_file)

    # Append total unique ORF count to stats file
    append_unique_orf_count(output_file, stats_file)


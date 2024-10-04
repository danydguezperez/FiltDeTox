import pandas as pd
from multiprocessing import Pool, cpu_count

def process_row(detox_row, blastp_df):
    contig_value = detox_row['contig']
    mask = blastp_df['qseqid'].str.contains(contig_value, na=False)
    matching_blastp_rows = blastp_df[mask]
    return [pd.concat([blastp_row, detox_row], axis=0) for _, blastp_row in matching_blastp_rows.iterrows()]

if __name__ == "__main__":
    # Use relative paths directly
    TransDecoder_output_path = 'TransDeTox/blastp.outfmt6.w_pct_hit_length'
    DeTox_output_path = 'TransDeTox/DeTox_output_Ss_SE_toxins.tsv'
    matched_content_path = 'TransDeTox/matched_content.tsv'  # Save in the current directory
    
    with open(TransDecoder_output_path, 'r') as original_file:
        lines = original_file.readlines()
    if lines[0].startswith('#'):
        lines[0] = lines[0][1:]
    with open(matched_content_path, 'w') as modified_file:
        modified_file.writelines(lines)

    print(f"TransDecoder BLASTp output file has been processed and saved as '{matched_content_path}'.")

    blastp_df = pd.read_csv(matched_content_path, sep='\t', header=0)
    detox_df = pd.read_csv(DeTox_output_path, sep='\t', header=0)
    
    with Pool(processes=cpu_count()) as pool:
        results = pool.starmap(process_row, [(row, blastp_df) for index, row in detox_df.iterrows()])
    
    matched_df = pd.DataFrame([item for sublist in results for item in sublist]).reset_index(drop=True)
    matched_df.to_csv(matched_content_path, sep='\t', index=False)
    print(f"Matched data saved to '{matched_content_path}'.")

    detox_df = pd.read_csv(DeTox_output_path, sep='\t', header=0)
    matched_content_df = pd.read_csv(matched_content_path, sep='\t', header=0)

    id_index = matched_content_df.columns.get_loc('ID')
    filtered_matched_content_df = matched_content_df.iloc[:, :id_index+1]

    combined_df = pd.merge(detox_df, filtered_matched_content_df, on='ID', how='left')
    
    combined_output_path = 'TransDeTox/combined_output.tsv'  # Save in the current directory
    combined_df.to_csv(combined_output_path, sep='\t', index=False)
    print(f"Combined data saved to '{combined_output_path}'.")


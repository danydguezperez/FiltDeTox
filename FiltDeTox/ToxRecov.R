.libPaths("~/R/x86_64-pc-linux-gnu-library/4.1")  # Adjust to the correct R version.
# Set the working directory to FiltDeTox
setwd("FiltDeTox")
# Dependencies required for running the script and their purpose:
# 1. dplyr: For data manipulation tasks such as filtering, grouping, and summarizing the data.
# 2. tidyr: For reshaping and tidying data, such as splitting or unnesting columns.
# 3. ggplot2: For creating visualizations, including the dot plot used in the final step.
# 4. stringdist: To calculate distances between ORFs or Gene_IDs based on binary or Levenshtein distances.
# 5. ggtree: To visualize hierarchical trees (dendrograms) generated from clustering.
# 6. ape: For working with phylogenetic trees and handling hierarchical clustering results.
# 7. cowplot: For combining multiple plots (the dendrogram and dot plot) into one final layout.
# 8. RColorBrewer: For defining custom color palettes used in the plots.

# To install the required packages, use the following commands:

# Install packages from CRAN () remove the "#" if package require installation:
# install.packages("dplyr")  # Data manipulation
# install.packages("tidyr")  # Data tidying
# install.packages("ggplot2")  # Visualization
# install.packages("stringdist")  # Distance calculations for ORFs or Gene_IDs
# install.packages("ape")  # Tree handling and clustering
# install.packages("cowplot")  # Combining plots
# install.packages("RColorBrewer")  # Custom color palettes

# Install Bioconductor manager to install ggtree:
# install.packages("BiocManager")
# BiocManager::install("ggtree")  # For tree visualization

# Alternatively, install ggtree using devtools (optional):
# install.packages("devtools")
# devtools::install_github("YuLab-SMU/ggtree")

# Load all required libraries at the beginning of your script:
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringdist)
library(ggtree)
library(ape)
library(cowplot)
library(RColorBrewer)

# Load necessary libraries to run this code
library(dplyr)  # For data manipulation
library(ggplot2)  # For plotting
library(RColorBrewer)  # For color palettes

# Step 1: Load the dataset
df <- read.delim('../ToxinKeyMatch/combined_output_keywords.tsv', sep = '\t', stringsAsFactors = FALSE)

# check duplicate rows based on the 'ID' column
df <- df %>%
  distinct(ID, .keep_all = TRUE)
nrow(df)
#4885

# Step 2: Define the classification logic using `case_when`
df <- df %>%
  mutate(FiltDeTox_Classification = case_when(
    # New Condition: If the wolfpsort_prediction does NOT contain 'extr' or 'E.R.', classify as 'Non-Toxins'
    !grepl("extr|E.R.", wolfpsort_prediction) ~ 'Non-Toxins',
    
    # Condition 1: Rating is SCD
    Rating == 'SCD' & pfam_ToxinKeywords == 'TRUE' & toxin_keywords == 'TRUE' ~ 'Toxins-Candidates',
    Rating == 'SCD' & pfam_ToxinKeywords == 'TRUE' & toxin_keywords == 'TRUE' & wolfpsort_prediction != 'extr' ~ 'Non-Toxins',
    Rating == 'SCD' & pfam_ToxinKeywords == 'TRUE' & toxin_keywords == 'FALSE' ~ 'Unlikely-Toxins',
    Rating == 'SCD' & pfam_ToxinKeywords == 'TRUE' & toxin_keywords == 'FALSE' & wolfpsort_prediction != 'extr' ~ 'Non-Toxins',
    Rating == 'SCD' & pfam_ToxinKeywords == 'FALSE' & toxin_keywords == 'TRUE' ~ 'Unlikely-Toxins',
    Rating == 'SCD' & pfam_ToxinKeywords == 'FALSE' & toxin_keywords == 'TRUE' & wolfpsort_prediction != 'extr' ~ 'Non-Toxins',
    Rating == 'SCD' & pfam_ToxinKeywords == 'FALSE' & toxin_keywords == 'FALSE' ~ 'Non-Toxins',
    
    # Condition 2: Rating is SCD! - both must be TRUE
    Rating == 'SCD!' & pfam_ToxinKeywords == 'TRUE' & toxin_keywords == 'TRUE' ~ 'Toxins-Candidates',
    Rating == 'SCD!' & pfam_ToxinKeywords == 'TRUE' & toxin_keywords == 'TRUE' & wolfpsort_prediction != 'extr' ~ 'Non-Toxins',
    Rating == 'SCD!' & pfam_ToxinKeywords == 'TRUE' & toxin_keywords == 'FALSE' ~ 'Unlikely-Toxins',
    Rating == 'SCD!' & pfam_ToxinKeywords == 'TRUE' & toxin_keywords == 'FALSE' & wolfpsort_prediction != 'extr' ~ 'Non-Toxins',
    Rating == 'SCD!' & pfam_ToxinKeywords == 'FALSE' & toxin_keywords == 'TRUE' ~ 'Unlikely-Toxins',
    Rating == 'SCD!' & pfam_ToxinKeywords == 'FALSE' & toxin_keywords == 'TRUE' & wolfpsort_prediction != 'extr' ~ 'Non-Toxins',
    Rating == 'SCD!' & pfam_ToxinKeywords == 'FALSE' & toxin_keywords == 'FALSE' ~ 'Non-Toxins',
    
    # Condition 3: Rating is SD, SD!, or D
    Rating %in% c('SD', 'SD!', 'D') & (toxin_keywords == 'TRUE' | pfam_ToxinKeywords == 'TRUE') ~ 'Unlikely-Toxins',
    Rating %in% c('SD', 'SD!', 'D') & (toxin_keywords == 'TRUE' | pfam_ToxinKeywords == 'TRUE') & wolfpsort_prediction != 'extr' ~ 'Non-Toxins',
    Rating %in% c('SD', 'SD!', 'D') & toxin_keywords == 'FALSE' & pfam_ToxinKeywords == 'FALSE' ~ 'Non-Toxins',
    
    # Condition 4: Rating is SC!
    Rating == 'SC!' & pfam_ToxinKeywords == 'TRUE' & toxin_keywords == 'TRUE' ~ 'Toxins-Candidates',
    Rating == 'SC!' & pfam_ToxinKeywords == 'TRUE' & toxin_keywords == 'FALSE' ~ 'Unlikely-Toxins',
    Rating == 'SC!' & pfam_ToxinKeywords == 'FALSE' & toxin_keywords == 'TRUE' ~ 'Unlikely-Toxins',
    Rating == 'SC!' & pfam_ToxinKeywords == 'FALSE' & toxin_keywords == 'FALSE' & grepl('uncharacterized', tolower(hit_descr)) ~ 'Unlikely-Toxins',
    Rating == 'SC!' & pfam_ToxinKeywords == 'FALSE' & toxin_keywords == 'FALSE' & !grepl('uncharacterized', tolower(hit_descr)) ~ 'Non-Toxins',
    
    # Condition 5: Rating is SC
    Rating == 'SC' & (toxin_keywords == 'TRUE' | pfam_ToxinKeywords == 'TRUE') ~ 'Toxins-Candidates',
    Rating == 'SC' & (toxin_keywords == 'TRUE' | pfam_ToxinKeywords == 'TRUE') & wolfpsort_prediction != 'extr' ~ 'Non-Toxins',
    Rating == 'SC' & (grepl('uncharacterized', tolower(hit_descr)) | hit_descr == "" | is.na(hit_descr)) ~ 'SCRs-WA',  # Modified condition to include empty hit_descr
    Rating == 'SC' & (grepl('uncharacterized', tolower(hit_descr)) | hit_descr == "" | is.na(hit_descr)) & wolfpsort_prediction != 'extr' ~ 'Non-Toxins',
    Rating == 'SC' & !grepl('uncharacterized', tolower(hit_descr)) & hit_descr != "" & !is.na(hit_descr) ~ 'Non-Toxins',
    
    # Condition 6: Rating is S or S!
    Rating %in% c('S', 'S!') & (toxin_keywords == 'TRUE' | pfam_ToxinKeywords == 'TRUE') & pct_hit_len_aligned > 60 ~ 'Unlikely-Toxins',
    Rating %in% c('S', 'S!') & (toxin_keywords == 'TRUE' | pfam_ToxinKeywords == 'TRUE') & pct_hit_len_aligned > 60 & wolfpsort_prediction != 'extr' ~ 'Non-Toxins',
    Rating %in% c('S', 'S!') & toxin_keywords == 'FALSE' & pfam_ToxinKeywords == 'FALSE' ~ 'Non-Toxins',
    
    # Additional conditions for other classifications...
    Rating %in% c('SBCD', 'SBCDT', 'SBC', 'SBCT') & wolfpsort_prediction == 'extr' ~ 'Toxins-Candidates',
    Rating %in% c('SBCD', 'SBCDT', 'SBC', 'SBCT') & wolfpsort_prediction != 'extr' ~ 'Unlikely-Toxins',
    
    Rating == 'SBDT' & wolfpsort_prediction == 'extr' ~ 'Toxins-Candidates',
    Rating == 'SBDT' & wolfpsort_prediction != 'extr' ~ 'Unlikely-Toxins',
    
    Rating == 'SBD' & pct_hit_len_aligned > 50 & wolfpsort_prediction == 'extr' ~ 'Toxins-Candidates',
    Rating == 'SBD' & pct_hit_len_aligned > 50 & wolfpsort_prediction != 'extr' ~ 'Unlikely-Toxins',
    Rating == 'SBD' & pct_hit_len_aligned <= 50 ~ 'Unlikely-Toxins',
    
    Rating == '*BD' & pct_hit_len_aligned > 60 & wolfpsort_prediction == 'extr' ~ 'Toxins-Candidates',
    Rating == '*BD' & pct_hit_len_aligned > 60 & wolfpsort_prediction != 'extr' ~ 'Unlikely-Toxins',
    Rating == '*BD' & pct_hit_len_aligned <= 60 ~ 'Unlikely-Toxins',
    
    Rating == '*B' & pct_hit_len_aligned > 60 & wolfpsort_prediction == 'extr' ~ 'Toxins-Candidates',
    Rating == '*B' & pct_hit_len_aligned > 60 & wolfpsort_prediction != 'extr' ~ 'Unlikely-Toxins',
    Rating == '*B' & pct_hit_len_aligned <= 60 ~ 'Unlikely-Toxins',
    
    Rating == 'BD' & pct_hit_len_aligned > 60 & wolfpsort_prediction == 'extr' ~ 'Toxins-Candidates',
    Rating == 'BD' & pct_hit_len_aligned > 60 & wolfpsort_prediction != 'extr' ~ 'Unlikely-Toxins',
    Rating == 'BD' & pct_hit_len_aligned <= 60 ~ 'Unlikely-Toxins',
    
    Rating == 'B' & pct_hit_len_aligned > 60 & wolfpsort_prediction == 'extr' ~ 'Toxins-Candidates',
    Rating == 'B' & pct_hit_len_aligned > 60 & wolfpsort_prediction != 'extr' ~ 'Unlikely-Toxins',
    Rating == 'B' & pct_hit_len_aligned <= 60 ~ 'Unlikely-Toxins',
    
    # Default case: Non-toxins
    TRUE ~ 'Non-Toxins'
  ))

# Step 3: Save each classification category to separate files
toxins_candidate <- df %>% filter(FiltDeTox_Classification == 'Toxins-Candidates')
write.table(toxins_candidate, file = 'Toxins_Candidates.tsv', sep = '\t', row.names = FALSE)

unlikely_toxin <- df %>% filter(FiltDeTox_Classification == 'Unlikely-Toxins')
write.table(unlikely_toxin, file = 'Unlikely_Toxins.tsv', sep = '\t', row.names = FALSE)

secreted_cysteine_rich <- df %>% filter(FiltDeTox_Classification == 'SCRs-WA')
write.table(secreted_cysteine_rich, file = 'SCRs_WA.tsv', sep = '\t', row.names = FALSE)

non_toxins <- df %>% filter(FiltDeTox_Classification == 'Non-Toxins')
write.table(non_toxins, file = 'Non_Toxins.tsv', sep = '\t', row.names = FALSE)

# Step 4: Generate and save classification statistics
stats <- data.frame(
  Category = c('Toxins-Candidates', 'Unlikely-Toxins', 'SCRs-WA', 'Non-Toxins'),
  Count = c(nrow(toxins_candidate), nrow(unlikely_toxin), nrow(secreted_cysteine_rich), nrow(non_toxins))
)
write.table(stats, file = 'FiltDeTox_Stats.tsv', sep = '\t', row.names = FALSE, col.names = TRUE)

# Step 5: Export the full dataframe with the new classification column
write.table(df, file = "Full_Classified_Data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Step 3: Save each classification category to separate files
# This saves each classification category to its own file for further analysis.

# Save Toxins Candidates
toxins_candidate <- df %>% filter(FiltDeTox_Classification == 'Toxins-Candidates')
write.table(toxins_candidate, file = 'Toxins_Candidates.tsv', sep = '\t', row.names = FALSE)

# Save Unlikely Toxins
unlikely_toxin <- df %>% filter(FiltDeTox_Classification == 'Unlikely-Toxins')
write.table(unlikely_toxin, file = 'Unlikely_Toxins.tsv', sep = '\t', row.names = FALSE)

# Save Secreted Cysteine-Rich Sequences (SCRs-WA)
secreted_cysteine_rich <- df %>% filter(FiltDeTox_Classification == 'SCRs-WA')
write.table(secreted_cysteine_rich, file = 'SCRs_WA.tsv', sep = '\t', row.names = FALSE)

# Save Non-Toxins
non_toxins <- df %>% filter(FiltDeTox_Classification == 'Non-Toxins')
write.table(non_toxins, file = 'Non_Toxins.tsv', sep = '\t', row.names = FALSE)

# Step 4: Generate and save classification statistics
# This generates a table of the counts for each classification and saves it to a file.
stats <- data.frame(
  Category = c('Toxins-Candidates', 'Unlikely-Toxins', 'SCRs-WA', 'Non-Toxins'),
  Count = c(nrow(toxins_candidate), nrow(unlikely_toxin), nrow(secreted_cysteine_rich), nrow(non_toxins))
)

write.table(stats, file = 'FiltDeTox_Stats.tsv', sep = '\t', row.names = FALSE, col.names = TRUE)

# Print a success message
print("Classification completed, files saved, and statistics generated.")

# Step 5: Export the full dataframe with the new classification column
write.table(df, file = "Full_Classified_Data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
print("Full dataframe with classifications exported to 'Full_Classified_Data.tsv'.")

#### Plot creation: Nested Pie Chart ####
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# 1. Load the classified dataset
df <- read.csv("Full_Classified_Data.tsv", sep = "\t")

# 2. Summarize the data by classification and rating
classification_summary <- df %>%
  group_by(FiltDeTox_Classification, Rating) %>%
  summarise(n = n()) %>%
  ungroup()

# 3. Calculate percentages for each classification and rating
classification_summary <- classification_summary %>%
  mutate(
    total_class = sum(n),  # Total count of all sequences
    total_per_class = ave(n, FiltDeTox_Classification, FUN = sum),  # Total per classification group
    perc_class = n / total_class * 100,  # Percent of the entire dataset
    perc_class_per_group = n / total_per_class * 100  # Percent within each classification group
  )

# 4. Reorder the classification and rating for better display
classification_summary$FiltDeTox_Classification <- factor(
  classification_summary$FiltDeTox_Classification, 
  levels = c("Toxins-Candidates", "Unlikely-Toxins", "SCRs-WA", "Non-Toxins")
)

classification_summary$Rating <- factor(classification_summary$Rating, levels = sort(unique(classification_summary$Rating)))

# 5. Define color palettes for classification and rating
num_classifications <- length(unique(classification_summary$FiltDeTox_Classification))
num_ratings <- length(unique(classification_summary$Rating))

classification_colors <- colorRampPalette(brewer.pal(12, "Set3"))(num_classifications)  # Colors for classification
rating_colors <- colorRampPalette(brewer.pal(9, "Set1"))(num_ratings)  # Colors for rating

# 6. Create labels with percentages for the legend
classification_labels <- classification_summary %>%
  group_by(FiltDeTox_Classification) %>%
  summarise(total = sum(n)) %>%
  mutate(label = paste0(FiltDeTox_Classification, " (", round(total / sum(total) * 100, 2), "%)"))

rating_labels <- classification_summary %>%
  group_by(Rating) %>%
  summarise(total = sum(n)) %>%
  mutate(label = paste0(Rating, " (", round(total / sum(total) * 100, 2), "%)"))

# 7. Create the Nested Pie-Donut Chart with percentages in the legend
# Save plot as PDF
pdf("NestedPie_FiltDeTox.pdf", width = 12, height = 8)
ggplot(classification_summary) +
  geom_bar(aes(x = 2, y = n, fill = FiltDeTox_Classification), stat = "identity", color = "white", width = 1) +  # Inner Pie (Classification)
  geom_bar(aes(x = 3, y = n, fill = Rating), stat = "identity", color = "white", width = 1) +  # Outer Donut (Rating)
  coord_polar(theta = "y") +  # Polar coordinates for pie chart
  scale_fill_manual(
    values = c(classification_colors, rating_colors), 
    labels = c(classification_labels$label, rating_labels$label)  # Apply labels with percentages
  ) +
  theme_void() +  # Remove axis and grid
  theme(
    legend.position = "right",           # Legend on the right
    legend.text = element_text(size = 12),   # Adjust legend text size
    legend.title = element_text(size = 14),  # Adjust legend title size
    plot.title = element_text(hjust = 0.5, size = 16),  # Center and size the title
    axis.text = element_blank(),    # Remove axis text
    axis.ticks = element_blank(),   # Remove axis ticks
    panel.grid = element_blank()    # Remove grid lines
  ) +
  labs(
    fill = "Classification / Rating",  # Legend title
    title = "FiltDeTox Classification with Rating (DeTox flags)"  # Plot title
  ) +
  xlim(0.5, 3.5)  # Adjust the x-axis limits to fit the nested rings

dev.off()

# 8. Save the plot as a high-resolution PNG file
# ggsave("NestedPie_FiltDeTox.png", width = 12, height = 8, dpi = 600)


# For PDF saving or other customized options, use RStudio's export options manually.

# Print message confirming save

# comment out the code below when running .sh script
# print("Nested Pie-Donut chart saved as 'NestedPie_FiltDeTox.png'.")

#### Plot creation: Nested Pie Chart (Flag_Sorted) ####
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# 1. Load the classified dataset
df <- read.csv("Full_Classified_Data.tsv", sep = "\t")

# 2. Summarize the data by classification and rating
classification_summary <- df %>%
  group_by(FiltDeTox_Classification, Rating) %>%
  summarise(n = n()) %>%
  ungroup()

# 3. Calculate percentages for each classification and rating
classification_summary <- classification_summary %>%
  mutate(
    total_class = sum(n),  # Total count of all sequences
    total_per_class = ave(n, FiltDeTox_Classification, FUN = sum),  # Total per classification group
    perc_class = n / total_class * 100,  # Percent of the entire dataset
    perc_class_per_group = n / total_per_class * 100  # Percent within each classification group
  )

# 4. Reorder the classification for better display
classification_summary$FiltDeTox_Classification <- factor(
  classification_summary$FiltDeTox_Classification, 
  levels = c("Toxins-Candidates", "Unlikely-Toxins", "SCRs-WA", "Non-Toxins")
)

# 4a. Sort the Rating factor by the total percentage contribution across all categories
rating_order <- classification_summary %>%
  group_by(Rating) %>%
  summarise(total_contribution = sum(n)) %>%
  arrange(total_contribution) %>%  # Sort in ascending order
  pull(Rating)

classification_summary$Rating <- factor(classification_summary$Rating, levels = rating_order)

# 5. Define color palettes for classification and rating
num_classifications <- length(unique(classification_summary$FiltDeTox_Classification))
num_ratings <- length(unique(classification_summary$Rating))

classification_colors <- colorRampPalette(brewer.pal(12, "Set3"))(num_classifications)  # Colors for classification
rating_colors <- colorRampPalette(brewer.pal(9, "Set1"))(num_ratings)  # Colors for rating

# 6. Create labels with percentages for the legend
classification_labels <- classification_summary %>%
  group_by(FiltDeTox_Classification) %>%
  summarise(total = sum(n)) %>%
  mutate(label = paste0(FiltDeTox_Classification, " (", round(total / sum(total) * 100, 2), "%)"))

rating_labels <- classification_summary %>%
  group_by(Rating) %>%
  summarise(total = sum(n)) %>%
  mutate(label = paste0(Rating, " (", round(total / sum(total) * 100, 2), "%)"))

# 7. Create the Nested Pie-Donut Chart with percentages in the legend (Flag_Sorted)
pdf("NestedPie_FiltDeTox_Flag_Sorted.pdf", width = 12, height = 8)
ggplot(classification_summary) +
  geom_bar(aes(x = 2, y = n, fill = FiltDeTox_Classification), stat = "identity", color = "white", width = 1) +  # Inner Pie (Classification)
  geom_bar(aes(x = 3, y = n, fill = Rating), stat = "identity", color = "white", width = 1) +  # Outer Donut (Rating)
  coord_polar(theta = "y") +  # Polar coordinates for pie chart
  scale_fill_manual(
    values = c(classification_colors, rating_colors), 
    labels = c(classification_labels$label, rating_labels$label)  # Apply labels with percentages
  ) +
  theme_void() +  # Remove axis and grid
  theme(
    legend.position = "right",           # Legend on the right
    legend.text = element_text(size = 12),   # Adjust legend text size
    legend.title = element_text(size = 14),  # Adjust legend title size
    plot.title = element_text(hjust = 0.5, size = 16),  # Center and size the title
    axis.text = element_blank(),    # Remove axis text
    axis.ticks = element_blank(),   # Remove axis ticks
    panel.grid = element_blank()    # Remove grid lines
  ) +
  labs(
    fill = "Classification / Rating",  # Legend title
    title = "FiltDeTox Classification with Rating (DeTox flags) - Flag_Sorted"  # Plot title
  ) +
  xlim(0.5, 3.5)  # Adjust the x-axis limits to fit the nested rings

dev.off()

# 8. Save the plot as a high-resolution PNG file
# ggsave("NestedPie_FiltDeTox_Flag_Sorted.png", width = 12, height = 8, dpi = 600)

# For PDF saving or other customized options, use RStudio's export options manually.

# Print message confirming save
# print("Nested Pie-Donut chart saved as 'NestedPie_FiltDeTox_Flag_Sorted.png'.")


### simple pie chart showing the relative percentage of each Rating (flag) within the "Toxins-Candidates" ###
# Load necessary libraries
library(dplyr)
library(ggplot2)

# 1. Load the classified dataset
df <- read.csv("Full_Classified_Data.tsv", sep = "\t")

# 2. Filter the data for the 'Toxins-Candidates' category
toxins_candidate_df <- df %>%
  filter(FiltDeTox_Classification == "Toxins-Candidates")

# 3. Summarize the data by Rating (flags) for the 'Toxins-Candidates' category
rating_summary <- toxins_candidate_df %>%
  group_by(Rating) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(percentage = (n / sum(n)) * 100)  # Calculate percentage

# 4. Sort the Rating by percentage in descending order for the legend
rating_summary <- rating_summary %>%
  arrange(desc(percentage)) %>%  # Sort in descending order by percentage
  mutate(Rating = factor(Rating, levels = Rating))  # Reorder the Rating factor based on percentage

# 5. Create labels with percentages for the legend
rating_summary <- rating_summary %>%
  mutate(label = paste0(Rating, " (", round(percentage, 2), "%)"))

# 6. Create the simple pie chart with percentages in the legend
pdf("Toxins_Candidate_Rating_PieChart.pdf", width = 8, height = 6)
ggplot(rating_summary, aes(x = "", y = percentage, fill = Rating)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +  # Convert to pie chart
  theme_void() +  # Remove axis and grid
  scale_fill_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(rating_summary$Rating))),
                    labels = rating_summary$label) +  # Apply labels with percentages
  labs(
    title = "Rating (flags) within Toxins-Candidates",
    fill = "Rating (Flags)"  # Legend title
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Center and size the title
    legend.title = element_text(size = 14),  # Adjust legend title size
    legend.text = element_text(size = 12)    # Adjust legend text size
  )

dev.off()

# 7. Save the pie chart as a PNG file
# ggsave("Toxins_Candidate_Rating_PieChart.png", width = 8, height = 6, dpi = 300)

# For PDF saving or other customized options, use RStudio's export options manually.

# Print message confirming save
# print("Pie chart for Toxins-Candidates Rating flags saved as 'Toxins_Candidate_Rating_PieChart.png'.")

### parsing and saving domain-wise statistics (number of ORFs, transcripts, genes, and total TPMs for each Pfam domain)
# Load necessary libraries
library(dplyr)
library(tidyr)

# Step 1: Load the Toxins_Candidates.tsv file into your R session
toxins_data <- read.table("Toxins_Candidates.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Step 2: Parse the `pfam.domains` column by semicolon and create new rows
toxins_data_parsed <- toxins_data %>%
  separate_rows(pfam.domains, sep = ";\\s*")  # Split by semicolon and optional spaces

# Step 3: Load the trusted Pfam domain list from 'ToxProt_domain_Keywords.tsv'
trusted_domains <- read.table("../ToxinKeyMatch/ToxProt_domain_Keywords.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(trusted_domains) <- c("pfam.domains")

# Step 4: Filter the parsed data to include only rows where `pfam.domains` matches trusted domains
toxins_data_filtered <- toxins_data_parsed %>%
  filter(pfam.domains %in% trusted_domains$pfam.domains)

# Step 5: Extract Gene, Transcript, and ORF Information from the 'ID' column
toxins_data_filtered <- toxins_data_filtered %>%
  mutate(
    gene = sub("_g.*", "", ID),  # Extract gene
    transcript = sub("_ORF.*", "", ID),  # Extract transcript
    orf = sub(".*_ORF\\.", "", ID)  # Extract ORF number
  )

# Step 6: Calculate domain-wise statistics (number of ORFs, transcripts, genes, and total TPMs for each Pfam domain)
domain_summary <- toxins_data_filtered %>%
  group_by(pfam.domains) %>%
  summarise(
    num_orfs = n_distinct(ID),                   # Number of unique ORFs
    num_transcripts = n_distinct(transcript),    # Number of unique transcripts
    num_genes = n_distinct(gene),                # Number of unique genes
    total_TPM = sum(as.numeric(TPM), na.rm = TRUE) # Total TPM (summing TPM values)
  ) %>%
  arrange(desc(num_orfs))  # Sort by the number of ORFs in descending order

# Step 7: Group by Pfam domain and collect ORF IDs and Gene IDs
toxins_data_filtered <- toxins_data_filtered %>%
  group_by(pfam.domains) %>%
  summarise(
    ORFs_ID = paste(unique(ID), collapse = "; "),    # Concatenate unique ORF IDs
    Gene_ID = paste(unique(gene), collapse = "; ")   # Concatenate unique Gene IDs
  )

# Step 8: Merge the ORF and Gene ID information with the domain summary
pfam_summary_with_ids <- domain_summary %>%
  left_join(toxins_data_filtered, by = "pfam.domains")  # Merge by pfam.domains

# Step 9: Save the final table with ORFs and Gene IDs
write.table(pfam_summary_with_ids, file = "Pfam_Domain_Summary_with_ORFs_Genes.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Optional: View the updated summary data
head(pfam_summary_with_ids)

# Print completion message
print("Pfam_Domain_Summary_with_ORFs_Genes.tsv has been created with domain statistics and ORFs_ID, Gene_ID.")


# This part shows the hierarchical clustering of Pfam domains, based on the similarity of the ORFs associated with each Pfam domain.
# The clustering was performed using binary distance to compare ORF presence/absence for each domain, highlighting relationships between domains based on shared ORFs.
# The resulting dendrogram and dot plot visualize both domain clustering (based on shared ORFs) and domain characteristics (ORF count and TPM).

# Install the required packages
# install.packages("stringdist")
# install.packages("ape")
# install.packages("cowplot")

# install.packages("BiocManager")
# BiocManager::install("ggtree")

# Optionally for installing ggtree
#install.packages("devtools")
#devtools::install_github("YuLab-SMU/ggtree")

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(stringdist)
library(ggtree)
library(ape)  # For tree handling
library(cowplot)

### Dendrogram reorder based on ORFs instead of genes
# Step 1: Load the Pfam_Domain_Summary_with_ORFs_Genes.tsv file into your R session
pfam_data <- read.delim("Pfam_Domain_Summary_with_ORFs_Genes.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)

# Step 2: Sort the data by pfam.domains
pfam_data_sorted <- pfam_data %>%
  arrange(pfam.domains)  # Sort alphabetically by pfam.domains

# Step 3: Parse ORFs_ID to calculate similarity between Pfam domains based on shared ORFs
# Create a binary matrix where rows are Pfam domains and columns are unique ORF IDs
orf_matrix <- pfam_data_sorted %>%
  separate_rows(ORFs_ID, sep = ";\\s*") %>%
  mutate(value = 1) %>%
  spread(ORFs_ID, value, fill = 0)

# Use only the ORF matrix without the pfam.domains column
row.names(orf_matrix) <- pfam_data_sorted$pfam.domains
orf_matrix <- orf_matrix[ , -1]  # Remove the pfam.domains column to keep only the ORF matrix

# Step 4: Perform hierarchical clustering on the ORF presence/absence matrix
distance_matrix <- dist(orf_matrix, method = "binary")  # Binary distance
hclust_result <- hclust(distance_matrix, method = "complete")  # Hierarchical clustering

# Convert hclust object to a phylo object for tree visualization
tree <- as.phylo(hclust_result)

# Step 5: Plot the tree using ggtree
p_tree <- ggtree(tree, layout = "rectangular") + 
  geom_tiplab(size = 3, hjust = -0.1) +  # Add Pfam domains as labels to the tips
  theme_tree2() +
  labs(title = "Hierarchical Clustering of Pfam Domains Based on ORFs")

# Step 6: Create the dot plot layer (size = ORFs, color = TPM)
p_dotplot_combined <- ggplot(pfam_data_sorted, aes(x = 1, y = pfam.domains)) +
  geom_point(aes(size = num_orfs, color = total_TPM), alpha = 0.7) +
  scale_size_continuous(range = c(2, 10)) +  # Adjust dot sizes based on ORF counts
  scale_color_viridis_c(option = "C", trans = "log")  # Log scale for color based on TPM
labs(size = "Number of ORFs", color = "Total TPM", y = "Pfam Domains") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Step 7: Combine the tree and dot plot into a single layout
combined_plot <- plot_grid(
  p_tree, p_dotplot_combined, 
  ncol = 2, align = "h", rel_widths = c(2, 1)  # Adjust the widths for better visibility
)

# Save the combined plot as a PDF
pdf("Dendrogram_and_DotPlot_ORFs_TPM_by_ORF.pdf", width = 14, height = 8)
print(combined_plot)
dev.off()  # Close the PDF device

# Display the combined plot
# print(combined_plot)

# Optional: Save the combined plot
# ggsave("Dendrogram_and_DotPlot_ORFs_TPM_by_ORF.png", plot = combined_plot, width = 14, height = 8)

# Print completion message
print("Combined dendrogram and dot plot saved as 'Dendrogram_and_DotPlot_ORFs_TPM_by_ORF.png'.")


### creating FASTA with mature peptides from SCRs-WA for in-silico bioactivity tests
# Load necessary library
library(dplyr)

# Step 1: Read the SCRs-WA.tsv file
# Assuming your SCRs-WA.tsv file is in the working directory
scrs_wa_data <- read.delim("SCRs_WA.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
nrow(scrs_wa_data)


# Step 2: Extract ORF names and mature peptide sequences
fasta_data <- scrs_wa_data %>%
  select(ID, mature_peptide) %>%
  filter(!is.na(mature_peptide))  # Ensure that there are no missing peptide sequences

# Step 3: Write to a FASTA formatted file
# Open a connection to write the file
fasta_file <- file("SCRs-WA.fasta", "w")

# Loop through each row of the data to write in FASTA format
for (i in 1:nrow(fasta_data)) {
  cat(">", fasta_data$ID[i], "\n", fasta_data$mature_peptide[i], "\n", file = fasta_file, sep = "")
}

# Close the file connection
close(fasta_file)
# Print confirmation
print("FASTA file 'SCRs-WA.fasta' created successfully.")

### creating FASTA with full precursor sequences from SCRs-WA
# Load necessary library
library(dplyr)

# Step 1: Read the SCRs-WA.tsv file
# Assuming your SCRs-WA.tsv file is in the working directory
scrs_wa_data <- read.delim("SCRs_WA.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)

# Step 2: Extract ORF names and precursor sequences (from the "Sequence" column)
fasta_data <- scrs_wa_data %>%
  select(ID, Sequence) %>%
  filter(!is.na(Sequence))  # Ensure that there are no missing sequences

# Step 3: Write to a FASTA formatted file
# Open a connection to write the file
fasta_file <- file("SCRs-WA_precursor.fasta", "w")

# Loop through each row of the data to write in FASTA format
for (i in 1:nrow(fasta_data)) {
  cat(">", fasta_data$ID[i], "\n", fasta_data$Sequence[i], "\n", file = fasta_file, sep = "")
}

# Close the file connection
close(fasta_file)
# Print confirmation
print("FASTA file 'SCRs-WA_precursor.fasta' created successfully.")

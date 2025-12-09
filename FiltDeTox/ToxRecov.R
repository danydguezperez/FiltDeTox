###############################################################################
# FiltDeTox: A Tool to Enhance and Filter Animal Toxins Identification
#            from RNAseq Analyses
#
# ToxRecov.R – classification, summaries, plots and exports
###############################################################################

############################
# Package installation     #
############################
# Uncomment the lines below if any package is missing.

# Core data manipulation and plotting
# install.packages("dplyr")        # data manipulation (filter, mutate, group_by, joins)
# install.packages("tidyr")        # reshaping tables (separate_rows, spread)
# install.packages("ggplot2")      # plotting (barplots, pie charts, dot plots)

# Tree, clustering and multi-panel plots
# install.packages("ape")          # hierarchical clustering as phylo objects
# install.packages("cowplot")      # combining multiple ggplot figures

# Colour and layout utilities
# install.packages("RColorBrewer") # qualitative colour palettes
# install.packages("grid")         # layout utilities (unit() for legend sizing)
# install.packages("ggnewscale")   # multiple fill/colour scales (split legends)

# Bioconductor packages
# install.packages("BiocManager")
# BiocManager::install("ggtree")   # tree visualisation for dendrograms
# BiocManager::install("Biostrings") # reading/writing FASTA and sequence handling


############################
# Library loading          #
############################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtree)
library(ape)
library(cowplot)
library(RColorBrewer)
library(grid)
library(ggnewscale)
library(Biostrings)   # used in Step 10 for full-length sequence handling


#######################################
### Step 1: Load and merge datasets ###
#######################################

# ToxinKeyMatch output (combined BLASTp + DeTox + keyword/domain information)
ToxinKeyMatch <- read.delim(
  "../ToxinKeyMatch/combined_output_keywords.tsv",
  sep = "\t",
  stringsAsFactors = FALSE
)

# HMMER-based mapping of ORFs to VenomZone toxin families
hhmer_orfs <- read.csv(
  "../hhmer_tx_VenomZone/hhmer_Tx_orf_mapping.csv",
  stringsAsFactors = FALSE
)

# Merge HMMER-based toxin family annotation into ToxinKeyMatch
# and ensure a single row per ORF (ID)
df <- ToxinKeyMatch %>%
  left_join(hhmer_orfs, by = c("ID" = "ORF")) %>%
  distinct(ID, .keep_all = TRUE) %>%
  mutate(ORF_precursor_length = nchar(Sequence))


#############################################
### Step 2: FiltDeTox classification rule ###
#############################################

# Assign each ORF to one of the four FiltDeTox categories:
#   - Toxins-Candidates
#   - Unlikely-Toxins
#   - SCRs-WA (secreted cysteine-rich with weak/absent annotation)
#   - Non-Toxins
df <- df %>%
  mutate(
    FiltDeTox_Classification = case_when(
      # 1) Not predicted as secreted
      !grepl("extr|E.R.", wolfpsort_prediction) ~ "Non-Toxins",
      
      # 2) Toxin family annotated by HMMER
      !is.na(ToxinFamily) & ORF_precursor_length > 500 ~ "Unlikely-Toxins",
      !is.na(ToxinFamily) & grepl("[*]", Rating)       ~ "Unlikely-Toxins",
      !is.na(ToxinFamily)                              ~ "Toxins-Candidates",
      
      # 3) Strong DeTox “toxin-like” flags
      Rating %in% c("SBCD", "SBCDT", "SBC", "SBCT") &
        ORF_precursor_length > 500                     ~ "Unlikely-Toxins",
      Rating %in% c("SBCD", "SBCDT", "SBC", "SBCT")    ~ "Toxins-Candidates",
      
      # 4) Strong keyword + Pfam toxin-domain evidence
      toxin_keywords == "TRUE" & pfam_ToxinKeywords == "TRUE" &
        ORF_precursor_length > 500                     ~ "Non-Toxins",
      toxin_keywords == "TRUE" & pfam_ToxinKeywords == "TRUE" ~ "Unlikely-Toxins",
      
      # 5) BD/B/SBD flags with long ORFs
      Rating %in% c("*BD", "*B", "BD", "B", "SBD") &
        ORF_precursor_length > 500                     ~ "Non-Toxins",
      Rating %in% c("*BD", "*B", "BD", "B", "SBD")     ~ "Unlikely-Toxins",
      
      # 6) SC flags – uncharacterised vs characterised sequences
      Rating == "SC" &
        (grepl("uncharacterized", tolower(hit_descr)) |
           hit_descr == "" | is.na(hit_descr))         ~ "SCRs-WA",
      Rating == "SC"                                   ~ "Non-Toxins",
      
      # 7) Default case
      TRUE                                            ~ "Non-Toxins"
    )
  )


#########################################################
### Step 3: Write per-class tables and summary stats  ###
#########################################################

# Split into four FiltDeTox categories
toxins_candidate       <- df %>% filter(FiltDeTox_Classification == "Toxins-Candidates")
unlikely_toxin         <- df %>% filter(FiltDeTox_Classification == "Unlikely-Toxins")
secreted_cysteine_rich <- df %>% filter(FiltDeTox_Classification == "SCRs-WA")
non_toxins             <- df %>% filter(FiltDeTox_Classification == "Non-Toxins")

# Per-category TSV exports
write.table(toxins_candidate,       "Toxins_Candidates.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(unlikely_toxin,         "Unlikely_Toxins.tsv",   sep = "\t", row.names = FALSE, quote = FALSE)
write.table(secreted_cysteine_rich, "SCRs_WA.tsv",           sep = "\t", row.names = FALSE, quote = FALSE)
write.table(non_toxins,             "Non_Toxins.tsv",        sep = "\t", row.names = FALSE, quote = FALSE)

# Category counts
stats <- data.frame(
  Category = c("Toxins-Candidates", "Unlikely-Toxins", "SCRs-WA", "Non-Toxins"),
  Count    = c(nrow(toxins_candidate),
               nrow(unlikely_toxin),
               nrow(secreted_cysteine_rich),
               nrow(non_toxins))
)

write.table(stats, "FiltDeTox_Stats.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# Full classified table
write.table(df, "Full_Classified_Data.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)


###########################################################
### Step 4: Sorted stacked barplots (classes and flags) ###
###########################################################

# These plots summarise the composition of:
#   - FiltDeTox classes (one stacked bar)
#   - DeTox Rating flags (one stacked bar)

df_bar <- read.csv("Full_Classified_Data.tsv", sep = "\t", stringsAsFactors = FALSE)

percent_labels <- function(x) paste0(round(x * 100, 1), "%")

#### 4A. Sorted stacked barplot: FiltDeTox classes ####

class_order_sorted <- df_bar %>%
  count(FiltDeTox_Classification) %>%
  arrange(n) %>%
  pull(FiltDeTox_Classification)

df_bar$FiltDeTox_Classification_sorted <- factor(
  df_bar$FiltDeTox_Classification,
  levels = class_order_sorted
)

class_cols <- c(
  "Toxins-Candidates" = "#F0E442",
  "Unlikely-Toxins"   = "#d95f02",
  "SCRs-WA"           = "#7570b3",
  "Non-Toxins"        = "#1b9e77"
)

p_class_sorted <- ggplot(
  df_bar,
  aes(x = "FiltDeTox classification (sorted)",
      fill = FiltDeTox_Classification_sorted)
) +
  geom_bar(position = "fill", colour = "black", width = 0.6) +
  coord_flip() +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(
    name   = "Percentage of sequences",
    labels = percent_labels,
    limits = c(0, 1)
  ) +
  scale_fill_manual(values = class_cols, drop = FALSE) +
  labs(
    fill  = "FiltDeTox classification",
    title = "StackedBar_FiltDeTox_Classification_Sorted"
  ) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid   = element_blank(),
    plot.title   = element_text(hjust = 0.5, size = 16),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10)
  )

ggsave("StackedBar_FiltDeTox_Classification_Sorted.png",
       p_class_sorted, width = 8, height = 3, dpi = 600)
ggsave("StackedBar_FiltDeTox_Classification_Sorted.pdf",
       p_class_sorted, width = 8, height = 3)

#### 4B. Sorted stacked barplot: DeTox Rating flags ####

df_flags_sorted <- df_bar %>% filter(!is.na(Rating))

rating_order_sorted <- df_flags_sorted %>%
  count(Rating) %>%
  arrange(n) %>%
  pull(Rating)

df_flags_sorted$Rating_sorted <- factor(
  df_flags_sorted$Rating,
  levels = rating_order_sorted
)

rating_palette_sorted <- colorRampPalette(
  brewer.pal(9, "Set1")
)(length(rating_order_sorted))
rating_cols_sorted <- setNames(rating_palette_sorted, rating_order_sorted)

p_flags_sorted <- ggplot(
  df_flags_sorted,
  aes(x = "DeTox flags (sorted)",
      fill = Rating_sorted)
) +
  geom_bar(position = "fill", colour = "black", width = 0.9) +
  coord_flip() +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(
    name   = "Percentage of sequences",
    labels = percent_labels,
    limits = c(0, 1)
  ) +
  scale_fill_manual(values = rating_cols_sorted, drop = FALSE) +
  labs(
    fill  = "DeTox flags (Rating)",
    title = "StackedBar_DeTox_Flags_Sorted"
  ) +
  guides(
    fill = guide_legend(
      reverse = TRUE,
      ncol    = 1
    )
  ) +
  theme_bw() +
  theme(
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    panel.grid       = element_blank(),
    plot.title       = element_text(hjust = 0.5, size = 16),
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(0.4, "cm"),
    legend.spacing.y = unit(0.02, "cm")
  )

ggsave("StackedBar_DeTox_Flags_Sorted.png",
       p_flags_sorted, width = 8, height = 3, dpi = 600)
ggsave("StackedBar_DeTox_Flags_Sorted.pdf",
       p_flags_sorted, width = 8, height = 3)


#######################################################
### Step 5: Nested pie charts (FiltDeTox + DeTox)   ###
#######################################################

# Inner ring: FiltDeTox classification (4 categories)
# Outer ring: DeTox Rating flags (one contiguous slice per flag)

# Helper to normalise Rating values to avoid artificial splits
clean_rating <- function(x) {
  x <- trimws(x)          # remove leading/trailing whitespace
  x <- toupper(x)         # harmonise case (s2 vs S2)
  x[x == "" | x == "NA"] <- NA
  x
}

##############################
## 5.1 Nested pie (unsorted) #
##############################

df_nested <- read.csv("Full_Classified_Data.tsv",
                      sep = "\t", stringsAsFactors = FALSE)

# Normalise Rating and drop missing
df_nested$Rating <- clean_rating(df_nested$Rating)
df_nested <- df_nested %>% filter(!is.na(Rating))

# Summary by FiltDeTox class and Rating (for inner proportions per class)
class_rating_summary <- df_nested %>%
  group_by(FiltDeTox_Classification, Rating) %>%
  summarise(n = n(), .groups = "drop")

# Inner ring: totals per FiltDeTox class
inner_summary <- class_rating_summary %>%
  group_by(FiltDeTox_Classification) %>%
  summarise(n = sum(n), .groups = "drop")

# Outer ring: totals per Rating across all classes
rating_summary_all <- df_nested %>%
  group_by(Rating) %>%
  summarise(n = n(), .groups = "drop")

# Factor order for classes (inner ring)
inner_summary$FiltDeTox_Classification <- factor(
  inner_summary$FiltDeTox_Classification,
  levels = c("Toxins-Candidates", "Unlikely-Toxins", "SCRs-WA", "Non-Toxins")
)

# Colour palettes (named vectors)
num_classifications <- length(levels(inner_summary$FiltDeTox_Classification))
num_ratings <- nrow(rating_summary_all)

classification_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Set3"))(num_classifications),
  levels(inner_summary$FiltDeTox_Classification)
)

rating_colors <- setNames(
  colorRampPalette(brewer.pal(9, "Set1"))(num_ratings),
  rating_summary_all$Rating
)

# Legend labels (with percentages)
classification_labels <- inner_summary %>%
  mutate(
    perc  = n / sum(n) * 100,
    label = paste0(FiltDeTox_Classification, " (", round(perc, 2), "%)")
  )

rating_labels <- rating_summary_all %>%
  mutate(
    perc  = n / sum(n) * 100,
    label = paste0(Rating, " (", round(perc, 2), "%)")
  )

class_breaks <- classification_labels$FiltDeTox_Classification
class_labs   <- classification_labels$label

rating_breaks <- rating_labels$Rating
rating_labs   <- rating_labels$label

# Nested pie with independent inner/outer rings
p_nested <- ggplot() +
  # Inner ring: FiltDeTox classes
  geom_bar(
    data = inner_summary,
    aes(x = 2, y = n, fill = FiltDeTox_Classification),
    stat  = "identity",
    colour = "white",
    width  = 1
  ) +
  scale_fill_manual(
    name   = "FiltDeTox classification",
    values = classification_colors,
    breaks = class_breaks,
    labels = class_labs
  ) +
  ggnewscale::new_scale_fill() +
  # Outer ring: global DeTox flags (one slice per Rating)
  geom_bar(
    data = rating_summary_all,
    aes(x = 3, y = n, fill = Rating),
    stat  = "identity",
    colour = "white",
    width  = 1
  ) +
  scale_fill_manual(
    name   = "DeTox rating",
    values = rating_colors,
    breaks = rating_breaks,
    labels = rating_labs
  ) +
  coord_polar(theta = "y") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.box      = "vertical",
    legend.text     = element_text(size = 12),
    legend.title    = element_text(size = 14),
    plot.title      = element_text(hjust = 0.5, size = 16),
    plot.subtitle   = element_text(hjust = 0.5, size = 12),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank()
  ) +
  labs(
    title    = "FiltDeTox classification with DeTox flags",
    subtitle = "Inner ring: FiltDeTox classification; outer ring: DeTox flags"
  ) +
  xlim(0.5, 3.5)

ggsave("NestedPie_FiltDeTox.png", p_nested, width = 12, height = 8, dpi = 600)
ggsave("NestedPie_FiltDeTox.pdf",  p_nested, width = 12, height = 8)


##############################################
## 5.2 Nested pie (Flag_Sorted – by rating) ##
##############################################

df_nested2 <- read.csv("Full_Classified_Data.tsv",
                       sep = "\t", stringsAsFactors = FALSE)

df_nested2$Rating <- clean_rating(df_nested2$Rating)
df_nested2 <- df_nested2 %>% filter(!is.na(Rating))

# Inner ring again from FiltDeTox class totals
inner_summary2 <- df_nested2 %>%
  group_by(FiltDeTox_Classification) %>%
  summarise(n = n(), .groups = "drop")

# Outer ring: totals per Rating, sorted by abundance
rating_summary_all2 <- df_nested2 %>%
  group_by(Rating) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(n)

rating_order <- rating_summary_all2$Rating

inner_summary2$FiltDeTox_Classification <- factor(
  inner_summary2$FiltDeTox_Classification,
  levels = c("Toxins-Candidates", "Unlikely-Toxins", "SCRs-WA", "Non-Toxins")
)

num_classifications2 <- length(levels(inner_summary2$FiltDeTox_Classification))
num_ratings2 <- nrow(rating_summary_all2)

classification_colors2 <- setNames(
  colorRampPalette(brewer.pal(12, "Set3"))(num_classifications2),
  levels(inner_summary2$FiltDeTox_Classification)
)

rating_colors2 <- setNames(
  colorRampPalette(brewer.pal(9, "Set1"))(num_ratings2),
  rating_order
)

classification_labels2 <- inner_summary2 %>%
  mutate(
    perc  = n / sum(n) * 100,
    label = paste0(FiltDeTox_Classification, " (", round(perc, 2), "%)")
  )

rating_labels2 <- rating_summary_all2 %>%
  mutate(
    perc  = n / sum(n) * 100,
    label = paste0(Rating, " (", round(perc, 2), "%)")
  )

class_breaks2 <- classification_labels2$FiltDeTox_Classification
class_labs2   <- classification_labels2$label

rating_breaks2 <- rating_labels2$Rating
rating_labs2   <- rating_labels2$label

p_nested_flag_sorted <- ggplot() +
  # Inner ring: FiltDeTox classes
  geom_bar(
    data = inner_summary2,
    aes(x = 2, y = n, fill = FiltDeTox_Classification),
    stat  = "identity",
    colour = "white",
    width  = 1
  ) +
  scale_fill_manual(
    name   = "FiltDeTox classification",
    values = classification_colors2,
    breaks = class_breaks2,
    labels = class_labs2
  ) +
  ggnewscale::new_scale_fill() +
  # Outer ring: DeTox flags sorted by abundance
  geom_bar(
    data = rating_summary_all2,
    aes(x = 3, y = n, fill = Rating),
    stat  = "identity",
    colour = "white",
    width  = 1
  ) +
  scale_fill_manual(
    name   = "DeTox rating",
    values = rating_colors2,
    breaks = rating_breaks2,
    labels = rating_labs2
  ) +
  coord_polar(theta = "y") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.box      = "vertical",
    legend.text     = element_text(size = 12),
    legend.title    = element_text(size = 14),
    plot.title      = element_text(hjust = 0.5, size = 16),
    plot.subtitle   = element_text(hjust = 0.5, size = 12),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank()
  ) +
  labs(
    title    = "FiltDeTox classification with DeTox flags",
    subtitle = "Inner ring: FiltDeTox classification; outer ring: DeTox flags"
  ) +
  xlim(0.5, 3.5)

ggsave("NestedPie_FiltDeTox_Flag_Sorted.png",
       p_nested_flag_sorted, width = 12, height = 8, dpi = 600)
ggsave("NestedPie_FiltDeTox_Flag_Sorted.pdf",
       p_nested_flag_sorted, width = 12, height = 8)



##############################################################
### Step 6: Simple pie chart for Rating within Toxins     ###
##############################################################

# Pie chart showing the proportion of DeTox flags restricted
# to the "Toxins-Candidates" category.

df_tc <- read.csv("Full_Classified_Data.tsv", sep = "\t", stringsAsFactors = FALSE)

df_tc <- df_tc %>%
  mutate(
    Rating = trimws(Rating),
    Rating = ifelse(Rating == "", NA, Rating)
  ) %>%
  filter(!is.na(Rating))

toxins_candidate_df <- df_tc %>%
  filter(FiltDeTox_Classification == "Toxins-Candidates")

rating_summary_tc <- toxins_candidate_df %>%
  group_by(Rating) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(percentage = (n / sum(n)) * 100) %>%
  arrange(desc(percentage)) %>%
  mutate(Rating = factor(Rating, levels = Rating)) %>%
  mutate(label = paste0(Rating, " (", round(percentage, 2), "%)"))

p_pie_tc <- ggplot(rating_summary_tc, aes(x = "", y = percentage, fill = Rating)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_manual(
    values = colorRampPalette(brewer.pal(9, "Set1"))(
      length(unique(rating_summary_tc$Rating))
    ),
    labels = rating_summary_tc$label
  ) +
  labs(
    title = "Rating (flags) within Toxins-Candidates",
    fill  = "DeTox flags"
  ) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  )

ggsave("Toxins_Candidate_Rating_PieChart.png",
       p_pie_tc, width = 8, height = 6, dpi = 300)
ggsave("Toxins_Candidate_Rating_PieChart.pdf",
       p_pie_tc, width = 8, height = 6)


########################################################################
### Step 7: Domain-wise statistics (ORFs, transcripts, genes, TPMs) ###
########################################################################

# Per-Pfam domain statistics restricted to Toxins-Candidates:
#   - number of ORFs / transcripts / genes
#   - total TPM
#   - list of ORF IDs and gene IDs per Pfam domain

toxins_data <- read.table("Toxins_Candidates.tsv",
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)

toxins_data_parsed <- toxins_data %>%
  separate_rows(pfam.domains, sep = ";\\s*")

trusted_domains <- read.table("../ToxinKeyMatch/ToxProt_domain_Keywords.tsv",
                              header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(trusted_domains) <- "pfam.domains"

toxins_data_filtered <- toxins_data_parsed %>%
  filter(pfam.domains %in% trusted_domains$pfam.domains) %>%
  mutate(
    gene       = sub("_g.*", "", ID),
    transcript = sub("_ORF.*", "", ID),
    orf        = sub(".*_ORF\\.", "", ID)
  )

domain_summary <- toxins_data_filtered %>%
  group_by(pfam.domains) %>%
  summarise(
    num_orfs        = n_distinct(ID),
    num_transcripts = n_distinct(transcript),
    num_genes       = n_distinct(gene),
    total_TPM       = sum(as.numeric(TPM), na.rm = TRUE),
    .groups         = "drop"
  ) %>%
  arrange(desc(num_orfs))

ids_per_domain <- toxins_data_filtered %>%
  group_by(pfam.domains) %>%
  summarise(
    ORFs_ID = paste(unique(ID),   collapse = "; "),
    Gene_ID = paste(unique(gene), collapse = "; "),
    .groups = "drop"
  )

pfam_summary_with_ids <- domain_summary %>%
  left_join(ids_per_domain, by = "pfam.domains")

write.table(pfam_summary_with_ids,
            file = "Pfam_Domain_Summary_with_ORFs_Genes.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)


#################################################################
### Step 8: Dendrogram and dot plot for Pfam domain profiles ###
#################################################################

pfam_data <- read.delim("Pfam_Domain_Summary_with_ORFs_Genes.tsv",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)

pfam_data_sorted <- pfam_data %>%
  arrange(pfam.domains)

orf_matrix <- pfam_data_sorted %>%
  separate_rows(ORFs_ID, sep = ";\\s*") %>%
  mutate(value = 1) %>%
  spread(ORFs_ID, value, fill = 0)

row.names(orf_matrix) <- pfam_data_sorted$pfam.domains
orf_matrix <- orf_matrix[, -1]

distance_matrix <- dist(orf_matrix, method = "binary")
hclust_result   <- hclust(distance_matrix, method = "complete")
tree            <- as.phylo(hclust_result)

tree$tip.label <- rep("", length(tree$tip.label))

p_tree <- ggtree(tree, layout = "rectangular") +
  geom_tiplab(size = 3, hjust = -0.1) +
  theme_tree2() +
  labs(title = "Hierarchical Clustering of Pfam Domains Based on ORF Presence")

p_dotplot_combined <- ggplot(pfam_data_sorted,
                             aes(x = 1, y = pfam.domains)) +
  geom_point(aes(size = num_orfs, color = total_TPM), alpha = 0.7) +
  scale_size_continuous(range = c(2, 10)) +
  scale_color_viridis_c(
    option = "C",
    trans  = "log",
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(
    size  = "Number of ORFs",
    color = "Total TPM",
    y     = "Pfam domains"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

combined_plot <- plot_grid(
  p_tree, p_dotplot_combined,
  ncol = 2, align = "h", rel_widths = c(2, 1)
)

ggsave("Dendrogram_and_DotPlot_ORFs_TPM_by_ORF.png",
       combined_plot, width = 14, height = 8)
ggsave("Dendrogram_and_DotPlot_ORFs_TPM_by_ORF.pdf",
       combined_plot, width = 14, height = 8)


############################################################
### Step 9: FASTA exports (mature and precursor)         ###
############################################################

# Helper functions to export mature peptide FASTA and precursor FASTA
# for each FiltDeTox category.

generate_mature_fasta <- function(tsv_file, fasta_file_name) {
  data <- read.delim(tsv_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  fasta_data <- data %>%
    select(ID, mature_peptide) %>%
    filter(!is.na(mature_peptide))
  
  fasta_file <- file(fasta_file_name, "w")
  for (i in seq_len(nrow(fasta_data))) {
    cat(">", fasta_data$ID[i], "\n",
        fasta_data$mature_peptide[i], "\n",
        file = fasta_file, sep = "")
  }
  close(fasta_file)
}

generate_precursor_fasta <- function(tsv_file, fasta_file_name) {
  data <- read.delim(tsv_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  fasta_data <- data %>%
    select(ID, Sequence) %>%
    filter(!is.na(Sequence))
  
  fasta_file <- file(fasta_file_name, "w")
  for (i in seq_len(nrow(fasta_data))) {
    cat(">", fasta_data$ID[i], "\n",
        fasta_data$Sequence[i], "\n",
        file = fasta_file, sep = "")
  }
  close(fasta_file)
}

# Export mature and precursor sequences for all four FiltDeTox categories
generate_mature_fasta("Non_Toxins.tsv",          "Non_Toxins_mature.fasta")
generate_precursor_fasta("Non_Toxins.tsv",       "Non_Toxins_precursor.fasta")

generate_mature_fasta("Toxins_Candidates.tsv",   "Toxins_Candidates_mature.fasta")
generate_precursor_fasta("Toxins_Candidates.tsv","Toxins_Candidates_precursor.fasta")

generate_mature_fasta("Unlikely_Toxins.tsv",     "Unlikely_Toxins_mature.fasta")
generate_precursor_fasta("Unlikely_Toxins.tsv",  "Unlikely_Toxins_precursor.fasta")

generate_mature_fasta("SCRs_WA.tsv",             "SCRs_WA_mature.fasta")
generate_precursor_fasta("SCRs_WA.tsv",          "SCRs_WA_precursor.fasta")


#########################################################################
### Step 10: Full-length matching using DeTox FASTA (Biostrings)     ###
#########################################################################

# This step identifies sequences that are labelled as full-length
# (type:complete) in the DeTox FASTA and propagates this information to all
# four FiltDeTox categories. For each category the script generates:
#   - <Category>_full-length.tsv (with a full.length column)
#   - <Category>_precursor_full-length-seqs.fasta (precursor sequences only)
#   - <Category>_full-length_stats.txt (short text summary)

# Locate DeTox FASTA in the parent directory (assumes a single FASTA file)
parent_dir  <- normalizePath("..")
fasta_files <- list.files(path = parent_dir, pattern = "\\.fasta$", full.names = TRUE)

if (length(fasta_files) == 0) {
  stop("No .fasta file found one level up from current working directory.")
}

detox_fasta_path  <- fasta_files[1]
detox_fasta_label <- basename(detox_fasta_path)

# Extract IDs for sequences annotated as type:complete
fasta_sequences <- readAAStringSet(detox_fasta_path)
headers         <- names(fasta_sequences)
complete_ids    <- headers[grepl("type:complete", headers)]
complete_ids    <- sapply(strsplit(complete_ids, " "), `[`, 1)

# Helper function to process each FiltDeTox category
process_full_length_category <- function(tsv_file, precursor_fasta_file, prefix,
                                         complete_ids, detox_fasta_label) {
  # Add full.length tag to TSV
  df_cat <- read.delim(tsv_file, sep = "\t", stringsAsFactors = FALSE)
  df_cat <- df_cat %>%
    mutate(full.length = ifelse(ID %in% complete_ids, "complete", "partial"))
  
  out_tsv <- paste0(prefix, "_full-length.tsv")
  write.table(df_cat, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Filter precursor FASTA to keep only complete sequences
  precursor_fasta     <- readAAStringSet(precursor_fasta_file)
  precursor_ids       <- names(precursor_fasta)
  precursor_clean_ids <- sapply(strsplit(precursor_ids, " "), `[`, 1)
  
  matching_ids   <- precursor_clean_ids %in% complete_ids
  filtered_fasta <- precursor_fasta[matching_ids]
  out_fasta      <- paste0(prefix, "_precursor_full-length-seqs.fasta")
  writeXStringSet(filtered_fasta, filepath = out_fasta)
  
  # Simple stats file for this category
  stats_file <- paste0(prefix, "_full-length_stats.txt")
  writeLines(
    c(
      paste("Number of complete sequences in", detox_fasta_label, ":",
            length(complete_ids)),
      paste("Number of complete sequences in", out_fasta, ":",
            length(filtered_fasta))
    ),
    con = stats_file
  )
}

# Apply full-length labelling to all four FiltDeTox categories
process_full_length_category(
  tsv_file            = "Toxins_Candidates.tsv",
  precursor_fasta_file= "Toxins_Candidates_precursor.fasta",
  prefix              = "Toxins_Candidates",
  complete_ids        = complete_ids,
  detox_fasta_label   = detox_fasta_label
)

process_full_length_category(
  tsv_file            = "Unlikely_Toxins.tsv",
  precursor_fasta_file= "Unlikely_Toxins_precursor.fasta",
  prefix              = "Unlikely_Toxins",
  complete_ids        = complete_ids,
  detox_fasta_label   = detox_fasta_label
)

process_full_length_category(
  tsv_file            = "SCRs_WA.tsv",
  precursor_fasta_file= "SCRs_WA_precursor.fasta",
  prefix              = "SCRs_WA",
  complete_ids        = complete_ids,
  detox_fasta_label   = detox_fasta_label
)

process_full_length_category(
  tsv_file            = "Non_Toxins.tsv",
  precursor_fasta_file= "Non_Toxins_precursor.fasta",
  prefix              = "Non_Toxins",
  complete_ids        = complete_ids,
  detox_fasta_label   = detox_fasta_label
)

#EDA Alicyclobacillaceae Genome MetaData
#BDS Capstone
#Lilly Hatch

#Load in libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(gt)
library(readr)

#Combining the NCBI genome metadata csvs
#Read in files
alicyclobacillaceaeGenomes1 <- read_csv("C:/Users/beaul/OneDrive/Desktop/Personal Documents/School/BDS Capstone/Data/Acidophiles/alicyclobacillaceaeGenomes1.csv")
alicyclobacillaceaeGenomes <- read_csv("C:/Users/beaul/OneDrive/Desktop/Personal Documents/School/BDS Capstone/Data/Acidophiles/AlicyclobacillacaeaGENOMES.csv")

#Assignment and Accession saved
assignmentLookup <- alicyclobacillaceaeGenomes1 %>%
  select(`Assembly Accession`, Assignment)

#Joined table with assignment and gene counts
allSpecies <- alicyclobacillaceaeGenomes %>%
  left_join(assignmentLookup, by = "Assembly Accession")


#Summary table
summary <- allSpecies %>%
  group_by(Assignment) %>%
  summarise(
    Total_Genomes = n(),
    Mean_Size_Mb = mean(Size) / 1000000,
    Median_Size_Mb = median(Size) / 1000000,
    Mean_Gene_Count = mean(`Gene Count`),
    Median_Gene_Count = median(`Gene Count`),
    Mean_N50_Mb = mean(`Contig N50`) / 1000000,
    Median_N50_Mb = median(`Contig N50`) / 1000000
    )
print(summary)

# Convert summary into fancy table
summary_table <- summary %>%
  gt() %>% 
  tab_header(title = "Summary of Selected Alicyclobacillaceae Species' Genomic Metadata",) %>%
  cols_label(
    Assignment = "Target Cohort",
    Total_Genomes = "Total Genomes",
    Mean_Size_Mb = "Mean Size (Mb)",
    Median_Size_Mb = "Median Size (Mb)",
    Mean_Gene_Count = "Mean Gene Count",
    Median_Gene_Count = "Median Gene Count",
    Mean_N50_Mb = "Mean N50 (Mb)",
    Median_N50_Mb = "Median N50 (Mb)") %>%
  # Format decimals for megabases 
  fmt_number(
    columns = c(Mean_Size_Mb, Median_Size_Mb, Mean_N50_Mb, Median_N50_Mb),
    decimals = 2) %>%
  tab_options(
    heading.title.font.size = px(18),
    heading.subtitle.font.size = px(13),
    column_labels.font.weight = "bold",
    table.border.top.color = "black",
    table.border.bottom.color = "black",
  )
summary_table

#Checking assembly level balance between groups
levelTable <- table(allSpecies$Level, allSpecies$Assignment, dnn = c("Level", "Assignment"))
addmargins(levelTable)
levelTable

# Boxplot of size
boxplot <- allSpecies %>% 
  ggplot(aes(x = Assignment, y = Size)) + 
  geom_boxplot() + 
  labs(
    title = "Size of Extremophile vs Mesophile Genomes", 
    x = "Assignment", 
    y = "Size"
  ) + 
  theme_bw()
boxplot

#Boxplot of gene count
geneCount_boxplot <- allSpecies %>% 
  ggplot(aes(x = Assignment, y = `Gene Count`)) + 
  geom_boxplot() + 
  labs(
    title = "Gene Count of Extremophile vs Mesophile Genomes", 
    x = "Assignment", 
    y = "Gene Count"
  ) + 
  theme_bw()
geneCount_boxplot

#scatterplot of size and gene count among both groups
scatterplot <- allSpecies %>%
  ggplot(aes(x = Size / 1000000, y = `Gene Count`, color = Assignment)) +
  geom_point() +
  labs(
    title = "Alicyclobacillaceae Scatter Plot of Size and Gene Count",
    x = "Size (Mb)",
    y = "Gene Count",
    color = "Assignment",
  ) +
  theme_bw()
scatterplot

pheno_mat <- allSpecies %>%
  select(`Assembly Accession`, Assignment) %>%
  rename(ids = `Assembly Accession`, pheno = Assignment)

# Assembly level distribution by Assignment
level_table_df <- allSpecies %>%
  count(Level, Assignment) %>%
  pivot_wider(names_from = Assignment, values_from = n, values_fill = 0)

level_gt_table <- level_table_df %>%
  gt() %>%
  tab_header(title = "Assembly Level Distribution of Selected Alicyclobacillaceae Genomes") %>%
  cols_label(
    Level = "Assembly Level"
  ) %>%
  tab_options(
    heading.title.font.size = px(18),
    column_labels.font.weight = "bold",
    table.border.top.color = "black",
    table.border.bottom.color = "black"
  )

level_gt_table

pangenome_data <- data.frame(
  Cohort = rep(c("Thermophile", "Acidophile"), each = 4),
  Category = rep(c("Core", "Soft Core", "Shell", "Cloud"), times = 2),
  Count = c(0, 0, 6619, 59504,      # Thermophile
            0, 0, 2454, 34918)      # Acidophile
)
pangenome_summary <- data.frame(
  Category = c("Core (99–100%)", "Soft Core (95–99%)", "Shell (15–95%)", "Cloud (0–15%)", "Total"),
  Thermophile = c(0, 0, 6619, 59504, 66123),
  Acidophile = c(0, 0, 2454, 34918, 37372)
)

pangenome_table <- pangenome_summary %>%
  gt() %>%
  tab_header(title = "Pangenome Gene Category Distribution by Cohort") %>%
  cols_label(
    Category = "Gene Category",
    Thermophile = "Thermophile Cohort",
    Acidophile = "Acidophile Cohort"
  ) %>%
  fmt_number(columns = c(Thermophile, Acidophile), decimals = 0, use_seps = TRUE) %>%
  tab_options(
    heading.title.font.size = px(18),
    column_labels.font.weight = "bold",
    table.border.top.color = "black",
    table.border.bottom.color = "black"
  )

pangenome_table


# Check the support/bootstrap values stored on the tree's internal nodes
tree$node.label

# Convert to numeric and check how many are below 75
support_vals <- as.numeric(tree$node.label)
sum(support_vals < 75, na.rm = TRUE)
length(support_vals)
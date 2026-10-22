#EDA Bacillaceae Genome MetaData
#BDS Capstone
#Lilly Hatch

#Load in libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(gt)

#Read in the Bacillaceae extremophile and mesophile dataset
extremos <- read.csv("C:/Users/beaul/OneDrive/Desktop/bacillaceaeExtremophiles.csv")
mesos <- read.csv("C:/Users/beaul/OneDrive/Desktop/bacillaceaeMesophiles.csv")

#Combining the two datasets into one
allSpecies <- bind_rows(extremos, mesos)

#Adding the Genus column data for all mesos because it was missed before and filling N/As
allSpecies <- allSpecies%>%
  mutate(across(where(is.character), ~na_if(., ""))) %>%
  mutate(Genus = tolower(str_extract(Organism.Scientific.Name, "^\\w+")))

#Summary table
summary <- allSpecies %>%
  group_by(Assignment) %>%
  summarise(
    Total_Genomes = n(),
    Distinct_Genera = n_distinct(Genus),
    Mean_Size_Mb = mean(Size) / 1000000,
    Median_Size_Mb = median(Size) / 1000000,
    Mean_Gene_Count = mean(Gene.Count),
    Median_Gene_Count = median(Gene.Count),
    Mean_N50_Mb = mean(Contig.N50) / 1000000,
    Median_N50_Mb = median(Contig.N50) / 1000000
  )
print(summary)

# Convert summary into fancy table
summary_table <- summary %>%
  gt() %>% 
  tab_header(title = "Summary of Selected Bacillaceae Species' Genomic Metadata",) %>%
  cols_label(
    Assignment = "Target Cohort",
    Total_Genomes = "Total Genomes",
    Distinct_Genera = "Distinct Genera",
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
  ggplot(aes(x = Assignment, y = Gene.Count)) + 
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
  ggplot(aes(x = Size / 1000000, y = Gene.Count, color = Assignment)) +
  geom_point() +
  labs(
    title = "Scatter Plot of Size and Gene Count",
    x = "Size (Mb)",
    y = "Gene Count",
    color = "Assignment",
  ) +
  theme_bw()
scatterplot

#Contingency table showing genera & extremophiles and mesophiles
ctable <- table(allSpecies$Genus, allSpecies$Assignment, dnn=c("Genus", "Assignment"))
addmargins(ctable)

ctable










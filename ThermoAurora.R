install.packages("reticulate")
install.packages(c("rlang", "devtools", "pryr"))
library(reticulate)
library(dplyr)
library(ggplot2)

reticulate::use_condaenv("auroraEnvironment", required = TRUE)

# Creating a conda environment with Python packages Aurora's ML functions need
conda_create("auroraEnvironment", packages = c("scikit-learn<1.6", "pandas<3.0", "numpy"))

devtools::install_github("hadley/pryr")
devtools::install_github("caitiecollins/treeWAS")
devtools::install_github("DalimilBujdos/aurora")
library(aurora)

save_dir_thermo <- "C:/Users/beaul/OneDrive/Documents/AuroraRun_Thermo/"
dir.create(save_dir_thermo, recursive = TRUE, showWarnings = FALSE)


bin_mat <- read.csv("//wsl$/Ubuntu/home/beaul/panaroo_output_thermo/gene_presence_absence.csv",
                    check.names = FALSE)
tree <- ape::read.tree("//wsl$/Ubuntu/home/beaul/panaroo_output_thermo/thermo_tree.newick")
pheno_mat <- allSpecies %>%
  select(Assembly.Accession, Assignment) %>%
  rename(ids = Assembly.Accession, pheno = Assignment)

colnames(bin_mat) <- sub("^(GC[AF]_[0-9]+\\.[0-9]+).*", "\\1", colnames(bin_mat))

#pre-filter check
cat("pheno_mat ids:", length(pheno_mat$ids), "\n")
cat("bin_mat genome columns:", ncol(bin_mat) - 3, "\n")
cat("tree tips:", length(tree$tip.label), "\n")

head(pheno_mat$ids)
head(colnames(bin_mat)[4:10])
head(tree$tip.label)

common_strains <- Reduce(intersect, list(pheno_mat$ids, colnames(bin_mat), tree$tip.label))
pheno_mat <- pheno_mat[pheno_mat$ids %in% common_strains, ]
tree <- ape::keep.tip(tree, common_strains)
bin_mat <- bin_mat[, c(1:3, which(colnames(bin_mat) %in% common_strains))]
gene_names <- bin_mat$Gene
bin_mat <- bin_mat[, -1:-3]          
bin_mat <- bin_mat %>%
  dplyr::mutate_all(~replace(., is.na(.) | . == "", 0)) %>%   
  dplyr::mutate_all(~replace(., . != 0, 1))
bin_mat <- sapply(bin_mat, as.numeric)
bin_mat <- t(bin_mat)
colnames(bin_mat) <- gene_names
bin_mat <- bin_mat[match(pheno_mat$ids, rownames(bin_mat)), ]

#shorten column names to avoid bug
gene_lookup_thermo <- data.frame(
  short_id      = as.character(seq_len(ncol(bin_mat))),
  original_name = colnames(bin_mat),
  stringsAsFactors = FALSE
)
colnames(bin_mat) <- gene_lookup_thermo$short_id

conda_path <- conda_list()[conda_list()$name == "auroraEnvironment", "python"]

#Checks before aurora_pheno
dim(bin_mat)
identical(rownames(bin_mat), pheno_mat$ids)
sum(is.na(bin_mat))
length(unique(gene_lookup_thermo$short_id)) == ncol(bin_mat)

thermo_results <- aurora_pheno(pheno_mat = pheno_mat,
                              bin_mat = bin_mat,
                              type_bin_mat = "custom",
                              tree = tree,
                              condaenv_path = conda_path,
                              save_dir = save_dir_thermo,
                              hamming_filter = FALSE)
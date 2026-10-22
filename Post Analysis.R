#Stats from aurora
total_variants_thermo <- ncol(bin_mat)
total_variants_thermo
thermo_results$abundance_filter_removed
total_variants_thermo - thermo_results$abundance_filter_removed
thermo_results$acr_filter_removed
sum(!is.na(gwas_table$`chisq p.val`))

#checking aurora derived bonferroni threshold against all the genes
thermo_results$bonferroni_threshold
n_significant_thermo <- sum(gwas_table$`chisq p.val` < thermo_results$bonferroni_threshold, na.rm = TRUE)
n_significant_thermo

# Top 10 highest F1 score
exploratory_genes_thermo <- gwas_table %>%
  filter(abs(`Extremophile F1` - 2/3) > 1e-6) %>%
  arrange(desc(`Extremophile F1`)) %>%
  slice_head(n = 10) %>%
  pull(Variant)

exploratory_genes_thermo


#top genes table creation
top_genes_gt <- top_genes_table %>%
  select(Variant, `Extremophile ratio`, `Mesophile ratio`, 
         `Extremophile std_residual`, `Extremophile F1`, Annotation, Functional_Category) %>%
  gt() %>%
  tab_header(title = "Top 10 Candidate Genes by F1 Score: Thermophile Cohort") %>%
  cols_label(
    Variant = "Gene",
    `Extremophile ratio` = "Extremophile Ratio",
    `Mesophile ratio` = "Mesophile Ratio",
    `Extremophile std_residual` = "Std. Residual",
    `Extremophile F1` = "F1 Score",
    Annotation = "Prokka Annotation",
    Functional_Category = "Functional Category"
  ) %>%
  fmt_number(columns = c(`Extremophile std_residual`, `Extremophile F1`), decimals = 3) %>%
  tab_options(
    heading.title.font.size = px(18),
    column_labels.font.weight = "bold",
    table.border.top.color = "black",
    table.border.bottom.color = "black"
  )

top_genes_gt

#Assigning functional category to genes
background_annotations <- raw_thermo$Annotation
topTenCategorized <- sapply(exploratory_genes_thermo, get_annotation)
topTenCategorized

keyword_categories <- list(
  "Transport"             = "transport|permease|ABC transporter|translocase",
  "DNA/RNA Processing"    = "polymerase|helicase|ligase|nuclease",
  "Stress Response"       = "heat shock|chaperone|universal stress|cold shock|repressor",
  "Metabolism"            = "dehydrogenase|synthase|kinase|reductase",
  "Cell Wall/Membrane"    = "membrane|peptidoglycan|cell wall",
  "Translation/Ribosomal" = "ribosomal|elongation factor|translation|tRNA"
)
enrichment_results <- lapply(names(keyword_categories), function(cat) {
  pattern <- keyword_categories[[cat]]
  sig_hits <- sum(str_detect(topTenCategorized, regex(pattern, ignore_case = TRUE)), na.rm = TRUE)
  bg_hits  <- sum(str_detect(background_annotations, regex(pattern, ignore_case = TRUE)), na.rm = TRUE)
  n_sig <- length(topTenCategorized)
  n_bg  <- length(background_annotations)
  p_val <- phyper(sig_hits - 1, bg_hits, n_bg - bg_hits, n_sig, lower.tail = FALSE)
  data.frame(Category = cat, Number_Genes_Categorized = sig_hits, Background_Count = bg_hits,
             Total_Genes = n_sig, Background_Total = n_bg, P_Value = p_val)
})
enrichment_table <- do.call(rbind, enrichment_results)
enrichment_table
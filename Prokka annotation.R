raw_thermo <- read.csv("//wsl$/Ubuntu/home/beaul/panaroo_output_thermo/gene_presence_absence.csv",
                       check.names = FALSE)

significant_genes <- c("ctsR", "mcsB", "dnaK")

# Match against the "Gene" column (handle the ~~~ paralog-merged naming too)
sig_annotations <- raw_thermo %>%
  filter(sapply(strsplit(as.character(Gene), "~~~"), function(x) any(x %in% significant_genes))) %>%
  pull(Annotation)

background_annotations <- raw_thermo$Annotation

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
  sig_hits <- sum(str_detect(sig_annotations, regex(pattern, ignore_case = TRUE)), na.rm = TRUE)
  bg_hits  <- sum(str_detect(background_annotations, regex(pattern, ignore_case = TRUE)), na.rm = TRUE)
  n_sig <- length(sig_annotations)
  n_bg  <- length(background_annotations)
  p_val <- phyper(sig_hits - 1, bg_hits, n_bg - bg_hits, n_sig, lower.tail = FALSE)
  data.frame(Category = cat, Sig_Count = sig_hits, Background_Count = bg_hits,
             Sig_Total = n_sig, Background_Total = n_bg, P_Value = p_val)
})

enrichment_table <- do.call(rbind, enrichment_results)
enrichment_table

enrichment_gt <- enrichment_table %>%
  gt() %>%
  tab_header(title = "Functional Category Analysis of Significant Genes: Thermophile Cohort") %>%
  cols_label(
    Category = "Functional Category",
    Sig_Count = "Significant Genes",
    Background_Count = "Background Genes",
    Sig_Total = "Total Significant",
    Background_Total = "Total Background",
    P_Value = "P-value"
  ) %>%
  fmt_number(columns = P_Value, decimals = 6) %>%
  tab_options(
    heading.title.font.size = px(18),
    column_labels.font.weight = "bold",
    table.border.top.color = "black",
    table.border.bottom.color = "black"
  )

enrichment_gt
# Genetic Plasticity in Extreme Environments: Comparative Genomic Analysis Using aurora

This repository contains the data selection, processing, and analysis pipelines for my Master of Science in Biological Data Science capstone project at Arizona State University. 

The project applies `aurora`, a machine learning-based microbial genome-wide association study (mGWAS) tool, to identify habitat-associated genomic regions in extremophiles while controlling for phylogenetic confounding.

### Overview
Extremophiles offer unique insights into evolutionary adaptation under strong selective pressure. This project shifts the focus from individual species adaptation to a broader pathway-convergent perspective. By analyzing two independent extremophile cohorts, the goal was to determine if distantly related lineages undergo convergent genomic change in similar functional pathways.

**Cohorts Analyzed:**
1. **Thermophiles:** Thermophilic and mesophilic *Bacillaceae* (n=32)
2. **Acidophiles:** Acidophilic and non-acidophilic *Alicyclobacillaceae* (n=41)

### Methodology & Pipeline
The analysis pipeline involves several stages of data cleaning, genome annotation, and machine learning models (Random Forest and Random Walk):

1. **Genome Selection (Python):** 
   - Automated querying and cross-referencing between ThermoBase, AciDB, and NCBI RefSeq.
   - Filtered for genome completeness (≥ 95%) and contamination (< 5%).
   - Programmatically selected closely related mesophilic controls within the same family/genus to isolate genomic differences reflecting habitat rather than broad evolutionary divergence.
2. **Pre-Processing & Annotation (Unix/Bash):**
   - **Prokka:** Annotated raw FASTA genome sequences to define gene regions.
   - **Panaroo:** Compiled sequence data into gene presence/absence matrices and alignments.
   - **FastTree:** Built phylogenetic trees to map lineage effects.
3. **mGWAS Analysis (R / `aurora`):**
   - Executed models to separate **Locus Effects** (genes functionally relevant to habitat adaptation) from **Lineage Effects** (genes inherited alongside adaptation but not causative).

### Results & Conclusion
While no single gene family reached Bonferroni-corrected significance, the top-ranked genes identified in the thermophile cohort (such as *ctsR* and *mcsB*) are known regulators of heat-shock stress response. 

More importantly, the results suggest that the habitat adaptation signal may be biologically real but is likely diffuse, polygenic, or masked by genome assembly quality and taxonomic breadth. 

A major conclusion of this project highlights the methodological challenges of conducting mGWAS on non-model extremophiles. Because published extremophile sequence data is so limited, there is a tricky balance between maximizing sample size and keeping cohorts closely related enough to ensure genomic differences reflect habitat rather than broad evolutionary divergence. 

*For a deeper dive into the methodology and challenges of building this pipeline, see my included Capstone Report and Presentation.*

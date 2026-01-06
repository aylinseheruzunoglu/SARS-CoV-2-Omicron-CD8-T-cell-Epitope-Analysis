# SARS-CoV-2 CD8+ T Cell Epitope Evolution & Conservation Analysis

## Project Overview
This repository contains a comprehensive bioinformatics suite designed to analyze the evolutionary trajectory of SARS-CoV-2 CD8+ T cell epitopes. The project evaluates the shift in the immunogenic landscape from the ancestral strain (Wuhan-Hu-1) through the emergence of recent Omicron sub-variants (BA.2.86, JN.1, LP.8.1, NB.1.8.1, and XEC).

The analysis is divided into five integrated pipelines that quantify epitope conservation, binding affinity shifts, and structural attrition across the viral Spike protein.

---

## Pipeline 1: Comparative Epitope & Mutational Landscape Analysis

### Description
This pipeline serves as the foundational comparative module of the project. It evaluates the shift in the **CD8+ T cell** immunogenic landscape from the ancestral strain to emergent Omicron lineages. The primary objective is to quantify how variant-defining mutations alter the physical and positional distribution of cytotoxic T lymphocyte (CTL) targets.

### Key Analytical Modules
* **Ancestral vs. Variant Preservation:** Quantifies the percentage of epitopes retained from the ancestral strain versus those lost due to amino acid substitutions.
* **Cross-Variant Overlap:** Identifies "hotspot" epitopes shared across multiple Omicron sub-variants to determine common evolutionary trajectories.
* **Mutation Impact Mapping:** Specifically analyzes how lineage-defining mutations disrupt the MHC Class I binding motifs of previously characterized ancestral epitopes.
* **De Novo Epitope Identification:** Isolates variant-specific epitopes created by new mutations that are absent in the ancestral sequence.
* **Topological Mapping:** Provides positional coordinates of conserved and mutated epitopes across the viral Spike protein.

---

## Pipeline 2: Quantitative Binding Affinity & Stability Analysis

### Description
This pipeline performs a comparative quantitative analysis of **HLA Class I** binding affinities. By calculating the Fold Change (FC) in elution percentile scores between ancestral and variant peptide–allele pairs, it categorizes the impact of mutations into three distinct classes: Increased, Decreased, or Unchanged binding strength.

### Key Analytical Modules
* **Peptide–Allele Pair Matching:** Specifically isolates unique pairs of peptide sequences and their corresponding HLA alleles to ensure a 1:1 comparison across variants.
* **Log2 Fold Change Calculation:** Normalizes the **NetMHCpan 4.1** percentile scores to visualize the magnitude of affinity shifts using a logarithmic scale.
* **Binding Strength Classification:**
    * *Increased:* Mutations that lower the percentile score (improving affinity).
    * *Decreased:* Mutations that raise the percentile score (weakening affinity).
    * *Unchanged:* Mutations that preserve the original binding properties.
* **Conservation Core Analysis:** Identifies the "Universal Preserved Core"—the subset of peptide–allele pairs that remain strong binders across every single variant analyzed (BA.2.86 through XEC).
* **Statistical Visualization:** Generates distribution histograms and comparative bar charts to provide a global view of the mutational impact on the immunome.

---

## Pipeline 3: Multi-Variant Intersection & Visual Analytics

### Description
This pipeline focuses on the collective immunological landscape of the Omicron sub-lineages. It integrates data from all analyzed variants to identify patterns of conservation and divergence using advanced set-theory visualizations. This is the primary engine for determining how much of the **CD8+ T cell** "target space" is shared across the evolutionary branch.

### Key Analytical Modules
* **Intersection Analysis (UpSet Plots):** Visualizes complex intersections of peptide–HLA pairs across five variants, identifying "universal" epitopes versus those unique to specific lineages.
* **Allele-Specific Distribution:** Generates per-variant histograms to reveal if certain **HLA Class I** alleles are more "mutation-resilient" than others.
* **Mutation Overlap Mapping:** Quantifies the exact percentage of the epitope–HLA repertoire that is physically impacted by variant-defining mutations.
* **Epitope Fate Tracking:** Categorizes the "life cycle" of an epitope into four sets: *Conserved, Unique, Lost,* and *Shared*.
* **Heatmap of Mutational Burden:** Provides a high-level matrix view of how mutations are distributed across peptide–HLA pairs.
* **Proportional Impact Assessment:** Uses pie charts to visualize the ratio of "Affected" vs. "Unaffected" epitopes.

---

## Pipeline 4: Lineage-Specific Epitope Attrition & Spike Domain Mapping

### Description
This pipeline shifts focus to a spatial and evolutionary analysis centered on the Spike (S) protein. It specifically investigates the "epitope decay" occurring as SARS-CoV-2 transitioned from BA.2.86/JN.1 into the most recent sub-variants. It maps lost epitopes to functional domains (e.g., RBD, NTD) to check for clustering in high-impact regions.

### Key Analytical Modules
* **Variant-Based Loss Quantification:** A comparative statistical analysis of peptide–HLA pairs viable in ancestral/early-Omicron strains but abolished in contemporary variants.
* **Evolutionary Attrition Tracking:** A heatmap analysis identifying "vulnerable" epitopes—those preserved in early transitions but systematically lost in subsequent lineages.
* **Structural Positional Mapping:** Uses JN.1 as a reference sequence to map the exact amino acid coordinates of the "disappearing immunome."
* **Spike Domain Localization:** Categorizes lost epitopes into structural sub-regions:
    * *S1 Subunit:* (NTD, RBD, SD1, SD2)
    * *S2 Subunit:* (Proteolytic cleavage sites, Fusion peptide)
* **Domain-Specific Summary Plots:** Visualizes which Spike domains are undergoing the fastest rate of **CD8+ T cell** epitope depletion.

---

## Pipeline 5: Conserved Core Immunome & HLA Locus Distribution

### Description
This final pipeline serves as the synthesis module, isolating the "Universal Conserved Core"—peptide–HLA pairs remaining functionally intact across the ancestral strain and all five Omicron sub-variants. It also provides a genetic perspective by grouping results by **HLA Locus (A, B, C)**.

### Key Analytical Modules
* **Strict Conservation Profiling:** Identifies epitope–HLA pairs that exist identically across all datasets regardless of binding score.
* **Affinity-Filtered Stability Analysis:** Implements stringent **NetMHCpan** percentile thresholds (e.g., <0.5%, <2%) to isolate only the most clinically relevant conserved epitopes.
* **Multi-Locus Distribution (Locus-wise):** Aggregates data at the HLA-locus level (**HLA-A, HLA-B, HLA-C**) to identify which loci drive the most stable **CD8+ T cell** responses.
* **Allelic Resilience Mapping:** Visualizes the count of conserved epitopes per specific allele.
* **Comparative Threshold Heatmaps:** Generates a matrix view of conservation levels across different binding strengths.

---

## Requirements & Dependencies
* Python 3.8+
* Pandas / NumPy
* Matplotlib / Seaborn (for static plots)
* UpSetPlot (for intersection analysis)
* NetMHCpan 4.1 (External tool or API)

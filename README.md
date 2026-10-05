# Human Microbiome Differential Abundance Workflows & Meta-Analyses

An R-based bioinformatics repository containing comprehensive statistical pipelines for identifying disease-associated microbial biomarkers and synthesizing cross-study consensus signatures. This portfolio leverages **curatedMetagenomicData** for primary shotgun metagenomic cohort profiling and **BugSigDB** for multi-center landscape evidence integration.

---

## 📂 Repository Architecture

```text
microbiome-differential-abundance/
├── README.md                          # Main repository landing page
├── cohort-level-metagenomics/          # curatedMetagenomicData workflows
│   ├── crc_gut_microbiome.Rmd         # Colorectal cancer progression staging
│   ├── delivery_mode_infants.Rmd      # Vaginal vs. C-section profiling
│   └── preterm_severity_analysis.Rmd  # Continuous gestational age modeling
└── cross_study_meta_analysis/        # BugSigDB Consensus Signatures workflows
    ├── general_womens_health.Rmd      # Broad landscape view across conditions
    ├── breast_cancer_vignette.Rmd     # Fecal/tissue profiling (gut-breast axis)
    ├── endometriosis_vignette.Rmd     # Focused condition-specific mapping
    └── pregnancy_maternal_child.Rmd   # Vaginal/gut shifts & vertical transmission
```

---

## 🔬 Workflow Summaries

### 1. Cohort-Level Metagenomics (`curatedMetagenomicData`)
These workflows evaluate microbial shifts directly from primary shotgun sequencing cohorts, handling both discrete experimental categories and continuous biological parameters:

*   **Colorectal Cancer Progression Staging (`ZellerG_2014`):**
    *   **Objective:** Tri-group differential testing comparing healthy controls (n=61), colorectal adenomas (n=42), and full-stage CRC (n=53).
    *   **Methods:** Runs count-based modeling using **ANCOM-BC2** alongside relative abundance-based biomarker discovery with **LEfSe**.
    *   **Key Finding:** Core disruptions overlap across multiple specific taxa (e.g., *Fusobacterium*, *Escherichia coli*, *Bacteroides* enriched in CRC), while whole-community alpha/beta diversity metrics remain broadly overlapping.
*   **Infant Delivery Mode Profiling (`WampachL_2018`):**
    *   **Objective:** Discrete abundance analysis comparing vaginal (n=7) and C-section (n=9) delivery pathways.
    *   **Methods:** Benchmarks five foundational differential ecology frameworks in parallel: **ANCOM-BC2**, **LEfSe**, **MAaSLin2**, **DESeq2**, and **Corncob**.
    *   **Statistical Control:** Restricts matrices to singular baseline observations per individual using early timeline filters to prevent longitudinal pseudo-replication.
*   **Continuous Gestational Age Severity (`BrooksB_2017`):**
    *   **Objective:** Modeling continuous preterm infant maturation scales against gestational age.
    *   **Methods:** Employs continuous linear predictors through **MAaSLin2** and **ANCOM-BC2**. 
    *   *Note:* LEfSe is explicitly omitted as it cannot evaluate non-discrete continuous predictors.

### 2. Cross-Study Meta-Analyses (`BugSigDB`)
These sub-analyses compile curated, peer-reviewed clinical signatures using **bugsigdbr** to assess global ecosystem concordance and calculate statistical enrichment thresholds:

*   **General Landscape View (n=580 signatures):**
    *   Integrates broad multi-condition datasets across reproductive, metabolic, and oncological profiles.
    *   Reveals a clear biological split between **vaginal-health communities** (*Lactobacillus*, *Gardnerella*, *Prevotella*) and **gut-health communities** (*Bacteroides*, *Faecalibacterium*, *Blautia*) responding to localized pathologies.
*   **Endometriosis (n=47 signatures):**
    *   A micro-scale condition mapping showing structural tissue profiles across anatomical body sites (feces, endometrium, cervical cavity).
*   **Breast Cancer (n=38 signatures):**
    *   Evaluates mucosal, salivary, and fecal biomarkers to map systemic modifications along the **gut-breast axis**.
*   **Pregnancy & Maternal-Child Health (n=248 signatures):**
    *   Surveys maternal gestational parameters (e.g., gestational diabetes), preterm risks, and neonatal transmission patterns stratified by delivery and feeding mechanics.

---

## 🛠️ Methodological & Statistical Checkpoints

To ensure exact reproducibility, the codebase enforces strict data structure rules:
1.  **Scale Alignment:** Enforces scale requirements per package—using raw estimated count matrices for count-based models (`ANCOM-BC2`, `DESeq2`) vs. parts-per-million scaling (1,000,000 column sums) via `lefser:::relativeAb()` to resolve execution crashes in relative abundance pipelines.
2.  **Temporal Filtering:** Leverages longitudinal tracking metrics (`days_from_first_collection`) to filter sample matrices down to unique single-subject entries, eliminating sample inflation.
3.  **Core Consensus Thresholding:** Sets standard filters requiring a taxon to be reported across ≥ 3 independent studies with ≥ 75% directional concordance before accepting it into a true "core consensus signature". For underpowered cohort horizons (e.g., Endometriosis), thresholds are transparently relaxed to ≥ 2 studies and ≥ 70% concordance.
4.  **Permutation Controls:** Utilizes Monte Carlo simulation engines via `getCriticalN()` to determine maximum background frequency distributions under empirical null models, checking them against real dataset frequencies at α = 0.05.

---

## 📦 Core Dependencies & Bioconductor Stack

To run these vignettes locally, install the required microbiome ecology stack within your R environment:

```r
# Initialize Bioconductor manager
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install core primary data and signature engines
BiocManager::install(c("curatedMetagenomicData", "bugsigdbr", "mia", "scater"))

# Install differential abundance frameworks
BiocManager::install(c("ANCOMBC", "lefser", "DESeq2"))
install.packages(c("Maaslin2", "corncob"))

# Install validation and visualization dependencies
BiocManager::install(c("waldronlab/bugSigSimple", "waldronlab/BugSigDBStats", "ComplexHeatmap"))
install.packages(c("tidyverse", "vegan", "kableExtra"))
```

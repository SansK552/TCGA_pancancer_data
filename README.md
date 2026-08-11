# TCGA Pan-Cancer Causal Gene Regulatory Network Atlas

**Identifying Master Regulatory Transcription Factors Across 33 Cancer Types using Machine Learning and Causal Inference**

This is a pan-cancer gene regulatory network (GRN) project built as part of my undergraduate research, where I applied machine learning and causal inference methods to bulk RNA-seq data from The Cancer Genome Atlas (TCGA). It covers 33 cancer types and around 10,000 patient samples with the goal of identifying which transcription factors act as "master regulators" behind gene expression in each cancer type causally, not just as a correlation.

The final results are packaged into a small interactive web portal so they're actually browsable, instead of sitting as raw output files across dozens of folders.

**🔗 Live Portal:** [sansk552.github.io/TCGA_pancancer_data](https://sansk552.github.io/TCGA_pancancer_data/)

---

### Why This Project

Cancer isn't driven only by individual mutated genes, a lot of it comes down to dysregulated **gene regulatory networks (GRNs)**, the transcription factor (TF) circuits that decide which genes get switched on or off and how that control breaks down in tumors. The question I wanted to explore was:

> *Which transcription factors act as "master regulators" driving gene expression changes in each cancer type and can this be shown causally, instead of just correlationally?*

This repository holds the consolidated results of that analysis across all major TCGA cancer types, along with the portal I built to explore them.

---

### What's in Here

- **33 TCGA cancer types analyzed** end-to-end (ACC, BLCA, BRCA, CESC, CHOL, COAD, DLBC, ESCA, GBM, HNSC, KICH, KIRC, KIRP, LAML, LGG, LIHC, LUAD, LUSC, MESO, OV, PAAD, PCPG, PRAD, READ, SARC, SKCM, STAD, TGCT, THCA, THYM, UCEC, UCS, UVM)
- **~10,000 tumor samples** processed through a reproducible RNA-seq preprocessing and network-inference pipeline
- **Regulatory network inference** using gradient boosting based methods (GRNBoost2) to build TF-target gene modules
- **Unsupervised clustering** (HDBSCAN, Louvain, K-Means, hierarchical) used to define patient/tumor states before causal analysis
- **Causal inference** via the MINER3 framework, used to separate genuine "master regulators" from the broader set of candidate TFs
- **A working portal** — a simple, self contained web page where you can search TFs, browse results by cancer type and compare across cancers

---

### How the Pipeline Works

#### 1. Data Acquisition & Preprocessing
- Raw RNA-seq count data pulled via recount3 (Bioconductor/R) for all 33 TCGA cancer cohorts
- Sample filtering: Primary Tumor samples prioritized, with a fallback to all non-normal tumor samples when a cohort had fewer than 50 samples; melanoma (SKCM) intentionally kept both Primary Tumor and Metastatic samples since that's biologically relevant for that cancer
- ID cleanup: Ensembl gene ID version-stripping → gene symbol mapping via `mygene` (in batches of 500) → duplicate collapsing by mean expression
- Normalization: log2(CPM + 1) transformation followed by gene wise Z-score standardization
- Final processed dataset: 10,507 samples across 33 cancer types

#### 2. Gene Regulatory Network Inference
- **pySCENIC** used for regulon inference, with **GRNBoost2** (a gradient-boosting regression method) inferring TF–target gene modules from the expression matrices
- For smaller cohorts, regulons were built directly from the GRNBoost2 modules skipping the motif-enrichment pruning step where sample size made it unreliable
- Regulon activity scored across all patient samples

#### 3. Causal Inference & Master Regulator Identification
- **MINER3** applied on top of the regulon activity matrices to move past correlation and identify actual causal drivers
- Patient states clustered using five different unsupervised methods (**HDBSCAN, Louvain, K-Means, hierarchical clustering, MINER F1**) before running causal analysis, to make sure the tumor states used weren't an artifact of one clustering choice
- Regulators are grouped into three tiers:
  - `master_attractors` — the broad candidate set
  - `selected_activators` — regulators with causal support ≥ 0.5 (fraction of impacted regulons aligned)
  - `master_regulators_activators` — the apex-node regulators with the strongest causal evidence
- The same causal logic is applied on the repressor side `master_repressors` (candidate set) and `selected_repressors` (causal support ≥ 0.5), capturing TFs that actively suppress transcriptional programs rather than drive them, which matters just as much for understanding what's holding a cancer state in place.
- The same threshold parameters were used across all 33 cancers so results are comparable cohort-to-cohort. Cancers with zero identified activators/repressors are reported as it is.

#### 4. Building the Portal
- All per-cancer results assembled into one dataset via a Python script (`assemble_portal_data.py`), producing `portal_data.json` plus per-cancer regulon-target files
- Four views: **TF Search**, **By Cancer**, **Pan-Cancer Overview** and **Cancer Summary**

---

### Repository Structure

```
TCGA_pancancer_data/
├── index.html                     # The interactive portal (single-page app)
├── portal_data.json               # Consolidated pan-cancer results (all 33 cancers)
├── portal_regulon_targets/        # Per-cancer TF → target gene lookup files
├── TCGA_ACC_files/                # Per-cancer-type results
├── TCGA_BLCA_files/
├── TCGA_BRCA_files/
│   ...                            # (33 per-cancer-type result directories)
├── TCGA_UVM_files/
└── .gitignore                     # Excludes large intermediate files (raw Z-score matrices, checkpoints)
```

---

### Getting Started

No installation needed to just look around:

```bash
# View live
open https://sansk552.github.io/TCGA_pancancer_data/
```

---

### What You Can Do With It

- **Search a transcription factor** (e.g. `MYC`, `TP53`, `ASCL1`) to see which cancer types it comes up as a master regulator in
- **Browse a specific cancer type** to see its full regulon network and top driver TFs
- **Compare across cancers** to spot which regulatory patterns are pan-cancer versus tissue-specific
- Use it as a **starting point for hypothesis generation** — flagging candidate TFs worth following up on experimentally

---

### Background

This project builds on the **gbmMINER** causal network framework and the **pySCENIC** regulon-inference toolkit, extending them to a pan-cancer scale across the full TCGA cohort instead of a single disease. It's part of ongoing research in the Department of Biochemical Engineering & Biotechnology at IIT Delhi, done under faculty supervision.

---

### Limitations

- This uses **bulk RNA-seq** not single-cell data, so cell-type-specific signals within a tumor are averaged out
- A few low-sample-size cohorts (e.g. DLBC, UVM) have weaker statistical power for causal inference
- Regulons with zero identified activators/repressors are kept in the results as genuine findings, not filtered out
- This is a research and hypothesis-generation tool, not a clinical or diagnostic one

---

### Author

**Sanskruti Krishnapurkar** — B.Tech Biotechnology, IIT Delhi (Dept. of Biochemical Engineering & Biotechnology)
Done as part of undergraduate research work in computational cancer genomics.

---

### Acknowledgments

- **The Cancer Genome Atlas (TCGA)** for the publicly available cancer data
- **pySCENIC** and **MINER3**, the open-source frameworks this pipeline is built on
- **recount3** (Bioconductor) for standardized RNA-seq count access

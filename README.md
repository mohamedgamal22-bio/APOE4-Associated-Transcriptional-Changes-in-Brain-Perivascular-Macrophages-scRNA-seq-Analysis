# Igf1-Regulated-BBB-through-PVMs

**Single-cell RNA-seq and pseudobulk differential expression reanalysis characterizing APOE genotype-associated transcriptional changes in mouse brain perivascular macrophages (PVMs).**

---

## Overview

Perivascular macrophages (PVMs) maintain cerebrovascular and blood–brain barrier (BBB) homeostasis through a cMAF-dependent transcriptional program that includes *Igf1*. In this project, we reanalyzed published mouse brain macrophage single-cell RNA-sequencing (scRNA-seq) data to evaluate *APOE4*- vs. *APOE3*-associated transcriptional alterations in PVM-enriched populations, focusing on *Igf1* and core PVM identity markers.

---

## Study Design & Workflow

* **Samples:** 6 biological samples total (`3 APOE3` vs. `3 APOE4` mouse brain macrophage datasets)
* **Cell Population:** **988** PVM-enriched cells (`400 APOE3` and `588 APOE4`)
* **Annotation Strategy:**
  * **Positive PVM markers:** *Folr2*, *Cd38*, *Lyve1*, *Mrc1*, *Cd163*
  * **Negative / low MHC-II markers:** *H2-Aa*, *H2-Ab1*, *H2-Eb1*, *Cd74*
* **Differential Expression Framework:** Sample-level raw count aggregation (pseudobulk) within the PVM-enriched cluster followed by **DESeq2** (`APOE4` vs. `APOE3`, significance threshold: `adjusted P < 0.05`).

---

## Repository Structure

The repository is organized sequentially according to the analytical pipeline:

```text
Igf1-Regulated-BBB-through-PVMs/
├── QC/                          # Quality control & doublet filtering scripts (e.g., doublets.filter.R)
├── Normalization and merging/   # Sample normalization & dataset merging workflows (e.g., Normalization.R)
├── Dimensionality reduction/    # Feature selection, PCA, Harmony batch integration & UMAP
├── Clustring/                   # Iterative clustering & PVM population marker annotation
├── Pseudobulk/                  # Sample-level count aggregation & DESeq2 differential expression
└── README.md                    # Project documentation

# FDCA: Single-Cell Metabolic Flux Landscapes in Breast Cancer Epithelium

## Overview
This project implements a **single-cell systems biology framework** combining:

- **scRNA-seq (10X Genomics)**
- **Automatic cell-type annotation (SingleR)**
- **Epithelial tumor vs normal stratification**
- **Flux-Driven Cellular Analysis (FDCA)**
- **Flux Balance Cellular Automata (FBCA)**

to decode **metabolic reprogramming in breast cancer epithelial cells** at single-cell and cluster levels.

The workflow reveals **heterogeneous metabolic states**, highlighting tumor-specific flux rewiring across glycolysis, oxidative phosphorylation, TCA cycle, and anabolic pathways.

---

## Biological Motivation

Cancer cells rewire metabolism to:
- Sustain uncontrolled proliferation
- Survive hypoxia and nutrient stress
- Support invasion and EMT

Traditional bulk RNA-seq masks this heterogeneity.  
**Single-cell metabolic flux inference** allows us to:

✔ Track metabolic plasticity  
✔ Separate tumor vs normal epithelium  
✔ Quantify pathway-level flux states  
✔ Model emergent tumor metabolic ecosystems  

---

## Key Concepts

### 🔹 FDCA (Flux-Driven Cellular Analysis)
A pathway-level scoring approach where **gene expression programs are translated into metabolic flux potentials** at the single-cell level.

---

## Dataset
- **Technology:** 10X Genomics scRNA-seq
- **Tissue:** Breast cancer (Invasive Ductal Carcinoma)
- **Input:** `filtered_feature_bc_matrix`
- **Species:** Human

---

## Software Requirements

### R Packages
```r
Seurat
tidyverse
patchwork
SingleR
celldex

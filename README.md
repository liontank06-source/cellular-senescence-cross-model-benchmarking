# Cross-Model Transcriptomic Benchmarking of Candidate Cellular Senescence Markers

Computational analysis of candidate cellular senescence markers across independent transcriptomic datasets, cell types, and senescence induction methods.

The project evaluates whether commonly used senescence-associated transcripts behave consistently across biological contexts or whether their direction and statistical support depend strongly on the experimental model.

Five pathway-prioritized candidate genes were evaluated:

**FOS, IL1A, CXCL8, CDKN2B, and ICAM1**

---

## Project Overview

Cellular senescence is commonly assessed using markers associated with cell-cycle arrest, inflammatory signaling, and the senescence-associated secretory phenotype (SASP). However, individual markers may behave differently depending on cell type, senescence trigger, and experimental design.

This project uses a multi-stage computational workflow to evaluate candidate-marker reproducibility across independent transcriptomic datasets.

The analysis includes:

1. Discovery of senescence-associated genes in BJ fibroblasts.
2. Pathway enrichment and candidate prioritization.
3. Independent validation in WI-38 replicative senescence.
4. Cross-model differential-expression analysis across eight senescence comparisons.
5. Random-effects meta-analysis.
6. Hartung-Knapp-Sidik-Jonkman uncertainty estimation.
7. Leave-one-model-out and shared-control sensitivity analyses.
8. Candidate reliability benchmarking.

---

## Research Question

**Do individual senescence-associated transcripts show reproducible directional behavior across different cellular and experimental models of senescence?**

The analysis specifically tests whether relatively stable growth-arrest-associated markers such as **CDKN2B** show greater directional robustness than inflammatory, secretory, adhesion-associated, and stress-response transcripts.

---

## Datasets

| Dataset / Stage | Biological Context | Repository Input | Purpose |
|---|---|---|---|
| BJ fibroblast discovery | Young vs. senescent BJ fibroblasts | GEO2R differential-expression table | Candidate discovery |
| GSE175533 | WI-38 replicative senescence | Processed endpoint and time-course candidate results | Independent validation |
| GSE130727 | Endothelial and fibroblast senescence models | Gene count matrix and matched comparison design | Cross-model benchmarking |

GSE130727 contains eight matched senescence comparisons spanning multiple cell types and induction methods, including ionizing radiation, replicative senescence, doxorubicin, and HRAS-induced senescence.

---

## Analysis Workflow

```text
BJ fibroblast discovery
        ↓
Differential-expression filtering
        ↓
GO / Reactome pathway prioritization
        ↓
Five candidate genes
        ↓
GSE175533 independent validation
        ↓
GSE130727 eight-model PyDESeq2 analysis
        ↓
Random-effects meta-analysis
        ↓
HKSJ confidence intervals
        ↓
Leave-one-model-out sensitivity analysis
        ↓
Candidate reliability benchmarking

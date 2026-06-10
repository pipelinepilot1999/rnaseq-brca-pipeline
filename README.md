# RNA-seq Differential Expression Pipeline
### Luminal A vs Triple Negative Breast Cancer (TNBC) | TCGA-BRCA

[![Snakemake](https://img.shields.io/badge/snakemake-9.22.0-brightgreen)](https://snakemake.readthedocs.io)
[![AWS](https://img.shields.io/badge/cloud-AWS%20EC2-orange)](https://aws.amazon.com)
[![Docker](https://img.shields.io/badge/container-Docker-blue)](https://docker.com)

## Overview

End-to-end reproducible RNA-seq pipeline for differential gene expression analysis, built with Snakemake, deployed on AWS EC2, and validated on mouse mammary gland data (GSE60450) before scaling to TCGA-BRCA human breast cancer data.

**Biological Question:** What genes are differentially expressed between Luminal A and Triple Negative Breast Cancer (TNBC) subtypes?

**Clinical Relevance:** TNBC has the worst prognosis among breast cancer subtypes and lacks targeted therapies. Identifying molecular differences between Luminal A and TNBC can inform therapeutic strategies.

---

## Pipeline Overview


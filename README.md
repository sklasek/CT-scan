---
title: "CT-scan"
author: "Scott Klasek"
date: "11/17/2020"
output: github_document
---

This repo shows codes for analyzing sediment microbial community responses to X-ray CT scanning.
Sequence data are 16S rRNA V4 hypervariable regions available at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA533633
For code showing the processing of ASVs from raw fastq sequences and the generation of phyloseq objects, see 01_sequence_processing.md
For code showing analysis of microbial community data starting from the phyloseq object ps4, see 02_sequence_analysis.md
ps4 is the phyloseq object that has been removed of unwanted taxa, decontaminated, and pruned of samples with very poor read depth

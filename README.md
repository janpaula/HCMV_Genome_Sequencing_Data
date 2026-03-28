# HCMV_Genome_Sequencing_Data

## Overview

This repository contains a bioinformatics workflow for the identification and analysis of genomic variants in **Human Cytomegalovirus (HCMV)** samples.
The analysis focuses on detecting SNPs and indels, characterizing their genomic distribution, and identifying genes with the highest mutation frequency.

The reference genome used in this study is **HCMV strain Merlin (NC_006273.2)**.

---

## Objectives

The main objectives of this study are:

* Identify genomic variants in HCMV sequencing data
* Annotate variants with genomic features
* Determine which viral genes accumulate the highest number of mutations
* Visualize variant distribution across the genome
* Identify potential mutation hotspots

---

## Dataset

The analysis was performed using sequencing data from **Human Cytomegalovirus (HCMV)**.

Reference genome:

* **NC_006273.2 – HCMV strain Merlin**

Annotation file:

* GFF3 annotation obtained from NCBI

---

## Workflow Overview

The analysis pipeline includes the following steps:

1. Quality control of sequencing reads
2. Alignment to the reference genome
3. Variant calling
4. Variant filtering
5. Functional annotation
6. Statistical analysis of mutation distribution
7. Visualization of genomic variants

---

## Tools and Software

The following tools were used in the analysis:

* BWA
* SAMtools
* BCFtools
* SnpEff
* R / Bioconductor

R packages used:

* VariantAnnotation
* GenomicRanges
* rtracklayer
* dplyr
* ggplot2
* Gviz

---

## Repository Structure

```
HCMV_variant_analysis/
│
├── data/                # Input data (reference genome, annotations)
├── scripts/             # Analysis scripts
├── results/             # Output files and figures
├── figures/             # Generated plots
└── README.md
```

---

## Example Results

### Variant Distribution Across Genes

A total of **2044 variants** overlapped coding regions of the HCMV genome.

Variant distribution across genes was heterogeneous, with several genes presenting higher mutation loads.

Genes with the highest mutation frequency included:

* RL1
* UL55
* UL83
* UL54

These genes may represent potential regions of increased evolutionary variability in the viral genome.

---

## Visualization

The analysis includes graphical representations such as:

* Variant density plots
* Gene mutation frequency barplots
* Genomic tracks of variants along the HCMV genome

---

## Reproducibility

To reproduce the analysis:

1. Clone the repository:

```
git clone https://github.com/your-username/HCMV_variant_analysis.git
```

2. Install required R packages.

3. Run the scripts in the `scripts/` directory following the workflow order.

---

## Author

Jânice Roberta de Paula
Bioinformatics Training Project

---

## License

This project is available for academic and educational purposes.

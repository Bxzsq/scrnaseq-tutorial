# snRNA-seq Analysis Pipeline

This project is a Nextflow-based pipeline for single-nucleus RNA sequencing (snRNA-seq) data analysis. It is designed to process raw sequencing data and perform downstream analysis and cell type annotation.

## Purpose

The main purpose of this project is to establish a standardized and reproducible workflow for snRNA-seq data analysis, from raw sequencing data to cell type annotation.

## Workflow

The main workflow is:

```text
Raw FASTQ
   ↓
FastQC
   ↓
Cell Ranger
   ↓
Gene expression matrix
   ↓
Seurat
   ├── Quality control
   ├── Normalization
   ├── Dimensionality reduction
   ├── Clustering
   └── Marker gene identification
   ↓
Cell type annotation
```

After Seurat analysis, cell type annotation can be performed using two approaches.

### Automatic Annotation

```text
Seurat results
      ↓
CellTypist
      ↓
Predicted cell types
      ↓
Annotation results
```

CellTypist is used to automatically assign cell types based on reference cell types.

### Manual Annotation

```text
Seurat clusters
      ↓
Marker genes
      ↓
Known cell-type markers
      ↓
Manual cell type annotation
```

Manual annotation is performed by examining cluster marker genes and assigning cell types based on known marker genes.

The results from automatic and manual annotation can be compared to evaluate the consistency of cell type assignments.

## Project Structure

```text
main.nf                 Main Nextflow workflow
nextflow.config         Pipeline configuration
modules/                Nextflow modules
bin/                    R and Python analysis scripts
docker/                 Dockerfiles
samplesheet.csv         Sample information
samplesheet_test.csv    Test sample information
results/                Analysis results
work/                   Nextflow working directory
```

## Running the Pipeline

Run the main workflow using Docker:

```bash
nextflow run main.nf \
    --input_mode fastq \
    -profile docker \
    -resume
```

The analysis results are saved in the `results/` directory.

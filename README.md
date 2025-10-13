# EGG: Epitope Generation through Gene-networks

A modular Snakemake pipeline for personalized neoantigen discovery and prioritization using patient-specific gene co-expression networks.

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0-brightgreen.svg)](https://snakemake.readthedocs.io)

## Overview

EGG is a comprehensive computational workflow for identifying and prioritizing tumor neoantigens from DNA and RNA sequencing data. Unlike traditional neoantigen prediction pipelines that focus solely on sequence-level features, EGG integrates **patient-specific gene co-expression networks** to identify neoantigens in biologically central and functionally relevant contexts.

### Key Features

- **Modular Architecture**: Seven independent Snakemake modules can be run separately or as a complete workflow
- **Multi-source Neoantigen Detection**: Identifies neoantigens from:
  - Somatic mutations (SNVs and indels)
  - Gene fusions
  - Alternative splicing events
- **Network-based Prioritization**: Uses LIONESS-derived personalized co-expression networks to assess functional importance
- **Automated Setup**: Automatically downloads references, installs dependencies, and configures environments
- **Reproducible**: Containerized tools (Docker) and versioned environments (Conda) ensure consistent results

### What Makes EGG Different?

EGG goes beyond traditional binding affinity predictions by incorporating:
- **Gene co-expression network topology** to identify functionally critical genes
- **Protein-protein interaction** data (HumanNet-XN) to prioritize biologically relevant candidates
- **Network centrality metrics** that capture a gene's importance to tumor cell function
- **Consensus scoring** that integrates network features with orthogonal evidence (DepMap essentiality, binding affinity, expression)

---

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Modules](#pipeline-modules)
- [Input Requirements](#input-requirements)
- [Execution Guide](#execution-guide)
- [Output](#output)
- [Resource Requirements](#resource-requirements)
- [Citation](#citation)
- [Contact](#contact)

---

## Installation

### Prerequisites

- **Conda/Mamba** (for environment management)
- **Snakemake** ≥6.0
- **Docker** (for containerized tools)

### Setup

1. **Clone the repository**
```bash
git clone https://github.com/yourusername/EGG.git
cd EGG
```

2. **Install Snakemake** (if not already installed)
```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
```

3. **That's it!** All other dependencies are managed automatically by the pipeline.

---

## Quick Start

### 1. Prepare your sample metadata

Create a CSV file (`samples.csv`) with the following structure:

```csv
sample_name,fastq_R1,fastq_R2,datatype,lane
SAMPLE_001,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,CancerRNA,L001
SAMPLE_001,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,CancerDNA,L001
SAMPLE_001,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,NormalDNA,L001
```

### 2. Run the complete pipeline

```bash
# Run all modules sequentially
snakemake --use-conda --cores 16 \
    --config csvfile="samples.csv" results_dir="results/"
```

### 3. Or run individual modules

```bash
# Example: Run only mutation analysis
snakemake --use-conda --snakefile scripts/final_scripts/mutation_analysis_module_1A.smk \
    --cores 16 --config csvfile="samples.csv" results_dir="results/"
```

---

## Pipeline Modules

EGG is organized into seven modular components, each handling a specific analysis step:

### Module 0: QC Filtering
**Purpose**: Quality control and preprocessing of raw sequencing data

**Tools**: FastQC, fastp

**Output**: Trimmed FASTQ files, QC reports

```bash
snakemake --use-conda --snakefile scripts/final_scripts/QC_filtering_module_0.smk \
    --cores {threads} --config csvfile="samples.csv" results_dir="output/"
```

---

### Module 1A: Mutation Analysis
**Purpose**: Somatic and germline variant calling from DNA sequencing

**Tools**: GATK (BQSR, Mutect2, HaplotypeCaller), dbSNP, Panel of Normals

**Output**: VCF files with somatic SNVs/indels and germline variants

```bash
snakemake --use-conda --snakefile scripts/final_scripts/mutation_analysis_module_1A.smk \
    --cores {threads} --config csvfile="samples.csv" results_dir="output/"
```

---

### Module 1B: RNA Analysis
**Purpose**: Multi-faceted RNA-seq analysis

**Tools**: 
- STAR-Fusion (gene fusions)
- ArcasHLA (HLA typing)
- FeatureCounts (gene expression quantification)

**Output**: Fusion calls, HLA types, gene expression matrices

```bash
snakemake --use-conda --snakefile scripts/final_scripts/RNA_analysis_module_1B.smk \
    --cores {threads} --config csvfile="samples.csv" results_dir="output/"
```

---

### Module 2A: Somatic Mutation Epitopes
**Purpose**: Predict MHC binding epitopes from somatic mutations

**Tools**: pVACseq

**Input**: Annotated VCF from Module 1A, HLA types from Module 1B

**Output**: Predicted neoepitopes from SNVs and indels

```bash
snakemake --use-conda --snakefile scripts/final_scripts/somatic_mutation_epitopes_2A.smk \
    --cores {threads} --config csvfile="samples.csv" results_dir="output/"
```

---

### Module 2B: Fusion Epitopes
**Purpose**: Predict neoepitopes from gene fusion breakpoints

**Tools**: pVACfuse

**Input**: Fusion transcripts from Module 1B

**Output**: Predicted neoepitopes spanning fusion junctions

```bash
snakemake --use-conda --snakefile scripts/final_scripts/fusion_epitopes_module_2B.smk \
    --cores {threads} --config csvfile="samples.csv" results_dir="output/"
```

---

### Module 2C: Splicing Epitopes
**Purpose**: Identify neoepitopes from alternative splicing events

**Tools**: RegTools, pVACsplice

**Input**: RNA-seq alignments (Module 1B), DNA variants (Module 1A)

**Output**: Predicted neoepitopes from tumor-specific splice junctions

```bash
snakemake --use-conda --snakefile scripts/final_scripts/splicing_epitopes_module_2C.smk \
    --cores {threads} --config csvfile="samples.csv" results_dir="output/"
```

---

### Co-expression Network Module
**Purpose**: Construct patient-specific gene co-expression networks and compute network centrality

**Tools**: LIONESS, HumanNet-XN PPI database

**Method**: 
- Generates personalized Pearson correlation networks
- Filters networks using PPI data and membrane protein enrichment
- Computes centrality metrics (degree, betweenness, strength, Weighted Connectivity Impact (WCI))
- Calculates Largest Component Impact (LCI) to identify structurally critical genes

**Output**: Network features for each gene per patient

```bash
snakemake --use-conda --snakefile scripts/final_scripts/sample_networks_generation.smk \
    --cores {threads} --config csvfile="samples.csv" results_dir="output/"
```

---

### Prioritization Module
**Purpose**: Integrate all evidence streams to rank neoepitope candidates

**Features Integrated**:
- Network topology (centrality, LCI)
- Gene essentiality (DepMap)
- HLA binding affinity
- Gene expression levels
- Subcellular localization
- Driver gene status (IntOGen)
- Cancer hallmarks (GO annotations)

**Method**: Borda consensus scoring across all features

**Output**: Ranked neoepitope table with comprehensive annotations

```bash
snakemake --use-conda --snakefile scripts/final_scripts/epitopes_prioritisation.smk \
    --cores {threads} --config csvfile="samples.csv" results_dir="output/"
```

---

## Input Requirements

### Sample Metadata CSV

The pipeline requires a CSV file with the following columns:

| Column | Description | Example |
|--------|-------------|---------|
| `sample_name` | Unique sample identifier | SAMPLE_001 |
| `fastq_R1` | Path to forward reads | /path/to/sample_R1.fastq.gz |
| `fastq_R2` | Path to reverse reads | /path/to/sample_R2.fastq.gz |
| `datatype` | Sample type | CancerRNA / CancerDNA / NormalDNA |
| `lane` | Sequencing lane identifier | L001 |

### Required Data Types

- **CancerRNA**: Tumor RNA-seq (for expression, fusions, splicing, HLA typing)
- **CancerDNA**: Tumor DNA-seq (for somatic mutations)
- **NormalDNA**: Matched normal DNA-seq (for filtering germline variants)

---

## Execution Guide

### Recommended Workflow

1. **Start with Module 0** (QC Filtering) on all samples
```bash
snakemake --use-conda --snakefile scripts/final_scripts/QC_filtering_module_0.smk \
    --cores 16 --config csvfile="samples.csv" results_dir="results/"
```

2. **Run Module 1A and 1B in parallel** (DNA and RNA analysis)
```bash
# Terminal 1: DNA analysis
snakemake --use-conda --snakefile scripts/final_scripts/mutation_analysis_module_1A.smk \
    --cores 8 --config csvfile="samples.csv" results_dir="results/"

# Terminal 2: RNA analysis
snakemake --use-conda --snakefile scripts/final_scripts/RNA_analysis_module_1B.smk \
    --cores 8 --config csvfile="samples.csv" results_dir="results/"
```

3. **Run epitope prediction modules** (2A, 2B, 2C)
   - Module 2A requires: Module 1A + HLA types (from 1B or provided externally)
   - Module 2B requires: Module 1B
   - Module 2C requires: Module 1A + Module 1B

4. **Generate co-expression networks**
```bash
snakemake --use-conda --snakefile scripts/final_scripts/sample_networks_generation.smk \
    --cores 16 --config csvfile="samples.csv" results_dir="results/"
```

5. **Run prioritization** (requires all previous modules)
```bash
snakemake --use-conda --snakefile scripts/final_scripts/epitopes_prioritisation.smk \
    --cores 16 --config csvfile="samples.csv" results_dir="results/"
```

### Module Dependencies

```
Module 0 (QC)
    ↓
    ├─→ Module 1A (DNA) ──→ Module 2A (SNV/Indel epitopes)
    │                   ↘
    │                    └→ Module 2C (Splice epitopes)
    │                       ↗
    └─→ Module 1B (RNA) ──→ Module 2B (Fusion epitopes)
                        ↓
                    Network Module
                        ↓
                All modules → Prioritization Module
```

---

## Output

### Final Prioritized Neoepitope Table

The prioritization module produces a comprehensive ranked table containing:

**Per-epitope information**:
- `sample_id`: Patient identifier
- `gene`: Gene harboring the alteration
- `variant_class`: SNV / indel / fusion / splice
- `peptide`: Amino acid sequence of the epitope
- `peptide_length`: Length of the peptide (typically 8-11 amino acids)
- `HLA_allele`: Predicted binding HLA allele

**Immunogenicity features**:
- `binding_affinity`: Predicted IC50 (nM)
- `binding_rank`: Percentile rank of binding affinity
- `expression`: Gene expression level (TPM/FPKM)

**Network features**:
- `degree`: Number of co-expressed genes
- `betweenness`: Centrality in co-expression network
- `betweennes_centrality`: Network information flow score
- `LCI`: Largest Component Impact (structural criticality)

**Functional annotations**:
- `DepMap_dependency`: Cancer cell line essentiality score
- `Gene_Ontology`: GO cellular component
- `driver_status`: IntOGen driver classification
- `cancer_hallmarks`: Associated cancer pathways

**Final scores**:
- `Borda_consensus_score`: Integrated multi-feature score
- `final_rank`: Overall epitope ranking

### Interactive Outputs

The pipeline also generates:
- Interactive HTML tables for exploring results
- Gene Ontology enrichment analysis
- Cancer hallmark enrichment per sample

---

## Resource Requirements

### Computational Resources

**Minimum**:
- 32 GB RAM
- 8 CPU cores
- 200 GB disk space

**Recommended**:
- 64+ GB RAM
- 16+ CPU cores
- 300 GB disk space (depends on how many samples are run)

### Reference Data

The pipeline automatically downloads required references (~117 GB):

```
resources/
├── genome_DNA/              # Human reference genome (DNA)
├── genome_RNA_fusion/       # STAR-Fusion reference library
├── GTF/                     # Gene annotations
├── Kallisto_ref/            # Kallisto indices
├── PON/                     # Panel of normals
├── VEP/                     # Variant Effect Predictor cache
├── dbSNP/                   # dbSNP database
├── gnomad/                  # gnomAD allele frequencies
├── arcasHLA/                # HLA reference sequences
├── agfusion/                # Gene fusion database
└── WGS_intervals/           # Genomic intervals for parallelization
```

These are downloaded once and cached locally for future runs.

---

## Configuration

### Adjusting pVACtools Parameters

Edit the Snakemake rule files directly to modify pVACseq/pVACfuse/pVACsplice parameters:
```python
# In scripts/final_scripts/somatic_mutation_epitopes_2A.smk
pvacseq run \
    --binding-threshold 500 \
    --percentile-threshold 2 \
    --minimum-fold-change 1
```

### Adjusting Prioritization Weights

Edit the Borda consensus configuration YAML:
```yaml
# config/prioritization_config.yaml
feature_weights:
  binding_affinity: 3
  expression: 2
  network_centrality: 2
  essentiality: 1
```

---

## Troubleshooting

### Common Issues

**Issue**: Docker permission errors
```bash
# Solution: Add your user to the docker group
sudo usermod -aG docker $USER
# Then log out and back in
```

**Issue**: Reference download failures
```bash
# Solution: Manually download and place in resources/ directory
# Check logs for the specific URL that failed
```

**Issue**: Out of memory errors
```bash
# Solution: Reduce parallel jobs or increase system RAM
snakemake --cores 8 --resources mem_mb=32000
```

---

## Citation

If you use EGG in your research, please cite:

```
[Your Citation Here]
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contact

For questions, issues, or feature requests:
- **GitHub Issues**: [github.com/yourusername/EGG/issues](github.com/yourusername/EGG/issues)
- **Email**: your.email@institution.edu

---

## Acknowledgments

EGG integrates and builds upon numerous open-source tools:
- **GATK** (Broad Institute)
- **pVACtools** (Griffith Lab)
- **STAR-Fusion** (Broad Institute)
- **HumanNet-XN** (Seoul National University)
- **LIONESS** (Glass Lab)

We thank the developers of these tools for making their software freely available.
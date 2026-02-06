# EGG: Epitope Generation Gateway

A modular Snakemake pipeline for personalized neoantigen discovery and prioritization using patient-specific gene co-expression networks.

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
- **Network centrality metrics** that capture a gene's importance to tumor cell function in particular related to resistance to cancer evolution and immune evasion
- **Consensus scoring** that integrates network features with orthogonal evidence (DepMap essentiality, binding affinity, expression)

---

## Table of Contents

- [Pipeline Implementation](#pipeline-implementation)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Modules](#pipeline-modules)
- [Input Requirements](#input-requirements)
- [Execution Guide](#execution-guide)
- [Output](#output)
- [Resource Requirements](#resource-requirements)
- [Citation](#citation)

---

## Pipeline Implementation

EGG is implemented as a Snakemake workflow cosisting of seven separate modules. where each module corresponds to a dedicated .smk file that can be executed independently or as part of an end-to-end run. 

Each module is encapsulated in its own Snakefile under scripts/final_scripts/ and follows a consistent interface (csvfile=..., results_dir=...) so outputs are organized under a single root directory and can be reused reliably by downstream modules. Sample inputs are provided through a standardized CSV metadata file (sample identity, FASTQ paths, datatype, lane), enabling uniform batch processing without per-sample manual configuration.

To minimize setup burden while preserving reproducibility, EGG provisions dependencies automatically at runtime: modules create and use local Conda environments via --use-conda, pull and run Docker images when required by specific tools, and download reference resources as needed as shown in more detail in the resource requirements section. Beyond installing Snakemake, Conda, and Docker, no additional manual installation should be necessary.

In practice, the workflow is commonly initialized with QC filtering (Module 0), after which DNA (Module 1A) and RNA (Module 1B) analyses can run in parallel. Downstream epitope generation modules consume these outputs: 2A uses somatic variants from 1A plus HLA alleles (from 1B or provided externally), 2B uses fusion calls from 1B, and 2C integrates RNA-derived splice events with matched DNA calls and HLA context. RNA-seq expression outputs feed the co-expression network module, which generates patient-specific network features used alongside orthogonal evidence (binding features, expression, essentiality, and biological annotations) in the prioritization module to produce a final per-sample ranked neoepitope table. This is intended to complement binding and expression-based ranking with tumor functional context, prioritizing epitopes arising from functionally indispensable / network-central genes that may be less likely to be lost under tumor evolution and selective pressure, thereby supporting identification of potentially more durable vaccine targets. The recommended execution order is explained in the Execution Guide.

## Installation

### Prerequisites

- **Conda/Mamba** (for environment management)
- **Snakemake** ≥6.0
- **Docker** (for containerized tools)

### Setup

1. **Clone the repository**
```bash
git clone https://gitlab.com/grecolab_group/epitope_generation_gateway
cd epitope_generation_gateway
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


### 2. Run each module separately with same result_dir name

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

2. **Run Module 1A and 1B in after Module 0** (DNA and RNA analysis)
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


### Output Organization

The pipeline generates outputs organized by module:

```
results/
├── 0_Filtering_and_QC/
│   ├── trimmed_fastq/                      # Quality-filtered FASTQ files
│   ├── qc_pre_filtering/                   # Initial FastQC reports
│   └── qc_post_filtering/                  # Post-trimming QC reports
│
├── 1A_mutation_analysis/
│   ├── VCF/                                # Raw variant calls
│   ├── VCF_filtered/                       # Filtered somatic variants
│   ├── VCF_germline/                       # Germline variant calls
│   ├── VCF_germline_filtered/              # Filtered germline variants
│   ├── bwa/                                # Alignment outputs
│   ├── dedup/                              # Duplicate-marked BAMs/metrics
│   ├── bqsr/                               # BQSR-processed BAMs
│   ├── merge/                              # Merge fastqs from same sample run across different lanes 
│   ├── metrics/                            # Alignment/QC metrics
│   └── filtering_tables/                   # Variant filtering tables/logs
│
├── 1B_RNA_fusion_HLA/
│   ├── ArcasHLA/                           # HLA typing results
│   ├── RNA_Counts/                         # Gene expression quantification
│   ├── StarFusionOut/                      # Fusion predictions
│   └── sorted_bam/                         # Coordinate-sorted RNA BAMs
│
├── 2A_somatic_mutation_epitopes/
│   ├── annotated_germline_VCF/             # Annotated germline variants
│   ├── annotated_germline_VCF_name_updated/# Germline VCFs with updated naming
│   ├── annotated_phased_VCF/               # Annotated phased variants
│   ├── annotated_somatic_VCF/              # Annotated somatic variants
│   ├── combined_VCF/                       # Combined VCFs (somatic/germline/etc.)
│   ├── kallisto_quantification/            # Expression quantification (kallisto)
│   ├── kallisto_somatic_VCF/               # Somatic VCFs used alongside expression
│   ├── phased_VCF/                         # Phased VCFs
│   ├── pvacSeq/                            # Neoepitope predictions
│   ├── regtools/                           # regtools outputs for splicing evidence
│   ├── sorted_VCF/                         # Sorted VCFs
│   └── tumor_only_VCF/                     # Tumor-only variant calls
│
├── 2B_fusion_epitopes/
│   ├── AGfusion/                           # Fusion annotations
│   └── pvacFuse/                           # Fusion neoepitope predictions
│
├── 2C_splicing_epitopes/
│   ├── annotated_somatic_VCF/              # Annotated somatic VCFs (splicing context)
│   ├── pvacSplice/                         # Splicing neoepitope predictions
│   └── regtools_genomic_VCF_genecode/      # Splice junction calls
│
├── sample_specific_networks/
│   ├── Network_Metrics_Betweenness/        # Betweenness centrality per gene
│   ├── Network_Metrics_Degree/             # Degree centrality per gene
│   ├── Network_Metrics_Full/               # Full metrics outputs
│   │   ├── Strength/                       # Strength metrics
│   │   └── WCI/                            # WCI metrics
│   ├── Network_Metrics_LargestComponentImpact/ # LCI scores per gene
│   ├── Network_Metrics_Strength/           # Strength centrality per gene
│   ├── Sample_Specific_Networks/           # Patient-specific networks
│   ├── Sample_Specific_Networks_Distances/ # Network distance computations
│   ├── Sample_Specific_Networks_PPI_filtered/
│   │   ├── filtered_networks_matrix/
│   │   ├── filtered_networks_rds/
│   │   ├── filtered_networks_rds_old/
│   │   └── qc_plots/
│   ├── Sample_Specific_Networks_PPI_filtered_copy/
│   │   └── filtered_networks_matrix/
│   ├── gene_lists/                         # Gene sets used for networks
│   └── normalised_counts/                  # Normalized expression for network construction
│
└── epitopes_prioritisation/
    ├── Borda_Epiotpes_Prioritisation/      # Consensus scoring results
    ├── annotated_epitopes/                 # Annotated epitopes
    ├── annotated_epitopes_with_network/    # Epitopes merged with network features
    ├── annotationhub_cache/                # AnnotationHub cache
    ├── combined_epitopes/                  # Merged epitopes from all sources
    └── final_epitopes/                     # Final ranked neoepitope outputs
        ├── HLA/
        ├── HTML/
        └── html_out/

```

### Final Prioritized Neoepitope Table

The prioritization module produces a comprehensive ranked table containing:

**Per-epitope information**:
- `sample_id`: Patient identifier
- `gene`: Gene harboring the alteration
- `variant_class`: SNV / indel / fusion / splice
- `peptide`: Amino acid sequence of the epitope
- `peptide_length`: Length of the epitope peptide
- `HLA_allele`: Predicted binding HLA allele

**Immunogenicity features**:
- `binding_affinity`: Predicted IC50 (nM)
- `binding_rank`: Percentile rank of binding affinity
- `expression`: Gene expression level (TPM/FPKM)

**Network features**:
- `degree`: Number of co-expressed genes
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
├── WGS_intervals/           # Genomic intervals
└── interval_list/           # WES genomic intervals
```

These are downloaded once and saved locally for future runs.



---

## Configuration
All editable YAML configs are found under:
Epitope_Generation_Gateway/scripts/final_scripts/config/

### Adjusting QC Filtering

Configure read-quality filtering (via fastp) by editing:
```
# Relative path: Epitope_Generation_Gateway/scripts/final_scripts/config/QC_filtering_config_0.yaml

fastp:
  detect_adapter: true          # Auto-detect adapters; set to false to use custom sequences
  custom_adapter_R1: ""         # Custom adapter for R1 (leave blank if auto-detecting)
  custom_adapter_R2: ""         # Custom adapter for R2 (leave blank if auto-detecting)
  min_length: 50                # Drop reads shorter than this after trimming
  quality: 20                   # Minimum Phred quality for a base to be considered "qualified"
  unqualified_base_limit: 30    # Max % of unqualified bases allowed per read
  cut_front: false              # Trim low-quality bases from the 5' end if true
  cut_tail: false               # Trim low-quality bases from the 3' end if true
  cut_window_size: 4            # Sliding window size for quality-based trimming
  cut_mean_quality: 20          # Mean quality threshold within the window to trigger trimming
```
### Adjusting Prioritization Weights

EGG uses a Borda consensus to combine feature-ranked scores. Edit the YAML to change which columns are used, their weights, and whether higher or lower values rank better:

```yaml
# config/borda_config.yaml
borda:
  # Choose columns, weight (any positive numbers; auto-renormalized),
  # and direction: "lower_better" or "higher_better".
  columns:
    Median.MT.IC50.Score:
      weight: 0.25
      direction: lower_better
    Depmap_survivability_score:
      weight: 0.25
      direction: lower_better
    Network_Betweenness:
      weight: 0.10
      direction: higher_better
    Network_Degree:
      weight: 0.10
      direction: higher_better
    Network_Impact:
      weight: 0.10
      direction: higher_better
    Network_Strength:
      weight: 0.10
      direction: higher_better
    Network_WCI:
      weight: 0.10
      direction: higher_better

  # How to break ties between equal ranks: average|min|max|first|random
  ties_method: average

  # NA handling policy:
  # - ignore: exclude NA for that feature and renormalize weights per row (no penalty)
  # - worst:  treat NA as worst rank (penalize missing values)
  # - drop:   drop rows with any NA across selected features (strict)
  na_policy: ignore

---


### Adjusting pVACtools Parameters

Edit the Snakemake rule files directly to modify pVACseq/pVACfuse/pVACsplice parameters:
```python
# In scripts/final_scripts/somatic_mutation_epitopes_2A.smk
pvacseq run \
    --binding-threshold 500 \
    --percentile-threshold 2 \
    --minimum-fold-change 1
```

## Troubleshooting

### Common Issues


**Issue**: Reference download failures
```bash
# Solution: Manually download and place in resources/ directory
# Check logs for the specific URL that failed
```
This troubleshooting secion will be updated as issues appear


---

## Citation

If you use EGG in your research, please cite:

```
[Citation]
```



## Acknowledgments

EGG integrates and builds upon numerous open-source tools:
- **GATK** (Broad Institute)
- **pVACtools** (Griffith Lab)
- **STAR-Fusion** (Broad Institute)
- **HumanNet-XN** (Seoul National University)
- **LIONESS** (Glass Lab)
- **DepMap** (Broad Institute)


We thank the developers of these tools for making their software freely available.

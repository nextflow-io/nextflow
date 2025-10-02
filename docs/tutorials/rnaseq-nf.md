(rnaseq-nf-page)=

# Getting started with rnaseq-nf

`rnaseq-nf` is a basic Nextflow pipeline for RNA-Seq analysis that performs quality control, transcript quantification, and result aggregation. The pipeline processes paired-end FASTQ files, generates quality control reports with FastQC, quantifies transcripts with Salmon, and produces a unified report with MultiQC.

This tutorial describes the architecture of the `rnaseq-nf` pipeline and provides instructions on how to run it.

## Pipeline architecture

The pipeline is organized into modular workflows and processes that coordinate data flow from input files through analysis steps to final outputs.

### Entry workflow

The entry workflow orchestrates the entire pipeline by coordinating input parameters and data flow:

```{mermaid}
flowchart TB
    subgraph " "
    subgraph params
    v0["transcriptome"]
    v1["reads"]
    v5["multiqc"]
    v2["outdir"]
    end
    v4([RNASEQ])
    v6([MULTIQC])
    v0 --> v4
    v1 --> v4
    v4 --> v6
    v5 --> v6
    end
```

<h3>Data flow:</h3>

1. The `transcriptome` and `reads` parameters are passed to the `RNASEQ` subworkflow, which performs indexing, quality control, and quantification
2. The outputs from `RNASEQ` along with the `multiqc` configuration are passed to the `MULTIQC` module, which aggregates results into a unified HTML report
3. The `outdir` parameter defines where all results are published

### `RNASEQ`

The `RNASEQ` subworkflow coordinates three processes that run in parallel and sequence:

```{mermaid}
flowchart TB
    subgraph RNASEQ
    subgraph take
    v1["read_pairs_ch"]
    v0["transcriptome"]
    end
    v2([INDEX])
    v3([FASTQC])
    v4([QUANT])
    subgraph emit
    v5["$out"]
    end
    v0 --> v2
    v1 --> v3
    v1 --> v4
    v2 --> v4
    v3 --> v5
    v4 --> v5
    end
```

<h3>Inputs (`take:`):</h3>

- `transcriptome`: Reference transcriptome file
- `read_pairs_ch`: Channel of paired-end read files

<h3>Process execution (`main:`):</h3>

- `INDEX` creates a Salmon index from the `transcriptome` (runs once)
- `FASTQC` analyzes the `read_pairs_ch` in parallel (runs independently for each sample)
- `QUANT` quantifies transcripts using both the index from **INDEX** and the `read_pairs_ch` (runs for each sample after INDEX completes)

<h3>Outputs (`emit:`):</h3>

- All outputs from `FASTQC` and `QUANT` are collected and emitted for downstream processing

### `MULTIQC`

The `MULTIQC` module aggregates all quality control and quantification outputs into a comprehensive HTML report.

<h3>Inputs:</h3>

- `RNASEQ` outputs: All collected outputs from the `RNASEQ` subworkflow (FastQC reports and Salmon quantification files)
- MultiQC config: Custom configuration files and branding (logo, styling)

<h3>Process execution:</h3>

- `MULTIQC` scans all input files, extracts metrics and statistics, and generates a unified report

<h3>Outputs:</h3>

- `multiqc_report.html`: A single consolidated HTML report providing an overview of:
  - General stats
  - Salmon fragment length distribution
  - FastQC quality control
  - Software versions

## Pipeline parameters
The pipeline behavior can be customized using command-line parameters to specify input data, output locations, and configuration files.

The pipeline accepts the following command-line parameters:

- `--reads`: Path to paired-end FASTQ files (default: `data/ggal/ggal_gut_{1,2}.fq`)
- `--transcriptome`: Path to reference transcriptome FASTA (default: `data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa`)
- `--outdir`: Output directory for results (default: `results`)
- `--multiqc`: Path to MultiQC configuration directory (default: `multiqc`)

## Execution profiles

Execution profiles allow you to customize how and where the pipeline runs by specifying the `-profile` flag. Multiple profiles can be combined by separating them with commas.

<h3>Container profiles</h3>

Container profiles specify which containerization technology to use for running the pipeline tools:

- `standard`: Use default Docker container
- `docker`: Explicitly use Docker
- `singularity`: Use Singularity containers
- `wave`: Use Wave container provisioning with Conda
- `wave-mirror`: Use Wave container mirroring strategy

:::{note}
The respective container tools must be installed to use these profiles.
:::

<h3>Environment profiles</h3>

Environment profiles manage software dependencies through package managers or specify architecture requirements:

- `conda`: Use Conda environment management
- `mamba`: Use Micromamba for faster dependency resolution
- `arm64`: Use ARM64 architecture support

:::{note}
The respective environment tools must be installed to use these profiles.
:::

<h3>Cloud and HPC profiles</h3>

Cloud and HPC profiles enable execution on distributed computing infrastructure and cloud storage:

- `slurm`: Run on SLURM-managed HPC clusters
- `batch`: Run on AWS Batch compute environments
- `google-batch`: Run on Google Cloud Batch
- `azure-batch`: Run on Azure Batch compute pools
- `s3-data`: Use input data stored in AWS S3
- `gs-data`: Use input data stored in Google Cloud Storage

:::{note}
To use the Cloud and HPC profiles, you must configure credentials, resource pools, and storage paths before execution.
:::

<h3>Other profiles</h3>

- `all-reads`: Process all FASTQ files matching `ggal_*_{1,2}.fq`

## Test data

The pipeline includes test data located in the `data/ggal/` directory for demonstration and validation purposes:

- Paired-end FASTQ files: Four tissue samples (gut, liver, lung, spleen) from *Gallus gallus* (chicken)
  - `ggal_gut_{1,2}.fq` - Default sample used when running with standard parameters
  - `ggal_liver_{1,2}.fq`
  - `ggal_lung_{1,2}.fq`
  - `ggal_spleen_{1,2}.fq`
- Reference transcriptome: `ggal_1_48850000_49020000.Ggal71.500bpflank.fa` - A subset of the chicken genome

:::{tip}
Use the `all-reads` profile to process all four tissue samples instead of just the default gut sample. See [Execution profiles](#execution-profiles) for more information.
:::

## Quick start

`rnaseq-nf` is a runnable pipeline. This section provides examples for running the pipeline with different configurations.

### Basic execution

Run the pipeline with default parameters using Docker:

```bash
nextflow run nextflow-io/rnaseq-nf -with-docker
```

### Configuring individual parameters

Override default parameters to use custom input files and output locations:

```bash
nextflow run nextflow-io/rnaseq-nf \
  --reads '/path/to/reads/*_{1,2}.fastq.gz' \
  --transcriptome '/path/to/transcriptome.fa' \
  --outdir 'my_results' \
  -with-docker
```

### Using profiles

Specify execution profiles to customize runtime environments and data sources:

```bash
# Use Conda for dependency management
nextflow run nextflow-io/rnaseq-nf -profile conda

# Run on a SLURM cluster
nextflow run nextflow-io/rnaseq-nf -profile slurm

# Combine multiple profiles: process all reads using Docker
nextflow run nextflow-io/rnaseq-nf -profile all-reads,docker
```

:::{tip}
See [Execution profiles](#execution-profiles) for more information about profiles.
:::

## Expected outputs

The `rnaseq-nf` pipeline generates the following outputs in the results directory:

```
results/
├── fastqc_<SAMPLE_ID>_logs/      # FastQC quality reports per sample
│   ├── <SAMPLE_ID>_1_fastqc.html
│   ├── <SAMPLE_ID>_1_fastqc.zip
│   ├── <SAMPLE_ID>_2_fastqc.html
│   └── <SAMPLE_ID>_2_fastqc.zip
└── multiqc_report.html           # Aggregated QC and Salmon report
```

The MultiQC report (multiqc_report.html) can be viewed in a web browser.

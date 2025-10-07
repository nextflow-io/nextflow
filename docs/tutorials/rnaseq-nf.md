(rnaseq-nf-page)=

# Getting started with rnaseq-nf

[`rnaseq-nf`](https://github.com/nextflow-io/rnaseq-nf) is a basic Nextflow pipeline for RNA-Seq analysis that performs quality control, transcript quantification, and result aggregation. The pipeline processes paired-end FASTQ files, generates quality control reports with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), quantifies transcripts with [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html), and produces a unified report with [MultiQC](https://seqera.io/multiqc/).

This tutorial describes the architecture of the [`rnaseq-nf`](https://github.com/nextflow-io/rnaseq-nf) pipeline and provides instructions for how to run it.

## Pipeline architecture

The pipeline is organized into modular workflows and processes that coordinate data flow from input files through analysis steps to final outputs.

### Entry workflow

The [entry workflow](https://github.com/nextflow-io/rnaseq-nf/blob/master/main.nf) orchestrates the entire pipeline by coordinating input parameters and data flow:

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

Data flow:

- The `transcriptome` and `reads` parameters are passed to the `RNASEQ` subworkflow, which performs indexing, quality control, and quantification.

- The outputs from `RNASEQ`, along with the MultiQC configuration (`multiqc`), are passed to the `MULTIQC` module, which aggregates results into a unified HTML report.

- The `outdir` parameter defines where all results are published.

### `RNASEQ`

The [`RNASEQ`](https://github.com/nextflow-io/rnaseq-nf/blob/master/modules/rnaseq.nf) subworkflow coordinates three processes that run in parallel and sequence:

```{mermaid}
flowchart TB
    subgraph RNASEQ
    subgraph take
    v0["read_pairs_ch"]
    v1["transcriptome"]
    end
    v2([INDEX])
    v4([FASTQC])
    v6([QUANT])
    subgraph emit
    v10["index"]
    v9["samples"]
    end
    v1 --> v2
    v0 --> v4
    v0 --> v6
    v2 --> v6
    v4 --> v9
    v6 --> v9
    v2 --> v10
    end
```

Inputs (`take:`):

- `read_pairs_ch`: Channel of paired-end read files
- `transcriptome`: Reference transcriptome file

Data flow (`main:`):

- [`INDEX`](https://github.com/nextflow-io/rnaseq-nf/blob/master/modules/index/main.nf) creates a Salmon index from the `transcriptome` input (runs once).

- [`FASTQC`](https://github.com/nextflow-io/rnaseq-nf/blob/master/modules/fastqc/main.nf) analyzes the samples in the `read_pairs_ch` channel in parallel (runs independently for each sample).

- [`QUANT`](https://github.com/nextflow-io/rnaseq-nf/blob/master/modules/quant/main.nf) quantifies transcripts using the index from `INDEX` and the samples in the `read_pairs_ch` channel (runs for each sample after `INDEX` completes).

Outputs (`emit:`):

- The results from `FASTQC` and `QUANT` are joined into a single channel for downstream processing.

- The index from `INDEX` is also emitted as a dataflow value.

### `MULTIQC`

The [`MULTIQC`](https://github.com/nextflow-io/rnaseq-nf/blob/master/modules/multiqc/main.nf) process aggregates all quality control and quantification outputs into a comprehensive HTML report.

Inputs:

- Input files: All collected outputs from the `RNASEQ` subworkflow (FastQC reports and Salmon quantification files).
- `config`: MultiQC configuration files and branding (logo, styling).

Process execution:

- `MULTIQC` scans all input files, extracts metrics and statistics, and generates a unified report.

Outputs:

- `multiqc_report.html`: A single consolidated HTML report providing an overview of:
  - General stats
  - Salmon fragment length distribution
  - FastQC quality control
  - Software versions

## Pipeline parameters

The pipeline behavior can be customized using command-line parameters to specify input data, output locations, and configuration files.

The pipeline accepts the following command-line parameters:

- `--reads`: Path to paired-end FASTQ files (default: `data/ggal/ggal_gut_{1,2}.fq`).

- `--transcriptome`: Path to reference transcriptome FASTA (default: `data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa`).

- `--outdir`: Output directory for results (default: `results`).

- `--multiqc`: Path to MultiQC configuration directory (default: `multiqc`).

## Configuration profiles

Configuration profiles allow you to customize how and where the pipeline runs by specifying the `-profile` flag. Multiple profiles can be specified as a comma-separated list. Profiles are defined in the [`nextflow.config`](https://github.com/nextflow-io/rnaseq-nf/blob/master/nextflow.config) file in the base directory.

<h3>Software profiles</h3>

Software profiles specify how software dependencies for processes should be provisioned:

- `conda`: Provision a Conda environment for each process based on its required Conda packages
- `docker`: Use a Docker container which contains all required dependencies
- `singularity`: Use a Singularity container which contains all required dependencies
- `wave`: Provision a Wave container for each process based on its required Conda packages

:::{note}
The respective container runtime or package manager must be installed to use these profiles.
:::

<h3>Execution profiles</h3>

Execution profiles specify the compute and storage environment used by the pipeline:

- `slurm`: Run on a SLURM HPC cluster
- `batch`: Run on AWS Batch
- `google-batch`: Run on Google Cloud Batch
- `azure-batch`: Run on Azure Batch

:::{note}
Depending on your environment, you may need to configure underlying infrastructure such as resource pools, storage, and credentials.
:::

## Test data

The pipeline includes test data in the [`data/ggal/`](https://github.com/nextflow-io/rnaseq-nf/tree/master/data/ggal) directory for demonstration and validation purposes:

- Paired-end FASTQ files from four tissue samples (gut, liver, lung, spleen):
  - `ggal_gut_{1,2}.fq`
  - `ggal_liver_{1,2}.fq`
  - `ggal_lung_{1,2}.fq`
  - `ggal_spleen_{1,2}.fq`

- Reference transcriptome:
  - `ggal_1_48850000_49020000.Ggal71.500bpflank.fa`

By default, only the `gut` sample is processed. You can use the `all-reads` profile to process all four tissue samples.

## Quick start

The [`rnaseq-nf`](https://github.com/nextflow-io/rnaseq-nf) pipeline is executable out-of-the-box. This section provides examples for running the pipeline with different configurations.

### Basic execution

Run the pipeline with default parameters using Docker:

```bash
nextflow run nextflow-io/rnaseq-nf -profile docker
```

### Configuring individual parameters

Override default parameters to use custom input files and output locations:

```bash
nextflow run nextflow-io/rnaseq-nf \
  --reads '/path/to/reads/*_{1,2}.fastq.gz' \
  --transcriptome '/path/to/transcriptome.fa' \
  --outdir 'my_results' \
  -profile docker
```

### Using profiles

Specify configuration profiles to customize runtime environments and data sources:

```bash
# Use Conda to provision software dependencies
nextflow run nextflow-io/rnaseq-nf -profile conda

# Run on a SLURM cluster
nextflow run nextflow-io/rnaseq-nf -profile slurm

# Combine multiple profiles: process all reads using Docker
nextflow run nextflow-io/rnaseq-nf -profile all-reads,docker
```

:::{tip}
See [Configuration profiles](#configuration-profiles) for more information about profiles.
:::

## Expected outputs

The [`rnaseq-nf`](https://github.com/nextflow-io/rnaseq-nf) pipeline produces the following outputs in the `results` directory:

```
results/
├── fastqc_<SAMPLE_ID>_logs/      # FastQC quality reports per sample
│   ├── <SAMPLE_ID>_1_fastqc.html
│   ├── <SAMPLE_ID>_1_fastqc.zip
│   ├── <SAMPLE_ID>_2_fastqc.html
│   └── <SAMPLE_ID>_2_fastqc.zip
└── multiqc_report.html           # Aggregated QC and Salmon report
```

The MultiQC report (`multiqc_report.html`) can be viewed in a web browser.

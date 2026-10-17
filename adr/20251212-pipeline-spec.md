# Pipeline spec

- Authors: Ben Sherman
- Status: accepted
- Deciders: Ben Sherman, Paolo Di Tommaso, Phil Ewels
- Date: 2025-12-12
- Tags: pipelines

## Summary

Provide a way for Nextflow to describe inherent properties of a pipeline that can be easily consumed by external systems.

## Problem Statement

A Nextflow pipeline is defined by Nextflow scripts (`.nf` files) and configuration (`.config` files). However, there are many aspects of a pipeline which are of interest to external systems, such as:

- Metadata (e.g. name, authors, license)
- Pipeline paramters and outputs
- Software dependency versions (e.g. modules)

Acquring this information directly from the source code requires parsing (or even executing) Nextflow code, which is generally not feasible for external systems. Additionally, it may be desirable to provide additional information that is not practical or otherwise does not belong in Nextflow code (e.g. display icons for pipeline parameters).

Primary use cases:

* **Viewing pipelines:** Display pipeline information (name, author, parameters, outputs) in an external user interface.

* **Form validation:** Validate pipeline parameters at launch time, prior to running the pipeline.

* **Pipeline chaining:** Validate a pipeline chain at launch time, allowing downstream pipeline inputs to reference upstream pipeline outputs that are compatible based on their respective pipeline specs.

- **Pipeline registry:** Enable pipelines to be published and executed as immutable software artifacts via the Nextflow registry, instead of cloning the source code repository.

## Solution

### Pipeline spec definition

The schema for pipeline specs is defined in [nextflow-io/schemas](https://github.com/nextflow-io/schemas/blob/main/pipeline/v1/schema.json). It was originally defined as the *meta-schema* for the [nf-core schema](https://nf-co.re/docs/nf-core-tools/pipelines/schema), a standard developed by the nf-core community to model pipeline parameters using JSON schema. The nf-core schema for a pipeline is typically defined as `nextflow_schema.json` in the project root.

Since the meta-schema was transferred to the `nextflow-io` GitHub organization, it is now considered an official Nextflow standard:

- The Nextflow language server uses the schema to provide code intelligence for pipeline parameters in Nextflow scripts.

- The Seqera Platform uses the schema to validate pipeline parameters at launch time.

- The `nf-schema` plugin, also under `nextflow-io`, uses the schema to validate pipeline parameters at runtime.

The pipeline spec adopts the structure of the nf-core schema, with only the following nominal changes:

- *nf-core schema* becomes *pipeline spec*
- *nf-core meta-schema* becomes *schema for pipeline specs*
- `nextflow_schema.json` becomes `nextflow_spec.json`

Preserving the structure of the original nf-core schema makes the migration process as easy as possible for users. At the same time, the nomenclature changes are needed to reduce confusion over different kinds of schemas and align with existing Nextflow standards (i.e. plugin specs, module specs).

The nf-core schema already defines the title, description, and parameters of a pipeline. The pipeline spec adds the following new properties:

- `version`: pipeline release version
- `contributors`: list of pipeline contributors (name, email, affiliation, etc)
- `documentation`: project documentation URL
- `homePage`: project home page
- `keywords`: relevant keywords
- `license`: project license
- `modules`: list of module versions used by the pipeline
- `requires`: runtime requirements
  - `nextflow`: Nextflow version constraint
  - `modules`: list of modules used by the pipeline
- `output`: list of pipeline outputs (name, type, description, etc)

Examples of these are shown in the following section on pipeline spec generation.

### Pipeline spec generation

Nextflow should be able to generate a pipeline spec from the pipeline source code:

- The parameter schema can be generated from the `params` block and associated record types.

- Samplesheet schemas (e.g. `schema_input.json`) can be generated from the record types used by corresponding parameters.

- The `output` section can be generated from the `output` block.

- Most of the other fields can be inferred from the `manifest` config scope in the main config file.

For example, given the following pipeline script and config:

**`main.nf`**

```groovy
params {
    // Samplesheet containing the input paired-end reads
    input: List<FastqPair>

    // The input transcriptome file
    transcriptome: Path

    // Directory containing multiqc configuration
    multiqc: Path = "${projectDir}/multiqc"
}

record FastqPair {
    id           : String
    fastq_1      : Path
    fastq_2      : Path?
    strandedness : Strandedness
}

enum Strandedness {
    FORWARD,
    REVERSE,
    UNSTRANDED,
    AUTO
}

workflow {
    // ...
}

output {
    // List of aligned samples
    samples: Channel<AlignedSample> {
        path { sample ->
            sample.fastqc >> 'fastqc/'
            sample.bam >> 'align/'
            sample.bai >> 'align/'
        }
        index {
            path 'samples.json'
        }
    }

    // MultiQC summary report
    multiqc_report: Path {
        path '.'
    }
}

record AlignedSample {
    id: String
    fastqc: Path
    bam: Path?
    bai: Path?
}
```

**`nextflow.config`**

```groovy
manifest {
    name = 'nf-core/rnaseq'
    contributors = [
        [
            name: 'Harshil Patel',
            affiliation: 'Seqera',
            github: '@drpatelh',
            contribution: ['author'],
            orcid: '0000-0003-2707-7940'
        ],
        [
            name: 'Phil Ewels',
            affiliation: 'Seqera',
            github: '@ewels',
            contribution: ['author'],
            orcid: '0000-0003-4101-2502'
        ],
    ]
    description = 'RNA sequencing analysis pipeline for gene/isoform quantification and extensive quality control.'
    nextflowVersion = '!>=25.04.3'
    version = '3.23.0'
}
```

The following pipeline spec should be produced:

**`nextflow_spec.json`**

```json
{
  // metadata
  "$schema": "https://raw.githubusercontent.com/nextflow/schemas/main/pipeline/v1/schema.json",
  "$id": "https://raw.githubusercontent.com/nf-core/rnaseq/refs/tags/3.23.0/nextflow_spec.json",
  "title": "nf-core/rnaseq",
  "description": "RNA sequencing analysis pipeline for gene/isoform quantification and extensive quality control.",
  "version": "3.23.0",
  "contributors": [
    {
      "name": "Harshil Patel",
      "affiliation": "Seqera",
      "github": "@drpatelh",
      "contribution": ["author"],
      "orcid": "0000-0003-2707-7940"
    },
    {
      "name": "Phil Ewels",
      "affiliation": "Seqera",
      "github": "@ewels",
      "contribution": ["author"],
      "orcid": "0000-0003-4101-2502"
    }
  ],

  // inputs
  "type": "object",
  "$defs": {
    "all_options": {
      "title": "Parameters",
      "type": "object",
      "properties": {
        "input": {
          "type": "string",
          "format": "file-path",
          "description": "Samplesheet containing the input paired-end reads",
          "schema": "assets/schema_input.json"
        },
        "transcriptome": {
          "type": "string",
          "format": "file-path",
          "description": "The input transcriptome file"
        },
        "multiqc": {
          "type": "string",
          "format": "directory-path",
          "description": "Directory containing multiqc configuration",
          "default": "${projectDir}/multiqc"
        }
      }
    }
  },
  "allOf": [
    {
      "$ref": "#/$defs/all_options"
    }
  ],

  // outputs
  "output": {
    "samples": {
      "description": "List of aligned samples",
      "schema": "assets/schema_samples.json",
      "path": "samples.json"
    },
    "multiqc_report": {
      "description": "MultiQC summary report",
      "type": "file",
      // (path)
    }
  },

  // software dependencies
  "requires": {
    "nextflow": "!>=25.04.3"
  }
}
```

**`assets/schema_input.json`**

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "id": {
        "type": "string",
      },
      "fastq_1": {
        "type": "string",
        "format": "file-path",
        "exists": true
      },
      "fastq_2": {
        "type": "string",
        "format": "file-path",
        "exists": true
      },
      "strandedness": {
        "type": "string",
        "enum": ["forward", "reverse", "unstranded", "auto"]
      },
    },
    "required": ["sample", "fastq_1", "strandedness"]
  }
}
```

**`assets/schema_samples.json`**

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "id": {
        "type": "string"
      },
      "fastqc": {
        "type": "string",
        "format": "directory-path"
      },
      "bam": {
        "type": "string",
        "format": "file-path"
      },
      "bai": {
        "type": "string",
        "format": "file-path"
      }
    },
    "required": ["id", "fastqc"]
  }
}
```

Notes:

- The `manifest` config options are effectively converted directly to JSON with only nominal changes, such as `manifest.name` -> `title` (preserve structure of original nf-core schema) and `nextflowVersion` -> `requires.nextflow` (leave space for module versions in the future).

- The parameter schema follows the structure of the nf-core schema, which defines *parameter groups* under `$defs` and combines them using JSON schema properties such as `allOf`. This section should be generated with sensible defaults since some properties (e.g. group name) can not be specified in pipeline code.

- Each output in the `output` section should specify either a type (e.g. `file`, `directory`) or a schema (e.g. if the output is a collection of records). Like parameters, the schema for an individual output should reference an external JSON schema file.

### Pipeline spec synchronization

The pipeline spec may contain additional fields that cannot be sourced from the pipeline code (e.g., the `fa_icon` property in the parameter schema). Such fields can be useful for external systems even if they aren't relevant to the pipeline execution.

As a result, the pipeline spec cannot be completely inferred from pipeline code. Instead, the generated pipeline spec should be treated as a "skeleton" that can be extended by the user with additional fields.

- When generating the pipeline spec, Nextflow should use any existing spec and preserve information that isn't inferred from pipeline code.

- Any inconsistencies between the existing spec and pipeline code (e.g. missing or extra parameters) should be reported as errors.

## Links

- [nextflow-io/schemas](https://github.com/nextflow-io/schemas)
- [nf-core schema](https://nf-co.re/docs/nf-core-tools/pipelines/schema)
- Examples: [nextflow_schema.json](https://github.com/nf-core/rnaseq/blob/3.23.0/nextflow_schema.json) and [schema_input.json](https://github.com/nf-core/rnaseq/blob/3.23.0/assets/schema_input.json)
- [JSON schema](https://json-schema.org/)

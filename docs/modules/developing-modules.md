(dev-modules-page)=

# Developing modules

Learn how to create modules for sharing through the Nextflow module registry.
For local module syntax (inclusion, aliases, templates, and binaries), see {ref}`module-page`.

(dev-modules-creating)=

## Creating a module

Use the `module create` command to scaffold a new module with the required files:

```console
$ nextflow module create myorg/my-module
```

If you omit the name, the command prompts you for details:

```console
$ nextflow module create
```

The command creates a module directory with the following files:

- `main.nf`: The module script containing your process definition.
- `meta.yml`: The module spec describing metadata, inputs, and outputs.
- `README.md`: Documentation for the module.

See {ref}`cli-module-create` for the full command reference.

(dev-modules-structure)=

## Module structure

Registry modules follow a standard directory structure:

```
modules/
└── myorg/
    └── my-module/
        ├── .module-info     # Integrity checksum (generated at install)
        ├── README.md        # Documentation (required for publishing)
        ├── main.nf          # Module script (required)
        ├── meta.yml         # Module spec (required for publishing)
        ├── resources/       # Optional: module binaries and resources
        └── templates/       # Optional: process templates
```

### main.nf

The `main.nf` file contains the process (or workflow/function) definition. For example, a simple module wrapping FastQC:

```nextflow
process FASTQC {
    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::fastqc=0.12.1'
    container "${ workflow.containerEngine == 'singularity'
        ? 'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0'
        : 'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip

    script:
    """
    fastqc $reads --threads $task.cpus
    """
}
```

### meta.yml

The `meta.yml` file describes the module's metadata, including its name, version, description, authors, and input/output specifications. Publishing requires this file, and the registry uses it to display module information and generate usage templates.

### README.md

The `README.md` file provides documentation for the module. It should describe what the module does, the tools it wraps, and any configuration requirements.

### templates

Include process script {ref}`templates <process-template>` alongside a module in the `templates/` directory.
See {ref}`module-templates` for details.

### resources

Modules can include binary scripts in the `resources/usr/bin/` directory that are locally scoped to the module's processes. See {ref}`module-binaries` for details.

(dev-modules-spec)=

## Generating a module spec

Use the `module spec` command to generate or update the `meta.yml` file from the module's `main.nf`:

```console
$ nextflow module spec myorg/my-module
```

Provide metadata fields directly to avoid `TODO` placeholders in the generated file:

```console
$ nextflow module spec \
    -namespace myorg \
    -version 1.0.0 \
    -description "Quality control of raw sequencing reads" \
    -license MIT \
    -author "@myname" \
    ./modules/myorg/my-module
```

Use `-dry-run` to preview the generated spec without writing to disk:

```console
$ nextflow module spec -dry-run myorg/my-module
```

If a `meta.yml` already exists, the command incorporates its content into the new file.

See {ref}`cli-module-spec` for the full command reference.

(dev-modules-validate)=

## Validating a module

Use the `module validate` command to check that a module is ready for publishing:

```console
$ nextflow module validate myorg/my-module
```

The command verifies that:

- All required files are present (`main.nf`, `meta.yml`, `README.md`).
- The module spec contains all required fields (name, version, description, license).

Validate by path:

```console
$ nextflow module validate ./modules/myorg/my-module
```

See {ref}`cli-module-validate` for the full command reference.

(dev-modules-testing)=

## Testing a module

Before publishing, test your module by running it directly:

```console
$ nextflow module run myorg/my-module --input 'test-data/*.fastq.gz'
```

The command executes the module as a standalone workflow, allowing you to verify that inputs are correctly declared, the process runs successfully, and outputs appear as expected.

For more thorough testing, create a small wrapper workflow that exercises the module:

```nextflow
include { MY_MODULE } from './modules/myorg/my-module'

workflow {
    input_ch = Channel.fromFilePairs('test-data/*_{1,2}.fastq.gz')
    MY_MODULE(input_ch)
    MY_MODULE.out.results.view()
}
```

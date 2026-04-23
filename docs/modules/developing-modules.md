(dev-modules-page)=

# Developing modules

Learn how to create modules and share them through the Nextflow module registry.

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

See {ref}`module create <cli-module-create>` for the full command reference.

(dev-modules-structure)=

## Module structure

Registry modules follow a standard directory structure:

```
modules/
└── myorg/
    └── my-module/
        ├── .module-info     # Integrity checksum
        ├── README.md        # Documentation
        ├── main.nf          # Module script
        ├── meta.yml         # Module spec
        ├── resources/       # Optional: Module resources
        └── templates/       # Optional: Process script templates
```

Local modules that are not intended for publishing do not need to follow this structure, although it is recommended as a best practice.

### main.nf

The `main.nf` file contains the process definition. For example, a simple module wrapping FastQC:

```nextflow
process FASTQC {
    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::fastqc=0.12.1'
    container 'biocontainers/fastqc:0.12.1--hdfd78af_0'

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

Registry modules are subject to the following constraints:

- The module must contain a single script called `main.nf`
- The module must define a single process
- The module must not define any named workflows
- The module may define an entry workflow to override the default behavior of `module run`.
- The module may define any number of functions

Local modules can define any number of processes, workflows, and functions. As a best practice, each process and named workflow should be defined in its own script.

### meta.yml

The `meta.yml` file contains the module's metadata, including its name, version, description, authors, and input/output specifications. The registry uses this file to display module information and generate usage templates.

### README.md

The `README.md` file provides documentation for the module. It should describe what the module does, the tools it wraps, and any configuration requirements.

(module-resources)=

### resources

Modules can include resource files in the `resources/` directory.

When running the module with {ref}`wave-page`, the contents of `resources/` are mounted into the root directory of the task container.

For example, given a module with the following structure:

```
my-module/
├── main.nf
└── resources/
    ├── data/
    |   └── file.txt
    └── usr/
        └── bin/
            └── hello.sh
```

The process script can use these files as follows:

```nextflow
process hello {
    container 'quay.io/nextflow/bash'

    script:
    """
    cat /data/file.txt
    hello.sh
    """
}
```

Module resources can be used without Wave or containerization, with the following limitations:

- The `nextflow.enable.moduleBinaries` feature flag must be enabled in the pipeline script.

- The pipeline work directory must be in a local or shared file system. Remote object storage is not supported without Wave.

- Only executable scripts in `resources/usr/bin/` are made accessible to the process script.

### templates

Modules can include process script {ref}`templates <process-template>` in the `templates/` directory.

For example, given a module with the following structure:

```
my-module/
|── main.nf
└── templates/
    └── hello.sh
```

The module's process can use the script template as follows:

```nextflow
process hello {
    input:
    val STR

    script:
    template 'hello.sh'
}
```

(dev-modules-spec)=

## Generating a module spec

Use the `module spec` command to generate or update the `meta.yml` file from the module's `main.nf`:

```console
$ nextflow module spec myorg/my-module
```

Use `-dry-run` to preview the generated spec without writing to disk:

```console
$ nextflow module spec -dry-run myorg/my-module
```

When generating the module spec for the first time, provide required fields directly to avoid `TODO` placeholders in the generated file:

```console
$ nextflow module spec \
    -namespace myorg \
    -version 1.0.0 \
    -description "Quality control of raw sequencing reads" \
    -license MIT \
    -author "@myname" \
    ./modules/myorg/my-module
```

When updating an existing module spec, it is incorporated into the new file.

See {ref}`module spec <cli-module-spec>` for the full command reference.

(dev-modules-validate)=

## Validating a module

Use the `module validate` command to check that a module is ready for publishing:

```console
$ nextflow module validate myorg/my-module
```

The command verifies that:

- All required files are present (`main.nf`, `meta.yml`, `README.md`).
- The module spec contains all required fields (name, version, description, and license).

See {ref}`module validate <cli-module-validate>` for the full command reference.

(dev-modules-testing)=

## Testing a module

Before publishing, test your module by running it directly:

```console
$ nextflow module run myorg/my-module --input 'test-data/*.fastq.gz'
```

The command executes the module as a standalone run, allowing you to verify that inputs are correctly declared, the process runs successfully, and the correct outputs are produced.

For more thorough testing, create a small wrapper workflow that exercises the module:

```nextflow
include { MY_MODULE } from './modules/myorg/my-module'

workflow {
    input_ch = channel.fromPath('test-data/*.fastq.gz')
    results_ch = MY_MODULE(input_ch)
    results_ch.view()
}
```

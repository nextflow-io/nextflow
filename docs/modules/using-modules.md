(using-modules-page)=

# Using modules

:::{versionadded} 26.04.0
:::

The Nextflow module system allows you to discover, install, and manage reusable modules from centralized registries.
This page describes how to use modules in your pipelines.

## Discovering modules

Search for available modules using the `module search` command:

```console
$ nextflow module search alignment
$ nextflow module search "quality control" -limit 10
```

Results include module names and descriptions.
Use `-output json` for machine-readable output.

See {ref}`cli-module-search` for the full command reference.

## Installing modules

Use the `module install` command to download modules from a registry into your project:

```console
$ nextflow module install nf-core/fastqc
$ nextflow module install nf-core/fastqc -version 0.0.0-0c7146d
```

:::{note}
Modules mirrored from nf-core do not follow standard semantic versioning.
Instead, they use the format `0.0.0-<hash>`, where the suffix is a short portion of the nf-core module's commit hash.
:::

Nextflow stores installed modules in the `modules/` directory and creates a `.module-info` file alongside the module to record installation metadata such as the module checksum and registry URL.

:::{tip}
Commit the `modules/` directory to your Git repository to ensure reproducibility.
:::

See {ref}`cli-module-install` for the full command reference.

## Listing installed modules

View all modules installed in your project with the `module list` command:

```console
$ nextflow module list
```

The output shows each module's name, installed version, and whether it has been modified locally.
Use `-output json` for machine-readable output.

See {ref}`cli-module-list` for the full command reference.

## Viewing module information

Use the `module view` command to view metadata and a usage template for a module:

```console
$ nextflow module view nf-core/fastqc
$ nextflow module view nf-core/fastqc -version 0.0.0-0c7146d
```

The output includes the module's version, URL, description, authors, maintainers, keywords, tools, input/output channels, and a generated usage template.
Use `-output json` for machine-readable output.

See {ref}`cli-module-view` for the full command reference.

## Including modules

Modules installed from a registry can be included by name:

```nextflow
include { FASTQC } from 'nf-core/fastqc'

workflow {
    reads = channel.fromPath('data/*.fastq')
    reads = reads.map { fastq -> tuple([id: fastq.baseName], fastq) }
    FASTQC(reads)
}
```

Local modules must be included by relative path:

```nextflow
include { FASTQC } from './modules/local/fastqc'
```

## Running modules directly

For ad-hoc tasks or testing, run a module directly without creating a wrapper workflow:

```console
$ nextflow module run nf-core/fastqc --meta.id test_sample --reads sample1_R1.fastq.gz
```

:::{tip}
Run `nextflow module view` to see the available inputs for a module.
:::

The command automatically downloads the module if it is not already installed.
It accepts all standard Nextflow run options (`-profile`, `-resume`, etc.):

```console
$ nextflow module run nf-core/fastqc \
    --meta.id test_sample \
    --reads sample1_R1.fastq.gz \
    -with-docker
```

Run a local module by specifying a path starting with `./` or `../`:

```console
$ nextflow module run ./modules/local/fastqc/main.nf \
    --meta.id test_sample \
    --reads sample1_R1.fastq.gz \
    -with-docker
```

:::{note}
The local module should define a single process, or else the command will fail.
:::

See {ref}`cli-module-run` for the full command reference.

## Updating modules

To update a module to a newer version, reinstall it with the desired version:

```console
$ nextflow module install nf-core/fastqc -version 0.0.0-c9h0bv4
```

Nextflow automatically verifies module integrity using a checksum stored in the `.module-info` file.
If the module has local modifications, it will not be updated.
Use the `-force` flag to overwrite local changes:

```console
$ nextflow module install nf-core/fastqc -version 0.0.0-c9h0bv4 -force
```

## Removing modules

Use the `module remove` command to uninstall a module from your project:

```console
$ nextflow module remove nf-core/fastqc
```

By default, Nextflow removes both the module files and the `.module-info` file.
Use flags to control this behavior:

- `-keep-files`: Remove the `.module-info` file but keep the module files in the `modules/` directory.
- `-force`: Remove the module directory even if it has no `.module-info` file or has local modifications.

See {ref}`cli-module-remove` for the full command reference.

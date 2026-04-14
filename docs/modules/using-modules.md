(using-modules-page)=

# Using modules

:::{versionadded} 26.04.0
:::

The Nextflow module system allows you to discover, install, and manage reusable modules from centralized registries.
This page describes how to use registry modules in your pipelines.
For local module syntax (inclusion, aliases, templates, and binaries), see {ref}`module-page`.

## Discovering modules

Search for available modules using the `module search` command:

```console
$ nextflow module search alignment
$ nextflow module search "quality control" -limit 10
```

Results include module names, versions, descriptions, and download statistics.
Use `-output json` for machine-readable output.

See {ref}`cli-module-search` for the full command reference.

## Installing modules

Use the `module install` command to download modules from a registry into your project:

```console
$ nextflow module install nf-core/fastqc
$ nextflow module install nf-core/fastqc -version 1.0.0
```

Nextflow stores installed modules in the `modules/` directory and creates a `.module-info` file alongside the module to record installation metadata such as the module checksum and registry URL.

Once installed, include a module by name rather than a relative path:

```nextflow
include { FASTQC } from 'nf-core/fastqc'

workflow {
    reads = Channel.fromFilePairs('data/*_{1,2}.fastq.gz')
    FASTQC(reads)
}
```

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

Use the `module info` command to view metadata and a usage template for a module:

```console
$ nextflow module info nf-core/fastqc
$ nextflow module info nf-core/fastqc -version 1.0.0
```

The output includes the module's description, authors, keywords, tools, input/output channels, and a generated usage template.
Use `-output json` for machine-readable output.

See {ref}`cli-module-info` for the full command reference.

## Running modules directly

For ad-hoc tasks or testing, run a module directly without creating a wrapper workflow:

```console
$ nextflow module run nf-core/fastqc --input 'data/*.fastq.gz'
```

The command automatically downloads the module if it is not already installed.
It accepts all standard Nextflow run options (`-profile`, `-resume`, etc.):

```console
$ nextflow module run nf-core/salmon \
    --reads reads.fq \
    --index salmon_index \
    -profile docker \
    -resume
```

Nextflow infers command-line params (e.g., `--reads`) from the module's declared inputs.
Run `nextflow module info` to see the available inputs for a module.

See {ref}`cli-module-run` for the full command reference.

## Updating modules

To update a module to a newer version, reinstall it with the desired version:

```console
$ nextflow module install nf-core/fastqc -version 1.1.0
```

If you have local modifications, Nextflow warns you and prevents the update.
Use the `-force` flag to override:

```console
$ nextflow module install nf-core/fastqc -version 1.1.0 -force
```

## Checksum verification

Nextflow automatically verifies module integrity using checksums stored in the `.module-info` file.
When you modify a module locally, Nextflow detects the change and prevents accidental overwrites during reinstallation.

Module integrity ensures that modules remain consistent with their registry versions unless you explicitly choose to override them.

## Removing modules

Use the `module remove` command to uninstall a module from your project:

```console
$ nextflow module remove nf-core/fastqc
```

By default, Nextflow removes both the module files and the `.module-info` file.
Use the following flags to control this behavior:

- `-keep-files`: Remove the `.module-info` file but keep the module files in the `modules/` directory.
- `-force`: Force removal even if the module has no `.module-info` file or has local modifications.

See {ref}`cli-module-remove` for the full command reference.

(modules-page)=

# Overview

Nextflow modules are reusable units of pipeline logic (i.e., processes, workflows, and functions) that you can share across projects.
By packaging common tasks as modules, you avoid duplicating code and benefit from community improvements.

There are two ways to use modules in Nextflow:

- [Local modules](#local-modules)
- [Registry modules](#registry-modules)

## Local modules

Local modules are Nextflow scripts stored directly in your project. You include definitions from local modules using the `include` keyword with a relative path:

```nextflow
include { FASTQC } from './modules/fastqc'
```

Local modules are well suited for project-specific logic that is not intended for reuse or sharing. See {ref}`module-page` for details on module inclusion syntax, aliases, templates, and binaries.

## Registry modules

:::{versionadded} 26.04.0
:::

The [Nextflow module registry](https://registry.nextflow.io) hosts registry modules in a centralized repository.
You manage them with the `nextflow module` command.
Registry modules follow a standard structure with metadata (`meta.yml`), documentation (`README.md`), and a module script (`main.nf`) for version management, integrity checking, and discoverability.

Key features of registry modules:

- **Discoverability**: Search for modules by keyword or name and browse available versions.
- **Version management**: Pin specific versions or use the latest release, with automatic integrity checking via checksums.
- **Direct execution**: Run modules as standalone workflows for ad-hoc tasks or testing without writing a wrapper script.
- **Standard structure**: Each module includes a script (`main.nf`), metadata (`meta.yml`), and documentation (`README.md`), enabling consistent tooling and automation.

:::{note}
Modules from the [nf-core](https://nf-co.re/) community are automatically mirrored to the Nextflow module registry under the `nf-core` namespace.
You can install and use them directly without any additional configuration.
:::

Install registry modules into your project and include them by name:

```console
$ nextflow module install nf-core/fastqc
```

```nextflow
include { FASTQC } from 'nf-core/fastqc'
```

For more information about registry modules, see:

- {ref}`using-modules-page` for details on discovering, installing, and managing registry modules.
- {ref}`dev-modules-page` for details on creating and publishing your own modules.
- {ref}`module-registry-page` for details on namespaces, access tokens, and registry configuration.

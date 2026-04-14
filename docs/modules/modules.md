(modules-page)=

# Overview

## What are modules

Nextflow modules are reusable units of pipeline logic -- processes, workflows, and functions -- that you can share across projects.
By packaging common tasks as modules, you avoid duplicating code and benefit from community improvements.

There are two ways to use modules in Nextflow: local modules and registry modules.

<h3>Local modules</h3>

Local modules are Nextflow scripts stored directly in your project. You include definitions from local modules using the `include` keyword with a relative path:

```nextflow
include { FASTQC } from './modules/fastqc'
```

Local modules are well suited for project-specific logic that is not intended for reuse across projects. See {ref}`module-page` for details on module inclusion syntax, aliases, templates, and binaries.

<h3>Registry modules</h3>

:::{versionadded} 26.04.0
:::

Centralized registries host registry modules, and you manage them with the `nextflow module` command.
They follow a standard structure with metadata (`meta.yml`), documentation (`README.md`), and a module script (`main.nf`) for version management, integrity checking, and discoverability.

Install registry modules into your project and include them by name:

```console
$ nextflow module install nf-core/fastqc
```

```nextflow
include { FASTQC } from 'nf-core/fastqc'
```

See {ref}`using-modules-page` for details on discovering, installing, and managing registry modules.

## Versioning

Registry modules use version strings to identify releases. Many modules use [Semantic Versioning](https://semver.org/):

* MAJOR: Increment for backward-incompatible changes.
* MINOR: Increment for backward-compatible feature additions.
* PATCH: Increment for backward-compatible bug fixes.

When installing a module, you can pin a specific version or let Nextflow use the latest available version:

```console
$ nextflow module install nf-core/fastqc -version 1.0.0
```

See {ref}`cli-module-install` for more information.

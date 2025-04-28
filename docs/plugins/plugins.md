(plugins-page)=

# Overview

## What are plugins

Nextflow plugins are extensions that enhance the functionality of the Nextflow workflow framework. They allow users to add new capabilities and integrate with external services without modifying core Nextflow code.

There are two types of Nextflow plugins: 

* Core plugins
* Third-party plugins

The main features of each plugin type are described below.

<h3>Core plugins</h3>

Core plugins do not require configuration. The latest versions of core plugins are automatically installed when a Nextflow pipeline requests them. Core plugins include:

* `nf-amazon`: Support for Amazon Web Services.
* `nf-azure`: Support for Microsoft Azure.
* `nf-cloudcache`: Support for the cloud cache.
* `nf-console`: Implement Nextflow [REPL console](https://seqera.io/blog/introducing-nextflow-console/).
* `nf-google`: Support for Google Cloud.
* `nf-tower`: Support for [Seqera Platform](https://seqera.io/platform/).
* `nf-wave`: Support for [Wave containers service](https://seqera.io/wave/).

Specific versions of core plugins can be declared in Nextflow configuration files or by using the `-plugins` option. See {ref}`using-plugins-page` for more information.

:::{note}
The automatic application of core plugins can be disabled by setting `NXF_PLUGINS_DEFAULT=false`. See {ref}`dev-plugins-env-var` for more information.
:::

<h3>Third-party plugins</h3>

Third-party plugins must be configured via Nextflow configuration files or at runtime. To configure a plugin via configuration files, use the `plugins` block. For example:

```
plugins {
    id 'nf-hello@0.5.0'
}
```

To configure plugins at runtime, use the `-plugins` option. For example:

```
nextflow run main.nf -plugins nf-hello@0.5.0
```

See {ref}`using-plugins-page` for more information.

## Nextflow plugin registry

The Nextflow plugin registry is a centralized repository of assembled plugins. It hosts an index of plugin metadata that supports plugin discovery, accessibility, and version tracking. The registry is integrated with the Nextflow runtime. Nextflow is able to automatically locate and download configured plugins. Developers can upload plugins to the registry, with built-in support from the Gradle plugin for Nextflow plugins.

## Versioning

Nextflow plugins follow the Semantic Versioning format: MAJOR.MINOR.PATCH (e.g., 0.5.0). This helps developers communicate the nature of changes in a release.

Version components:

* MAJOR: Increment for incompatible changes.
* MINOR: Increment for backward-compatible feature additions.
* PATCH: Increment for backward-compatible bug fixes.

Optional pre-release and build metadata can be added (e.g., 1.2.1-alpha+001) as extensions to the base version format.
See {ref}`dev-plugins-versioning` for more information about the specification.

(plugins-page)=

# Overview

## What are plugins

Nextflow plugins are extensions that enhance the functionality of Nextflow. They allow users to add new capabilities and integrate with external services without bloating their pipeline code or modifying Nextflow itself.

There are two types of Nextflow plugins: core plugins and third-party plugins. The main features of each plugin type are described below.

<h3>Core plugins</h3>

Core plugins do not require configuration. The latest versions of core plugins are automatically installed when a Nextflow pipeline requests them. Core plugins include:

* `nf-amazon`: Support for Amazon Web Services.
* `nf-azure`: Support for Microsoft Azure.
* `nf-cloudcache`: Support for the cloud cache.
* `nf-console`: Implementation of the Nextflow [REPL console](https://seqera.io/blog/introducing-nextflow-console/).
* `nf-google`: Support for Google Cloud.
* `nf-tower`: Support for [Seqera Platform](https://seqera.io/platform/).
* `nf-wave`: Support for [Wave containers service](https://seqera.io/wave/).

Specific versions of core plugins can be declared in Nextflow configuration files or by using the `-plugins` option. See {ref}`using-plugins-page` for more information.

:::{note}
The automatic retrieval of core plugins can be disabled by setting `NXF_PLUGINS_DEFAULT=false`. See {ref}`dev-plugins-env-vars` for more information.
:::

<h3>Third-party plugins</h3>

Third-party plugins must be configured via Nextflow configuration files or at runtime. To configure a plugin via configuration files, use the `plugins` block. For example:

```groovy
plugins {
    id 'nf-hello@0.5.0'
}
```

To configure plugins at runtime, use the `-plugins` option. For example:

```bash
nextflow run main.nf -plugins nf-hello@0.5.0
```

See {ref}`using-plugins-page` for more information.

## Versioning

Nextflow plugins are free to use any versioning convention. Many plugins, including all core plugins, use [Semantic Versioning](https://semver.org/), which helps developers communicate the kind of changes in a release.

Semantic versions have the form MAJOR.MINOR.PATCH (e.g., 0.5.0):

* MAJOR: Increment for backward-incompatible changes.
* MINOR: Increment for backward-compatible feature additions.
* PATCH: Increment for backward-compatible bug fixes.

Optional pre-release and build metadata can be added (e.g., 1.2.1-alpha+001) as extensions to the base version format.

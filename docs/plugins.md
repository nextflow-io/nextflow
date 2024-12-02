(plugins-page)=

# Plugins

Nextflow has a plugin system that allows the use of extensible components that are downloaded and installed at runtime.

(plugins-core)=

## Core plugins

The following functionalities are provided via plugin components, and they make part of the Nextflow *core* plugins:

- `nf-amazon`: Support for Amazon Web Services.
- `nf-azure`: Support for Microsoft Azure.
- `nf-cloudcache`: Support for the cloud cache (see `NXF_CLOUDCACHE_PATH` under {ref}`config-env-vars`).
- `nf-console`: Implement Nextflow [REPL console](https://www.nextflow.io/blog/2015/introducing-nextflow-console.html).
- `nf-google`: Support for Google Cloud.
- `nf-tower`: Support for [Seqera Platform](https://seqera.io) (formerly Tower Cloud).
- `nf-wave`: Support for [Wave containers](https://seqera.io/wave/) service.

## Using plugins

The core plugins do not require any configuration. They are automatically installed when the corresponding feature is requested by a Nextflow pipeline. You can still specify them as described below, e.g. if you want to pin the version of a plugin, however if you try to use a plugin version that isn't compatible with your Nextflow version, Nextflow will fail.

You can enable a plugin by declaring it in your Nextflow configuration:

```groovy
plugins {
    id 'nf-hello@0.1.0'
}
```

Or you can use the `-plugins` command line option:

```bash
nextflow run <pipeline> -plugins nf-hello@0.1.0
```

The plugin identifier consists of the plugin name and plugin version separated by a `@`. Multiple plugins can be specified in the configuration with multiple `id` declarations, or on the command line as a comma-separated list. When specifying plugins via the command line, any plugin declarations in the configuration file are ignored.

The plugin version is optional. If it is not specified, Nextflow will download the latest plugin version that is compatible with your Nextflow version. In general, it recommended that you not specify the plugin version unless you actually want to stick to that version, such as for [offline usage](#offline-usage).

The core plugins are documented in this documentation. For all other plugins, please refer to the plugin's code repository for documentation and support.

## Offline usage

To use Nextflow plugins in an offline environment:

1. {ref}`Install Nextflow <install-nextflow>` on a system with an internet connection. Do not use the "all" package, as this does not allow the use of custom plugins.

2. Download any additional plugins by running `nextflow plugin install <pluginId,..>`. Alternatively, simply run your pipeline once and Nextflow will download all of the plugins that it needs.

3. Copy the `nextflow` binary and `$HOME/.nextflow` folder to your offline environment.

4. In your Nextflow configuration file, specify each plugin that you downloaded, both name and version, including default plugins. This will prevent Nextflow from trying to download newer versions of plugins.

Nextflow caches the plugins that it downloads, so as long as you keep using the same Nextflow version and pin your plugin versions in your config file, Nextflow will use the locally installed plugins and won't try to download them from the Internet.

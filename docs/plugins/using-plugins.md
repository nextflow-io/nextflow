(using-plugins-page)=

# Using plugins

Nextflow core plugins require no additional configuration. When a pipeline uses a core plugin, Nextflow automatically downloads and uses the latest compatible plugin version. In contrast, third-party plugins must be explicitly declared. When a pipeline uses one or more third-party plugins, Nextflow must be configured to download and use the plugin.

(using-plugins-identifiers)=

## Identifiers

A plugin identifier consists of the plugin name and version, separated by an `@` symbol:

```
nf-hello@0.5.0
```

The plugin version is optional. If it is not specified, Nextflow will download the latest version of the plugin that meets the minimum Nextflow version requirements specified by the plugin.

:::{note}
Plugin versions are required for {ref}`offline usage <using-plugins-offline>`.
:::

:::{versionadded} 25.04.0
:::

The plugin version can be prefixed with `~` to pin the major and minor versions and allow the latest patch release to be used. For example, `nf-amazon@~2.9.0` will resolve to the latest version matching `2.9.x`. When working offline, Nextflow will resolve version ranges against the local plugin cache defined by `NXF_PLUGINS_DIR`.

:::{tip}
It is recommended to pin the major and minor version of each plugin to minimize the risk of breaking changes while allowing patch updates to be used automatically.
:::

(using-plugins-config)=

## Configuration

Plugins can be configured via Nextflow configuration files or at runtime.

To configure a plugin via configuration files, use the `plugins` block. For example:

```nextflow
plugins {
    id 'nf-hello@0.5.0'
    id 'nf-amazon@2.9.0'
}
```

To configure plugins at runtime, use the `-plugins` option. For example:

```bash
nextflow run main.nf -plugins nf-hello@0.5.0,nf-amazon@2.9.0
```

:::{note}
Plugin declarations in Nextflow configuration files are ignored when specifying plugins via the `-plugins` option.
:::

## Caching

When Nextflow downloads plugins, it caches them in the directory specified by `NXF_PLUGINS_DIR` (`$HOME/.nextflow/plugins` by default).

(using-plugins-offline)=

## Offline usage

When running Nextflow in an offline environment, any required plugins must be downloaded and moved into the offline environment prior to any runs.

To use Nextflow plugins in an offline environment:

1. Install a self-contained version of Nextflow in an environment with an internet connection. See {ref}`install-standalone` for more information.

2. Run `nextflow plugin install <PLUGIN_NAME>@<VERSION>` for each required plugin to download it. Alternatively, run the pipeline once, which will automatically download all plugins required by the pipeline.

3. Copy the `nextflow` binary and `$HOME/.nextflow` directory to the offline environment.

4. Specify each plugin and its version in Nextflow configuration files or at runtime. See {ref}`using-plugins-config` for more information.

    :::{warning}
    Nextflow will attempt to download newer versions of plugins if their versions are not set. See {ref}`using-plugins-identifiers` for more information.
    :::

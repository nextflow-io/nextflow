(using-plugins-page)=

# Using plugins

Nextflow core plugins require no additional configuration. When a pipeline uses a core plugin, Nextflow automatically downloads and uses the latest compatible plugin version. In contrast, third-party plugins must be explicitly declared. When a pipeline uses one or more third-party plugins, Nextflow must be configured to download and use the plugin.

(using-plugins-identifiers)=

## Identifiers

A plugin identifier consists of the plugin name and version, separated by an `@` symbol:

```console
nf-hello@0.5.0
```

The plugin version is optional. If it is not specified, Nextflow will download the latest version of the plugin that meets the minimum Nextflow version requirements specified by the plugin.

:::{note}
Plugin versions are required for {ref}`offline usage <using-plugins-offline>`.
:::

:::{versionadded} 25.02.0-edge.
:::

The plugin version can be prefixed with `~` to pin the major and minor versions and allow the latest patch release to be used. For example, `nf-amazon@~2.9.0` will resolve to the latest version matching `2.9.x`. When working offline, Nextflow will resolve version ranges against the local plugin cache defined by `NXF_PLUGINS_DIR`.

:::{tip}
It is recommended to pin plugin versions to the major and minor versions and allow the latest patch update.
:::

(using-plugins-config)=

## Configuration

Plugins can be configured via nextflow configuration files or at runtime.

To configure a plugin via configuration files, use the `plugins` block. For example:

```nextflow
plugins {
    id 'nf-hello@0.5.0'
}
```

To configure plugins at runtime, use the `-plugins` option. For example:

```bash
nextflow run main.nf -plugins nf-hello@0.5.0
```

To configure multiple plugins at runtime, use the `-plugins` option and a comma-separated list. For example:

```bash
nextflow run main.nf -plugins nf-hello@0.5.0,nf-amazon@2.9.0
```

:::{note}
Plugin declarations in nextflow configuration files are ignored when specifying plugins via the `-plugins` option.
:::

## Caching

When Nextflow downloads plugins, it caches them in the directory specified by `NXF_PLUGINS_DIR` (`$HOME/.nextflow/plugins` by default). It tracks the plugins and their versions in a local cache located at `.nextflow/cache/<SESSION_ID>`. If the Nextflow version and plugin versions match those from a previous run, the cached plugins are reused.

(using-plugins-offline)=

## Offline usage

Nextflow plugins may be required by some pipelines in an offline environment. Plugins must be manually downloaded and moved to the offline environment.

To use Nextflow plugins in an offline environment:

1. Install a self-contained version of Nextflow in an environment with an internet connection. See {ref}`install-standalone` for more information.
2. Run `nextflow plugin install <PLUGIN_NAME>@<VERSION>` to download the plugins.

    :::{tip}
    Running a pipeline in an environment with an internet connection will also download the plugins used by that pipeline.
    :::

3. Copy the `nextflow` binary and `$HOME/.nextflow` directory to the offline environment.
4. Specify each plugin and its version in Nextflow configuration files or at runtime. See {ref}`using-plugins-config` for more information.

    :::{warning}
    Nextflow will attempt to download newer versions of plugins if their versions are not set. See {ref}`using-plugins-identifiers` for more information.
    :::
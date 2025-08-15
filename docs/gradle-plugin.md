(gradle-plugin-page)=

# Using the Nextflow Gradle plugin

The [Nextflow Gradle plugin](https://github.com/nextflow-io/nextflow-plugin-gradle) simplifies plugin development by configuring default dependencies needed for Nextflow integration and defining Gradle tasks for building, testing, and publishing Nextflow plugins.

:::{note}
Nextflow plugins can be developed without the Gradle plugin. However, this approach is only suggested if you are an advanced developer and your project is incompatible with the Gradle plugin.
:::

(gradle-plugin-create)=

## Creating a plugin

:::{versionadded} 25.04.0
:::

The easiest way to get started with the Nextflow Gradle plugin is to use the `nextflow plugin create` sub-command, which creates a plugin project based on the [Nextflow plugin template](https://github.com/nextflow-io/nf-plugin-template/), which in turn uses the Gradle plugin.

To create a Nextflow plugin with the Gradle plugin, simply run `nextflow plugin create` on the command line. It will prompt you for your plugin name, organization name, and project path.

See {ref}`dev-plugins-extension` for more information about available plugin extension points.

## Building a plugin

To build a plugin, simply run `make assemble`.

Plugins can also be installed locally without being published. To install a plugin locally:

1. In the plugin root directory, run `make install`.

    :::{note}
    Running `make install` will add your plugin to your `$HOME/.nextflow/plugins` directory.
    :::

2. Run your pipeline:

    ```bash
    nextflow run main.nf -plugins <PLUGIN_NAME>@<VERSION>
    ```

    :::{note}
    Plugins can also be configured via nextflow configuration files. See {ref}`using-plugins-page` for more information.
    :::

(gradle-plugin-test)=

## Testing a plugin

<h3> Unit tests </h3>

Unit tests are small, focused tests designed to verify the behavior of individual plugin components.

To run unit tests:

1. Develop your unit tests. See [MyObserverTest.groovy](https://github.com/nextflow-io/nf-plugin-template/blob/main/src/test/groovy/acme/plugin/MyObserverTest.groovy) in [nf-plugin-template](https://github.com/nextflow-io/nf-plugin-template) for unit test examples.

2. In the plugin root directory, run `make test`.

<h3> End-to-end tests </h3>

End-to-end tests are comprehensive tests that verify the behavior of an entire plugin as it would be used in a Nextflow pipeline. See [nf-hello](https://github.com/nextflow-io/nf-hello) for an example of an end-to-end test for a Nextflow plugin.

(gradle-plugin-publish)=

## Publishing a plugin

The Nextflow Gradle plugin allows you to publish your plugin to the Nextflow plugin registry from the command line.

:::{note}
The Nextflow plugin registry is currently available as a private beta. Contact [info@nextflow.io](mailto:info@nextflow.io) for more information.
:::

To publish your plugin:

1. Create a file named `$HOME/.gradle/gradle.properties`, where `$HOME` is your home directory.

2. Add the following properties:

    ```bash
    pluginRegistry.accessToken=<REGISTRY_ACCESS_TOKEN>
    ```

    Replace `<REGISTRY_ACCESS_TOKEN>` with your plugin registry access token.

3. Run `make release`.

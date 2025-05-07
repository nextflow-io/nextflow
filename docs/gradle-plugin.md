(gradle-plugin-page)=

# Using the Gradle plugin

The [Gradle plugin for Nextflow plugins](https://github.com/nextflow-io/nextflow-plugin-gradle) simplifies plugin development by configuring default dependencies needed for Nextflow integration and incorporates custom Gradle tasks that streamline building, testing, and publishing Nextflow plugins. The [`nf-plugin-template`](https://github.com/nextflow-io/nf-plugin-template/) is a scaffold for plugin development and incorporates the Gradle plugin by default. You can use the `nextflow plugin create` sub-command to create plugins scaffolds with the [`nf-plugin-template`](https://github.com/nextflow-io/nf-plugin-template/) scaffold.

:::{note}
Nextflow Plugins can be developed without the Gradle plugin. However, this approach is only suggested if you are an advanced developer and your project is incompatible with the Gradle plugin.
:::

(gradle-plugin-create)=

## Creating a plugin

:::{versionadded} 25.04.0
:::

To create a Nextflow plugin with the Gradle plugin:

1. Run `nextflow plugin create`.
    - When prompted `Enter plugin name:`, enter your plugin name.
    - When prompted `Enter organization:`, enter your organization.
    - When prompted `Enter project path:`, enter the project path in your local file system.
    - When prompted `All good, are you OK to continue [y/N]?`, enter `y`.
2. Develop your plugin extension points. See {ref}`dev-plugins-extension` for descriptions and examples.
3. In the plugin root directory, run `make assemble`.

(gradle-plugin-install)=

## Installing a plugin

Plugins can be installed locally without being packaged, uploaded, and published.

To install a plugin locally:

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


(gradle-plugin-unit-test)=

## Unit testing a plugin

Unit tests are small, focused tests designed to verify the behavior of individual plugin components and are an important part of software development.

To run unit tests:

1. Develop your unit tests. See [MyObserverTest.groovy](https://github.com/nextflow-io/nf-plugin-template/blob/main/src/test/groovy/acme/plugin/MyObserverTest.groovy) in [nf-plugin-template](https://github.com/nextflow-io/nf-plugin-template/tree/main) for unit test examples.
2. In the plugin root directory, run `make test`.

(gradle-plugin-package)=

## Packaging, uploading, and publishing a plugin

The Gradle plugin for Nextflow plugins simplifies publishing your plugin to the Nextflow Plugin Registry.

:::{note}
The Nextflow Plugin Registry is currently available as private beta technology. Contact [info@nextflow.io](mailto:info@nextflow.io) to learn how to get access.
:::

To package, upload, and publish your plugin to the Nextflow Plugin Registry:

1. Create a file named `$HOME/.gradle/gradle.properties`, where `$HOME` is your home directory.
2. Add the following properties:

    ```bash
    pluginRegistry.accessToken=<REGISTRY_ACCESS_TOKEN>
    ```

    Replace <REGISTRY_ACCESS_TOKEN> with your plugin registry access token.

3. Run `make release`.
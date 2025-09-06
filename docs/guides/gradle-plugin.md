(gradle-plugin-page)=

# Using the Nextflow Gradle plugin

Nextflow provides a new `plugin create` command that simplifies the creation of Nextflow plugins. This command leverages the [nf-plugin-template](https://github.com/nextflow-io/nf-plugin-template) project, which uses the [Nextflow Gradle plugin](https://github.com/nextflow-io/nextflow-plugin-gradle) to streamline plugin development.

The [Nextflow Gradle plugin](https://github.com/nextflow-io/nextflow-plugin-gradle) configures default dependencies needed for Nextflow integration and defines Gradle tasks for building, testing, and publishing Nextflow plugins. The Gradle plugin is versioned and published to the [Gradle Plugin Portal](https://plugins.gradle.org/), allowing developers to manage it like any other dependency. As the plugin ecosystem evolves, the Gradle plugin will enable easier maintenance and adoption of improvements. This page introduces [Nextflow Gradle plugin](https://github.com/nextflow-io/nextflow-plugin-gradle) and how to use it.

:::{note}
The Nextflow Gradle plugin and plugin registry are currently available as a public preview. See the {ref}`Migrating to the Nextflow plugin registry <plugin-registry-page>` for more information.
:::

(gradle-plugin-create)=

## Creating a plugin

:::{versionadded} 25.04.0
:::

The best way to create a plugin with the Nextflow Gradle plugin is to use the `nextflow plugin create` sub-command and create a plugin project based on the [Nextflow plugin template](https://github.com/nextflow-io/nf-plugin-template/). See {ref}`dev-plugins-template` for more information about the Nextflow plugin template.

To create Nextflow plugins with the Gradle plugin:

1. Run `nextflow plugin create`.

2. Follow the prompts to add your plugin name, organization name, and project path.

3. Develop your plugin extension points. See {ref}`dev-plugins-extension-points` for more information.

4. Develop your tests. See {ref}`gradle-plugin-test` for more information.

5. In the plugin root directory, run `make assemble`.

## Installing a plugin

Plugins can be installed locally without being published.

To install plugins locally:

1. In the plugin root directory, run `make install`.

    :::{note}
    Running `make install` will add your plugin to your `$HOME/.nextflow/plugins` directory.
    :::

2. Run your pipeline:

    ```bash
    nextflow run main.nf -plugins <PLUGIN_NAME>@<VERSION>
    ```

    Replace `<PLUGIN_NAME>@<VERSION>` with your plugin name and version.

    :::{note}
    Plugins can also be configured via Nextflow configuration files. See {ref}`using-plugins-page` for more information.
    :::

(gradle-plugin-test)=

## Testing a plugin

Testing your Nextflow plugin requires unit tests and end-to-end tests.

<h3>Unit tests</h3>

Unit tests are small, focused tests designed to verify the behavior of individual plugin components.

To run unit tests:

1. Develop your unit tests. See [MyObserverTest.groovy](https://github.com/nextflow-io/nf-plugin-template/blob/main/src/test/groovy/acme/plugin/MyObserverTest.groovy) in the [plugin template](https://github.com/nextflow-io/nf-plugin-template) for an example unit test.

2. In the plugin root directory, run `make test`.

<h3>End-to-end tests</h3>

End-to-end tests are comprehensive tests that verify the behavior of an entire plugin as it would be used in Nextflow pipelines.

End-to-end tests should be tailored to the needs of your plugin, but generally take the form of a small Nextflow pipeline. See the `validation` directory in the [plugin template](https://github.com/nextflow-io/nf-plugin-template) for an example end-to-end test.

(gradle-plugin-publish)=

## Publishing a plugin

The Nextflow Gradle plugin allows you to publish plugins to the [Nextflow plugin registry](https://registry.nextflow.io/) from the command line.

To publish plugins to the [Nextflow plugin registry](https://registry.nextflow.io/):

1. Create a file named `$HOME/.gradle/gradle.properties`, where `$HOME` is your home directory.

2. Add the following property:

    ```
    npr.apiKey=<API_KEY>
    ```

    Replace `<API_KEY>` with your plugin registry API key. See {ref}`plugin-registry-access-token` for more information about creating an API key.

3. Run `make release`.

## Additional resources

For additional Nextflow Gradle plugin configuration options, see the [Nextflow Gradle plugin](https://github.com/nextflow-io/nextflow-plugin-gradle) repository.

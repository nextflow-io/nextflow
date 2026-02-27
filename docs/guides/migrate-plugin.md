(migrate-plugin-page)=

# Migrating to the Nextflow plugin registry

The [Nextflow plugin registry](https://registry.nextflow.io/) is a central repository for Nextflow plugins. It hosts an index of plugin metadata that supports plugin discovery, accessibility, and version tracking.

The [legacy plugin index](https://github.com/nextflow-io/plugins) is deprecated in favor of the Nextflow plugin registry. Starting with Nextflow 25.10, the plugin registry is the only method for downloading plugins, replacing the legacy plugin index.

This page describes the registry's impact on users and developers, and how to migrate existing plugins.

## Migration impact on plugin users

No immediate actions are required for plugin users. The plugin configuration has not changed.

## Impact on plugin developers

Plugin developers will need to update their plugin to publish to the Nextflow plugin registry instead of the legacy plugin index. The best way to do this is to migrate to the {ref}`Nextflow Gradle plugin<gradle-plugin-page>`, which simplifies the development process and supports publishing to the plugin registry from the command line.

## Migrating to the Nextflow Gradle plugin

To migrate an existing Nextflow plugin to the {ref}`Nextflow Gradle plugin<gradle-plugin-page>`:

1. Remove the following files and folders:
    - `buildSrc/`
    - `launch.sh`
    - `plugins/build.gradle`

2. If your plugin has a `plugins` directory, move the `src` directory to the project root.

    :::{note}
    Plugin sources should be in `src/main/groovy` or `src/main/java`.
    :::

3. Replace the contents of `settings.gradle` with the following:

    ```groovy
    rootProject.name = '<PLUGIN_NAME>'
    ```

    Replace `<PLUGIN_NAME>` with your plugin name.

4. In the project root, create a new `build.gradle` file with the following configuration:

    ```groovy
    // Plugins
    plugins {
        id 'io.nextflow.nextflow-plugin' version '1.0.0-beta.6'
    }

    // Dependencies (optional)
    dependencies {
        <DEPENDENCY>
    }

    // Plugin version
    version = '<PLUGIN_VERSION>'

    nextflowPlugin {
        // Minimum Nextflow version
        nextflowVersion = '<MINIMUM_NEXTFLOW_VERSION>'

        // Plugin metadata
        provider = '<PROVIDER>'
        className = '<CLASS_NAME>'
        extensionPoints = [
            '<EXTENSION_POINT>'
        ]

    }
    ```

    Replace the following:

    - `<DEPENDENCY>`: (Optional) Your plugins dependency libraries. For example, `commons-io:commons-io:2.18.0`.
    - `<PLUGIN_VERSION>:` Your plugin version. For example, `0.5.0`.
    - `<MINIMUM_NEXTFLOW_VERSION>`: The minimum Nextflow version required to run your plugin. For example, `25.04.0`.
    - `<PROVIDER>`: Your name or organization. For example, `acme`.
    - `<CLASS_NAME>`: Your plugin class name. For example, `acme.plugin.MyPlugin`.
    - `<EXTENSION_POINT>`: Your extension point identifiers that the plugin will implement or expose. For example, `acme.plugin.MyFactory`.

5. Replace the contents of `Makefile` with the following:

    ```Makefile
    # Build the plugin
    assemble:
        ./gradlew assemble

    clean:
        rm -rf .nextflow*
        rm -rf work
        rm -rf build
        ./gradlew clean

    # Run plugin unit tests
    test:
        ./gradlew test

    # Install the plugin into local nextflow plugins dir
    install:
        ./gradlew install

    # Publish the plugin
    release:
        ./gradlew releasePlugin
    ```

6. Update `README.md` with information about the structure of your plugin.

7. In the plugin root directory, run `make assemble`.

Alternatively, use the `nextflow plugin create` command to re-create your plugin with the plugin template and add your existing plugin code. See {ref}`dev-plugins-template` for more information about the plugin template.

## Publishing to the Nextflow plugin registry

The Nextflow Gradle plugin supports publishing plugins from the command line. See {ref}`gradle-plugin-publish` for more information.

The legacy plugin index is **closed to new pull requests**. It will be kept up-to-date with the plugin registry and will remain available indefinitely to ensure continued support for older versions of Nextflow.

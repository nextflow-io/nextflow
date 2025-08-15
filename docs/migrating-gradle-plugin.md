(migrating-plugin-page)=

# Migrating to the Nextflow plugin registry

The Nextflow plugin ecosystem is evolving to support a more robust and user-friendly experience by simplifying the development, publishing, and discovery of Nextflow plugins. This page introduces the Nextflow plugin registry, the Nextflow Gradle plugin, and how to migrate to them.

## Overview

### Nextflow plugin registry

The Nextflow plugin registry is a central repository for Nextflow plugins. It hosts an index of plugin metadata that supports plugin discovery, accessibility, and version tracking. The registry is integrated with the Nextflow runtime. It is intended as a replacement for the [plugins index](https://github.com/nextflow-io/plugins) hosted on GitHub.

:::{note}
The Nextflow plugin registry is currently available as a private beta. Contact [info@nextflow.io](mailto:info@nextflow.io) for more information.
:::

### Nextflow Gradle plugin

The [Nextflow Gradle plugin](https://github.com/nextflow-io/nextflow-plugin-gradle) simplifies the development of Nextflow plugins. It provides default configuration required for Nextflow integration, as well as custom Gradle tasks for building, testing, and publishing plugins.

The Gradle plugin is versioned and published to the [Gradle Plugin Portal](https://plugins.gradle.org/), allowing developers to manage it like any other dependency. As the plugin ecosystem evolves, this Gradle plugin will enable easier maintenance and adoption of ongoing improvements to the Nextflow plugin framework.

## Impact on plugin users

If you are a plugin user, no immediate actions are required. The plugin configuration has not changed.

## Impact on plugin developers

Developers are encouraged to migrate to the Nextflow Gradle plugin in order to take advantage of the simplified development and publishing process.

To migrate an existing Nextflow plugin:

1. Remove the following files and folders:
    - `buildSrc/`
    - `launch.sh`
    - `plugins/build.gradle`

2. If your plugin uses a `plugins` directory, move the `src` directory to the project root.

    :::{note}
    Plugin sources should be in `src/main/groovy` or `src/main/java`.
    :::

3. Replace the contents of `settings.gradle` with the following:

    ```groovy
    rootProject.name = '<PLUGIN_NAME>'
    ```

    Replace `PLUGIN_NAME` with your plugin name.

4. In the project root, create a new `build.gradle` file with the following configuration:

    ```groovy
    // Plugins
    plugins {
        id 'io.nextflow.nextflow-plugin' version '0.0.1-alpha4'
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

        publishing {
            registry {
                authToken = project.findProperty('pluginRegistry.accessToken')
            }
        }
    }
    ```

    Replace the following:

    - `DEPENDENCY`: (Optional) Your plugins dependency libraries—for example, `commons-io:commons-io:2.18.0`.
    - `PLUGIN_VERSION:` Your plugin version—for example, `0.5.0`.
    - `MINIMUM_NEXTFLOW_VERSION`: The minimum Nextflow version required to run your plugin—for example, `25.04.0`.
    - `PROVIDER`: Your name or organization—for example, `acme`.
    - `CLASS_NAME`: Your plugin class name—for example, `acme.plugin.MyPlugin`.
    - `EXTENSION_POINT`: Your extension point identifiers that the plugin will implement or expose—for example, `acme.plugin.MyFactory`.

5. Replace the contents of `Makefile` with the following:

    ```
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

The Nextflow Gradle plugin also supports publishing plugins. See {ref}`gradle-plugin-publish` for more information.

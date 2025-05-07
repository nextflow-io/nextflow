(migrating-plugin-page)=

# Migrating to the Nextflow Gradle plugin

This page introduces the Nextflow Gradle plugin, the Nextflow plugin registry, and how to migrate to the new plugin framework.


## Improvements to the plugin framework

The Nextflow plugin ecosystem is evolving to support a more robust and user-friendly experience by simplifying plugin development, streamlining publishing and discovery, and improving how plugins are loaded into workflows. These improvements make plugins more accessible, maintainable, and interoperable with Nextflow.

### Nextflow Gradle plugin 

The Nextflow Gradle plugin simplifies and standardizes the development of Nextflow plugins. It configures default dependencies required for Nextflow integration and introduces custom Gradle tasks to streamline building, testing, packaging, and publishing plugins.

The Gradle plugin is versioned and published to the [Gradle Plugin Portal](https://plugins.gradle.org/), allowing developers to manage it like any other dependency. As the plugin ecosystem evolves, this Gradle plugin will enable easier maintenance and adoption of ongoing improvements to the Nextflow plugin framework.

### Nextflow Plugin Registry

The Nextflow plugin registry is a centralized repository of assembled plugins. It hosts an index of plugin metadata that supports plugin discovery, accessibility, and version tracking. The registry is integrated with the Nextflow runtime. Nextflow will automatically locate and download configured plugins.

:::{note}
The Nextflow Plugin Registry is currently available as private beta technology. Contact [info@nextflow.io](mailto:info@nextflow.io) to learn how to get access.
:::

## Impact on users and developers

The impact of the Nextflow Gradle plugin differs for plugin users and developers.

### Plugin Users

If you are a plugin user, no immediate actions are required. The plugin configuration has not changed.

### Plugin developers

Developers are encouraged to migrate to the Nextflow Gradle plugin and benefit from features that simplify plugin development and integration with the wider plugin ecosystem.

To migrate an existing Nextflow plugin:

1. Remove the following files and folders:
    - `buildSrc/`
    - `nextflow.config`
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
                url = 'https://nf-plugins-registry.dev-tower.net/api'
                authToken = project.findProperty('pluginRegistry.accessToken')
            }
        }
    }
    ```

    Replace the following:

    - `DEPENDENCY`: (Optional) Your plugins dependency libraries—for example, `commons-io:commons-io:2.18.0`.
    - `PLUGIN_VERSION:` Your plugin version—for example, `0.5.0`.
    - `MINIMUM_NEXTFLOW_VERSION`: The minimum Nextflow version required to run your plugin—for example, `25.03.0-edge`.
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

(migrating-plugin-page)=

# Migrating to the Gradle plugin for Nextflow plugins

This page introduces the Gradle plugin for Nextflow plugins, the Nextflow plugin registry, and how to migrate to the new plugin framework.


## Improvements to the plugin framework

The Nextflow plugin ecosystem is evolving to support a more robust and user-friendly experience by simplifying plugin development, streamlining publishing and discovery, and improving how plugins are loaded into workflows. These improvements make plugins more accessible, maintainable, and interoperable with Nextflow.

### Gradle plugin for Nextflow plugins 

The Gradle plugin for Nextflow plugins simplifies and standardizes the development of Nextflow plugins. It configures default dependencies required for Nextflow integration and introduces custom Gradle tasks to streamline building, testing, packaging, and publishing plugins.

The Gradle plugin is versioned and published to the [Gradle Plugin Portal](https://plugins.gradle.org/), allowing developers to manage it like any other dependency. As the plugin ecosystem evolves, this Gradle plugin will enable easier maintenance and adoption of ongoing improvements to the Nextflow plugin framework.

### Nextflow plugin registry

The Nextflow plugin registry is a centralized repository of assembled plugins. It hosts an index of plugin metadata that supports plugin discovery, accessibility, and version tracking. The registry is integrated with the Nextflow runtime. Nextflow will automatically locate and download configured plugins.

## Impact on users and developers

The impact of the Gradle plugin for Nextflow plugins differs for plugin users and developers.

### Plugin Users

If you are a plugin user, no immediate actions are required. The plugin configuration has not changed.

### Plugin developers

Developers are encouraged to migrate to the Gradle plugin for Nextflow plugins and benefit from features that simplify plugin development and integration with the wider plugin ecosystem.

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
        id 'io.nextflow.nextflow-plugin' version '0.0.1-alpha3'
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
            github {
                repository = '<GITHUB_REPOSITORY>'
                userName = project.findProperty('github_username')
                authToken = project.findProperty('github_access_token')
                email = project.findProperty('github_commit_email')
                indexUrl = '<GITHUB_INDEX_URL>'
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
    - `GITHUB_REPOSITORY`: Your GitHub plugin repository name—for example, `nextflow-io/nf-plugin-template`.
    - `GITHUB_INDEX_URL`: The URL of your fork of the plugins index repository—for example, [`https://github.com/nextflow-io/plugins/blob/main/plugins.json`](https://github.com/nextflow-io/plugins/blob/main/plugins.json).

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

The Gradle plugin for Nextflow plugins also supports publishing plugins. See {ref}`gradle-plugin-package` for more information.

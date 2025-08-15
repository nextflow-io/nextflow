(plugin-registry-page)=

# Migrating to the Nextflow plugin registry

The Nextflow plugin ecosystem is evolving to support a more robust and user-friendly experience by simplifying the development, publishing, and discovery of Nextflow plugins. This page introduces the Nextflow plugin registry, the Nextflow Gradle plugin, and how to migrate to them.

:::{note}
The Nextflow plugin registry and Gradle plugin are currently available as a public preview. Plugin developers can access the registry by contacting [info@nextflow.io](mailto:info@nextflow.io) to obtain registry access tokens.
:::

## Overview

### Nextflow plugin registry

The Nextflow plugin registry is a central repository for Nextflow plugins. It hosts an index of plugin metadata that supports plugin discovery, accessibility, and version tracking. Nextflow 25.04 and later can use the plugin registry as a drop-in replacement for the [legacy plugin index](https://github.com/nextflow-io/plugins) hosted on GitHub.

### Nextflow Gradle plugin

The [Nextflow Gradle plugin](https://github.com/nextflow-io/nextflow-plugin-gradle) simplifies the development of Nextflow plugins. It provides default configuration required for Nextflow integration, as well as custom Gradle tasks for building, testing, and publishing plugins.

The Gradle plugin is versioned and published to the [Gradle Plugin Portal](https://plugins.gradle.org/), allowing developers to manage it like any other dependency. As the plugin ecosystem evolves, the Gradle plugin will enable easier maintenance and adoption of ongoing improvements to the Nextflow plugin framework.

## Timeline

The [legacy plugin index](https://github.com/nextflow-io/plugins) will be deprecated in favor of the Nextflow plugin registry.

:::{note}
The following timeline is tentative and subject to modification.
:::

<h3>Nextflow 25.04</h3>

The Nextflow plugin registry is available as a private beta. Nextflow 25.04 can use the Nextflow plugin registry as an opt-in feature. The Nextflow plugin registry will be automatically kept up-to-date with the [legacy plugin index](https://github.com/nextflow-io/plugins).

During this time, plugin developers are encouraged to experiment with the Gradle plugin and plugin registry.

<h3>Nextflow 25.10</h3>

The Nextflow plugin registry will be generally available. Nextflow 25.10 will use the plugin registry by default. The legacy plugin index will be **closed to new pull requests**.

Developers will be required to publish to the Nextflow plugin registry. To ensure continued support for older versions of Nextflow, the legacy plugin index will be automatically kept up-to-date with the Nextflow plugin registry.

<h3>Nextflow 26.04</h3>

Nextflow 26.04 will only be able to use the Nextflow plugin registry.

At some point in the future, the legacy plugin index will be **frozen** -- it will no longer receives updates from the Nextflow plugin registry. To ensure continued support for older versions of Nextflow, the legacy plugin index will remain available indefinitely.

## Impact on plugin users

No immediate actions are required for plugin users. The plugin configuration has not changed.

## Impact on plugin developers

Plugin developers will need to update their plugin to publish to the Nextflow plugin registry instead of the legacy plugin index. The best way to do this is to migrate to the Nextflow Gradle plugin, which simplifies the development process and supports publishing to the plugin registry from the command line.

### Migrating to the Nextflow Gradle plugin

To migrate an existing Nextflow plugin:

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

    Replace `PLUGIN_NAME` with your plugin name.

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

        publishing {
            registry {
                authToken = project.findProperty('npr.apiKey')
            }
        }
    }
    ```

    Replace the following:

    - `DEPENDENCY`: (Optional) Your plugins dependency libraries. For example, `commons-io:commons-io:2.18.0`.
    - `PLUGIN_VERSION:` Your plugin version. For example, `0.5.0`.
    - `MINIMUM_NEXTFLOW_VERSION`: The minimum Nextflow version required to run your plugin. For example, `25.04.0`.
    - `PROVIDER`: Your name or organization. For example, `acme`.
    - `CLASS_NAME`: Your plugin class name. For example, `acme.plugin.MyPlugin`.
    - `EXTENSION_POINT`: Your extension point identifiers that the plugin will implement or expose. For example, `acme.plugin.MyFactory`.

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

### Publishing to the Nextflow plugin registry

The Nextflow Gradle plugin supports publishing plugins from the command line. See {ref}`gradle-plugin-publish` for more information.

Once you migrate to the Gradle plugin, you will no longer be able to publish to the legacy plugin index. See the [transition timeline](#timeline) for more information.

## Plugin Creation and Publishing Workflow

### Components

* **Nextflow** - https://github.com/nextflow-io/nextflow/
* **Plugin template** - https://github.com/nextflow-io/nf-plugin-template  
* **Gradle plugin for Nextflow** - https://github.com/nextflow-io/nextflow-plugin-gradle
* **Plugin registry** - https://registry.nextflow.io/

### Workflow Steps

#### 1. Plugin Name Registration
- User claims a unique plugin name via the Plugin registry (e.g., `nf-foo`)
- Registry admin approves the plugin name claim (one-time process)

#### 2. Plugin Creation  
- User creates a new plugin using: `nextflow plugin create <plugin-name> <provider-name>`
- Nextflow generates plugin scaffold using the Plugin template
- Interactive mode available: `nextflow plugin create` (prompts for inputs)

#### 3. Plugin Development
- User completes plugin implementation based on the generated scaffold
- Plugin structure follows standardized template with:
  - Gradle build configuration using Nextflow Gradle plugin
  - Extension points for custom functionality
  - Unit and end-to-end tests
  - Documentation templates

#### 4. Local Testing
- Install plugin locally: `make install` 
- Test plugin: `make test`
- Build plugin: `make assemble`

#### 5. Plugin Publishing
- Configure registry API key in `$HOME/.gradle/gradle.properties`:
  ```
  npr.apiKey=<YOUR_API_KEY>
  ```
- Release plugin: `make release`
- The Gradle plugin handles upload to the Plugin registry
- Plugin becomes accessible for public distribution
- For additional configuration options, see the [Nextflow Gradle plugin documentation](https://github.com/nextflow-io/nextflow-plugin-gradle)

### Key Features
- **Registry-based distribution**: Centralized plugin discovery and versioning
- **Standardized structure**: Consistent development experience across plugins  
- **Automated publishing**: Command-line publishing with proper validation
- **Backward compatibility**: Gradual migration from legacy GitHub-based system

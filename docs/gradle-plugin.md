(gradle-plugin-page)=

# Using the Gradle plugin

The [Gradle plugin for Nextflow plugins](https://github.com/nextflow-io/nextflow-plugin-gradle) simplifies plugin development by configuring default dependencies needed for Nextflow integration and incorporates custom Gradle tasks that streamline building, testing, and publishing Nextflow plugins. This guide describes how to use the Gradle plugin for plugin development.

:::{note}
Nextflow Plugins can be developed without the Gradle plugin. However, this approach is only suggested if you are an advanced developer and your project is incompatible with the Gradle plugin.
:::

(gradle-plugin-create)=

## Creating a plugin

The [nf-hello](https://github.com/nextflow-io/nf-hello/tree/gradle-plugin-example) plugin uses the Gradle plugin and is a valuable starting point for developers.

To create a Nextflow plugin with the Gradle plugin:

1. Fork the [nf-hello](https://github.com/nextflow-io/nf-hello/tree/gradle-plugin-example) plugin. See {ref}`nf-hello-page` for more information.
2. Rename the forked `nf-hello` directory with your plugin name.
3. Replace the contents of `settings.gradle` with the following:

    ```groovy
    rootProject.name = '<PLUGIN_NAME>'
    ```

    Replace `PLUGIN_NAME` with your plugin name.

4. Replace the contents of `build.gradle` with the following:

    ```groovy
    // Plugins
    plugins {
        id 'io.nextflow.nextflow-plugin' version '0.0.1-alpha'
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

    - `DEPENDENCY`: (optional) your plugins dependency libraries—for example, `commons-io:commons-io:2.18.0`.
    - `PLUGIN_VERSION:` your plugin version—for example, `0.5.0`.
    - `MINIMUM_NEXTFLOW_VERSION`: the minimum Nextflow version required to run your plugin—for example, `24.11.0-edge`.
    - `PROVIDER`: your name or organization—for example, `nextflow`.
    - `CLASS_NAME`: your plugin class name—for example, `nextflow.hello.HelloPlugin`.
    - `EXTENSION_POINT`: your extension point identifiers that the plugin will implement or expose—for example, `nextflow.hello.HelloFactory`.
    - `GITHUB_REPOSITORY`: your GitHub plugin repository name—for example, `nextflow-io/nf-hello`.
    - `GITHUB_INDEX_URL`: the URL of your fork of the plugins index repository—for example, [`plugins.json`](https://github.com/username/plugins/blob/main/plugins.json)</code>.
5. Develop your plugin extension points. See {ref}`dev-plugins-extension` for descriptions and examples.
6. In the plugin root directory, run `make assemble`.

(gradle-plugin-install)=

## Installing a plugin

Plugins can be installed locally without being packaged, uploaded, and published.

To install a plugin locally:

1. In the plugin root directory, run `make install`.

    :::{note}
    Running `make install` will add your plugin to your `$HOME/.nextflow/plugins` directory.
    :::

2. Configure your plugin. See {ref}`using-plugins-page` for more information.
3. Run your pipeline:

    ```bash
    nextflow run main.nf
    ```

(gradle-plugin-unit-test)=

## Unit testing a plugin

Unit tests are small, focused tests designed to verify the behavior of individual plugin components and are an important part of software development.

To run unit tests:

1. Develop your unit tests. See [HelloDslTest.groovy](https://github.com/nextflow-io/nf-hello/blob/gradle-plugin-example/src/test/groovy/nextflow/hello/HelloDslTest.groovy) in the [nf-hello](https://github.com/nextflow-io/nf-hello/tree/gradle-plugin-example) plugin for unit test examples.
2. In the plugin root directory, run `make test`.

(gradle-plugin-package)=

## Packaging, uploading, and publishing a plugin

The Gradle plugin for Nextflow plugins simplifies publishing your plugin.

To package, upload, and publish your plugin:

1. Fork the [Nextfow plugins index repository](https://github.com/nextflow-io/plugins).
2. In the plugin root directory, open `build.gradle` and ensure that:
    * `github.repository` matches the plugin repository.
    * `github.indexUrl` matches your fork of the plugins index repository.
3. Create a file named `$HOME/.gradle/gradle.properties` and add the following:

    ```bash
    github_username=<GITHUB_USERNAME>
    github_access_token=<GITHUB_ACCESS_TOKEN>
    github_commit_email=<GITHUB_EMAIL>
    ```

    Replace the following:
    * `GITHUB_USERNAME`: your GitHub username granting access to the plugin repository.
    * `GITHUB_ACCESS_TOKEN`: your GitHub access token with permission to upload and commit changes to the plugin repository.
    * `GITHUB_EMAIL`: your email address associated with your GitHub account.
4. Run `make release`.
5. Create a pull request against the [Nextfow plugins index repository](https://github.com/nextflow-io/plugins) from your fork.
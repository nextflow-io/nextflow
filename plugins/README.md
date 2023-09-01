# Nextflow plugins

This directory contains the subprojects for Nextflow core plugins.

## Plugin structure

Each plugin subproject has its own `build.gradle` which defines the main build actions and the required dependencies.

Minimal dependencies shown below:

```groovy
dependencies {
    compileOnly project(':nextflow')
    compileOnly 'org.slf4j:slf4j-api:1.7.10'
    compileOnly 'org.pf4j:pf4j:3.4.1'

    testImplementation project(':nextflow')
    testImplementation "org.codehaus.groovy:groovy:3.0.5"
    testImplementation "org.codehaus.groovy:groovy-nio:3.0.5"
}
```

Each plugin subproject directory name must begin with the prefix `nf-` and must include a file named `src/resources/META-INF/MANIFEST.MF` which contains the plugin metadata. The manifest content looks like the following:

```
Manifest-Version: 1.0
Plugin-Class: the.plugin.ClassName
Plugin-Id: the-plugin-id
Plugin-Provider: Your Name or Organization
Plugin-Version: 0.0.0
```

## Environment variables

* `NXF_PLUGINS_MODE`: the plugin execution mode, either `prod` for production or `dev` for development (see below for details).

* `NXF_PLUGINS_DIR`: the path where the plugin archives are loaded and stored (default: `$NXF_HOME/plugins` in production, `./plugins` in development).

* `NXF_PLUGINS_DEFAULT`: whether to use the default plugins when no plugins are specified in the Nextflow configuration (default: `true`).

* `NXF_PLUGINS_DEV`: comma-separated list of development plugin root directories.

## Development environment

When running in development the plugin system uses the `DevPluginClasspath` to load plugins classes from each plugin project build path e.g. `plugins/nf-amazon/build/classes` and `plugins/nf-amazon/build/target/libs` (for deps libraries).

## The plugins repository

The plugins meta-info are published to the [nextflow-io/plugins](https://github.com/nextflow-io/plugins) GitHub repository, specifically [this file](https://raw.githubusercontent.com/nextflow-io/plugins/main/plugins.json).

The repository index has the following structure:

```json
[
  {
    "id": "nf-amazon",
    "releases": [
      {
        "version": "0.2.0",
        "url": "https://github.com/nextflow-io/nf-amazon/releases/download/0.2.0/nf-amazon-0.2.0.zip",
        "date": "2020-10-12T10:05:44.28+02:00",
        "sha512sum": "9e9e33695c1a7c051271..."
      }
    ]
  },
  /* ... */
]
```

## Plugins

A plugin is a ZIP file containing the plugin classes and the required dependencies.

Nextflow core plugins are stored in the corresponding GitHub project release page. However, this is not a
strict requirement, it has been chosen to simplify the build deployment process and provide a more consistent
download experience keeping all of them with the GitHub [nextflow-io](https://github.com/nextflow-io) organization.

### Installation

Plugins must be declared in the `nextflow.config` file using the plugins scope, for example:

```groovy
plugins {
    id 'nf-amazon@0.2.0'
}
```

If a plugin is not locally available, Nextflow checks the repository index for the download URL, downloads and extracts the plugin archive, and installs the plugin into the directory specified by `NXF_PLUGINS_DIR` (default: `$NXF_HOME/plugins`).

Since each Nextflow run can have a different set of plugins (and versions), each Nextflow run keeps a local plugins directory called `.nextflow/plr/<session-id>` which links the exact set of plugins required for the given Nextflow run.

If no plugins are specified in the Nextflow configuration, the default plugins are automatically added. The set of default plugins is defined in the Nextflow resources file included in the distribution runtime `modules/nextflow/src/main/resources/META-INF/plugins-info.txt`.

To disable the use of default plugins, set the environment variable `NXF_PLUGINS_DEFAULT=false`.

## Gradle Tasks

### `makeZip`

Create the plugin archive and JSON meta file in the subproject `build/libs` directory.

```console
$ ls -l1 plugins/nf-tower/build/libs/
nf-tower-0.1.0.jar
nf-tower-0.1.0.json
nf-tower-0.1.0.zip
```

### `copyPluginLibs`

Copy plugin dependencies JAR files into the subproject `build/target/libs` directory.
This is needed only when launching the plugin in development mode.

### `copyPluginZip`

Copy the plugin archive to the root project build directory, i.e. `build/plugins`.

### `uploadPlugin`

Upload the plugin archive and meta file to the corresponding GitHub repository. Options:

* `release`: the plugin version e.g. `1.0.1`
* `repo`: the GitHub repository name e.g. `nf-amazon`
* `owner`: the GitHub owning organization e.g. `nextflow-io`
* `skipExisting`: do not upload a file already existing (only if the checksum is the same, default: `true`).
* `dryRun`: execute the tasks without uploading file (default: `false`).
* `overwrite`: prevent to overwrite a remote file already existing (default: `false`).
* `userName`: the user name for authenticate GitHub API requests
* `authToken`: the personal token to authenticate GitHub API requests

### `upload`

Upload the plugin archive and meta file.

### `publishIndex`

Upload the plugins index to the repository hosted at [nextflow-io/plugins](https://github.com/nextflow-io/plugins), which makes them publicly accessible at [this URL](https://raw.githubusercontent.com/nextflow-io/plugins/main/plugins.json).

## Additional Resources

* https://pf4j.org/
* https://proandroiddev.com/understanding-gradle-the-build-lifecycle-5118c1da613f

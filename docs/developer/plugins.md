
# Plugins

This page describes how to create, test, and publish third-party plugins.

## Plugin structure

The best way to get started with your own plugin is to refer to the [nf-hello](https://github.com/nextflow-io/nf-hello) repository. This repository provides a minimal plugin implementation with several examples of different extension points and instructions for building, testing, and publishing.

Plugins can be written in Java or Groovy.

The minimal dependencies are as follows:

```groovy
dependencies {
    compileOnly project(':nextflow')
    compileOnly 'org.slf4j:slf4j-api:1.7.10'
    compileOnly 'org.pf4j:pf4j:3.4.1'

    testImplementation project(':nextflow')
    testImplementation "org.codehaus.groovy:groovy:4.0.24"
    testImplementation "org.codehaus.groovy:groovy-nio:4.0.23"
}
```

The plugin subproject directory name must begin with the prefix `nf-` and must include a file named `src/resources/META-INF/MANIFEST.MF` which contains the plugin metadata. The manifest content looks like the following:

```
Manifest-Version: 1.0
Plugin-Class: the.plugin.ClassName
Plugin-Id: the-plugin-id
Plugin-Provider: Your Name or Organization
Plugin-Version: 0.0.0
```

## Extension points

Nextflow's plugin system exposes a variety of extension points for plugins. This section describes how to use these extension points when writing a plugin, as well as how they are used in a pipeline.

:::{note}
If you would like to implement something in a plugin that isn't covered by any of the following sections, feel free to create an issue on GitHub and describe your use case. In general, any class in the Nextflow codebase that implements `ExtensionPoint` can be extended by a plugin, and existing plugins are a great source of examples when writing new plugins.
:::

:::{note}
Plugin extension points must be added to `extensions.idx` in the plugin repository to make them discoverable. See the [`nf-hello`](https://github.com/nextflow-io/nf-hello) plugin for an example.
:::

### Commands

Plugins can define custom CLI commands that can be executed with the `nextflow plugin` command.

To implement a plugin-specific command, implement the `PluginExecAware` interface in your plugin entrypoint (the class that extends `BasePlugin`). Alternatively, you can implement the `PluginAbstractExec` trait, which provides an abstract implementation with some boilerplate code. This trait requires you to implement two methods, `getCommands()` and `exec()`:

```groovy
import nextflow.cli.PluginAbstractExec
import nextflow.plugin.BasePlugin

class MyPlugin extends BasePlugin implements PluginAbstractExec {
    @Override
    List<String> getCommands() {
        [ 'hello' ]
    }

    @Override
    int exec(String cmd, List<String> args) {
        if( cmd == 'hello' ) {
            println "Hello! You gave me these arguments: ${args.join(' ')}"
            return 0
        }
        else {
            System.err.println "Invalid command: ${cmd}"
            return 1
        }
    }
}
```

You can then execute this command using the `nextflow plugin` command:

```bash
nextflow plugin my-plugin:hello --foo --bar
```

See the {ref}`cli-plugin` CLI command for usage information.

### Configuration

Plugins can access the resolved Nextflow configuration through the session object using `session.config.navigate()`. Several extension points provide the session object for this reason. This method allows you to query any configuration option in a safe manner -- if the option isn't defined, it will return `null`. A common practice is to define any configuration for your plugin in a custom config scope.

For example, you can query a config option in a trace observer hook:

```groovy
import nextflow.Session
import nextflow.trace.TraceObserver

class MyObserver implements TraceObserver {

    @Override
    void onFlowCreate(Session session) {
        final message = session.config.navigate('myplugin.create.message')
        println message
    }
}
```

You can then set this option in your config file:

```groovy
// dot syntax
myplugin.create.message = "I'm alive!"

// closure syntax
myplugin {
    create {
        message = "I'm alive!"
    }
}
```

### Executors

Plugins can define custom executors that can then be used with the `executor` process directive.

To implement an executor, create a class in your plugin that extends the [`Executor`](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/executor/Executor.groovy) class and implements the `ExtensionPoint` interface. Add the `@ServiceName` annotation to your class with the name of your executor:

```groovy
import nextflow.executor.Executor
import nextflow.util.ServiceName
import org.pf4j.ExtensionPoint

@ServiceName('my-executor')
class MyExecutor extends Executor implements ExtensionPoint {

    // ...

}
```

You can then use this executor in your pipeline:

```nextflow
process foo {
    executor 'my-executor'
}
```

:::{tip}
Refer to the source code of Nextflow's built-in executors to see how to implement the various components of an executor. You might be able to implement most of your executor by simply reusing existing code.
:::

### Functions

:::{versionadded} 22.09.0-edge
:::

Plugins can define custom functions, which can then be included in Nextflow pipelines.

To implement a custom function, create a class in your plugin that extends the `PluginExtensionPoint` class, and implement your function with the `Function` annotation:

```groovy
import nextflow.Session
import nextflow.plugin.extension.Function
import nextflow.plugin.extension.PluginExtensionPoint

class MyExtension extends PluginExtensionPoint {

    @Override
    void init(Session session) {}

    @Function
    String reverseString(String origin) {
        origin.reverse()
    }

}
```

You can then use this function in your pipeline:

```nextflow
include { reverseString } from 'plugin/my-plugin'

channel.of( reverseString('hi') )
```

You can also use an alias:

```nextflow
include { reverseString as anotherReverseMethod } from 'plugin/my-plugin'
```

### Operators

:::{versionadded} 22.04.0
:::

Plugins can define custom channel factories and operators, which can then be included into Nextflow pipelines.

To implement a custom factory or operator, create a class in your plugin that extends the `PluginExtensionPoint` class, and implement your function with the `Factory` or `Operator` annotation:

```groovy
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Session
import nextflow.plugin.extension.Factory
import nextflow.plugin.extension.Operator
import nextflow.plugin.extension.PluginExtensionPoint

class MyExtension extends PluginExtensionPoint {

    @Override
    void init(Session session) {}

    @Factory
    DataflowWriteChannel fromQuery(Map opts, String query) {
        // ...
    }

    @Operator
    DataflowWriteChannel sqlInsert(DataflowReadChannel source, Map opts) {
        // ...
    }

}
```

You can then use them in your pipeline:

```nextflow
include { sqlInsert; fromQuery as fromTable } from 'plugin/nf-sqldb'

def sql = 'select * from FOO'
channel
    .fromTable(sql, db: 'test', emitColumns: true)
    .sqlInsert(into: 'BAR', columns: 'id', db: 'test')
```

The above snippet is based on the [nf-sqldb](https://github.com/nextflow-io/nf-sqldb) plugin. The `fromQuery` factory 
is included under the alias `fromTable`.

### Process directives

Plugins that implement a [custom executor](#executors) will likely need to access {ref}`process directives <process-directives>` that affect the task execution. When an executor receives a task, the process directives can be accessed through that task's configuration. As a best practice, custom executors should try to support all process directives that have executor-specific behavior and are relevant to your executor.

Nextflow does not provide the ability to define custom process directives in a plugin. Instead, you can use the {ref}`process-ext` directive to provide custom process settings to your executor. Try to use specific names that are not likely to conflict with other plugins or existing pipelines.

For example, a custom executor can use existing process directives as well as a custom setting through the `ext` directive:

```groovy
class MyExecutor extends Executor {

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        final cpus = task.config.cpus
        final memory = task.config.memory
        final myOption = task.config.ext.myOption

        println "This task is configured with cpus=${cpus}, memory=${memory}, myOption=${myOption}"

        // ...
    }

    // ...

}
```

### Trace observers

A *trace observer* in Nextflow is an entity that can listen and react to workflow events, such as when a workflow starts, a task completes, a file is published, etc. Several components in Nextflow, such as the execution report and DAG visualization, are implemented as trace observers.

Plugins can define custom trace observers that react to workflow events with custom behavior. To implement a trace observer, create a class that implements the `TraceObserver` trait and another class that implements the `TraceObserverFactory` interface. Implement any of the hooks defined in `TraceObserver`, and implement the `create()` method in your observer factory:

```groovy
// MyObserverFactory.groovy
import nextflow.Session
import nextflow.trace.TraceObserver
import nextflow.trace.TraceObserverFactory

class MyObserverFactory implements TraceObserverFactory {

    @Override
    Collection<TraceObserver> create(Session session) {
        final enabled = session.config.navigate('myplugin.enabled')
        return enabled ? [ new MyObserver() ] : []
    }
}

// MyObserver.groovy
import java.nio.file.Path

import nextflow.processor.TaskHandler
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

class MyObserver implements TraceObserver {
    
    @Override
    void onFlowBegin() {
        println "Okay, let's begin!"
    }

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        println "I completed a task! It's name is '${handler.task.name}'"
    }

    @Override
    void onProcessCached(TaskHandler handler, TraceRecord trace) {
        println "I found a task in the cache! It's name is '${handler.task.name}'"
    }

    @Override
    void onFilePublish(Path destination, Path source) {
        println "I published a file! It's located at ${path.toUriString()}"
    }
    
    @Override
    void onFlowError(TaskHandler handler, TraceRecord trace) {
        println "Uh oh, something went wrong..."
    }

    @Override
    void onFlowComplete() {
        println 'All done!'
    }
}
```

You can then use your trace observer by simply enabling the plugin in your pipeline. In the above example, the observer must also be enabled with a config option:

```groovy
myplugin.enabled = true
```

Refer to the `TraceObserver` [source code](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/trace/TraceObserver.groovy) for descriptions of the available workflow events.

## Publishing plugins

Nextflow resolves plugins through a plugin registry, which stores metadata for each plugin version, including the publishing date, checksum, and download URL for the plugin binary. The default registry is located on GitHub at [nextflow-io/plugins](https://github.com/nextflow-io/plugins/), specifically [this file](https://raw.githubusercontent.com/nextflow-io/plugins/main/plugins.json).

To publish a plugin, you must create the plugin release, publish it to GitHub, and make a pull request against the main plugin registry with the metadata for your plugin release.

A plugin release is a ZIP archive containing the compiled plugin classes, the required dependencies, and a JSON file with the plugin metadata. See the [`nf-hello`](https://github.com/nextflow-io/nf-hello) example to see how to create a plugin release.

Here is an example meta file for a plugin release:

```json
{
    "version": "0.2.0",
    "url": "https://github.com/nextflow-io/nf-amazon/releases/download/0.2.0/nf-amazon-0.2.0.zip",
    "date": "2020-10-12T10:05:44.28+02:00",
    "sha512sum": "9e9e33695c1a7c051271..."
}
```

(testing-plugins)=

## Testing plugins

Plugins must be declared in the `nextflow.config` file using the `plugins` scope, for example:

```groovy
plugins {
    id 'nf-amazon@0.2.0'
}
```

If a plugin is not locally available, Nextflow checks the repository index for the download URL, downloads and extracts the plugin archive, and installs the plugin into the directory specified by `NXF_PLUGINS_DIR` (default: `${NXF_HOME}/plugins`).

Since each Nextflow run can have a different set of plugins (and versions), each run keeps a local plugins directory called `.nextflow/plr/<session-id>` which links the exact set of plugins required for the given run.

Additionally, the "default plugins" (defined in the Nextflow resources file `modules/nextflow/src/main/resources/META-INF/plugins-info.txt`) are always made available for use. To disable the use of default plugins, set the environment variable `NXF_PLUGINS_DEFAULT=false`.

When running in development mode, the plugin system uses the `DevPluginClasspath` to load plugin classes from each plugin project build path, e.g. `plugins/nf-amazon/build/classes` and `plugins/nf-amazon/build/target/libs` (for deps libraries).

## Environment variables

The following environment variables are available when developing and testing plugins:

`NXF_PLUGINS_MODE`
: The plugin execution mode, either `prod` for production or `dev` for development (see below for details).

`NXF_PLUGINS_DIR`
: The path where the plugin archives are loaded and stored (default: `$NXF_HOME/plugins` in production, `./plugins` in development).

`NXF_PLUGINS_DEFAULT`
: Whether to use the default plugins when no plugins are specified in the Nextflow configuration (default: `true`).

`NXF_PLUGINS_DEV`
: Comma-separated list of development plugin root directories.

`NXF_PLUGINS_TEST_REPOSITORY`
: :::{versionadded} 23.04.0
  :::
: Comma-separated list of URIs for additional plugin registries or meta files, which will be used in addition to the default registry.
: The URI should refer to a plugins repository JSON file or a specific plugin JSON meta file. In the latter case it should match the pattern `https://host.name/some/path/<plugin id>-X.Y.Z-meta.json`.

: For example:

    ```bash
    # custom plugin repository at https://github.com/my-org/plugins
    export NXF_PLUGINS_TEST_REPOSITORY="https://raw.githubusercontent.com/my-org/plugins/main/plugins.json"

    # custom plugin release
    export NXF_PLUGINS_TEST_REPOSITORY="https://github.com/nextflow-io/nf-hello/releases/download/0.3.0/nf-hello-0.3.0-meta.json"

    nextflow run <pipeline> -plugins nf-hello
    ```

: This variable is useful for testing a plugin release before publishing it to the main registry.

## Gradle Tasks

The following build tasks are defined in the `build.gradle` of each core plugin as well as the `nf-hello` example plugin.

### `makeZip`

Create the plugin archive and JSON meta file in the subproject `build/libs` directory.

```console
$ ls -l1 plugins/nf-tower/build/libs/
nf-tower-0.1.0.jar
nf-tower-0.1.0.json
nf-tower-0.1.0.zip
```

### `copyPluginLibs`

Copy plugin dependencies JAR files into the subproject `build/target/libs` directory. This is needed only when launching the plugin in development mode.

### `copyPluginZip`

Copy the plugin archive to the root project build directory, i.e. `build/plugins`.

### `uploadPlugin`

Upload the plugin archive and meta file to the corresponding GitHub repository. Options:

`release`
: The plugin version, e.g. `1.0.1`.

`repo`
: The GitHub repository name, e.g. `nf-amazon`.

`owner`
: The GitHub owning organization, e.g. `nextflow-io`.

`skipExisting`
: Do not upload a file if it already exists, i.e. checksum is the same (default: `true`).

`dryRun`
: Execute the tasks without uploading file (default: `false`).

`overwrite`
: Prevent to overwrite a remote file already existing (default: `false`).

`userName`
: The user name used to authenticate GitHub API requests.

`authToken`
: The personal token used to authenticate GitHub API requests.

### `upload`

Upload the plugin archive and meta file.

### `publishIndex`

Upload the plugins index to the repository hosted at [nextflow-io/plugins](https://github.com/nextflow-io/plugins), which makes them publicly accessible at [this URL](https://raw.githubusercontent.com/nextflow-io/plugins/main/plugins.json).

## Additional Resources

* [PF4J](https://pf4j.org/)
* [Understanding Gradle: The Build Lifecycle](https://proandroiddev.com/understanding-gradle-the-build-lifecycle-5118c1da613f)

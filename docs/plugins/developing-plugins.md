(dev-plugins-page)=

# Developing plugins

This page describes how to develop plugins for Nextflow.

(dev-plugins-template)=

## Nextflow plugin template

The [Nextflow plugin template](https://github.com/nextflow-io/nf-plugin-template/) is a scaffold for plugin development. It uses [Gradle](https://gradle.org/), a build automation tool optimized for Java and Groovy projects, as well as the [Nextflow Gradle plugin](https://github.com/nextflow-io/nextflow-plugin-gradle).

You can use the `nextflow plugin create` sub-command to create plugins using the plugin template. See {ref}`gradle-plugin-create` for more information.

:::{note}
The Nextflow Gradle plugin is currently available as a public preview. See {ref}`Migrating to the Nextflow plugin registry <plugin-registry-page>` for more information.
:::

(dev-plugins-structure)=

### Structure

The plugin template includes the source directories, build configuration files, and metadata required for development, testing, and publishing. Depending on the developer’s preference, plugins can be written in Java or Groovy.

For example, a plugin created from the plugin template with the name `nf-hello` and organization `nextflow` will have the following structure:

```console
nf-hello
├── .github/workflows/build.yml
├── COPYING
├── Makefile
├── README.md
├── build.gradle
├── gradle
│   └── wrapper
│       ├── gradle-wrapper.jar
│       └── gradle-wrapper.properties
├── gradlew
├── settings.gradle
├── src
│   ├── main
│   │   └── groovy
│   │       └── nextflow
│   │           └── hello
│   │               ├── HelloExtension.groovy
│   │               ├── HelloFactory.groovy
│   │               ├── HelloObserver.groovy
│   │               └── HelloPlugin.groovy
│   └── test
│       └── groovy
│           └── nextflow
│               └── hello
│                   └── HelloObserverTest.groovy
└── validation
    └── main.nf
    └── nextflow.config
```

This structure contains the following key files and folders:

- `.github/workflows/build.yml`: GitHub Action which implements continuous integration for the plugin.

- `build.gradle`: The Gradle build script.

- `COPYING`: The project license, detailing the terms under which the code can be used and distributed.​

- `gradle/wrapper/`: Helper files for the Gradle Wrapper.

- `gradlew`: The Gradle Wrapper script, which allows you to use Gradle without installing it into your environment.

- `Makefile`: Defines common tasks for building, testing, and publishing the plugin with Make.​

- `README.md`: The project README, which provides an overview of the project, including its purpose, features, and instructions for usage and development.​

- `settings.gradle`: The Gradle project configuration, which specifies project-specific settings such as the project name and included modules.

- `src/main/groovy/<ORGANIZATION>/`: The main source directory, which contains the plugin source code and resources.​

- `src/test/groovy/<ORGANIZATION>/`: The test source directory, which contains the plugin unit tests.

- `validation`: A small Nextflow pipeline which serves as an end-to-end test for the plugin.

The plugin template also implements the following example features:

- A custom trace observer that prints a message when the workflow starts and when the workflow completes (see `HelloObserver`).

- A custom function called `sayHello` (see `HelloExtension`).

### Nextflow Gradle plugin

The [Nextflow Gradle plugin](https://github.com/nextflow-io/nextflow-plugin-gradle) simplifies the development of Nextflow plugins. It provides default configuration required for Nextflow integration, as well as custom Gradle tasks for building, testing, and publishing plugins.

It is versioned and published to the [Gradle Plugin Portal](https://plugins.gradle.org/), and can be declared and managed like any other dependency in the `build.gradle` file:

```nextflow
plugins {
    id 'io.nextflow.nextflow-plugin' version '1.0.0-beta.6'
}
```

:::{note}
You can develop Nextflow plugins without the Gradle plugin. However, this approach is only suggested if you are an advanced developer and your project is incompatible with the Gradle plugin.
:::

### Make commands

The plugin template includes a Makefile which wraps the most important Gradle tasks provided by the Nextflow Gradle plugin.

These tasks can be executed with [Make](https://www.gnu.org/software/make/). For example:

```bash
make assemble
```

The following `make` commands are available:

`assemble`
: Compiles the Nextflow plugin code and assembles it into a zip archive.

`install`
: Installs the plugin into the local Nextflow plugins directory.

`release`
: Publishes the plugin. See {ref}`gradle-plugin-publish` for more information.

`test`
: Runs plugin unit tests. See {ref}`gradle-plugin-test` for more information.

(dev-plugins-extension-points)=

## Extension points

Nextflow’s plugin system exposes various extension points. This section gives examples of typical extension points and how to use them.

### Commands

Plugins can define custom CLI commands that are executable with the `nextflow plugin` command.

To implement a plugin-specific command, implement the `PluginExecAware` interface in your plugin entry point (the class that extends `BasePlugin`). Alternatively, implement the `PluginAbstractExec` trait, which provides an abstract implementation with some boilerplate code. This trait requires you to implement the `getCommands()` and `exec()` methods. For example:

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

The command can be run using the `nextflow plugin` command:

```bash
nextflow plugin my-plugin:hello --alpha --beta
```

See the {ref}`cli-plugin` for usage information.

### Configuration

Plugins can access the resolved Nextflow configuration through the session object using `session.config.navigate()`. Several extension points provide the session object for this reason. This method allows you to query any configuration option safely. If the option isn’t defined, it will return null.

A common practice is to use a custom config scope to define any configuration for your plugin. For example:

```groovy
import nextflow.Session
import nextflow.trace.TraceObserver

class MyObserver implements TraceObserver {

    @Override
    void onFlowCreate(Session session) {
        final message = session.config.navigate('myplugin.createMessage')
        println message
    }
}
```

This option can then be set in your configuration file:

```groovy
// dot syntax
myplugin.createMessage = "I'm alive!"

// block syntax
myplugin {
    createMessage = "I'm alive!"
}
```

:::{versionadded} 25.04.0
:::

Plugins can declare their configuration options by implementing the `ConfigScope` interface and declaring each config option as a field with the `@ConfigOption` annotation. For example:

```groovy
import nextflow.config.schema.ConfigOption
import nextflow.config.schema.ConfigScope
import nextflow.config.schema.ScopeName
import nextflow.script.dsl.Description

@ScopeName('myplugin')
@Description('''
    The `myplugin` scope allows you to configure the `nf-myplugin` plugin.
''')
class MyPluginConfig implements ConfigScope {

    @ConfigOption
    @Description('''
        Message to print to standard output when a run is initialized.
    ''')
    final String createMessage

    // no-arg constructor is required to enable validation of config options
    MyPluginConfig() {
    }

    MyPluginConfig(Map opts) {
        this.createMessage = opts.createMessage
    }
}
```

This approach is not required to support plugin config options. However, it allows Nextflow to recognize plugin definitions when validating the configuration. See {ref}`config-scopes-page` for more information.

### Executors

Plugins can define custom executors that can be used with the `executor` process directive.

To implement an executor, create a class in your plugin that extends the [`Executor`](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/executor/Executor.groovy) class and implements the `ExtensionPoint` interface. Add the `@ServiceName` annotation to your class with the name of your executor. For example:

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

```groovy
process hello {
    executor 'my-executor'

    // ...
}
```

:::{tip}
See the source code of Nextflow's built-in executors for examples of how to implement various executor components.
:::

### Filesystems

Plugins can define custom filesystems that Nextflow can use to interact with external storage systems using a single interface. For more information about accessing remote files, see {ref}`remote-files`.

To implement a custom filesystem, create a class in your plugin that extends [`FileSystemProvider`](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/nio/file/spi/FileSystemProvider.html). Implement the `getScheme()` method to define the URI scheme for your filesystem. For example:

```groovy
import java.nio.file.spi.FileSystemProvider

class MyFileSystemProvider extends FileSystemProvider {

    @Override
    String getScheme() {
        return 'myfs'
    }

    // ...
}
```

You can then use this filesystem in your pipeline:

```nextflow
input = file('myfs://<PATH_TO_INPUT_FILE>')
```

See [Developing a Custom File System Provider](https://docs.oracle.com/javase/8/docs/technotes/guides/io/fsp/filesystemprovider.html) for more information and the `nf-amazon` plugin (`S3FileSystemProvider`) for an example of a custom filesystem.

:::{tip}
Custom filesystems are an advanced plugin extension. Before creating a new filesystem, check that your use case cannot already be supported by an existing filesystem such as HTTP or S3.
:::

### Functions

:::{versionadded} 22.10.0
:::

Plugins can define custom functions that can be included in Nextflow pipelines.

To implement a custom function, create a plugin class that extends the `PluginExtensionPoint` class and implement your function with the `Function` annotation. For example:

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

You can then add this function to your pipeline:

```nextflow
include { reverseString } from 'plugin/my-plugin'

channel.of( reverseString('hi') )
```

Alternatively, you can use an alias:

```nextflow
include { reverseString as anotherReverseMethod } from 'plugin/my-plugin'
```

### Operators

:::{versionadded} 22.04.0
:::

Plugins can define custom channel factories and operators that can then be included in pipelines.

To implement a custom channel factory or operator, create a class in your plugin that extends the `PluginExtensionPoint` class and implement your function with the `Factory` or `Operator` annotation. For example:

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

You can then use the custom channel factories or operators in your pipeline:

```nextflow
include { sqlInsert; fromQuery as fromTable } from 'plugin/nf-sqldb'

def sql = 'select * from FOO'
channel
    .fromTable(sql, db: 'test', emitColumns: true)
    .sqlInsert(into: 'BAR', columns: 'id', db: 'test')
```

:::{note}
The above snippet is based on the [nf-sqldb](https://github.com/nextflow-io/nf-sqldb) plugin. The `fromQuery` factory is included under the alias `fromTable`.
:::

:::{tip}
Before creating a custom operator, consider whether the operator can be defined as a [function](#functions) that can be composed with existing operators such as `map` or `subscribe`. Functions are easier to implement and can be used anywhere in your pipeline, not just channel logic.
:::

### Process directives

Plugins that implement a custom executor will likely need to access {ref}`process directives <process-directives>` that affect the task execution. When an executor receives a task, the process directives can be accessed through that task’s configuration. Custom executors should try to support all process directives that have executor-specific behavior and are relevant to the executor.

Nextflow does not provide the ability to define custom process directives in a plugin. Instead, use the {ref}`process-ext` directive to provide custom process settings to your executor. Use specific names that are not likely to conflict with other plugins or existing pipelines.

For example, a custom executor can use existing process directives and a custom setting through the `ext` directive:

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

(plugins-trace-observers)=

### Trace observers

:::{versionchanged} 25.04.0
The `TraceObserver` interface is now deprecated. Use [TraceObserverV2](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/trace/TraceObserverV2.groovy) and [TraceObserverFactoryV2](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/trace/TraceObserverFactoryV2.groovy) instead.
:::

A *trace observer* is an entity that can listen and react to workflow events, such as when a workflow starts, a task is completed, or a file is published. Several components in Nextflow, such as the execution report and DAG visualization, are implemented as trace observers.

Plugins can define custom trace observers that react to workflow events with custom behavior. To implement a trace observer, create a class that implements the `TraceObserver` trait and another class that implements the `TraceObserverFactory` interface. Implement any of the hooks defined in `TraceObserver` and implement the `create()` method in your observer factory. For example:

```groovy
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

```nextflow
myplugin.enabled = true
```

See the [`TraceObserver` source code](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/trace/TraceObserver.groovy) for descriptions of the available workflow events.

(dev-plugins-env-vars)=

## Environment variables

The following environment variables are available to develop and test plugins:

`NXF_PLUGINS_MODE`
: The plugin execution mode. Either `prod` for production or `dev` for development.

`NXF_PLUGINS_DIR`
: The path where the plugin archives are loaded and stored (default: `$NXF_HOME/plugins` in production and `./plugins` in development).

`NXF_PLUGINS_DEFAULT`
: Whether to use the default plugins when no plugins are specified in the Nextflow configuration (default: true).

`NXF_PLUGINS_DEV`
: Comma-separated list of development plugin root directories.

`NXF_PLUGINS_TEST_REPOSITORY`
: :::{versionadded} 23.04.0.
  :::
: Comma-separated list of URIs for additional plugin registries or meta files, which will be used in addition to the default registry.

: The URI should refer to a plugin repository JSON file or a specific plugin JSON meta file. In the latter case, it should match the pattern `https://host.name/some/path/<PLUGIN_NAME>-X.Y.Z-meta.json`. For example:

    ```bash
    # custom plugin repository at https://github.com/my-org/plugins
    export NXF_PLUGINS_TEST_REPOSITORY="https://raw.githubusercontent.com/my-org/plugins/main/plugins.json"

    # custom plugin release
    export NXF_PLUGINS_TEST_REPOSITORY="https://github.com/nextflow-io/nf-hello/releases/download/0.3.0/nf-hello-0.3.0-meta.json"

    nextflow run main.nf -plugins nf-hello
    ```

: This variable is useful for testing a plugin release before publishing it to the main registry.

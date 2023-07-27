(plugins-page)=

# Plugins

Nextflow has a plugin system that allows the use of extensible components that are downloaded and installed at runtime.

## Core plugins

The following functionalities are provided via plugin components, and they make part of the Nextflow *core* plugins:

- `nf-amazon`: Support for Amazon cloud.
- `nf-azure`: Support for Azure cloud.
- `nf-console`: Implement Nextflow [REPL console](https://www.nextflow.io/blog/2015/introducing-nextflow-console.html).
- `nf-ga4gh`: Support [GA4GH APIs](https://www.ga4gh.org/).
- `nf-google`: Support for Google cloud.
- `nf-tower`: Support for [Tower](https://tower.nf) cloud platform.
- `nf-wave`: Support for [Wave containers](https://seqera.io/wave/) service.


## Using plugins

The core plugins do not require any configuration. They are automatically installed when the corresponding feature is requested by a Nextflow pipeline. You can still specify them as described below, e.g. if you want to pin the version of a plugin, however if you try to use a plugin version that isn't compatible with your Nextflow version, Nextflow will fail.

You can enable a plugin by declaring it in your Nextflow configuration:

```groovy
plugins {
    id 'nf-hello@0.1.0'
}
```

Or you can use the `-plugins` command line option:

```bash
nextflow run <pipeline> -plugins nf-hello@0.1.0
```

The plugin identifier consists of the plugin name and plugin version separated by a `@`. Multiple plugins can be specified 
in the configuration with multiple `id` declarations, or on the command line as a comma-separated list. When specifying 
plugins via the command line, any plugin declarations in the configuration file are ignored.

The default plugins are documented in this documentation. For all other plugins, please refer to the plugin's code repository 
for documentation and support.

## Writing plugins

To get started with your own plugin, refer to the [nf-hello](https://github.com/nextflow-io/nf-hello) repository, 
which provides a minimal plugin implementation with several examples of different extension points, as well as instructions 
for building, testing, and publishing.

Nextflow's plugin system exposes a variety of extension points for plugins. The following sections describe how to use 
these extension points when writing a plugin, as well as how they are used in a pipeline.

:::{note}
If you would like to implement something in a plugin that isn't covered by any of the following sections, feel free to 
create an issue on GitHub and describe your use case. In general, any class in the Nextflow codebase that implements 
`ExtensionPoint` can be extended by a plugin, and existing plugins are a great source of examples when writing new plugins.
:::

:::{note}
Plugin extension points must be added to `extensions.idx` in the plugin repository to make them discoverable. 
See the `nf-hello` plugin for an example.
:::

### Commands

Plugins can define custom CLI commands that can be executed with the `nextflow plugin` command.

To implement a plugin-specific command, implement the `PluginExecAware` interface in your plugin entrypoint 
(the class that extends `BasePlugin`). Alternatively, you can implement the `PluginAbstractExec` trait, which 
provides an abstract implementation with some boilerplate code. This trait requires you to implement two methods, 
`getCommands()` and `exec()`:

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

Plugins can access the resolved Nextflow configuration through the session object using `session.config.navigate()`. 
Several extension points provide the session object for this reason. This method allows you to query any configuration 
option in a safe manner -- if the option isn't defined, it will return `null`. A common practice is to define any 
configuration for your plugin in a custom config scope.

Here is an example of querying a config option in a trace observer hook:

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

```groovy
process foo {
    executor 'my-executor'
}
```

:::{tip}
Refer to the source code of Nextflow's built-in executors to see how to implement the various components of an executor. 
You might be able to implement most of your executor by simply reusing existing code.
:::

### Functions

:::{versionadded} 22.09.0-edge
:::

Plugins can define custom Groovy functions, which can then be included into Nextflow pipelines.

To implement a custom function, create a class in your plugin that extends the `PluginExtensionPoint` class, and implement 
your function with the `Function` annotation:

```groovy
import nextflow.Session
import nextflow.plugin.extension.Function
import nextflow.plugin.extension.PluginExtensionPoint

class MyExtension implements PluginExtensionPoint {

    @Override
    void init(Session session) {}

    @Function
    String reverseString(String origin) {
        origin.reverse()
    }

}
```

You can then use this function in your pipeline:

```groovy
include { reverseString } from 'plugin/my-plugin'

channel.of( reverseString('hi') )
```

You can also use an alias:

```groovy
include { reverseString as anotherReverseMethod } from 'plugin/my-plugin'
```

### Operators

:::{versionadded} 22.04.0
:::

Plugins can define custom channel factories and operators, which can then be included into Nextflow pipelines.

To implement a custom factory or operator, create a class in your plugin that extends the `PluginExtensionPoint` class, 
and implement your function with the `Factory` or `Operator` annotation:

```groovy
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Session
import nextflow.plugin.extension.Factory
import nextflow.plugin.extension.Operator
import nextflow.plugin.extension.PluginExtensionPoint

class MyExtension implements PluginExtensionPoint {

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

```groovy
include { sqlInsert; fromQuery as fromTable } from 'plugin/nf-sqldb'

def sql = 'select * from FOO'
channel
    .fromTable(sql, db: 'test', emitColumns: true)
    .sqlInsert(into: 'BAR', columns: 'id', db: 'test')
```

The above snippet is based on the [nf-sqldb](https://github.com/nextflow-io/nf-sqldb) plugin. The `fromQuery` factory 
is included under the alias `fromTable`.

### Trace observers

A *trace observer* in Nextflow is an entity that can listen and react to workflow events, such as when a workflow starts, 
a task completes, a file is published, etc. Several components in Nextflow, such as the execution report and DAG visualization, 
are implemented as trace observers.

Plugins can define custom trace observers that react to workflow events with custom behavior. To implement a trace observer, 
create a class that implements the `TraceObserver` trait and another class that implements the `TraceObserverFactory` interface. 
Implement any of the hooks defined in `TraceObserver`, and implement the `create()` method in your observer factory:

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

You can then use your trace observer by simply enabling the plugin in your pipeline. In the above example, the observer 
must also be enabled with a config option:

```groovy
myplugin.enabled = true
```

Refer to the `TraceObserver` [source code](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/trace/TraceObserver.groovy) for descriptions of the available workflow events.

## Plugin registry

Nextflow resolves plugins through a plugin registry, which stores metadata for each plugin version, including the publishing date, 
checksum, and download URL for the plugin binary. The default registry is located on GitHub at [nextflow-io/plugins](https://github.com/nextflow-io/plugins/).

To publish a plugin release to the main registry, simply create a pull request with the requested plugin metadata.

(testing-plugins)=

### Testing plugins

:::{versionadded} 23.04.0
:::

You can also use a different plugin registry with the `NXF_PLUGINS_TEST_REPOSITORY` environment variable. This setting 
is useful for testing a plugin release before publishing it to the main registry. It can refer to the JSON file for a 
custom registry or a plugin release.

For example:

```bash
# custom registry at https://github.com/my-org/plugins
export NXF_PLUGINS_TEST_REPOSITORY="https://raw.githubusercontent.com/my-org/plugins/main/plugins.json"

# custom plugin release
export NXF_PLUGINS_TEST_REPOSITORY="https://github.com/nextflow-io/nf-hello/releases/download/0.3.0/nf-hello-0.3.0-meta.json"

nextflow run <pipeline> -plugins nf-hello
```

## Offline usage

To use Nextflow plugins in an offline environment:

1. Download the {ref}`"all" release <getstarted-install>` of Nextflow, which comes with the following default plugins: `nf-amazon`, `nf-google`, `nf-tower`.

2. Download any additional plugins by running `nextflow plugin install <pluginId,..>`. Alternatively, simply run your pipeline once and Nextflow will download all of the plugins that it needs.

3. Copy the `nextflow` binary and `$HOME/.nextflow` folder to your offline environment.

4. In your Nextflow configuration file, specify each plugin that you downloaded, both name and version, including default plugins. This will prevent Nextflow from trying to download newer versions of plugins.

Nextflow caches the plugins that it downloads, so as long as you keep using the same Nextflow version and pin your plugin versions in your config file, Nextflow will use the locally installed plugins and won't try to download them from the Internet.

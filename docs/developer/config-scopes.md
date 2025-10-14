(config-scopes-page)=

# Configuration scopes

This page provides guidance on defining configuration scopes in the Nextflow runtime.

## Overview

The Nextflow configuration is defined as a collection of *scope classes*. Each scope class defines the set of available options, including their name, type, and an optional description for a specific configuration scope.

Scope classes are used to generate a configuration spec, which is in turn used for several purposes:

- Validating config options at runtime (`nextflow run` and `nextflow config`)

- Providing code intelligence in the language server (validation, hover hints, code completion)

- Generating reference documentation (in progress)

Scope classes are also used by the runtime itself as type-safe domain objects. This way, the construciton of domain objects from the configuration map is isolated from the rest of the runtime.

## Definition

### Config scopes

A *config scope* is defined as a class that implements the `ConfigScope` interface. Top-level scope classes must have the `@ScopeName` annotation, which defines the name of the config scope.

For example:

```groovy
package nextflow.hello

import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName

@ScopeName('hello')
class HelloConfig implements ConfigScope {
}
```

A scope class must provide a default constructor, so that it can be instantiated as an extension point. If no such constructor is defined, the config scope will not be detected by Nextflow. In the above example, this constructor is implicitly defined because no constructors were declared.

The fully-qualified class name (in this case, `nextflow.hello.HelloConfig`) must be included in the list of extension points.

### Config options

A *config option* is defined as a field with the `@ConfigOption` annotation. The field name determines the name of the config option.

For example:

```groovy
    @ConfigOption
    String createMessage
```

The `@ConfigOption` annotation can specify an optional set of types that are valid in addition to the field type. For example, the `fusion.tags` option, which accepts either a String or Boolean, is declared as follows:

```groovy
    @ConfigOption(types=[Boolean])
    String tags
```

The field type and any additional types are included in the config spec, allowing them to be used for validation.

The field type can be any Java or Groovy class, but in practice it should be a class that can be constructed from primitive values (numbers, booleans, strings). For example, `Duration` and `MemoryUnit` are standard Nextflow types that can each be constructed from an integer or string.

### Nested scopes

A *nested scope* is defined as a field whose type is an implementation of `ConfigScope`. The field name determines the name of the nested scope.

The scope class referenced by the field type defines config options and scopes in the same manner as top-level scope classes. Unlike top-level scopes, nested scope classes do not need to use the `@ScopeName` annotation or provide a default constructor.

See `ExecutorConfig` and `ExecutorRetryConfig` for an example of how a nested scope is defined and constructed.

### Placeholder scopes

A *placeholder scope* is a config scope that applies to a collection of user-defined names.

For example, the `azure.batch.pools` scope allows the user to define a set of named pools, where each pool is configured with a standard set of options such as `autoScale`, `lowPriority`, `maxVmCount`, etc. These options are defined in a placeholder scope with a placeholder name of `<name>`. Thus, the generic name for the `autoScale` option is `azure.batch.pools.<name>.autoScale`.

A placeholder scope is defined as a field with type `Map<String, P>`, where `P` is a nested scope class which defines the scope options. The field should have the `@PlaceholderName` annotation which defines the placeholder name (e.g. `<name>`).

See `AzBatchOpts` and `AzPoolOpts` for an example of how placeholder scopes are defined and constructed.

### Descriptions

Top-level scope classes and config options should use the `@Description` annotation to provide a description of the scope or option. This description is included in the config spec, which is in turn used by the language server to provide hover hints.

For example:

```groovy
@ScopeName('hello')
@Description('''
    The `hello` scope controls the behavior of the `nf-hello` plugin.
''')
class HelloConfig implements ConfigScope {

    @ConfigOption
    @Description('''
        Message to print to standard output when a run is initialized.
    ''')
    String createMessage
}
```

Nested scopes and placeholder scopes may also use this annotation, but will inherit the description of the top-level scope by default.

### Best practices

The Nextflow runtime adheres the following best practices where appropriate:

- Config options should be declared as public and final, so that the scope class can be used as an immutable domain object.

- Scope classes should define a constructor that initializes each field from a map, casting each map property to the required type and providing default values as needed.

- In cases where an option defaults to an environment variable, the environment map should be provided as an additional constructor argument rather than accessing the system environment directly.

- In cases where an option with a primitive type (e.g., `int`, `float`, `boolean`) can be unspecified without a default value, it should be declared with the equivalent reference type (e.g. `Integer`, `Float`, `Boolean`), otherwise it should use the primitive type.

- In cases where an option represents a path, it should be declared as a `String` and allow clients to construct paths as needed, since path construction may depend on plugins which aren't yet loaded.

For example:

```groovy
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName

@ScopeName('hello')
class HelloConfig implements ConfigScope {

    @ConfigOption
    final String createMessage

    @ConfigOption
    final boolean verbose

    HelloConfig() {}

    HelloConfig(Map opts, Map env) {
        this.createMessage = opts.createMessage ?: env.get('NXF_HELLO_CREATE_MESSAGE')
        this.verbose = opts.verbose as boolean
    }
}
```

## Usage

### Runtime

Nextflow validates the config map after it is loaded. Top-level config scopes are loaded by the plugin system as extension points and converted into a config spec, which is used to validate the config map.

Plugins are loaded after the config is loaded and before it is validated, since plugins can also define config scopes. If a third-party plugin declares a config scope, it must be explicitly enabled in order to validate config options from the plugin. Otherwise, Nextflow will report these options as unrecognized.

Core plugins are loaded automatically based on other config options. Therefore, Nextflow only validates config from a core plugin when that plugin is loaded. Otherwise, any config options from the plugin are ignored -- they are neither validated nor reported as unrecognized.

For example, when the `process.executor` config option is set to `'awsbatch'`, the `nf-amazon` is automatically loaded. In this case, all options in the `aws` config scope will be validated. If the executor is not set to `'awsbatch'`, all `aws` options will be ignored. This way, config files can be validated appropriately without loading additional core plugins that won't be used by the run.

The scope classes themselves can be used to construct domain objects on-demand from the config map. For example, an `ExecutorConfig` can be constructed from the `executor` config scope as follows:

```groovy
new ExecutorConfig( Global.session.config.executor as Map ?: Collections.emptyMap() )
```

:::{note}
In practice, it is better to avoid the use of `Global` and provide an instance of `Session` to the client class instead.
:::

### JSON Schema

Config scope classes can be converted into a config spec with the `SpecNode` class, which uses reflection to extract metadata such as scope names, option names, types, and descriptions. This spec is rendered to JSON and used by the language server at build-time to provide code intelligence such as code completion and hover hints.

### Documentation

The config spec described above can be rendered to Markdown using the `MarkdownRenderer` class. It produces a Markdown document approximating the {ref}`config-options` page.

This approach to docs generation is not yet complete, and has not been incorporated into the build process yet. However, it can be used to check for discrepancies between the source code and docs when making changes. The documentation should match the `@Description` annotations as closely as possible, but may contain additional details such as version notes and extra paragraphs.

(updating-syntax-page)=

# Updating Nextflow syntax

This page describes how to update Nextflow scripts and config files to adhere to the {ref}`Nextflow language specification <syntax-page>`, also known as the "strict syntax".

:::{note}
If you are still using DSL1, see {ref}`dsl1-page` to learn how to migrate your Nextflow pipelines to DSL2 first, before consulting this guide.
:::

## Preparing for strict syntax

The strict syntax is a subset of DSL2. While DSL2 in practice allows any Groovy syntax, the strict syntax allows only a subset of Groovy syntax that is appropriate for Nextflow scripts and config files. This new specification enables more specific error reporting and more consistent code, and it will allow the Nextflow language to evolve independently of Groovy.

The strict syntax is currently only enforced by the Nextflow language server, which is provided as part of the {ref}`vscode-page` for Nextflow. However, the strict syntax will be gradually adopted by the Nextflow CLI in future versions, and will eventually be the only way to write Nextflow code.

New language features will be generally implemented as part of the strict syntax, and not the current "lenient" DSL2 parser, with few exceptions. Therefore, it will be important to prepare for the strict syntax in order to use new language features in the future.

This section highlights the most common errors encountered when updating to the strict syntax, and shows how to resolve them. In general, the amount of required changes depends on the amount of custom Groovy code that you have.

### Removed syntax

<h3>Import declarations</h3>

In Groovy, the `import` declaration can be used to import external classes:

```groovy
import groovy.json.JsonSlurper

def json = new JsonSlurper().parseText(json_file.text)
```

In the strict syntax, use the fully qualified name instead:

```nextflow
def json = new groovy.json.JsonSlurper().parseText(json_file.text)
```

<h3>Class declarations</h3>

Some users use custom classes in Nextflow to define helper functions or custom record types. Instead, define helper functions as standalone functions in a script and move custom record classes to the `lib` directory. Enums, a special type of class, are supported. However, they cannot be included across modules at this time.

:::{note}
Record types will be addressed in a future version of the Nextflow language specification.
:::

<h3>Mixing script declarations and statements</h3>

A script may contain any of the following top-level declarations:

- Feature flags
- Include declarations
- Parameter declarations
- Workflows
- Processes
- Functions
- Output block

Alternatively, a script may contain only statements, also known as a _code snippet_:

```nextflow
println 'Hello world!'
```

Code snippets are treated as an implicit entry workflow:

```nextflow
workflow {
    println 'Hello world!'
}
```

Script declarations and statements cannot be mixed at the same level. All statements must reside within script declarations unless the script is a code snippet:

```nextflow
process foo {
    // ...
}

// incorrect -- move into entry workflow
// println 'Hello world!'

// correct
workflow {
    println 'Hello world!'
}
```

:::{note}
Mixing statements and script declarations was necessary in DSL1 and allowed in DSL2. However, it is no longer supported in the strict syntax in order to simplify the language and to ensure that top-level statements are only executed when the script is executed directly and not when it is included as a module.
:::

<h3>Assignment expressions</h3>

In Groovy, variables can be assigned as part of an expression:

```groovy
foo(x = 1, y = 2)
```

In the strict syntax, assign variables as statements instead:

```nextflow
x = 1
y = 2
foo(x, y)
```

In Groovy, the increment (`++`) and decrement (`--`) arithmetic operators are supported.

```nextflow
x++
x--
```

In the strict syntax, use `+=` and `-=` instead:

```nextflow
x += 1
x -= 1
```

<h3>For and while loops</h3>

In Groovy, loop statements, such as `for` and `while`, are supported:

```groovy
for (rseqc_module in ['read_distribution', 'inner_distance', 'tin']) {
    if (rseqc_modules.contains(rseqc_module))
        rseqc_modules.remove(rseqc_module)
}
```

In the strict syntax, use higher-order functions, such as the `each` method, instead:

```nextflow
['read_distribution', 'inner_distance', 'tin'].each { rseqc_module ->
    if (rseqc_modules.contains(rseqc_module))
        rseqc_modules.remove(rseqc_module)
}
```

Lists, maps, and sets provide several functions (e.g., `collect`, `find`, `findAll`, `inject`) for iteration. See [Groovy standard library](https://docs.groovy-lang.org/latest/html/groovy-jdk/overview-summary.html) for more information.

<h3>Switch statements</h3>

In Groovy, switch statements can be used for pattern matching on a value:

```groovy
switch (aligner) {
case 'bowtie2':
    // ...
    break
case 'bwamem':
    // ...
    break
case 'dragmap':
    // ...
    break
case 'snap':
    // ...
    break
default:
    // ...
}
```

In the strict syntax, use if-else statements instead:

```nextflow
if (aligner == 'bowtie2') {
    // ...
} else if (aligner == 'bwamem') {
    // ...
} else if (aligner == 'dragmap') {
    // ...
} else if (aligner == 'snap') {
    // ...
} else {
    // ...
}
```

<h3>Spread operator</h3>

In Groovy, the _spread_ operator can be used to flatten a nested list:

```groovy
ch.map { meta, bambai -> [meta, *bambai] }
```

In the strict syntax, enumerate the list elements explicitly instead:

```groovy
// alternative 1
ch.map { meta, bambai -> [meta, bambai[0], bambai[1]] }

// alternative 2
ch.map { meta, bambai ->
    def (bam, bai) = bambai
    [meta, bam, bai]
}
```

<h3>Implicit environment variables</h3>

In Nextflow DSL1 and DSL2, you can reference environment variables directly in strings:

```nextflow
println "PWD = ${PWD}"
```

In the strict syntax, use `System.getenv()` instead:

```nextflow
println "PWD = ${System.getenv('PWD')}"
```

:::{versionadded} 24.11.0-edge
The `env()` function can be used instead of `System.getenv()`. For example,
`println "PWD = ${env('PWD')}"`.
:::

### Restricted syntax

The following patterns are still supported but have been restricted. That is, some syntax variants have been removed.

<h3>Variable declarations</h3>

In Groovy, variables can be declared in many different ways:

```groovy
def a = 1
final b = 2
def c = 3, d = 4
def (e, f) = [5, 6]
String str = 'foo'
def Map meta = [:]
```

In the strict syntax, declare variables with `def` and do not specify a type:

```nextflow
def a = 1
def b = 2
def (c, d) = [3, 4]
def (e, f) = [5, 6]
def str = 'foo'
def meta = [:]
```

Similarly, for the strict syntax, declare functions with `def` and do not specify a return type or parameter type:

```nextflow
/**
 * You can use comments to denote types, for example:
 *
 * @param x: Map
 * @param y: String
 * @param z: Integer
 * @return List
 */
def foo(x, y, z) {
    // ...
}
```

:::{note}
Type annotations are useful for providing type checking at runtime. The language server will not report errors or warnings for Groovy-style type annotations at this time.

Instead, type annotations will be addressed in a future version of the Nextflow language specification, at which point the language server will provide a way to automatically migrate Groovy-style type annotations to the new syntax.
:::

<h3>Strings</h3>

Groovy supports a wide variety of strings, including multi-line strings, dynamic strings, slashy strings, multi-line dynamic slashy strings, and more.

The strict syntax supports single- and double-quoted strings, multi-line strings, and slashy strings.

Slashy strings cannot be interpolated:

```nextflow
def id = 'SRA001'
assert 'SRA001.fastq' ~= /${id}\.f(?:ast)?q/
```

Use a double-quoted string instead:

```nextflow
def id = 'SRA001'
assert 'SRA001.fastq' ~= "${id}\\.f(?:ast)?q"
```

Slashy strings cannot span multiple lines:

```groovy
/
Patterns in the code,
Symbols dance to match and find,
Logic unconfined.
/
```

Use a multi-line string instead:

```nextflow
"""
Patterns in the code,
Symbols dance to match and find,
Logic unconfined.
"""
```

Dollar slashy strings are not supported:

```groovy
$/
echo "Hello world!"
/$
```

Use a multi-line string instead:

```nextflow
"""
echo "Hello world!"
"""
```

<h3>Type conversions</h3>

In Groovy, soft and hard casts are supported:

```groovy
def map = (Map) readJson(json)  // soft cast
def map = readJson(json) as Map // hard cast
```

In the strict syntax, only hard casts are supported. However, hard casts are discouraged because they can cause unexpected behavior if used improperly. Use a Groovy-style type annotation instead:

```groovy
def Map map = readJson(json)
```

Nextflow will raise an error at runtime if the `readJson()` function does not return a `Map`.

In cases where you want to explicitly convert a value to a different type, use an explicit method. For example, to parse a string as a number:

```groovy
def x = '42' as Integer
def x = '42'.toInteger()    // preferred
```

<h3>Process env inputs and outputs</h3>

In Nextflow DSL1 and DSL2, the name of a process `env` input/output can be specified with or without quotes:

```nextflow
process PROC {
    input:
    env FOO
    env 'BAR'

    // ...
}
```

In the strict syntax, specify the name with quotes:

```nextflow
process PROC {
    input:
    env 'FOO'
    env 'BAR'

    // ...
}
```

<h3>Implicit process script section</h3>

In Nextflow DSL1 and DSL2, the process `script:` section label can almost always be omitted:

```nextflow
process greet {
    input:
    val greeting

    """
    echo '${greeting}!'
    """
}
```

In the strict syntax, the `script:` label can be omitted only if there are no other sections:

```nextflow
process sayHello {
    """
    echo 'Hello world!'
    """
}

process greet {
    input:
    val greeting

    script:
    """
    echo '${greeting}!'
    """
}
```

<h3>Workflow onComplete/onError handlers</h3>

{ref}`Workflow handlers <workflow-handlers>` (i.e. `workflow.onComplete` and `workflow.onError`) can be defined in several different ways in a script, but are typically defined as top-level statements and without an equals sign:

```nextflow
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
```

Because the strict syntax does not allow statements to be mixed with script declarations, workflow handlers must be defined in the entry workflow:

```nextflow
workflow {
    // ...

    workflow.onComplete = {
        println "Pipeline completed at: $workflow.complete"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    }
}
```

:::{note}
A more concise syntax for workflow handlers will be addressed in a future version of the Nextflow language specification.
:::

### Deprecated syntax

The following patterns are deprecated. The language server reports _paranoid warnings_ for these patterns. Paranoid warnings are disabled by default. Enable them by selecting **Nextflow > Paranoid Warnings** in the {ref}`extension settings <vscode-settings>`. These warnings may become errors in the future.

<h3>Implicit closure parameter</h3>

In Groovy, a closure with no parameters is assumed to have a single parameter named `it`:

```nextflow
ch | map { it * 2 }
```

In the strict syntax, implicit closure parameters are not supported. Declare the parameter explicitly instead:

```nextflow
ch | map { v -> v * 2 }   // correct
ch | map { it -> it * 2 } // also correct
```

<h3>Using params outside the entry workflow</h3>

While params can be used anywhere in the pipeline code, they are only intended to be used in the entry workflow.

In the strict syntax, processes and workflows should receive params as explicit inputs:

```nextflow
process foo {
    input:
    val foo_args

    // ...
}

workflow bar {
    take:
    bar_args

    // ...
}

workflow {
    foo(params.foo_args)
    bar(params.bar_args)
}
```

<h3>Process each input</h3>

The `each` process input is deprecated. Use the `combine` or `cross` operator to explicitly repeat over inputs in the calling workflow.

<h3>Process when section</h3>

The process `when` section is deprecated. Use conditional logic, such as an `if` statement or the `filter` operator, to control the process invocation in the calling workflow.

<h3>Process shell section</h3>

The process `shell` section is deprecated. Use the `script` block instead. The VS Code extension provides syntax highlighting and error checking to help distinguish between Nextflow variables and Bash variables.

### Configuration syntax

See {ref}`config-syntax` for a comprehensive description of the configuration language.

Currently, Nextflow parses config files as Groovy scripts, allowing the use of scripting constructs like variables, helper functions, try-catch blocks, and conditional logic for dynamic configuration:

```groovy
def getHostname() {
    // ...
}

def hostname = getHostname()
if (hostname == 'small') {
    params.max_memory = 32.GB
    params.max_cpus = 8
}
else if (hostname == 'large') {
    params.max_memory = 128.GB
    params.max_cpus = 32
}
```

The strict config syntax does not support functions, and only allows statements (e.g., variables and if statements) within closures. You can achieve the same dynamic configuration by using a dynamic include:

```groovy
includeConfig ({
    def hostname = // ...
    if (hostname == 'small')
        return 'small.config'
    else if (hostname == 'large')
        return 'large.config'
    else
        return '/dev/null'
}())
```

The include source is a closure that is immediately invoked. It includes a different config file based on the return value of the closure. Including `/dev/null` is equivalent to including nothing.

Each conditional configuration is defined in a separate config file:

```groovy
// small.config
params.max_memory = 32.GB
params.max_cpus = 8

// large.config
params.max_memory = 128.GB
params.max_cpus = 32
```

## Preserving Groovy code

There are two ways to preserve Groovy code, for example, if the code is difficult to update: using the `lib` directory or writing a plugin.

Any Groovy code can be moved into the `lib` directory, which supports the full Groovy language. This approach is useful for temporarily preserving some Groovy code until it can be updated later and incorporated into a Nextflow script. See the {ref}`lib directory <lib-directory>` documentation for more information.

For Groovy code that is more complicated or depends on third-party libraries, it may be better to create a plugin. Plugins can define custom functions that can be included by Nextflow scripts like a module. Furthermore, plugins can be easily re-used across different pipelines. See {ref}`plugins-dev-page` for more information on how to develop plugins.

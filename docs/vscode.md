
# VS Code integration

The [Nextflow VS Code extension](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow) provides IDE support for Nextflow pipelines.

## Features

### Syntax highlighting

The VS Code extension highlights Nextflow scripts and config files for better readability.

### Diagnostics

The extension highlights source code in red for errors and yellow for warnings.

To view all diagnostics for the workspace, open the **Problems** tab. Here, you can search for diagnostics by diagnostic message, filename, and so on.

### Hover hints

When you hover over certain source code elements, such as variable names and function calls, the extension provides a tooltip with related information, such as the definition and/or documentation for the element.

### Code navigation

The **Outline** section in the Explorer panel lists top-level definitions when you view a script. Include declarations in scripts and config files act as links, and ctrl-clicking them opens the corresponding script or config file.

To view the definition of a symbol (e.g., a workflow, process, function, or variable), right-click the symbol and select **Go to Definition**. Ctrl-click the symbol to view the definition. Ctrl-click the definition to show all references.

### Code completion

The extension suggests auto-completions for variable names, function names, config settings, and other symbols as you type. The extension also provides several snippets for common script declarations, such as processes and workflows.

### Formatting

The extension can format your scripts and config files based on a standard set of formatting rules. Rules can be customized using the **Nextflow > Formatting** [extension settings](#settings).

Use the **Format Document** command in the command palette to format the current file.

### Renaming symbols

The extension can rename all references of a symbol (e.g., a workflow, process, function, or variable) throughout the code.

To rename a symbol, right-click the symbol, select **Rename Symbol**, and enter a new name.

### Parameter schema

If a `nextflow_schema.json` file exists in the same directory as a script with an entry workflow, the extension uses the schema to provide validation, hover hints, and code completion for params in the entry workflow.

### DAG preview for workflows

The extension can generate a workflow DAG that includes the workflow inputs, outputs, and any processes or workflows that are called by the selected workflow. The workflow DAG is displayed in a new panel to the side.

To preview the DAG of a workflow, select the **Preview DAG** CodeLens above the workflow definition.

:::{note}
The **Preview DAG** CodeLens is only available when the script does not contain any errors.
:::

## Syntax guide

The language server parses scripts and config files according to the {ref}`Nextflow language specification <syntax-page>`. It enforces a stricter syntax compared to the Nextflow CLI. As a result, the language server is able to perform more extensive error checking and provide more specific error messages. However, unlike the Nextflow CLI which allows all Groovy syntax, the Nextflow language specification is not a superset of Groovy. You may need to adjust your code to adhere to the strict syntax, especially if you use more advanced Groovy syntax.

This section describes some of the most common unsupported features and how to address them. For a full description of the strict syntax, refer to the Nextflow language specification.

:::{note}
The "Nextflow language specification" is a strict specification of DSL2, not a new DSL version. The Nextflow language will be defined by this specification moving forward, rather than new DSL versions.
:::

:::{note}
You can move unsupported code into the `lib` directory or a plugin, both of which support the full Groovy language.
:::

### Excluded syntax

**Import declarations**

In Groovy, the `import` declaration can be used to import external classes:

```groovy
import groovy.json.JsonSlurper

def json = new JsonSlurper().parseText(json_file.text)
```

Instead, use the fully qualified name directly:

```nextflow
def json = new groovy.json.JsonSlurper().parseText(json_file.text)
```

**Class declarations**

Some users use custom classes in Nextflow to define helper functions or custom record types. Instead, helper functions can be defined as standalone functions in a script. Custom record classes must be moved to the `lib` directory.

:::{note}
Record types will be addressed in a future version of the Nextflow language specification.
:::

:::{note}
Enums, a special type of class, are supported, but cannot be included across modules at this time.
:::

**Mixing script declarations and statements**

A script may contain any of the following top-level declarations:

- Feature flags
- Includes
- Parameter declarations
- Workflows
- Processes
- Functions
- Output block

Alternatively, a script may contain only statements, also known as a "code snippet":

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
Mixing statements and script declarations was necessary in DSL1 and allowed in DSL2. However, it is no longer supported in the Nextflow language specification in order to simplify the language and to ensure that top-level statements are only executed when the script is executed directly and not when it is included as a module.
:::

**Assignment expressions**

In Groovy, you can assign a variable as part of an expression:

```groovy
foo(x = 1, y = 2)
```

The Nextflow language specification only allows assignment as statements:

```nextflow
x = 1
y = 2
foo(x, y)
```

The increment (`++`) and decrement (`--`) operators are no longer supported. Use `+=` and `-=` instead:

```nextflow
x += 1 // x++
x -= 1 // x--
```

**For and while loops**

Groovy supports loop statements such as `for` and `while`:

```groovy
for (rseqc_module in ['read_distribution', 'inner_distance', 'tin']) {
    if (rseqc_modules.contains(rseqc_module))
        rseqc_modules.remove(rseqc_module)
}
```

The Nextflow language specification does not support loop statements. Use higher-order functions like the `each` method instead:

```nextflow
['read_distribution', 'inner_distance', 'tin'].each { rseqc_module ->
    if (rseqc_modules.contains(rseqc_module))
        rseqc_modules.remove(rseqc_module)
}
```

Lists, maps, and sets provide several functions (e.g., `collect`, `find`, `findAll`, `inject`) for iteration. See [Groovy standard library](https://docs.groovy-lang.org/latest/html/groovy-jdk/overview-summary.html) for more information.

**Switch statements**

Groovy supports switch statements for pattern matching on a value:

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

The Nextflow language specification does not support switch statements. Use if-else statements instead:

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

**Slashy dollar strings**

Groovy supports a wide variety of strings, including multi-line strings, dynamic strings, slashy strings, multi-line dynamic slashy strings, and more.

The Nextflow language specification supports single- and double-quoted strings, multi-line strings, and slashy strings. Dynamic slashy strings are not supported:

```groovy
def logo = /--cl-config 'custom_logo: "${multiqc_logo}"'/
```

Use a double-quoted string instead:

```nextflow
def logo = "--cl-config 'custom_logo: \"${multiqc_logo}\"'"
```

Slashy dollar strings are not supported:

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

**Implicit environment variables**

In Nextflow DSL1 and DSL2, you can reference environment variables directly in strings:

```nextflow
println "PWD = ${PWD}"
```

The Nextflow language specification does not support implicit environment variables. Use `System.getenv()` instead:

```nextflow
println "PWD = ${System.getenv('PWD')}"
```

### Restricted syntax

The following patterns are still supported but have been restricted, i.e. some syntax variants have been removed.

**Variable declarations**

In Groovy, variables can be declared in many different ways:

```groovy
def a = 1
final b = 2
def c = 3, d = 4
def (e, f) = [5, 6]
String str = 'foo'
def Map meta = [:]
```

In Nextflow, variables should be declared with `def` and should not specify a type:

```nextflow
def a = 1
def b = 2
def (c, d) = [3, 4]
def (e, f) = [5, 6]
def str = 'foo'
def meta = [:]
```

Similarly, functions should be declared with `def` and should not specify a return type or parameter types:

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

To ease the migration of existing scripts, the language server only reports warnings for Groovy-style type annotations and implicit variable declarations. These warnings will become errors in the future.

:::{note}
Type annotations and static type checking will be addressed in a future version of the Nextflow language specification.
:::

**Process env inputs/outputs**

In Nextflow DSL1 and DSL2, the name of a process `env` input/output can be specified with or without quotes:

```nextflow
process PROC {
    input:
    env FOO
    env 'BAR'
}
```

The Nextflow language specification requires the name to be specified with quotes:

```nextflow
process PROC {
    input:
    env 'FOO'
    env 'BAR'
}
```

**Implicit process script section**

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

The Nextflow language specification allows the `script:` label to be omitted only if there are no other sections:

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

### Deprecated syntax

The following patterns are deprecated. The language server reports "future warnings" for these patterns. Future warnings are disabled by default. Enable them by deselecting **Nextflow > Suppress Future Warnings** in the [extension settings](#settings). These warnings may become errors in the future.

**Implicit closure parameter**

In Groovy, a closure with no parameters is assumed to have a single parameter named `it`. The Nextflow language specification does not support implicit closure parameters. Declare the parameter explicitly instead:

```nextflow
ch | map { it * 2 }       // deprecated
ch | map { v -> v * 2 }   // correct
ch | map { it -> it * 2 } // also correct
```

**Using params outside the entry workflow**

While params can be used anywhere in the pipeline code, they are only intended to be used in the entry workflow.

Processes and workflows should receive params as explicit inputs:

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

**Process each input**

The `each` process input is deprecated. Use the `combine` or `cross` operator to explicitly repeat over inputs in the calling workflow.

**Process when section**

The process `when` section is deprecated. Use conditional logic, such as an `if` statement or the `filter` operator, to control the process invocation in the calling workflow.

### Configuration syntax

See {ref}`config-syntax` for a comprehensive description of the configuration language.

Currently, Nextflow parses config files as Groovy scripts, allowing the use of scripting constructs like variables, helper functions, and conditional logic for dynamic configuration. For example:

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
})()
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

## Troubleshooting

In the event of an error, you can stop or restart the language server from the command palette. See [Commands](#commands) for the set of available commands.

Report issues at [nextflow-io/vscode-language-nextflow](https://github.com/nextflow-io/vscode-language-nextflow) or [nextflow-io/language-server](https://github.com/nextflow-io/language-server). When reporting, include a minimal code snippet that reproduces the issue and any error logs from the server. To view logs, open the **Output** tab and select **Nextflow Language Server** from the dropdown. Enable **Nextflow > Debug** in the [extension settings](#settings) to show additional log messages while debugging.

## Limitations

- The language server does not detect certain filesystem changes, such as changing the current Git branch. Restart the language server from the command palette to sync it with your workspace.

- The language server does not recognize configuration options from third-party plugins and will report "Unrecognized config option" warnings for them.

- The language server provides limited support for Groovy scripts in the `lib` directory. Errors in Groovy scripts are not reported as diagnostics, and changing a Groovy script does not automatically re-compile the Nextflow scripts that reference it. Edit the Nextflow script or close and re-open it to refresh the diagnostics.

## Commands

The following commands are available from the command palette:

- Restart language server
- Stop language server

## Settings

The following settings are available:

`nextflow.debug`
: Enable debug logging and debug information in hover hints.

`nextflow.files.exclude`
: Configure glob patterns for excluding folders from being searched for Nextflow scripts and configuration files.

`nextflow.formatting.harshilAlignment`
: Use the [Harshil Alignment™️](https://nf-co.re/docs/contributing/code_editors_and_styling/harshil_alignment) when formatting Nextflow scripts and config files.

`nextflow.java.home`
: Specifies the folder path to the JDK. Use this setting if the extension cannot find Java automatically.

`nextflow.suppressFutureWarnings`
: Hide warnings for future changes, deprecations, and removals.

## Language server

Most of the functionality of the VS Code extension is provided by the [Nextflow language server](https://github.com/nextflow-io/language-server), which implements the [Language Server Protocol (LSP)](https://microsoft.github.io/language-server-protocol/) for Nextflow scripts and config files.

The language server is distributed as a standalone Java application. It can be integrated with any editor that functions as an LSP client. Currently, only the VS Code integration is officially supported, but community contributions for other editors are welcome.

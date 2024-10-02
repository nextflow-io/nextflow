
# VS Code integration

The [Nextflow VS Code extension](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow) provides IDE support for Nextflow pipelines.

## Features

### Syntax highlighting

Nextflow scripts and config files are highlighted for better readability.

### Diagnostics

Source code is highlighted red or yellow when there are errors or warnings, respectively.

To view all diagnostics for the workspace, open the "Problems" tab. From there, you can search for diagnostics based on the diagnostic message, filename, etc.

An additional class of warnings, known as "future warnings", will report deprecated features or other issues that are not errors but may become errors in the future. These warnings are disabled by default, and can be enabled by unselecting "Nextflow > Suppress Future Warnings" in the [extension settings](#settings).

### Hover hints

Some source code elements such as definitions, variable names, and function calls will provide a tooltip when they are hovered over with the cursor, with related information such as the definition of a variable or function and documentation.

### Code navigation

When viewing a script, the "Outline" section in the Explorer panel will list the top-level script definitions, including workflows, processes, and functions.

Include declarations in both scripts and config files are links, and will open the corresponding script or config file when selected via ctrl-click.

To view the definition of a symbol (e.g. a workflow, process, function, or variable), right-click on the symbol and select "Go to Definition". You can also view all references of a symbol, as well as the call hierarchy of a function call.

### Code completion

The VS Code extension can make suggestions while typing, such as for a variable name, function name, or config setting, depending on the context.

Additionally, there are several prepared snippets for common script declarations such as processes and workflows.

### Formatting

The "Format Document" command will format the current script or config file based on a standard set of whitespace rules.

The "Nextflow > Formatting" [extension settings](#settings) can be used to customize the formatting.

### Renaming symbols

To rename a symbol, right-click on the symbol, select "Rename Symbol", and enter the new name. The symbol definition and all references will be renamed.

### Parameter schema

If there is a `nextflow_schema.json` in the same directory as a script with an entry workflow, the schema will be used to validate params and provide hover hints and code completions in the entry workflow.

### DAG preview for workflows

To preview the DAG of a workflow, select the "Preview DAG" CodeLens above the workflow definition. The workflow DAG will be opened in a new panel to the side. The DAG includes the workflow inputs and outputs, as well as any processes and workflows that are called in the selected workflow.

:::{note}
The "Preview DAG" CodeLens will not be displayed if there are any errors in the script.
:::

## Syntax guide

The language server parses scripts and config files according to the [Nextflow language specification](https://deploy-preview-5336--nextflow-docs-staging.netlify.app/reference/syntax). It is more strict than the Nextflow CLI, which allows all Groovy syntax. As a result, the language server can perform more extensive error checking and provide more specific error messages. However, some Groovy syntax is no longer supported.

This section describes some of the most common features that are no longer supported and how to work around them. Refer to the Nextflow language specification for a comprehensive description of what is supported.

In general, any code can be moved into the `lib` directory or into a plugin, both of which support the full Groovy language.

:::{note}
The "Nextflow language specification" is essentially a strict specification of DSL2. It is not a new DSL version. Moving forward, the Nextflow language will be defined by this specification rather than new DSL versions.
:::

### Excluded syntax

**Import declarations**

Groovy has an `import` declaration for importing external classes:

```groovy
import groovy.json.JsonSlurper

def json = new JsonSlurper().parseText(json_file.text)
```

Instead, you can use the fully qualified name directly:

```nextflow
def json = new groovy.json.JsonSlurper().parseText(json_file.text)
```

**Class declarations**

Class declarations are sometimes used in Nextflow to define helper functions or custom record types. Helper functions can be defined as standalone functions in a script. Class declarations should be moved to the `lib` directory.

:::{note}
Record types will be addressed in a future version of the Nextflow language specification.
:::

:::{note}
Enums, a special type of class, are supported, but cannot be included across modules at this time.
:::

**Mixing script declarations and statements**

A script may contain any of the following "top-level" declarations:

- feature flags
- includes
- parameter declarations
- workflows
- processes
- functions
- output block

Alternatively, a script may contain only statements, e.g. as a code snippet:

```nextflow
println 'Hello world!'
```

Code snippets are treated as an "implicit" entry workflow:

```nextflow
workflow {
    println 'Hello world!'
}
```

However, script declarations and statements cannot be mixed at the same level. This practice was necessary in DSL1 and allowed in DSL2, but is not supported by the Nextflow language specification. Unless the script is a code snippet, all statements must reside within script declarations. Typically, such statements should be placed in the entry workflow:

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

**Variable declarations**

Groovy provides many different ways to declare variables, including multi-declarations and optional type annotations:

```groovy
def a = 1
final b = 2
def c = 3, d = 4
def (e, f) = [5, 6]
String str = 'foo'
def Map meta = [:]
```

The Nextflow language specification requires that variables be declared with `def`:

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

To ease the migration of existing scripts, the language server will only report warnings for Groovy-style type annotations and implicit variable declarations (in processes and workflows). These warnings will become errors in the future.

:::{note}
Type annotations and static type checking will be addressed in a future version of the Nextflow language specification.
:::

**Variable assignments**

In Groovy, a variable can be assigned as part of an expression:

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

Loop statements are not supported in the Nextflow language specification. Use higher-order functions such as the `each` method instead:

```nextflow
['read_distribution', 'inner_distance', 'tin'].each { rseqc_module ->
    if (rseqc_modules.contains(rseqc_module))
        rseqc_modules.remove(rseqc_module)
}
```

Lists, maps, and sets provide several such functions (e.g. `collect`, `find`, `findAll`, `inject`) for iteration. Refer to the [Groovy standard library](https://docs.groovy-lang.org/latest/html/groovy-jdk/overview-summary.html) for more information.

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

Groovy supports a wide variety of strings, including multi-line strings, dynamic strings, slashy strings, multi-line dynamic slashy strings, and so on.

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

In Nextflow DSL1 and DSL2, environment variables can be referenced directly in strings:

```nextflow
println "PWD = ${PWD}"
```

The Nextflow language specification does not support implicit environment variables. Use `System.getenv()` instead:

```nextflow
println "PWD = ${System.getenv('PWD')}"
```

### Deprecated syntax

**Implicit closure parameter**

In Groovy, a closure that doesn’t define any parameters is assumed to have a single parameter named `it`. The Nextflow language specification does not support implicit closure parameters. The language server will report a warning for implicit closure parameters. Declare the parameter explicitly instead:

```nextflow
ch | map { it * 2 }       // deprecated
ch | map { it -> it * 2 } // correct
ch | map { v -> v * 2 }   // also correct
```

**Using params outside the entry workflow**

While params can be used anywhere in the pipeline code, they are only intended to be used in the entry workflow. The language server will report a warning for params used outside the entry workflow.

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

The `each` process input is deprecated, and the language server will report a warning for it. Use the `combine` or `cross` operator to explicitly repeat over inputs in the calling workflow.

**Process when section**

The process `when` section is deprecated, and the language server will report a warning for it. Use conditional logic, such as an `if` statement or `filter` operator, to control the process execution in the calling workflow.

### Configuration syntax

See [Configuration](https://deploy-preview-5336--nextflow-docs-staging.netlify.app/config#syntax) for a comprehensive description of what is supported in the configuration language.

Historically, config files were parsed as Groovy scripts, which allowed the use of scripting constructs like variables, helper functions, and conditional logic to perform dynamic configuration. For example:

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

The strict config syntax does not support functions, and only allows statements like variables and if statements in closure expressions. This kind of configuration can be achieved with a dynamic include:

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

The include source is a closure that is immediately invoked. It will include a different config file based on some condition. Including `/dev/null` is equivalent to including nothing.

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

When debugging, it can be useful to enable "Nextflow > Debug" in the [extension settings](#settings).

Issues can be reported at [nextflow-io/vscode-language-nextflow](https://github.com/nextflow-io/vscode-language-nextflow) or [nextflow-io/language-server](https://github.com/nextflow-io/language-server). When reporting an issue, please include a minimal code snippet that triggers the issue as well as any error traces from the server logs. To view the server logs, open the "Output" tab and select "Nextflow Language Server" in the dropdown list.

## Limitations

- The language server cannot currently detect certain filesystem changes, such as switching to a different Git branch. Restart the language server from the command palette to ensure that it is in sync with your workspace.

- The language server does not currently support configuration options from third-party plugins, so it will report "Unrecognized config option" warnings for those settings.

- The language server does not currently support the `lib` directory, so it will report an error for any reference to the `lib` directory in a script.

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

The language server is distributed as a standalone Java application. It can be integrated with any editor that functions as an LSP client. Currently only the VS Code integration is officially supported, but we welcome community contributions for other editors.

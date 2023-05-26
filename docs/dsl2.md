(dsl2-page)=

# DSL 2

Nextflow provides a syntax extension that allows the definition of module libraries and simplifies the writing of complex data analysis pipelines.

To enable this feature you need to define the following directive at the beginning of your workflow script:

```groovy
nextflow.enable.dsl=2
```

:::{versionchanged} 22.03.0-edge
Nextflow uses DSL2 if no version is specified explicitly. You can restore the previous behavior by setting the following environment variable:

```bash
export NXF_DEFAULT_DSL=1
```
:::

:::{versionchanged} 22.03.0-edge
The DSL version specification (either 1 or 2) can also be specified in the Nextflow configuration file using the same notation shown above.
:::

:::{versionchanged} 22.11.0-edge
Support for DSL1 was removed from Nextflow.
:::

## Function

Nextflow allows the definition of custom functions in the workflow script using the following syntax:

```groovy
def <function name> ( arg1, arg, .. ) {
    <function body>
}
```

For example:

```groovy
def foo() {
    'Hello world'
}

def bar(alpha, omega) {
    alpha + omega
}
```

The above snippet defines two simple functions, that can be invoked in the workflow script as `foo()` which returns the `Hello world` string and `bar(10,20)` which returns the sum of two parameters (`30` in this case).

:::{note}
Functions implicitly return the result of the last evaluated statement.
:::

The keyword `return` can be used to explicitly exit from a function and return the specified value. For example:

```groovy
def fib( x ) {
    if( x <= 1 )
        return x
    else
        fib(x-1) + fib(x-2)
}
```

## Process

### Process definition

The new DSL separates the definition of a process from its invocation. The process definition follows the usual syntax as described in the {ref}`process documentation <process-page>`. The only difference is that the `from` and `into` channel declarations have to be omitted.

Then a process can be invoked as a function in the `workflow` scope, passing the expected input channels as parameters as if it were a custom function. For example:

```groovy
nextflow.enable.dsl=2

process foo {
    output:
      path 'foo.txt'

    script:
      """
      your_command > foo.txt
      """
}

process bar {
    input:
      path x

    output:
      path 'bar.txt'

    script:
      """
      another_command $x > bar.txt
      """
}

workflow {
    data = channel.fromPath('/some/path/*.txt')
    foo()
    bar(data)
}
```

:::{warning}
A process component can be invoked only once in the same workflow context.
:::

### Process composition

Processes having matching *input-output* declaration can be composed so that the output of the first process is passed as input to the next process. Taking in consideration the previous example, it's possible to write the following:

```groovy
workflow {
    bar(foo())
}
```

### Process output

A process output can also be accessed using the `out` attribute on the corresponding process object. For example:

```groovy
workflow {
    foo()
    bar(foo.out)
    bar.out.view()
}
```

When a process defines two or more output channels, each of them can be accessed using the array element operator e.g. `out[0]`, `out[1]`, etc. or using *named outputs* (see below).

### Process named output

The `emit` option can be added to the process output definition to assign a name identifier. This name can be used to reference the channel within the caller scope. For example:

```groovy
process foo {
  output:
    path '*.bam', emit: samples_bam

  '''
  your_command --here
  '''
}

workflow {
    foo()
    foo.out.samples_bam.view()
}
```

### Process named stdout

The `emit` option can be used also to name the stdout:

```groovy
process sayHello {
    input:
        val cheers

    output:
        stdout emit: verbiage

    script:
    """
    echo -n $cheers
    """
}

workflow {
    things = channel.of('Hello world!', 'Yo, dude!', 'Duck!')
    sayHello(things)
    sayHello.out.verbiage.view()
}
```

:::{note}
Optional params for a process input/output are always prefixed with a comma, except for `stdout`. Because `stdout` does not have an associated name or value like other types, the first param should not be prefixed.
:::

## Workflow

### Workflow definition

The `workflow` keyword allows the definition of sub-workflow components that enclose the invocation of one or more processes and operators:

```groovy
workflow my_pipeline {
    foo()
    bar( foo.out.collect() )
}
```

For example, the above snippet defines a workflow component, named `my_pipeline`, that can be invoked from another workflow component definition as any other function or process with `my_pipeline()`.

### Workflow parameters

A workflow component can access any variable and parameter defined in the outer scope:

```groovy
params.data = '/some/data/file'

workflow my_pipeline {
    if( params.data )
        bar(params.data)
    else
        bar(foo())
}
```

### Workflow input

A workflow component can declare one or more input channels using the `take` keyword. For example:

```groovy
workflow my_pipeline {
    take: data
    main:
        foo(data)
        bar(foo.out)
}
```

:::{warning}
When the `take` keyword is used, the beginning of the workflow body must be identified with the `main` keyword.
:::

Then, the input can be specified as an argument in the workflow invocation statement:

```groovy
workflow {
    my_pipeline( channel.from('/some/data') )
}
```

:::{note}
Workflow inputs are always channels by definition. If a basic data type is provided instead, such as a number, string, list, etc, it is implicitly converted to a {ref}`value channel <channel-type-value>`.
:::

### Workflow output

A workflow component can declare one or more output channels using the `emit` keyword. For example:

```groovy
workflow my_pipeline {
    main:
      foo(data)
      bar(foo.out)
    emit:
      bar.out
}
```

Then, the result of the `my_pipeline` execution can be accessed using the `out` property, i.e. `my_pipeline.out`. When multiple output channels are declared, use the array bracket notation to access each output channel as described for the [Process output](#process-output) definition.

### Workflow named output

If the output channel is assigned to an identifier in the `emit` declaration, such identifier can be used to reference the channel within the caller scope. For example:

```groovy
workflow my_pipeline {
   main:
     foo(data)
     bar(foo.out)
   emit:
     my_data = bar.out
}
```

Then, the result of the above snippet can accessed using `my_pipeline.out.my_data`.

### Workflow entrypoint

A workflow definition which does not declare any name (also known as *implicit workflow*) is the entry point of execution for the workflow application.

:::{note}
Implicit workflow definition is ignored when a script is included as a module. This allows the writing of a workflow script that can be used either as a library module or as an application script.
:::

:::{tip}
A different workflow entrypoint can be specified using the `-entry` command line option.
:::

### Workflow composition

Workflows defined in your script or imported with [Module inclusion](#module-inclusion) can be invoked and composed as any other process in your application.

```groovy
workflow flow1 {
    take: data
    main:
        foo(data)
        bar(foo.out)
    emit:
        bar.out
}

workflow flow2 {
    take: data
    main:
        foo(data)
        baz(foo.out)
    emit:
        baz.out
}

workflow {
    take: data
    main:
        flow1(data)
        flow2(flow1.out)
}
```

:::{note}
Nested workflow execution determines an implicit scope. Therefore the same process can be invoked in two different workflow scopes, like for example `foo` in the above snippet that is used both in `flow1` and `flow2`. The workflow execution path, along with the process names, determines the *fully qualified process name* that is used to distinguish the two different process invocations, i.e. `flow1:foo` and `flow2:foo` in the above example.
:::

:::{tip}
The fully qualified process name can be used as a valid {ref}`process selector <config-process-selectors>` in the `nextflow.config` file and it has priority over the simple process name.
:::

## Modules

The new DSL allows the definition of *module scripts* that can be included and shared across workflow applications.

A module script (or simply, module) can contain the definition of functions, processes and workflows as described in the previous sections.

:::{note}
Functions, processes and workflows are globally referred to as *components*.
:::

### Module inclusion

A component defined in a module script can be imported into another Nextflow script using the `include` keyword.

For example:

```groovy
include { foo } from './some/module'

workflow {
    data = channel.fromPath('/some/data/*.txt')
    foo(data)
}
```

The above snippet includes a process with name `foo` defined in the module script in the main execution context. This way, `foo` can be invoked in the `workflow` scope.

Nextflow implicitly looks for the script file `./some/module.nf` resolving the path against the *including* script location.

:::{note}
Relative paths must begin with the `./` prefix. Also, the `include` statement must be defined **outside** of the workflow definition.
:::

(dsl2-module-directory)=

### Module directory

:::{versionadded} 22.10.0
:::

The module can be defined as a directory whose name matches the module name and contains a script named `main.nf`. For example:

```
some
└-module
  └-main.nf
```

When defined as a directory the module needs to be included specifying the module directory path:

```groovy
include { foo } from './some/module'
```

Module directories allows the use of module scoped binaries scripts. See [Module binaries](#module-binaries) for details.

### Multiple inclusions

A Nextflow script allows the inclusion of an arbitrary number of modules and components. When multiple components need to be included from the same module script, the component names can be specified in the same inclusion using the curly brackets notation as shown below:

```groovy
include { foo; bar } from './some/module'

workflow {
    data = channel.fromPath('/some/data/*.txt')
    foo(data)
    bar(data)
}
```

### Module aliases

When including a module component, it's possible to specify an *alias* with the `as` keyword. This allows the inclusion and the invocation of components with the same name in your script using different names. For example:

```groovy
include { foo } from './some/module'
include { foo as bar } from './other/module'

workflow {
    foo(some_data)
    bar(other_data)
}
```

The same is possible when including the same component multiple times from the same module script as shown below:

```groovy
include { foo; foo as bar } from './some/module'

workflow {
    foo(some_data)
    bar(other_data)
}
```

### Module parameters

A module script can define one or more parameters using the same syntax of a Nextflow workflow script:

```groovy
params.foo = 'Hello'
params.bar = 'world!'

def sayHello() {
    println "$params.foo $params.bar"
}
```

Then, parameters are inherited from the including context. For example:

```groovy
params.foo = 'Hola'
params.bar = 'Mundo'

include {sayHello} from './some/module'

workflow {
    sayHello()
}
```

The above snippet prints:

```
Hola Mundo
```

:::{note}
The module inherits the parameters defined *before* the `include` statement, therefore any further parameter set later is ignored.
:::

:::{tip}
Define all pipeline parameters at the beginning of the script *before* any `include` declaration.
:::

The option `addParams` can be used to extend the module parameters without affecting the external scope. For example:

```groovy
include {sayHello} from './some/module' addParams(foo: 'Ciao')

workflow {
    sayHello()
}
```

The above snippet prints:

```
Ciao world!
```

Finally, the include option `params` allows the specification of one or more parameters without inheriting any value from the external environment.

(module-templates)=

### Module templates

The module script can be defined in an external {ref}`template <process-template>` file. With DSL2 the template file can be placed under the `templates` directory where the module script is located.

For example, let's suppose to have a project L with a module script defining 2 processes (P1 and P2) and both use templates. The template files can be made available under the local `templates` directory:

```
Project L
|─myModules.nf
└─templates
  |─P1-template.sh
  └─P2-template.sh
```

Then, we have a second project A with a workflow that includes P1 and P2:

```
Pipeline A
└-main.nf
```

Finally, we have a third project B with a workflow that includes again P1 and P2:

```
Pipeline B
└-main.nf
```

With the possibility to keep the template files inside the project L, A and B can use the modules defined in L without any changes. A future project C would do the same, just cloning L (if not available on the system) and including its module script.

Beside promoting sharing modules across pipelines, there are several advantages in keeping the module template under the script path:

1. module components are *self-contained*,
2. module components can be tested independently from the pipeline(s) importing them,
3. it is possible to create libraries of module components.

Ultimately, having multiple template locations allows a more structured organization within the same project. If a project has several module components, and all them use templates, the project could group module scripts and their templates as needed. For example:

```
baseDir
|─main.nf
|─Phase0-Modules
  |─mymodules1.nf
  |─mymodules2.nf
  └─templates
    |─P1-template.sh
    └─P2-template.sh
|─Phase1-Modules
  |─mymodules3.nf
  |─mymodules4.nf
  └─templates
    |─P3-template.sh
    └─P4-template.sh
└─Phase2-Modules
  |─mymodules5.nf
  |─mymodules6.nf
  └─templates
    |─P5-template.sh
    |─P6-template.sh
    └─P7-template.sh
```

(module-binaries)=

### Module binaries

:::{versionadded} 22.10.0
:::

Modules can define binary scripts that are locally scoped to the processes defined by the tasks.

To enable this feature add the following setting in pipeline configuration file:

```groovy
nextflow.enable.moduleBinaries = true
```

The binary scripts must be placed in the module directory names `<module-dir>/resources/usr/bin`:

```
<module-dir>
|─main.nf
└─resources
  └─usr
    └─bin
      |─your-module-script1.sh
      └─another-module-script2.py
```

Those scripts will be accessible as any other command in the tasks environment, provided they have been granted the Linux execute permissions.

:::{note}
This feature requires the use of a local or shared file system as the pipeline work directory or {ref}`wave-page` when using cloud based executors.
:::

## Channel forking

Using the new DSL, Nextflow channels are automatically forked when connecting two or more consumers.

For example:

```groovy
channel
    .from('Hello','Hola','Ciao')
    .set{ cheers }

cheers
    .map{ it.toUpperCase() }
    .view()

cheers
    .map{ it.reverse() }
    .view()
```

The same is valid for the result (channel) of a process execution. Therefore a process output can be consumed by two or more processes without the need to fork it using the `into` operator, making the writing of workflow scripts more fluent and readable.

## Pipes

### The *pipe* operator

Nextflow processes and operators can be composed using the `|` *pipe* operator. For example:

```groovy
process foo {
    input:
    val data

    output:
    val result

    exec:
    result = "$data world"
}

workflow {
   channel.from('Hello','Hola','Ciao') | foo | map { it.toUpperCase() } | view
}
```

The above snippet defines a process named `foo` and invokes it passing the content of the `data` channel. The result is then piped to the {ref}`operator-map` operator which converts each string to uppercase and finally, the last {ref}`operator-view` operator prints it.

### The *and* operator

The `&` *and* operator allows feeding of two or more processes with the content of the same channel(s). For example:

```groovy
process foo {
    input:
    val data

    output:
    val result

    exec:
    result = "$data world"
}

process bar {
    input:
    val data

    output:
    val result

    exec:
    result = data.toUpperCase()
}

workflow {
    channel.from('Hello') | map { it.reverse() } | (foo & bar) | mix | view
}
```

In the above snippet the channel emitting the `Hello` string is piped with the {ref}`operator-map` which reverses the string value. Then, the result is passed to both `foo` and `bar` processes which are executed in parallel. Each process outputs a channel, and the two channels are merged into a single channel using the {ref}`operator-mix` operator. Finally the result is printed using the {ref}`operator-view` operator.

:::{tip}
The break-line operator `\` can be used to split long statements over multiple lines. The above snippet can also be written as:

```groovy
workflow {
    channel.from('Hello') \
      | map { it.reverse() } \
      | (foo & bar) \
      | mix \
      | view
}
```
:::

## DSL2 migration notes

- DSL2 final version is activated using the declaration `nextflow.enable.dsl=2` in place of `nextflow.preview.dsl=2`.

- Process inputs of type `set` have to be replaced with {ref}`tuple <process-input-tuple>`.

- Process outputs of type `set` have to be replaced with {ref}`tuple <process-out-tuple>`.

- Process output option `mode flatten` is no longer available. Replace it using the {ref}`operator-flatten` operator on the corresponding output channel.

- Anonymous and unwrapped includes are not supported anymore. Replace them with an explicit module inclusion. For example:

  ```groovy
  include './some/library'
  include bar from './other/library'

  workflow {
    foo()
    bar()
  }
  ```

  Should be replaced with:

  ```groovy
  include { foo } from './some/library'
  include { bar } from './other/library'

  workflow {
    foo()
    bar()
  }
  ```

- The use of unqualified value and file elements into input tuples is not allowed anymore. Replace them with a corresponding
  `val` or `path` qualifier:

  ```groovy
  process foo {
  input:
    tuple X, 'some-file.bam'

  script:
    '''
    your_command --in $X some-file.bam
    '''
  }
  ```

  Use:

  ```groovy
  process foo {
  input:
    tuple val(X), path('some-file.bam')

  script:
    '''
    your_command --in $X some-file.bam
    '''
  }
  ```

- The use of unqualified value and file elements into output tuples is not allowed anymore. Replace them with a corresponding
  `val` or `path` qualifier:

  ```groovy
  process foo {
  output:
    tuple X, 'some-file.bam'

  script:
    X = 'some value'
    '''
    your_command > some-file.bam
    '''
  }
  ```

  Use:

  ```groovy
  process foo {
  output:
    tuple val(X), path('some-file.bam')

  script:
    X = 'some value'
    '''
    your_command > some-file.bam
    '''
  }
  ```

- The `bind` channel method (and corresponding `<<` operator) has been deprecated in DSL2.

- The `choice` operator has been deprecated in DSL2. Use {ref}`operator-branch` instead.

- The `close` operator has been deprecated in DSL2.

- The `countBy` operator has been deprecated in DSL2.

- The `create` channel factory has been deprecated in DSL2.

- The `fork` operator has been renamed to {ref}`operator-multimap`.

- The `groupBy` operator has been deprecated in DSL2. Replace it with {ref}`operator-grouptuple`

- The `into` operator has been deprecated in DSL2 since it's not needed anymore.

- The `print` and `println` operators have been deprecated in DSL2. Use {ref}`operator-view` instead.

- The `route` operator has been deprecated in DSL2.

- The `separate` operator has been deprecated in DSL2.

- The `spread` operator has been deprecated in DSL2. Replace it with {ref}`operator-combine`.

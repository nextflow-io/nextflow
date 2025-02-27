(module-page)=

# Modules

Nextflow scripts can include **definitions** (workflows, processes, and functions) from other scripts. When a script is included in this way, it is referred to as a **module**. Modules can be included by other modules or pipeline scripts and can even be shared across workflows.

:::{note}
Modules were introduced in DSL2. If you are still using DSL1, see the {ref}`dsl1-page` page to learn how to migrate your Nextflow pipelines to DSL2.
:::

## Module inclusion

You can include any definition from a module into a Nextflow script using the `include` keyword.

For example:

```nextflow
include { foo } from './some/module'

workflow {
    data = channel.fromPath('/some/data/*.txt')
    foo(data)
}
```

The above snippet imports a process named `foo`, defined in the module, into the main execution context. This way, `foo` can be invoked in the `workflow` scope.

Nextflow implicitly looks for the script file `./some/module.nf`, resolving the path against the *including* script location.

Module includes are subject to the following rules:

- Relative paths must begin with the `./` prefix.
- Include statements are not allowed from within a workflow. They must occur at the script level.

(module-directory)=

## Module directory

:::{versionadded} 22.10.0
:::

A module can be defined as a directory with the same name as the module and with a script named `main.nf`. For example:

```
some
└── module
    └── main.nf
```

When defined as a directory, the module must be included by specifying the module directory path:

```nextflow
include { foo } from './some/module'
```

Module directories allow the use of module scoped binaries scripts. See [Module binaries] for details.

## Multiple inclusions

A Nextflow script can include any number of modules, and an `include` statement can import any number of definitions from a module. Multiple definitions can be included from the same module by using the syntax shown below:

```nextflow
include { foo; bar } from './some/module'

workflow {
    data = channel.fromPath('/some/data/*.txt')
    foo(data)
    bar(data)
}
```

(module-aliases)=

## Module aliases

When including definition from a module, it's possible to specify an *alias* with the `as` keyword. Aliasing allows you to avoid module name clashes, by assigning them different names in the including context. For example:

```nextflow
include { foo } from './some/module'
include { foo as bar } from './other/module'

workflow {
    foo(some_data)
    bar(other_data)
}
```

You can also include the same definition multiple times under different names:

```nextflow
include { foo; foo as bar } from './some/module'

workflow {
    foo(some_data)
    bar(other_data)
}
```

(module-params)=

## Module parameters

:::{deprecated} 24.07.0-edge
As a best practice, parameters should be used in the entry workflow and passed to workflows, processes, and functions as explicit inputs.
:::

A module can define parameters using the same syntax as a Nextflow workflow script:

```nextflow
params.foo = 'Hello'
params.bar = 'world!'

def sayHello() {
    println "$params.foo $params.bar"
}
```

When including a module, the module will first use parameters from the including context. For example:

```nextflow
params.foo = 'Hola'
params.bar = 'Mundo'

include { sayHello } from './some/module'

workflow {
    sayHello()
}
```

The above snippet prints:

```
Hola Mundo
```

:::{note}
The module inherits the parameters defined *before* the `include` statement, therefore any parameters set afterwards will not be used by the module.
:::

:::{tip}
It is best to define all pipeline parameters *before* any `include` statements.
:::

The `addParams` option can be used to pass parameters to the module without adding them to the including scope.

```nextflow
params.foo = 'Hola'
params.bar = 'Mundo'

include { sayHello } from './some/module' addParams(foo: 'Ciao')

workflow {
    sayHello()
}
```

The above snippet prints:

```
Ciao Mundo
```

Alternatively, the `params` option can be used to pass parameters to module without adding them to the including scope, *and* without inheriting any parameters from the including scope.

```nextflow
params.foo = 'Hola'
params.bar = 'Mundo'

include { sayHello } from './some/module' params(foo: 'Ciao')

workflow {
    sayHello()
}
```

The above snippet prints:

```
Ciao world!
```

(module-templates)=

## Module templates

Process script {ref}`templates <process-template>` can be included alongside a module in the `templates` directory.

For example, suppose we have a project L with a module that defines two processes, P1 and P2, both of which use templates. The template files can be made available in the local `templates` directory:

```
Project L
|── myModules.nf
└── templates
    |── P1-template.sh
    └── P2-template.sh
```

Then, we have a second project A with a workflow that includes P1 and P2:

```
Pipeline A
└── main.nf
```

Finally, we have a third project B with a workflow that also includes P1 and P2:

```
Pipeline B
└── main.nf
```

With the possibility to keep the template files inside the project L, A and B can use the modules defined in L without any changes. A future project C would do the same, just cloning L (if not available on the system) and including its module.

Beside promoting the sharing of modules across pipelines, there are several advantages to keeping the module template under the script path:

1. Modules are self-contained
2. Modules can be tested independently from the pipeline(s) that import them
3. Modules can be made into libraries

Having multiple template locations enables a structured project organization. If a project has several modules, and they all use templates, the project could group module scripts and their templates as needed. For example:

```
baseDir
|── main.nf
|── Phase0-Modules
    |── mymodules1.nf
    |── mymodules2.nf
    └── templates
        |── P1-template.sh
        |── P2-template.sh
|── Phase1-Modules
    |── mymodules3.nf
    |── mymodules4.nf
    └── templates
        |── P3-template.sh
        └── P4-template.sh
└── Phase2-Modules
    |── mymodules5.nf
    |── mymodules6.nf
    └── templates
        |── P5-template.sh
        |── P6-template.sh
        └── P7-template.sh
```

(module-binaries)=

## Module binaries

:::{versionadded} 22.10.0
:::

Modules can define binary scripts that are locally scoped to the processes defined by the tasks.

To enable this feature, set the following flag in your pipeline script or configuration file:

```nextflow
nextflow.enable.moduleBinaries = true
```

The binary scripts must be placed in the module directory names `<module-dir>/resources/usr/bin`:

```
<module-dir>
|── main.nf
└── resources
    └── usr
        └── bin
            |── your-module-script1.sh
            └── another-module-script2.py
```

Those scripts will be made accessible like any other command in the task environment, provided they have been granted the Linux execute permissions.

:::{note}
This feature requires the use of a local or shared file system for the pipeline work directory, or {ref}`wave-page` when using cloud-based executors.
:::

## Sharing modules

Modules are designed to be easy to share and re-use across different pipelines, which helps eliminate duplicate work and spread improvements throughout the community. While Nextflow does not provide an explicit mechanism for sharing modules, there are several ways to do it:

- Simply copy the module files into your pipeline repository
- Use [Git submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules) to fetch modules from other Git repositories without maintaining a separate copy
- Use the [nf-core](https://nf-co.re/tools#modules) CLI to install and update modules with a standard approach used by the nf-core community

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
include { cat } from './some/module'

workflow {
    data = channel.fromPath('/some/data/*.txt')
    cat(data)
}
```

The above snippet imports a process named `cat`, defined in the module, into the main execution context. This way, `cat` can be invoked in the `workflow` scope.

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
include { hello } from './some/module'
```

Module directories allow the use of module scoped binary scripts. See [Module binaries] for details.

## Multiple inclusions

A Nextflow script can include any number of modules, and an `include` statement can import any number of definitions from a module. Multiple definitions can be included from the same module by using the syntax shown below:

```nextflow
include { cat; wc } from './some/module'

workflow {
    data = channel.fromPath('/some/data/*.txt')
    cat(data)
    wc(data)
}
```

(module-aliases)=

## Module aliases

When including definition from a module, it's possible to specify an *alias* with the `as` keyword. Aliasing allows you to avoid module name clashes, by assigning them different names in the including context. For example:

```nextflow
include { cat as cat_alpha } from './some/module'
include { cat as cat_beta } from './other/module'

workflow {
    cat_alpha(some_data)
    cat_beta(other_data)
}
```

You can also include the same definition multiple times under different names:

```nextflow
include { cat as cat_alpha; cat as cat_beta } from './some/module'

workflow {
    cat_alpha(some_data)
    cat_beta(other_data)
}
```

(module-params)=

## Module parameters

:::{deprecated} 24.07.0-edge
As a best practice, parameters should be used in the entry workflow and passed to workflows, processes, and functions as explicit inputs.
:::

A module can define parameters using the same syntax as a Nextflow workflow script:

```nextflow
params.message = 'Hello'
params.target = 'world!'

def sayHello() {
    println "$params.message $params.target"
}
```

When including a module, the module will first use parameters from the including context. For example:

```nextflow
params.message = 'Hola'
params.target = 'Mundo'

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
params.message = 'Hola'
params.target = 'Mundo'

include { sayHello } from './some/module' addParams(message: 'Ciao')

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
params.message = 'Hola'
params.target = 'Mundo'

include { sayHello } from './some/module' params(message: 'Ciao')

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

The binary scripts must be placed in the module directory named `<module-dir>/resources/usr/bin`:

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

Modules are designed to be easy to share and re-use across different pipelines, which helps eliminate duplicate work and spread improvements throughout the community. There are several ways to share modules:

- Use the Nextflow module registry system (recommended, see below)
- Simply copy the module files into your pipeline repository
- Use [Git submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules) to fetch modules from other Git repositories without maintaining a separate copy
- Use the [nf-core](https://nf-co.re/tools#modules) CLI to install and update modules with a standard approach used by the nf-core community

(module-registry)=

## Registry-based modules

:::{versionadded} 26.04.0
:::

Nextflow provides a module registry system that enables you to install, share, and manage reusable modules from centralized registries. This system provides version management, integrity checking, and seamless integration with the Nextflow DSL.

### Installing modules from a registry

Use the `module install` command to download modules from a registry:

```console
$ nextflow module install nf-core/fastqc
$ nextflow module install nf-core/fastqc -version 1.0.0
```

Installed modules are stored in the `modules/` directory and can be included using the registry syntax with the `@` prefix:

```nextflow
include { FASTQC } from '@nf-core/fastqc'

workflow {
    reads = Channel.fromFilePairs('data/*_{1,2}.fastq.gz')
    FASTQC(reads)
}
```

### Running modules directly

For ad-hoc tasks or testing, you can run a module directly without creating a workflow:

```console
$ nextflow module run nf-core/fastqc --input 'data/*.fastq.gz'
```

This command accepts all standard Nextflow options (`-profile`, `-resume`, etc.) and automatically downloads the module if not already installed.

### Managing module versions

Module versions are tracked in `nextflow_spec.json` in your project directory:

```json
{
  "modules": {
    "@nf-core/fastqc": "1.0.0",
    "@nf-core/bwa-align": "1.2.0"
  }
}
```

You can also configure versions in `nextflow.config`:

```nextflow
modules {
    '@nf-core/fastqc' = '1.0.0'
    '@nf-core/bwa-align' = '1.2.0'
}
```

When you run your workflow, Nextflow automatically installs or updates modules to match the specified versions.

### Discovering modules

Search for available modules using the `module search` command:

```console
$ nextflow module search alignment
$ nextflow module search "quality control" -limit 10
```

List installed modules in your project:

```console
$ nextflow module list
```

### Module integrity protection

Nextflow automatically verifies module integrity using checksums. If you modify a module locally, Nextflow will detect the change and prevent accidental overwrites:

```console
$ nextflow module install nf-core/fastqc -version 1.1.0
Warning: Module @nf-core/fastqc has local modifications. Use -force to override.
```

Use the `-force` flag to override local modifications when needed.

### Removing modules

Use the `module remove` command to uninstall a module:

```console
$ nextflow module remove nf-core/fastqc
```

By default, both the local module files and the entry in `nextflow_spec.json` are removed. Use the flags below to control this behaviour:

- `-keep-files` — Remove the entry from `nextflow_spec.json` but keep the local module files
- `-keep-config` — Remove the local module files but keep the entry in `nextflow_spec.json`

### Viewing module information

Use the `module info` command to display metadata and a usage template for a module:

```console
$ nextflow module info nf-core/fastqc
$ nextflow module info nf-core/fastqc -version 1.0.0
```

The output includes the module description, authors, keywords, tools, inputs, outputs, and a ready-to-use command-line template. Use `-json` to get machine-readable output.

### Publishing modules

To share your own modules, use the `module publish` command:

```console
$ nextflow module publish myorg/my-module
```

The argument can be either a `scope/name` reference (for an already-installed module) or a local directory path containing the module files.

Your module directory must include:

- `main.nf` - The module entry point
- `meta.yaml` - Module metadata (name, description, version, etc.)
- `README.md` - Module documentation

Authentication is required for publishing and can be provided via the `NXF_REGISTRY_TOKEN` environment variable or in your configuration:

```nextflow
registry {
    apiKey = 'YOUR_REGISTRY_TOKEN'
}
```

Use `-dry-run` to validate your module structure without uploading:

```console
$ nextflow module publish myorg/my-module -dry-run
```

### Registry configuration

By default, Nextflow uses the public registry at `https://registry.nextflow.io`. You can configure alternative or additional registries:

```nextflow
registry {
    url = [
        'https://private.registry.myorg.com',
        'https://registry.nextflow.io'
    ]
    apiKey = '${MYORG_TOKEN}'
}
```

Registries are queried in the order specified until a module is found. The `apiKey` is used only for the primary (first) registry.

### Module directory structure

Registry modules follow a standard directory structure:

```
modules/
└── @scope/
    └── module-name/
        ├── main.nf          # Module entry point (required)
        ├── meta.yaml        # Module metadata (required for publishing)
        ├── README.md        # Documentation (required for publishing)
        ├── .checksum        # Integrity checksum (generated automatically)
        ├── templates/       # Optional: process templates
        └── resources/       # Optional: module binaries and resources
```

The `modules/` directory should be committed to your Git repository to ensure reproducibility.

See the {ref}`cli-page` documentation for complete details on all module commands.

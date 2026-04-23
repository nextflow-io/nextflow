(process-reference)=

# Process reference

This page lists the task properties, input/output methods, and directives available in {ref}`process <process-page>` definitions.

## Task properties

The following task properties are defined in the process body:

`task.attempt`
: The current task attempt.

`task.exitStatus`
: *Available only in `script:` and `shell:` blocks*
: The exit code returned by the task script.
: The exit code is only available after the task has been executed (e.g., the [errorStrategy](#errorstrategy) directive).

`task.hash`
: *Available only in `exec:` blocks*
: The task hash.

`task.index`
: The process-level task index.

`task.name`
: *Available only in `exec:` blocks*
: The task name.

`task.previousException`
: :::{versionadded} 24.10.0
  :::
: The exception reported by the previous task attempt.
: Since the exception is available after a failed task attempt, it can only be accessed when retrying a failed task execution, i.e., when `task.attempt` is greater than 1.

`task.previousTrace`
: :::{versionadded} 24.10.0
  :::
: The trace record associated with the previous task attempt.
: Since the trace record is available after a failed task attempt, it can only be accessed when retrying a failed task execution, i.e., when `task.attempt` is greater than 1. See {ref}`trace-report` for a list of available fields.

: :::{note}
  The trace fields `%cpu` and `%mem` can be accessed as `pcpu` and `pmem`, respectively.
  :::


`task.process`
: The name of the process that spawned the task.

`task.workDir`
: *Available only in `exec:` blocks*
: The unique directory path for the task.

:::{note}
[Directive values](#directives) for a task can be accessed via `task.<directive>`. See {ref}`task-directive-values` for more information.
:::

(process-reference-typed)=

## Inputs and outputs (typed)

:::{versionadded} 25.10.0
:::

:::{note}
Typed processes require the `nextflow.preview.types` feature flag to be enabled in every script that uses them. The syntax and behavior may change in future releases.
:::

### Stage directives

The following directives can be used in the `stage:` section of a typed process:

`env( name: String, String value )`
: Declares an environment variable with the specified name and value in the task environment.

`stageAs( value: Path, filePattern: String )`
: Stages a file into the task directory under the given alias.

`stageAs( value: Iterable<Path>, filePattern: String )`
: Stages a collection of files into the task directory under the given alias.

`stdin( value: String )`
: Stages the given value as the standard input (i.e., `stdin`) to the task script.

### Outputs

The following functions are available in the `output:` and `topic:` sections of a typed process:

`env( name: String ) -> String`
: Returns the value of an environment variable from the task environment.

`eval( command: String ) -> String`
: Returns the standard output of the specified command, which is executed in the task environment after the task script completes.

`file( pattern: String, [options] ) -> Path`
: Returns a file from the task environment that matches the specified pattern.

: Available options:

  `followLinks: Boolean`
  : When `true`, target files are returned in place of any matching symlink (default: `true`).

  `glob: Boolean`
  : When `true`, the file name is interpreted as a glob pattern (default: `true`).

  `hidden: Boolean`
  : When `true`, hidden files are included in the matching output files (default: `false`).

  `includeInputs: Boolean`
  : When `true` and the file name is a glob pattern, any input files matching the pattern are also included in the output (default: `false`).

  `maxDepth: Integer`
  : Maximum number of directory levels to visit (default: no limit).

  `optional: Boolean`
  : When `true`, the task will not fail if the given file is missing (default: `false`).

  `type: String`
  : Type of paths returned, either `file`, `dir` or `any` (default: `any`, or `file` if the given file name contains a double star (`**`)).

`files( pattern: String, [options] ) -> Set<Path>`
: Returns files from the task environment that match the given pattern.

: Supports the same options as `file()` (except for `optional`).

`stdout() -> String`
: Returns the standard output of the task script.

(process-reference-legacy)=

## Inputs and outputs (legacy)

(process-reference-inputs)=

### Inputs

`val( identifier )`

: Declare a variable input. The received value can be any type, and it will be made available to the process body (i.e. `script`, `shell`, `exec`) as a variable given by `identifier`.

`file( identifier | stageName )`

: :::{deprecated} 19.10.0
  Use `path` instead.
  :::

: Declare a file input. The received value can be any type, and it will be staged into the task directory. If the received value is not a file or collection of files, it is implicitly converted to a string and written to a file.

: The argument can be an identifier or string. If an identifier, the received value will be made available to the process body as a variable. If a string, the received value will be staged into the task directory under the given alias.

`path( identifier | stageName )`

: Declare a file input. The received value should be a file or collection of files and will be staged into the task directory.

: :::{tip}
  See {ref}`process-multiple-input-files` for more information about accepting collections of files.
  :::

: The argument can be an identifier or string. If an identifier, the received value will be made available to the process body as a variable. If a string, the received value will be staged into the task directory under the given alias.

: Available options:

  `arity`
  : :::{versionadded} 23.10.0
    :::
  : Specify the number of expected files. Can be a number, e.g. `'1'`, or a range, e.g. `'1..*'`. If a task receives an invalid number of files for this `path` input, it will fail.

  `name`
  : Specify how the file should be named in the task work directory. Can be a name or a pattern.

  `stageAs`
  : Alias of `name`.

`env( name )`

: Declare an environment variable input. The received value should be a string, and it will be exported to the task environment as an environment variable given by `name`.

`stdin`

: Declare a `stdin` input. The received value should be a string, and it will be provided as the standard input (i.e. `stdin`) to the task script. It should be declared only once for a process.

`tuple( arg1, arg2, ... )`

: Declare a tuple input. Each argument should be an input declaration such as `val`, `path`, `env`, or `stdin`.

: The received value should be a tuple with the same number of elements as the `tuple` declaration, and each received element should be compatible with the corresponding `tuple` argument. Each tuple element is treated the same way as if it were a standalone input.

(process-reference-outputs)=

### Outputs

`val( value )`

: Declare a variable output. The argument can be any value, and it can reference any output variables defined in the process body (i.e. variables declared without the `def` keyword).

`file( pattern )`

: :::{deprecated} 19.10.0
  Use `path` instead.
  :::

: Declare a file output. It receives the output files from the task environment that match the given pattern.

: Multiple patterns can be specified using the colon separator (`:`). The union of all files matched by each pattern will be collected.

`path( pattern, [options] )`

: Declare a file output. It receives the output files from the task environment that match the given pattern.

: Available options:

  `arity`
  : :::{versionadded} 23.10.0
    :::
  : Specify the number of expected files. Can be a number or a range. If a task produces an invalid number of files for this `path` output, it will fail.

  : If the arity is `1`, a single file will be emitted. Otherwise, a list will always be emitted, even if only one file is produced.

  : :::{warning}
    If the arity is not specified, a single file or list will be emitted based on whether a single file or multiple files are produced at runtime, resulting potentially in an output channel with a mixture of files and file collections.
    :::

  `followLinks`
  : When `true`, target files are returned in place of any matching symlink (default: `true`).

  `glob`
  : When `true`, the specified name is interpreted as a glob pattern (default: `true`).

  `hidden`
  : When `true`, hidden files are included in the matching output files (default: `false`).

  `includeInputs`
  : When `true` and the output path is a glob pattern, any input files matching the pattern are also included in the output (default: `false`).

  `maxDepth`
  : Maximum number of directory levels to visit (default: no limit).

  `type`
  : Type of paths returned, either `file`, `dir` or `any` (default: `any`, or `file` if the specified file name pattern contains a double star (`**`)).

`env( name )`

: Declare an environment variable output. It receives the value of the environment variable (given by `name`) from the task environment.

: :::{versionchanged} 24.04.0
  Prior to this version, if the environment variable contained multiple lines of output, the output would be compressed to a single line by converting newlines to spaces.
  :::

`stdout`

: Declare a `stdout` output. It receives the standard output of the task script.

`eval( command )`

: :::{versionadded} 24.04.0
  :::

: Declare an `eval` output. It receives the standard output of the given command, which is executed in the task environment after the task script.

: If the command fails, the task will also fail.

`tuple( arg1, arg2, ... )`

: Declare a tuple output. Each argument should be an output declaration such as `val`, `path`, `env`, `stdin`, or `eval`. Each tuple element is treated the same way as if it were a standalone output.

(process-additional-options)=

### Generic options

The following options are available for all process outputs:

`emit: <name>`

: Defines the name of the output channel.

`optional: true | false`

: When `true`, the task will not fail if the specified output is missing (default: `false`).

`topic: <name>`

: :::{versionadded} 25.04.0
  :::

: Send the output to a {ref}`topic channel <channel-topic>` with the given name.

(process-reference-directives)=

## Directives

(process-accelerator)=

### accelerator

The `accelerator` directive defines the number of hardware accelerators (e.g. GPUs) required by each task execution. For example:

```nextflow
process hello {
    accelerator 4, type: 'nvidia-tesla-k80'

    script:
    """
    your_gpu_enabled --command --line
    """
}
```

The above example requests 4 GPUs of type `nvidia-tesla-k80` for each task.

:::{note}
This directive is only used by certain executors. Refer to the {ref}`executor-page` page to see which executors support this directive.
:::

:::{note}
Additional options may be required to fully enable the use of accelerators. When using containers with GPUs, you must pass the GPU drivers through to the container. For Docker, this requires the option `--gpus all` in the `docker run` command. For Apptainer/Singularity, this requires the option `--nv`. The specific implementation details depend on the accelerator and container type being used.
:::

The following options are available:

`request: Integer`
: The number of requested accelerators.
: Specifying this directive with a number (e.g., `accelerator 4`) is equivalent to the `request` option (e.g., `accelerator request: 4`).

`type: String`
: The accelerator type.
: The meaning of this option depends on the target execution platform. See the platform-specific documentation for more information about the available accelerators:

  - [Google Cloud](https://cloud.google.com/compute/docs/gpus/)
  - [Kubernetes](https://kubernetes.io/docs/tasks/manage-gpus/scheduling-gpus/#clusters-containing-different-types-of-gpus)

: This option is not supported for AWS Batch. You can control the accelerator type indirectly through the allowed instance types in your Compute Environment. See the [AWS Batch FAQs](https://aws.amazon.com/batch/faqs/?#GPU_Scheduling_) for more information.

(process-afterscript)=

### afterScript

The `afterScript` directive executes a custom (Bash) snippet immediately *after* the main process has run. This may be useful to clean up your staging area.

When combined with the [container](#container) directive, the `afterScript` is executed outside the specified container. In other words, the `afterScript` is always executed in the host environment.

(process-arch)=

### arch

The `arch` directive defines the CPU architecture to build the software in use by the process' task. For example:

```nextflow
process blast {
    spack 'blast-plus@2.13.0'
    arch 'linux/x86_64', target: 'cascadelake'

    script:
    """
    blastp -query input_sequence -num_threads ${task.cpus}
    """
}
```

The example above declares that the CPU generic architecture is `linux/x86_64` (X86 64 bit), and more specifically that the microarchitecture is `cascadelake` (a specific generation of Intel CPUs).

This directive is currently used by the following Nextflow functionalities:

- by the [spack](#spack) directive, to build microarchitecture-optimized applications;
- by the {ref}`wave-page` service, to build containers for one of the generic families of CPU architectures (see below);
- by the `spack` strategy within {ref}`wave-page`, to optimize the container builds for specific CPU microarchitectures.

Allowed values for the `arch` directive are as follows, grouped by equivalent family (choices available for the sake of compatibility):
- X86 64 bit: `linux/x86_64`, `x86_64`, `linux/amd64`, `amd64`
- ARM 64 bit: `linux/aarch64`, `aarch64`, `linux/arm64`, `arm64`, `linux/arm64/v8`
- ARM 64 bit, older generation: `linux/arm64/v7`

Examples of values for the architecture `target` option are `cascadelake`, `icelake`, `zen2` and `zen3`. See the [Spack documentation](https://spack.readthedocs.io/en/latest/basic_usage.html#support-for-specific-microarchitectures) for the full and up-to-date list of meaningful targets.

(process-array)=

### array

:::{versionadded} 24.04.0
:::

The `array` directive submits tasks as *job arrays* for executors that support it.

A job array is a collection of jobs with the same resource requirements and the same script (parameterized by an index). Job arrays incur significantly less scheduling overhead compared to individual jobs, and as a result they are preferred by HPC schedulers where possible.

The directive should be specified with a given array size, along with an executor that supports job arrays. For example:

```nextflow
process hello {
    executor 'slurm'
    array 100

    script:
    """
    your_command --here
    """
}
```

Nextflow currently supports job arrays for the following executors:

- {ref}`awsbatch-executor`
- {ref}`google-batch-executor`
- {ref}`lsf-executor`
- {ref}`pbs-executor`
- {ref}`pbspro-executor`
- {ref}`sge-executor`
- {ref}`slurm-executor`

A process using job arrays collects tasks and submits each batch as a job array when it is ready. Any "leftover" tasks are submitted as a partial job array.

Once a job array is submitted, each "child" task is executed as an independent job. Any tasks that fail (and can be retried) are retried without interfering with the tasks that succeeded. Retried tasks are submitted individually rather than through a job array, in order to allow for the use of {ref}`dynamic resources <dynamic-task-resources>`.

The following directives must be uniform across all tasks in a process that uses job arrays, because these directives are specified once for the entire job array:

- {ref}`process-accelerator`
- {ref}`process-clusterOptions`
- {ref}`process-cpus`
- {ref}`process-disk`
- {ref}`process-machineType`
- {ref}`process-memory`
- {ref}`process-queue`
- {ref}`process-resourcelabels`
- {ref}`process-resourcelimits`
- {ref}`process-time`

For cloud-based executors like AWS Batch, or when using Fusion with any executor, the following additional directives must be uniform:

- {ref}`process-container`
- {ref}`process-containerOptions`

When using Wave, the following additional directives must be uniform:

- {ref}`process-conda`

(process-beforescript)=

### beforeScript

The `beforeScript` directive executes a custom (Bash) snippet *before* the main process script is run. This may be useful to initialize the underlying cluster environment or for other custom initialization.

For example:

```nextflow
process hello {
    beforeScript 'source /cluster/bin/setup'

    script:
    """
    echo 'hello'
    """
}
```

When the process is containerized (using the [container](#container) directive), the `beforeScript` is executed in the container only if the executor is *container-native* (e.g. cloud batch executors, Kubernetes). Otherwise, the `beforeScript` is executed outside the container.

(process-cache)=

### cache

The `cache` directive controls whether and how task executions are cached.

By default, cached task executions are re-used when the pipeline is launched with the {ref}`resume <getstarted-resume>` option. The `cache` directive can be used to disable caching for a specific process:

```nextflow
process hello {
    cache false

    // ...
}
```

See {ref}`cache-resume-page` for more information.

The following options are available:

`false`
: Disable caching.

`true` (default)
: Enable caching. Input file metadata (name, size, last updated timestamp) are included in the cache keys.

`'deep'`
: Enable caching. Input file content is included in the cache keys.

`'lenient'`
: Enable caching. Minimal input file metadata (name and size only) are included in the cache keys.
: This strategy provides a workaround for incorrect caching invalidation observed on shared file systems due to inconsistent file timestamps.

(process-clusteroptions)=

### clusterOptions

The `clusterOptions` directive specifies additional submission options for grid executors. You can use it to specify options for your cluster that are not supported directly by other process directives.

The cluster options can be a string:

```nextflow
process hello {
    clusterOptions '-x 1 -y 2'

    // ...
}
```

:::{versionchanged} 24.04.0
Prior to this version, grid executors that require each option to be on a separate line in the job script would attempt to split multiple options using a variety of different conventions. Multiple options can now be specified more clearly using a string list as shown below.
:::

The cluster options can also be a string list:

```nextflow
process hello {
    clusterOptions '-x 1', '-y 2', '--flag'

    // ...
}
```

Grid executors that require one option per line will write each option to a separate line, while grid executors that allow multiple options per line will write all options to a single line, the same as with a string. This form is useful to control how the options are split across lines when it is required by the scheduler.

:::{note}
This directive is only used by grid executors. Refer to the {ref}`executor-page` page to see which executors support this directive.
:::

:::{warning}
While you can use the `clusterOptions` directive to specify options that are supported as process directives (`queue`, `memory`, `time`, etc), you should not use both at the same time, as it will cause undefined behavior. Most HPC schedulers will either fail or simply ignore one or the other.
:::

(process-conda)=

### conda

The `conda` directive defines the set of [Conda](https://conda.io) packages required by each task. For example:

```nextflow
process hello {
    conda 'bwa=0.7.15'

    script:
    """
    your_command --here
    """
}
```

Nextflow automatically creates an environment for each unique set of Conda packages.

The name of the desired channel for a specific package can be specified using the standard Conda notation, e.g. `bioconda::bwa=0.7.15`. Multiple packages can be specified separating them with a blank space, e.g. `bwa=0.7.15 fastqc=0.11.5`.

The `conda` directive can also accept a Conda environment file path or the path of an existing Conda environment. See {ref}`conda-page` for more information.

(process-container)=

### container

The `container` directive defines the container required by each task. For example:

```nextflow
process hello_docker {
    container 'busybox:latest'

    script:
    """
    your_command --here
    """
}
```

The corresponding container runtime (e.g. Docker, Singularity) should be running on the compute nodes where tasks are executed. See {ref}`container-page` for the container runtimes supported by Nextflow.

:::{note}
This directive is ignored by {ref}`native processes <process-native>` (i.e. `exec` processes).
:::

(process-containeroptions)=

### containerOptions

The `containerOptions` directive specifies additional container options for the underlying container runtime (e.g. Docker, Singularity). For example:

```nextflow
process hello_docker {
    container 'busybox:latest'
    containerOptions '--volume /data/db:/db'

    output:
    path 'output.txt'

    script:
    """
    your_command --data /db > output.txt
    """
}
```

The above example provides a custom volume mount for a specific process.

:::{warning}
This directive is not supported by the {ref}`k8s-executor` executor.
:::

(process-cpus)=

### cpus

The `cpus` directive defines the number of CPUs required by each task execution. For example:

```nextflow
process blast {
    cpus 8

    script:
    """
    blastp -query input_sequence -num_threads ${task.cpus}
    """
}
```

This directive is required for tasks that execute multi-process or multi-threaded commands/tools and it is meant to reserve enough CPUs when a pipeline task is executed through a cluster resource manager.

See also: [disk](#disk), [memory](#memory), [time](#time), [queue](#queue), {ref}`dynamic-task-resources`

(process-debug)=

### debug

The `debug` directive prints the standard output of each task to the pipeline standard output.

For example:

```nextflow
process hello {
    debug true

    script:
    """
    echo Hello
    """
}
```

Prints:

```
Hello
```

Removing the `debug` directive or setting it to `false` in the above example will cause `Hello` to not be printed.

(process-disk)=

### disk

The `disk` directive defines the amount of disk storage required by each task execution. For example:

```nextflow
process hello {
    disk 2.GB

    script:
    """
    your_command --here
    """
}
```

The following suffixes can be used to specify disk values:

- `B`: Bytes
- `KB`: Kilobytes
- `MB`: Megabytes
- `GB`: Gigabytes
- `TB`: Terabytes

See {ref}`stdlib-types-memoryunit` for more information.

:::{note}
The `disk` directive is only used by certain executors. Refer to the {ref}`executor-page` page to see which executors support this directive.
:::

See also: [cpus](#cpus), [memory](#memory), [time](#time), [queue](#queue), {ref}`dynamic-task-resources`

(process-error-strategy)=

### errorStrategy

The `errorStrategy` directive defines how to handle task failures.

A task failure occurs when the executed script returns a non-zero exit code. By default, the pipeline run is aborted.

The following error strategies are available:

`'terminate'` (default)
: When a task fails, terminate the pipeline immediately and report an error. Pending and running jobs are killed.

`'finish'`
: When a task fails, wait for submitted and running tasks to finish and then terminate the pipeline, reporting an error.

`'ignore'`
: When a task fails, ignore it and continue the pipeline execution. If the `workflow.failOnIgnore` config option is set to `true`, the pipeline will report an error (i.e. return a non-zero exit code) upon completion. Otherwise, the pipeline will complete successfully.
: See the {ref}`stdlib-namespaces-workflow` namespace for more information.

`'retry'`
: When a task fails, retry it.

When setting the `errorStrategy` directive to `ignore` the process doesn't stop on an error condition, it just reports a message notifying you of the error event.

For example:

```nextflow
process hello {
    errorStrategy 'ignore'

    // ...
}
```

In this case, the workflow will complete successfully and return an exit status of 0. However, if you set `workflow.failOnIgnore = true` in your Nextflow configuration, the workflow will return a non-zero exit status and report the failed tasks as an error.

The `retry` error strategy retries failed tasks. For example:

```nextflow
process hello {
    errorStrategy 'retry'

    // ...
}
```

The number of times a failing process is re-executed is defined by the [maxRetries](#maxretries) and [maxErrors](#maxerrors) directives.

:::{tip}
More complex strategies depending on the task exit status or other parametric values can be defined using a dynamic `errorStrategy`. See {ref}`dynamic-directives` for details.
:::

See also: [maxErrors](#maxerrors), [maxRetries](#maxretries), {ref}`dynamic-task-resources`

(process-executor)=

### executor

The `executor` directive defines the underlying system where tasks are executed. For example:

```nextflow
process hello {
    executor 'slurm'

    // ...
}
```

Commonly used executors include:

- `awsbatch`: [AWS Batch](https://aws.amazon.com/batch/)
- `azurebatch`: [Azure Batch](https://azure.microsoft.com/en-us/services/batch/)
- `google-batch`: [Google Cloud Batch](https://cloud.google.com/batch)
- `k8s`: [Kubernetes](https://kubernetes.io/) cluster
- `local`: The local machine where the pipeline is launched
- `lsf`: [Platform LSF](http://en.wikipedia.org/wiki/Platform_LSF) job scheduler
- `slurm`: [SLURM](https://en.wikipedia.org/wiki/Slurm_Workload_Manager) workload manager

Each executor supports additional configuration options under the `executor` config scope. See {ref}`executor-page` for more information.

(process-ext)=

### ext

The `ext` is a generic directive for user-defined properties. For example:

```nextflow
process star {
    container "biocontainers/star:${task.ext.version}"

    input:
    path genome
    tuple val(sampleId), path(reads)

    script:
    """
    STAR --genomeDir $genome --readFilesIn $reads ${task.ext.args ?: ''}
    """
}
```

In the above example, the process container version is controlled by `ext.version`, and the script supports additional command line arguments through `ext.args`.

The `ext` directive can be set in the process definition:

```nextflow
process hello {
    ext version: '2.5.3', args: '--alpha --beta'

    // ...
}
```

Or in the Nextflow configuration:

```groovy
process.ext.version = '2.5.3'
process.ext.args = '--alpha --beta'
```

(process-fair)=

### fair

:::{versionadded} 23.04.0
:::

The `fair` directive, when enabled, guarantees that process outputs will be emitted in the order in which they were received. For example:

```nextflow
process hello {
    fair true

    input:
    val x

    output:
    tuple val(task.index), val(x)

    script:
    """
    sleep \$((RANDOM % 3))
    """
}

workflow {
    channel.of('A','B','C','D') | hello | view
}
```

The above example produces:

```
[1, A]
[2, B]
[3, C]
[4, D]
```

(process-hints)=

### hints

The `hints` directive specifies executor-specific hints as key-value pairs. Each executor uses the hints it recognizes and ignores the rest. Hint values must be strings.

Unprefixed keys are available to **every** executor -- any executor that recognizes the key consumes it. Prefixing a key with an executor name (e.g. `awsbatch/...`) restricts the hint to that executor only. For example:

```nextflow
process hello {
    hints consumableResources: 'my-license=1'

    script:
    """
    your_command --here
    """
}
```

To restrict a hint to a single executor, prefix the key with the executor name:

```nextflow
    hints 'awsbatch/consumableResources': 'my-license=1'
```

When the same hint is provided both unprefixed and with a matching executor prefix, the prefixed form takes precedence for that executor.

Calling `hints` multiple times in a process definition accumulates entries, with later calls overwriting entries for the same key. Setting `hints` via configuration (e.g. in `nextflow.config`) replaces the entire map.

See {ref}`executor-page` to see which hints are recognized by each executor.

(process-label)=

### label

The `label` directive attaches a custom label to the process. For example:

```nextflow
process hello {
    label 'big_mem'

    script:
    """
    your_command --here
    """
}
```

A label may contain alphanumeric characters or `_`. It must start and end with an alphabetic character.

The same label can be applied to multiple processes. Multiple labels can be applied to the same process by using the `label` directive multiple times.

Process labels are used to apply shared process configuration via `withLabel` selectors. They are not recorded in execution logs, trace reports, or lineage metadata. See {ref}`config-process-selectors` for more information.

:::{note}
To tag individual task executions for logging and debugging, use [tag](#tag). To tag cloud computing resources for cost tracking, use [resourceLabels](#resourcelabels). To attach metadata labels to output files for lineage tracking, use the `label` {ref}`output directive <workflow-output-def>` in the `output` block.
:::

(process-machinetype)=

### machineType

The `machineType` can be used to specify a predefined Google Compute Platform [machine type](https://cloud.google.com/compute/docs/machine-types) when running using the {ref}`Google Batch <google-batch-executor>`, or when using auto-pools with {ref}`Azure Batch <azurebatch-executor>`.

For example:

```nextflow
process hello {
    machineType 'n1-highmem-8'

    script:
    """
    your_command --here
    """
}
```

See also: [cpus](#cpus), [memory](#memory)

(process-maxerrors)=

### maxErrors

The `maxErrors` directive defines the maximum number of task failures allowed for a process when using the `retry` error strategy. For example:

```nextflow
process hello {
    errorStrategy 'retry'
    maxErrors 5

    script:
    """
    echo 'do this as that .. '
    """
}
```

In the above example, the run will fail if the `hello` process accrues more than 5 failures across all of its task executions.

By default, there is no limit. However, the run can still fail if an individual task exceeds the number of retries allowed by the `maxRetries` directive.

See also: [errorStrategy](#errorstrategy), [maxRetries](#maxretries)

(process-maxforks)=

### maxForks

The `maxForks` directive defines the maximum number of concurrent task executions for a process. For example:

```nextflow
process hello {
    maxForks 1

    script:
    """
    your_command --here
    """
}
```

The above example forces the `hello` process to execute tasks sequentially.

By default, there is no limit. However, the number of concurrent tasks can still be limited globally by the number of CPUs (for local tasks) and the `executor.queueSize` config option.

(process-maxretries)=

### maxRetries

The `maxRetries` directive defines the maximum number of times a task can be retried when using the `retry` error strategy. For example:

```nextflow
process hello {
    errorStrategy 'retry'
    maxRetries 3

    script:
    """
    echo 'do this as that .. '
    """
}
```

In the above example, the run will fail if any task executed by `hello` fails more than three times.

By default, only one retry per task is allowed. However, the run can still fail if the total number of failures for the process exceeds the number allowed by the `maxErrors` directive.

See also: [errorStrategy](#errorstrategy), [maxErrors](#maxerrors)

(process-maxsubmitawait)=

### maxSubmitAwait

The `maxSubmitAwait` directive defines how long a task can remain in submission queue without being executed. Tasks that exceed this duration in the queue will fail.

It can be used with the `retry` error strategy to re-submit tasks to a different queue or with different resource requirements. For example:

```nextflow
process hello {
    errorStrategy 'retry'
    maxSubmitAwait 10.m
    maxRetries 3
    queue "${task.submitAttempt==1 ? 'spot-compute' : 'on-demand-compute'}"

    script:
    """
    your_command --here
    """
}
```

In the above example, each task is submitted to the `spot-compute` queue on the first attempt (`task.submitAttempt==1`). If a task remains in the queue for more than 10 minutes, it fails and is re-submitted to the `on-demand-compute` queue.

(process-memory)=

### memory

The `memory` directive defines how much memory is required by each task execution. For example:

```nextflow
process hello {
    memory 2.GB

    script:
    """
    your_command --here
    """
}
```

The following suffixes can be used to specify memory values:

- `B`: Bytes
- `KB`: Kilobytes
- `MB`: Megabytes
- `GB`: Gigabytes
- `TB`: Terabytes

See {ref}`stdlib-types-memoryunit` for more information.

See also: [cpus](#cpus), [disk](#disk), [time](#time), [queue](#queue), {ref}`dynamic-task-resources`

(process-module)=

### module

The `module` directive defines the set of [Environment Modules](http://modules.sourceforge.net/) required by each task, if supported by your compute environment. For example:

```nextflow
process blast {
    module 'ncbi-blast/2.2.27'

    script:
    """
    blastp -query <etc..>
    """
}
```

Multiple modules can be specified using the `:` separator:

```nextflow
process blast {
    module 'ncbi-blast/2.2.27:t_coffee/10.0:clustalw/2.1'

    script:
    """
    blastp -query <etc..>
  """
}
```

(process-penv)=

### penv

The `penv` directive defines the parallel environment to use when submitting tasks to the {ref}`SGE <sge-executor>` resource manager. For example:

```nextflow
process blast {
    cpus 4
    penv 'smp'
    executor 'sge'

    script:
    """
    blastp -query input_sequence -num_threads ${task.cpus}
    """
}
```

Refer to your cluster documentation or your system administrator to determine whether this feature is supported in your environment.

(process-pod)=

### pod

The `pod` directive defines pod-specific settings, such as environment variables, secrets, and config maps, when using the {ref}`k8s-executor` executor.

For example:

```nextflow
process echo {
    pod env: 'MESSAGE', value: 'hello world'

    script:
    """
    echo $MESSAGE
    """
}
```

The above snippet defines an environment variable named `MESSAGE` whose value is `'hello world'`.

Pod settings can be specified in Nextflow configuration:

```groovy
// single setting
process.pod = [env: 'MESSAGE', value: 'hello world']

// multiple settings
process.pod = [
    [env: 'MESSAGE', value: 'hello world'],
    [secret: 'my-secret/key1', mountPath: '/etc/file.txt']
]
```

The following options are available:

`affinity: <config>`
: Specifies the pod [affinity](https://kubernetes.io/docs/concepts/scheduling-eviction/assign-pod-node/#affinity-and-anti-affinity) with the given configuration.

`annotation: '<name>', value: '<value>'`
: *Can be specified multiple times*
: Defines a pod [annotation](https://kubernetes.io/docs/concepts/overview/working-with-objects/annotations/) with the given name and value.

`automountServiceAccountToken: true | false`
: Specifies whether to [automount service account token](https://kubernetes.io/docs/tasks/configure-pod-container/configure-service-account/#opt-out-of-api-credential-automounting) into the pod (default: `true`).

`config: '<configMap>/<key>', mountPath: '</absolute/path>'`
: *Can be specified multiple times*
: Mounts a [ConfigMap](https://kubernetes.io/docs/tasks/configure-pod-container/configure-pod-configmap/) with name and optional key to the given path. If the key is omitted, the path is interpreted as a directory and all entries in the `ConfigMap` are exposed in that path.

`csi: '<config>', mountPath: '</absolute/path>'`
: *Can be specified multiple times*
: Mounts a [CSI ephemeral volume](https://kubernetes.io/docs/concepts/storage/ephemeral-volumes/#csi-ephemeral-volumes) with the given configuration to the given path.

`emptyDir: <config>, mountPath: '</absolute/path>'`
: *Can be specified multiple times*
: Mounts an [emptyDir](https://kubernetes.io/docs/concepts/storage/volumes/#emptydir) with the given configuration to the given path.

`env: '<name>', config: '<configMap>/<key>'`
: *Can be specified multiple times*
: Defines an environment variable whose value is defined by the given [ConfigMap](https://kubernetes.io/docs/tasks/configure-pod-container/configure-pod-configmap/) and key.

`env: '<name>', fieldPath: '<fieldPath>'`
: *Can be specified multiple times*
: Defines an environment variable whose value is defined by the given [field path](https://kubernetes.io/docs/tasks/inject-data-application/environment-variable-expose-pod-information/#use-pod-fields-as-values-for-environment-variables) value.

: For example, the following pod option:

  ```groovy
  process.pod = [env: 'MY_NODE_NAME', fieldPath: 'spec.nodeName']
  ```

  Maps to the following pod spec:

  ```yaml
  env:
    - name: MY_NODE_NAME
      valueFrom:
        fieldRef:
          fieldPath: spec.nodeName
  ```

`env: '<name>', secret: '<secret>/<key>'`
: *Can be specified multiple times*
: Defines an environment variable whose value is defined by the given [Secret](https://kubernetes.io/docs/concepts/configuration/secret/) and key.

`env: '<name>', value: '<value>'`
: *Can be specified multiple times*
: Defines an environment variable with the given name and value.

`hostPath: '/host/absolute/path', mountPath: '</pod/absolute/path>'`
: :::{versionadded} 23.10.0
  :::
: *Can be specified multiple times*
: Allows creating [hostPath](https://kubernetes.io/docs/concepts/storage/volumes/#hostpath) volume and access it with the specified `mountPath` in the pod.

`imagePullPolicy: 'IfNotPresent' | 'Always' | 'Never'`
: Specifies the [image pull policy](https://kubernetes.io/docs/concepts/containers/images/#image-pull-policy) used by the pod to pull the container image.

`imagePullSecret: '<name>'`
: Specifies the [image pull secret](https://kubernetes.io/docs/concepts/containers/images/#specifying-imagepullsecrets-on-a-pod) used to access a private container image registry.

`label: '<name>', value: '<value>'`
: *Can be specified multiple times*
: Defines a pod [label](https://kubernetes.io/docs/concepts/overview/working-with-objects/labels/) with the given name and value.

`nodeSelector: <config>`
: Specifies the [node selector](https://kubernetes.io/docs/concepts/scheduling-eviction/assign-pod-node/#nodeselector) with the given configuration.

: The configuration can be a map or a string:

  ```groovy
  // map
  process.pod = [nodeSelector: [disktype: 'ssd', cpu: 'intel']]

  // string
  process.pod = [nodeSelector: 'disktype=ssd,cpu=intel']
  ```

`priorityClassName: '<name>'`
: Specifies the [priority class name](https://kubernetes.io/docs/concepts/scheduling-eviction/pod-priority-preemption/) for pods.

`privileged: true | false`
: Specifies whether the pod should run as a *privileged* container (default: `false`).

`runAsUser: '<uid>'`
: Specifies the user ID with which to run the container. Shortcut for the `securityContext` option.

`runtimeClassName: '<name>'`
: Specifies the [runtime class](https://kubernetes.io/docs/concepts/containers/runtime-class/).

`schedulerName: '<name>'`
: Specifies which [scheduler](https://kubernetes.io/docs/tasks/extend-kubernetes/configure-multiple-schedulers/#specify-schedulers-for-pods) is used to schedule the container.

`secret: '<secret>/<key>', mountPath: '</absolute/path>'`
: *Can be specified multiple times*
: Mounts a [Secret](https://kubernetes.io/docs/concepts/configuration/secret/) with name and optional key to the given path. If the key is omitted, the path is interpreted as a directory and all entries in the `Secret` are exposed in that path.

`securityContext: <config>`
: Specifies the pod [security context](https://kubernetes.io/docs/tasks/configure-pod-container/security-context/) with the given configuration.

`toleration: <config>`
: *Can be specified multiple times*
: Specifies the pod [toleration](https://kubernetes.io/docs/concepts/scheduling-eviction/taint-and-toleration/) with the given configuration.

: The configuration should be a map corresponding to a single toleration rule. For example, the following pod options:

  ```groovy
  process.pod = [
      [toleration: [key: 'key1', operator: 'Equal', value: 'value1', effect: 'NoSchedule']],
      [toleration: [key: 'key1', operator: 'Exists', effect: 'NoSchedule']],
  ]
  ```

  Maps to the following pod spec:

  ```yaml
  tolerations:
    - key: "key1"
      operator: "Equal"
      value: "value1"
      effect: "NoSchedule"
    - key: "key1"
      operator: "Exists"
      effect: "NoSchedule"
  ```

`ttlSecondsAfterFinished`
: :::{versionadded} 24.04.0
  :::
: Specifies the [TTL mechanism](https://kubernetes.io/docs/concepts/workloads/controllers/job/#ttl-mechanism-for-finished-jobs) for finished jobs in seconds. Applies to both successful and failed jobs.

`volumeClaim: '<name>', mountPath: '</absolute/path>' [, subPath: '<path>', readOnly: true | false]`
: *Can be specified multiple times*
: Mounts a [Persistent volume claim](https://kubernetes.io/docs/concepts/storage/persistent-volumes/) with the given name to the given path.
: The `subPath` option can be used to mount a sub-directory of the volume instead of its root.
: The `readOnly` option can be used to mount the volume as read-only (default: `false`)

(process-publishdir)=

### publishDir

:::{note}
{ref}`Workflow outputs <workflow-output-def>` can be used instead of `publishDir`. See {ref}`migrating-workflow-outputs` to learn how to migrate existing code.
:::

The `publishDir` directive publishes matching process output files to a target directory. For example:

```nextflow
process hello {
    publishDir '/data/chunks'

    output:
    path 'chunk_*'

    script:
    """
    printf 'Hola' | split -b 1 - chunk_
    """
}
```

The above example publishes the `chunk_*` output files into the `/data/chunks` directory.

Only files that match the declaration in the `output` block are published, not all the outputs of the process.

The `publishDir` directive can be specified more than once in order to publish output files to different target directories based on different rules.

By default, files are published via *symbolic link* from the task directory to the target directory. Use the `mode` option to control this behavior:

```nextflow
process hello {
    publishDir '/data/chunks', mode: 'copy', overwrite: false

    output:
    path 'chunk_*'

    script:
    """
    printf 'Hola' | split -b 1 - chunk_
    """
}
```

:::{warning}
Output files are published *asynchronously* after the task execution, so they may not be immediately available in the publish directory during the pipeline run. Downstream processes should access output files through the declared process outputs, not the publish directory.
:::

Available options:

`contentType`
: :::{versionadded} 22.10.0
  :::
: *Experimental: currently only supported for S3.*
: Allow specifying the media content type of the published file a.k.a. [MIME type](https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/MIME_Types). If set to `true`, the content type is inferred from the file extension (default: `false`).

`enabled`
: Enable or disable the publish rule depending on the boolean value specified (default: `true`).

`failOnError`
: :::{versionchanged} 24.04.0
  The default value was changed from `false` to `true`
  :::
: When `true` abort the execution if some file can't be published to the specified target directory or bucket for any cause (default: `true`)

`mode`
: The file publishing method. Can be one of the following values:

  - `'copy'`: Copies the output files into the publish directory.
  - `'copyNoFollow'`: Copies the output files into the publish directory without following symlinks ie. copies the links themselves.
  - `'link'`: Creates a hard link in the publish directory for each output file.
  - `'move'`: Moves the output files into the publish directory. **Note**: this is only supposed to be used for a *terminal* process i.e. a process whose output is not consumed by any other downstream process.
  - `'rellink'`: Creates a relative symbolic link in the publish directory for each output file.
  - `'symlink'`: Creates an absolute symbolic link in the publish directory for each output file (default).

`overwrite`
: When `true` any existing file in the target directory will be overridden (default: `true` during normal pipeline execution and `false` when pipeline execution is `resumed`).

`path`
: Specifies the directory where files need to be published. **Note**: the syntax `publishDir '/some/dir'` is a shortcut for `publishDir path: '/some/dir'`.

`pattern`
: Specifies a [glob][glob] file pattern that selects which files to publish from the overall set of output files.

`saveAs`
: A closure which, given the name of the file being published, returns the actual file name or a full path where the file is required to be stored. This can be used to rename or change the destination directory of the published files dynamically by using a custom strategy. Return the value `null` from the closure to *not* publish a file. This is useful when the process has multiple output files, but you want to publish only some of them.

`storageClass`
: :::{versionadded} 23.04.0
  :::
: *Experimental: currently only supported for S3.*
: Allow specifying the storage class to be used for the published file.

`tags`
: *Experimental: currently only supported for S3.*
: Allow the association of arbitrary tags with the published file e.g. `tags: [MESSAGE: 'Hello world']`.

(process-queue)=

### queue

The `queue` directive defines the queue to which tasks should be submitted, for executors that support queues. For example:

```nextflow
process hello {
    queue 'long'
    executor 'slurm'

    script:
    """
    your_command --here
    """
}
```

Some executors can accept multiple queue names as a comma-separated string:

```nextflow
queue 'short,long,cn-el6'
```

However, this is not generally supported by cloud executors such as AWS Batch, Azure Batch, and Google Batch.

See {ref}`executor-page` to see which executors support this directive.

(process-resourcelabels)=

### resourceLabels

The `resourceLabels` directive attaches custom name-value pairs to task executions, for executors that support it. For example:

```nextflow
process hello {
    resourceLabels region: 'some-region', user: 'some-username'

    script:
    """
    your_command --here
    """
}
```

Resource labels are attached to underlying resources such as cloud VMs, and are intended for operational purposes such as cost tracking. They are not recorded in lineage metadata.

Resource labels are currently supported by the following executors:

- {ref}`awsbatch-executor`
- {ref}`azurebatch-executor`
- {ref}`google-batch-executor`
- {ref}`k8s-executor`
- {ref}`seqera-executor`

:::{note}
The limits and the syntax of the corresponding executor should be taken into consideration when using resource labels.
:::

:::{versionadded} 23.10.0
Resource labels in Azure are added to auto-pools, rather than jobs, in order to facilitate cost analysis. A new pool will be created for each new set of resource labels. Therefore, it is recommended to also set `azure.batch.deletePoolsOnCompletion = true` when using process-specific resource labels.
:::

See also: [label](#label) (for shared process configuration), [tag](#tag) (for per-task identification)

(process-resourcelimits)=

### resourceLimits

:::{versionadded} 24.04.0
:::

The `resourceLimits` directive defines environment-specific limits for task resource requests.

Resource limits can be specified in a process:

```nextflow
process hello {
    resourceLimits cpus: 24, memory: 768.GB, time: 72.h

    script:
    """
    your_command --here
    """
}
```

Or in the Nextflow configuration:

```nextflow
process.resourceLimits = [
    cpus: 24,
    memory: 768.GB,
    time: 72.h
]
```

Resource limits can be defined for the following directives:

- [cpus](#cpus)
- [disk](#disk)
- [memory](#memory)
- [time](#time)

When a task resource request exceeds the corresponding limit, the task resources are automatically reduced to comply with these limits before the job is submitted.

Resource limits are a useful way to prevent tasks with {ref}`dynamic resources <dynamic-task-resources>` from requesting more resources than can be provided by an executor (e.g. a task requests 32 cores but the largest node in the cluster has 24).

(process-scratch)=

### scratch

The `scratch` directive executes each task in a temporary directory that is local to the compute node.

This is useful when executing tasks on an executor with a shared filesystem, because it decreases the network overhead of reading and writing files. Only the files declared as process outputs are copied to the pipeline work directory.

For example:

```nextflow
process hello {
    scratch true

    output:
    path 'data_out'

    script:
    """
    your_command --here
    """
}
```

It can also be specified in the Nextflow configuration:

```groovy
process.scratch = true
```

By default, the `scratch` directive uses the `$TMPDIR` environment variable in the underlying node as the base scratch directory. If `$TMPDIR` is not defined, then it creates a scratch directory using the `mktemp` command.

Each task creates a subdirectory within the base scratch directory and automatically deletes it upon completion.

:::{note}
Cloud-based executors enable `scratch` by default since the pipeline work directory resides in object storage.
:::

The following values are supported:

`false`
: Do not use a scratch directory.

`true`
: Create a scratch directory in the directory defined by the `$TMPDIR` environment variable, or `$(mktemp /tmp)` if `$TMPDIR` is not set.

`'$YOUR_VAR'`
: Create a scratch directory in the directory defined by the given environment variable, or `$(mktemp /tmp)` if that variable is not set. The value must use single quotes, otherwise the environment variable will be evaluated in the pipeline script context.

`'/my/tmp/path'`
: Create a scratch directory in the specified directory.

`'ram-disk'`
: Create a scratch directory in the RAM disk `/dev/shm/`.

(process-secret)=

### secret

The `secret` directive allows a process to access secrets.

For example:

```nextflow
process hello_secret {
    secret 'MY_ACCESS_KEY'
    secret 'MY_SECRET_KEY'

    script:
    """
    your_command --access \$MY_ACCESS_KEY --secret \$MY_SECRET_KEY
    """
}
```

Each secret is provided to the task as an environment variable.

See {ref}`secrets-page` for more information.

:::{note}
Secrets can only be used with the local or grid executors (e.g., Slurm or Grid Engine). Secrets can be used with AWS Batch and Google Batch when launched from Seqera Platform.
:::

(process-directive-shell)=

### shell

The `shell` directive defines a custom shell command for process scripts. By default, script blocks are executed with `/bin/bash -ue`.

```nextflow
process hello {
    shell '/bin/bash', '-euo', 'pipefail'

    script:
    """
    your_command --here
    """
}
```

It can also be specified in the Nextflow configuration:

```groovy
process.shell = ['/bin/bash', '-euo', 'pipefail']
```

(process-spack)=

### spack

The `spack` directive defines the set of [Spack](https://spack.io) packages required by each task. For example:

```nextflow
process hello {
    spack 'bwa@0.7.15'

    script:
    """
    your_command --here
    """
}
```

Nextflow automatically creates a Spack environment for each unique set of packages.

Multiple packages can be specified separating them with a blank space, e.g. `bwa@0.7.15 fastqc@0.11.5`.

The `spack` directive also accepts a Spack environment file path or the path of an existing Spack environment. See {ref}`spack-page` for more information.

(process-stageinmode)=

### stageInMode

The `stageInMode` directive defines how input files are staged into the task work directory.

The following modes are supported:

`'copy'`
: Input files are staged in the task work directory by creating a copy.

`'link'`
: Input files are staged in the task work directory by creating a hard link for each of them.

`'rellink'`
: Input files are staged in the task work directory by creating a symbolic link with a relative path for each of them.

`'symlink'`
: Input files are staged in the task work directory by creating a symbolic link with an absolute path for each of them (default).

(process-stageoutmode)=

### stageOutMode

The `stageOutMode` directive defines how output files are staged out from the scratch directory to the task work directory.

The following modes are supported:

`'copy'`
: Output files are copied from the scratch directory to the work directory.

`'fcp'`
: :::{versionadded} 23.04.0
  :::
: Output files are copied from the scratch directory to the work directory by using the [fcp](https://github.com/Svetlitski/fcp) utility (note: it must be available in the task environment).

`'move'`
: Output files are moved from the scratch directory to the work directory.

`'rclone'`
: :::{versionadded} 23.04.0
  :::
: Output files are copied from the scratch directory to the work directory by using the [rclone](https://rclone.org) utility (note: it must be available in the task environment).

`'rsync'`
: Output files are copied from the scratch directory to the work directory by using the `rsync` utility.

See also: [scratch](#scratch)

(process-storedir)=

### storeDir

The `storeDir` directive stores task outputs in a permanent *store directory* instead of the work directory.

On subsequent runs, each task is executed only if the declared output files do not exist in the store directory. When the files are present, the task is skipped and these files are used as the task outputs.

The following example shows how to use the `storeDir` directive to create a directory containing a BLAST database for each species specified by an input parameter:

```nextflow
process make_blast_db {
    storeDir '/db/genomes'

    input:
    path species

    output:
    path "${dbName}.*"

    script:
    dbName = species.baseName
    """
    makeblastdb -dbtype nucl -in ${species} -out ${dbName}
    """
}
```

:::{warning}
If a process uses `storeDir` and all of its outputs are optional, the process will always be skipped, even if the store directory is empty. This issue can be avoided by specifying at least one required file output.
:::

:::{warning}
The `storeDir` directive should not be used to publish outputs. Use the [publishDir](#publishdir) directive or {ref}`workflow outputs <workflow-output-def>` instead.
:::

(process-tag)=

### tag

The `tag` directive defines a custom identifier for each task execution. For example:

```nextflow
process hello {
    tag "$code"

    input:
    val code

    script:
    """
    echo $code
    """
}

workflow {
    ch_codes = channel.of('alpha', 'gamma', 'omega')
    hello(ch_codes)
}
```

The above example logs each task with its corresponding tag:

```
[6e/28919b] Submitted process > hello (alpha)
[d2/1c6175] Submitted process > hello (gamma)
[1c/3ef220] Submitted process > hello (omega)
```

Tags are a useful way to track related tasks in a pipeline run. Tasks can be identified by tag in the {ref}`execution log <execution-log>` and the {ref}`trace report <trace-report>`.

:::{note}
The name of a task in both {ref}`reports <tracing-page>` and {ref}`lineage <data-lineage-page>` is defined as `<process> (<tag>)`.
:::

:::{note}
The `tag` directive is not related to the [label](#label) directive. Process labels are only used for shared process configuration, not for tracking.
:::

(process-time)=

### time

The `time` directive defines the maximum runtime for each task. For example:

```nextflow
process hello {
    time 1.h

    script:
    """
    your_command --here
    """
}
```

The following suffixes can be used to specify duration values:

- `ms`: milliseconds
- `s`: seconds
- `m`: minutes
- `h`: hours
- `d`: days

See {ref}`stdlib-types-duration` for more information.

:::{note}
This directive is only used by certain executors. Refer to the {ref}`executor-page` page to see which executors support this directive.
:::

See also: [cpus](#cpus), [disk](#disk), [memory](#memory), [queue](#queue), {ref}`dynamic-task-resources`

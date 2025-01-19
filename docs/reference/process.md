(process-reference)=

# Process reference

This page lists the task properties, input/output methods, and directives available in {ref}`process <process-page>` definitions.

## Task properties

The following task properties are defined in the process body:

`task.attempt`
: The current task attempt.

`task.exitStatus`
: The exit code of the task script. Only applicable for processes with a `script:` or `shell:` block.
: Since the exit code is only available after the task has been executed, it can only be used by certain process directives such as [errorStrategy](#errorstrategy).

`task.hash`
: *Available only in `exec:` blocks*
: The task unique hash ID.

`task.id`
: The pipeline-level task index. Corresponds to `task_id` in the {ref}`execution trace <trace-report>`.

`task.index`
: The process-level task index.

`task.name`
: *Available only in `exec:` blocks*
: The current task name.

`task.previousException`
: :::{versionadded} 24.10.0
  :::
: The exception reported by the previous task attempt.
: Since the exception is available after a failed task attempt,
  it can only be accessed when retrying a failed task execution, and therefore when `task.attempt` is greater than 1.

`task.previousTrace`
: :::{versionadded} 24.10.0
  :::
: The trace record associated with the previous task attempt.
: Since the trace record is available after a failed task attempt,
  it can only be accessed when retrying a failed task execution, and therefore when `task.attempt` is greater than 1.
: This is useful when retrying a task execution to access the previous task attempt runtime metrics e.g. used memory and CPUs.

`task.process`
: The current process name.

`task.workDir`
: *Available only in `exec:` blocks*
: The task unique directory.

Additionally, the [directive values](#directives) for the given task can be accessed via `task.<directive>`.

(process-reference-inputs)=

## Inputs

`val( identifier )`

: Declare a variable input. The received value can be any type, and it will be made available to the process body (i.e. `script`, `shell`, `exec`) as a variable given by `identifier`.

`file( identifier | stageName )`

: :::{deprecated} 19.10.0
  Use `path` instead.
  :::

: Declare a file input. The received value can be any type, and it will be staged into the task directory. If the received value is not a file or collection of files, it is implicitly converted to a string and written to a file.

: The argument can be an identifier or string. If an identifier, the received value will be made available to the process body as a variable. If a string, the received value will be staged into the task directory under the given alias.

`path( identifier | stageName )`

: Declare a file input. The received value can be any type, and it will be staged into the task directory. If the received value is not a file or collection of files, it is implicitly converted to a string and written to a file.

: The argument can be an identifier or string. If an identifier, the received value will be made available to the process body as a variable. If a string, the received value will be staged into the task directory under the given alias.

: Available options:

  `arity`
  : :::{versionadded} 23.09.0-edge
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

## Outputs

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
  : :::{versionadded} 23.09.0-edge
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

: :::{versionchanged} 23.12.0-edge
  Prior to this version, if the environment variable contained multiple lines of output, the output would be compressed to a single line by converting newlines to spaces.
  :::

`stdout`

: Declare a `stdout` output. It receives the standard output of the task script.

`eval( command )`

: :::{versionadded} 24.02.0-edge
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

: :::{versionadded} 23.11.0-edge
  :::

: *Experimental: may change in a future release.*

: Defines the {ref}`channel topic <channel-topic>` to which the output will be sent.

(process-reference-directives)=

## Directives

(process-accelerator)=

### accelerator

:::{versionadded} 19.09.0-edge
:::

The `accelerator` directive allows you to request hardware accelerators (e.g. GPUs) for the task execution. For example:

```nextflow
process foo {
    accelerator 4, type: 'nvidia-tesla-k80'

    script:
    """
    your_gpu_enabled --command --line
    """
}
```

The above examples will request 4 GPUs of type `nvidia-tesla-k80`.

:::{note}
This directive is only used by certain executors. Refer to the {ref}`executor-page` page to see which executors support this directive.
:::

:::{note}
Additional options may be required to fully enable the use of accelerators. When using containers with GPUs, you must pass the GPU drivers through to the container. For Docker, this requires the option `--gpus all` in the docker run command. For Apptainer/Singularity, this requires the option `--nv`. The specific implementation details depend on the accelerator and container type being used.
:::

:::{note}
The accelerator `type` option depends on the target execution platform. Refer to the platform-specific documentation for details on the available accelerators:

- [Google Cloud](https://cloud.google.com/compute/docs/gpus/)
- [Kubernetes](https://kubernetes.io/docs/tasks/manage-gpus/scheduling-gpus/#clusters-containing-different-types-of-gpus)

The accelerator `type` option is not supported for AWS Batch. You can control the accelerator type indirectly through the allowed instance types in your Compute Environment. See the [AWS Batch FAQs](https://aws.amazon.com/batch/faqs/?#GPU_Scheduling_) for more information.
:::

(process-afterscript)=

### afterScript

The `afterScript` directive allows you to execute a custom (Bash) snippet immediately *after* the main process has run. This may be useful to clean up your staging area.

:::{note}
When combined with the {ref}`container directive <process-container>`, the `afterScript` will be executed outside the specified container. In other words, the `afterScript` is always executed in the host environment.
:::

(process-arch)=

### arch

The `arch` directive allows you to define the CPU architecture to build the software in use by the process' task. For example:

```nextflow
process cpu_task {
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

- by the [spack](#spack) directive, to build microarchitecture-optimised applications;
- by the {ref}`wave-page` service, to build containers for one of the generic families of CPU architectures (see below);
- by the `spack` strategy within {ref}`wave-page`, to optimise the container builds for specific CPU microarchitectures.

Allowed values for the `arch` directive are as follows, grouped by equivalent family (choices available for the sake of compatibility):
- X86 64 bit: `linux/x86_64`, `x86_64`, `linux/amd64`, `amd64`
- ARM 64 bit: `linux/aarch64`, `aarch64`, `linux/arm64`, `arm64`, `linux/arm64/v8`
- ARM 64 bit, older generation: `linux/arm64/v7`

Examples of values for the architecture `target` option are `cascadelake`, `icelake`, `zen2` and `zen3`. See the [Spack documentation](https://spack.readthedocs.io/en/latest/basic_usage.html#support-for-specific-microarchitectures) for the full and up-to-date list of meaningful targets.

(process-array)=

### array

:::{versionadded} 24.04.0
:::

:::{warning} *Experimental: may change in a future release.*
:::

The `array` directive allows you to submit tasks as *job arrays* for executors that support it.

A job array is a collection of jobs with the same resource requirements and the same script (parameterized by an index). Job arrays incur significantly less scheduling overhead compared to individual jobs, and as a result they are preferred by HPC schedulers where possible.

The directive should be specified with a given array size, along with an executor that supports job arrays. For example:

```nextflow
process cpu_task {
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

A process using job arrays will collect tasks and submit each batch as a job array when it is ready. Any "leftover" tasks will be submitted as a partial job array.

Once a job array is submitted, each "child" task is executed as an independent job. Any tasks that fail (and can be retried) will be retried without interfering with the tasks that succeeded. Retried tasks are submitted individually rather than through a job array, in order to allow for the use of {ref}`dynamic resources <dynamic-task-resources>`.

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

(process-beforescript)=

### beforeScript

The `beforeScript` directive allows you to execute a custom (Bash) snippet *before* the main process script is run. This may be useful to initialise the underlying cluster environment or for other custom initialisation.

For example:

```nextflow
process foo {
  beforeScript 'source /cluster/bin/setup'

  script:
  """
  echo bar
  """
}
```

When the process is containerized (using the {ref}`process-container` directive), the `beforeScript` will be executed in the container only if the executor is *container-native* (e.g. cloud batch executors, Kubernetes). Otherwise, the `beforeScript` will be executed outside the container.

(process-cache)=

### cache

The `cache` directive allows you to store the process results to a local cache. When the cache is enabled *and* the pipeline is launched with the {ref}`resume <getstarted-resume>` option, any task executions that are already cached will be re-used. See the {ref}`cache-resume-page` page for more information about how the cache works.

The cache is enabled by default, but you can disable it for a specific process by setting the `cache` directive to `false`. For example:

```nextflow
process noCacheThis {
  cache false
  // ...
}
```

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

The `clusterOptions` directive allows the usage of any native configuration option accepted by your cluster submit command. You can use it to request non-standard resources or use settings that are specific to your cluster and not supported out of the box by Nextflow.

The cluster options can be a string:

```nextflow
process foo {
  clusterOptions '-x 1 -y 2'
  // ...
}
```

:::{versionchanged} 24.04.0
Prior to this version, grid executors that require each option to be on a separate line in the job script would attempt to split multiple options using a variety of different conventions. Multiple options can now be specified more clearly using a string list as shown below.
:::

The cluster options can also be a string list:

```nextflow
process foo {
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

The `conda` directive allows for the definition of the process dependencies using the [Conda](https://conda.io) package manager.

Nextflow automatically sets up an environment for the given package names listed by in the `conda` directive. For example:

```nextflow
process foo {
  conda 'bwa=0.7.15'

  script:
  """
  your_command --here
  """
}
```

Multiple packages can be specified separating them with a blank space e.g. `bwa=0.7.15 fastqc=0.11.5`. The name of the channel from where a specific package needs to be downloaded can be specified using the usual Conda notation i.e. prefixing the package with the channel name as shown here `bioconda::bwa=0.7.15`.

The `conda` directive also allows the specification of a Conda environment file path or the path of an existing environment directory. See the {ref}`conda-page` page for further details.

(process-container)=

### container

The `container` directive allows you to execute the process script in a [Docker](http://docker.io) container.

It requires the Docker daemon to be running in machine where the pipeline is executed, i.e. the local machine when using the *local* executor or the cluster nodes when the pipeline is deployed through a *grid* executor.

For example:

```nextflow
process runThisInDocker {
  container 'dockerbox:tag'

  script:
  """
  your_command --here
  """
}
```

Simply replace in the above script `dockerbox:tag` with the name of the Docker image you want to use.

:::{tip}
Containers are a very useful way to execute your scripts in a reproducible self-contained environment or to run your pipeline in the cloud.
:::

:::{note}
This directive is ignored for processes that are {ref}`executed natively <process-native>`.
:::

(process-containeroptions)=

### containerOptions

The `containerOptions` directive allows you to specify any container execution option supported by the underlying container engine (ie. Docker, Singularity, etc). This can be useful to provide container settings only for a specific process e.g. mount a custom path:

```nextflow
process runThisWithDocker {
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

:::{warning}
This feature is not supported by the {ref}`k8s-executor` and {ref}`google-lifesciences-executor` executors.
:::

(process-cpus)=

### cpus

The `cpus` directive allows you to define the number of (logical) CPU required by the process' task. For example:

```nextflow
process big_job {
  cpus 8
  executor 'sge'

  script:
  """
  blastp -query input_sequence -num_threads ${task.cpus}
  """
}
```

This directive is required for tasks that execute multi-process or multi-threaded commands/tools and it is meant to reserve enough CPUs when a pipeline task is executed through a cluster resource manager.

See also: [penv](#penv), [memory](#memory), [time](#time), [queue](#queue), [maxForks](#maxforks)

(process-debug)=

### debug

By default the `stdout` produced by the commands executed in all processes is ignored. By setting the `debug` directive to `true`, you can forward the process `stdout` to the current top running process `stdout` file, showing it in the shell terminal.

For example:

```nextflow
process sayHello {
  debug true

  script:
  """
  echo Hello
  """
}
```

```
Hello
```

Without specifying `debug true`, you won't see the `Hello` string printed out when executing the above example.

(process-disk)=

### disk

The `disk` directive allows you to define how much local disk storage the process is allowed to use. For example:

```nextflow
process big_job {
    disk '2 GB'
    executor 'cirrus'

    script:
    """
    your_command --here
    """
}
```

The following memory unit suffix can be used when specifying the disk value:

| Unit | Description |
| ---- | ----------- |
| B    | Bytes       |
| KB   | Kilobytes   |
| MB   | Megabytes   |
| GB   | Gigabytes   |
| TB   | Terabytes   |

See {ref}`stdlib-types-memoryunit` for more information.

:::{note}
This directive is only used by certain executors. Refer to the {ref}`executor-page` page to see which executors support this directive.
:::

See also: [cpus](#cpus), [memory](#memory) [time](#time), [queue](#queue) and {ref}`dynamic-task-resources`.

(process-echo)=

### echo

:::{deprecated} 22.04.0
Use `debug` instead
:::

(process-error-strategy)=

### errorStrategy

The `errorStrategy` directive allows you to define how the process manages an error condition. By default, when an error status is returned by the executed script (i.e. when it ends with a non-zero exit status), the process stops immediately, forcing the entire pipeline to terminate.

The following error strategies are available:

`terminate` (default)
: When a task fails, terminate the pipeline immediately and report an error. Pending and running jobs are killed.

`finish`
: When a task fails, wait for submitted and running tasks to finish and then terminate the pipeline, reporting an error.

`ignore`
: When a task fails, ignore it and continue the pipeline execution. If the `workflow.failOnIgnore` config option is set to `true`, the pipeline will report an error (i.e. return a non-zero exit code) upon completion. Otherwise, the pipeline will complete successfully.
: See {ref}`stdlib-constants` for more information on `workflow.failOnIgnore`.

`retry`
: When a task fails, retry it.

When setting the `errorStrategy` directive to `ignore` the process doesn't stop on an error condition, it just reports a message notifying you of the error event.

For example:

```nextflow
process ignoreAnyError {
  errorStrategy 'ignore'

  // ...
}
```

In this case, the workflow will complete successfully and return an exit status of 0. However, if you set `workflow.failOnIgnore = true` in your Nextflow configuration, the workflow will return a non-zero exit status and report the failed tasks as an error.

The `retry` error strategy allows you to re-submit for execution a process returning an error condition. For example:

```nextflow
process retryIfFail {
  errorStrategy 'retry'

  // ...
}
```

The number of times a failing process is re-executed is defined by the [maxRetries](#maxretries) and [maxErrors](#maxerrors) directives.

:::{tip}
More complex strategies depending on the task exit status or other parametric values can be defined using a dynamic `errorStrategy`. See {ref}`dynamic-directives` for details.
:::

See also: [maxErrors](#maxerrors), [maxRetries](#maxretries) and {ref}`dynamic-task-resources`.

(process-executor)=

### executor

The `executor` defines the underlying system where processes are executed. By default a process uses the executor defined globally in the `nextflow.config` file.

The `executor` directive allows you to configure what executor has to be used by the process, overriding the default configuration.

The following executors are available:

| Name                  | Executor                                                                                    |
| --------------------- | ------------------------------------------------------------------------------------------- |
| `awsbatch`            | [AWS Batch](https://aws.amazon.com/batch/) service                                          |
| `azurebatch`          | [Azure Batch](https://azure.microsoft.com/en-us/services/batch/) service                    |
| `condor`              | [HTCondor](https://research.cs.wisc.edu/htcondor/) job scheduler                            |
| `google-lifesciences` | [Google Genomics Pipelines](https://cloud.google.com/life-sciences) service                 |
| `k8s`                 | [Kubernetes](https://kubernetes.io/) cluster                                                |
| `local`               | The computer where `Nextflow` is launched                                                   |
| `lsf`                 | [Platform LSF](http://en.wikipedia.org/wiki/Platform_LSF) job scheduler                     |
| `moab`                | [Moab](http://www.adaptivecomputing.com/moab-hpc-basic-edition/) job scheduler              |
| `nqsii`               | [NQSII](https://www.rz.uni-kiel.de/en/our-portfolio/hiperf/nec-linux-cluster) job scheduler |
| `oge`                 | Alias for the `sge` executor                                                                |
| `pbs`                 | [PBS/Torque](http://en.wikipedia.org/wiki/Portable_Batch_System) job scheduler              |
| `pbspro`              | [PBS Pro](https://www.pbsworks.com/) job scheduler                                          |
| `sge`                 | Sun Grid Engine / [Open Grid Engine](http://gridscheduler.sourceforge.net/)                 |
| `slurm`               | [SLURM](https://en.wikipedia.org/wiki/Slurm_Workload_Manager) workload manager              |
| `uge`                 | Alias for the `sge` executor                                                                |

The following example shows how to set the process's executor:

```nextflow
process doSomething {
  executor 'sge'

  // ...
}
```

:::{note}
Each executor supports additional directives and `executor` configuration options. See {ref}`executor-page` for more information.
:::

(process-ext)=

### ext

The `ext` is a special directive used as *namespace* for user custom process directives. This can be useful for advanced configuration options. For example:

```nextflow
process mapping {
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
process mapping {
  ext version: '2.5.3', args: '--foo --bar'

  // ...
}
```

Or in the Nextflow configuration:

```groovy
process.ext.version = '2.5.3'
process.ext.args = '--foo --bar'
```

(process-fair)=

### fair

:::{versionadded} 22.12.0-edge
:::

The `fair` directive, when enabled, guarantees that process outputs will be emitted in the order in which they were received. For example:

```nextflow
process foo {
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
    channel.of('A','B','C','D') | foo | view
}
```

The above example produces:

```
[1, A]
[2, B]
[3, C]
[4, D]
```

(process-label)=

### label

The `label` directive allows the annotation of processes with mnemonic identifier of your choice. For example:

```nextflow
process bigTask {
  label 'big_mem'

  script:
  """
  your_command --here
  """
}
```

The same label can be applied to more than a process and multiple labels can be applied to the same process using the `label` directive more than one time.

:::{note}
A label must consist of alphanumeric characters or `_`, must start with an alphabetic character and must end with an alphanumeric character.
:::

Labels are useful to organize workflow processes in separate groups which can be referenced in the configuration file to select and configure a subset of processes having similar computing requirements. See {ref}`config-process-selectors` for more information.

See also: [resourceLabels](#resourcelabels)

(process-machinetype)=

### machineType

:::{versionadded} 19.07.0
:::

The `machineType` can be used to specify a predefined Google Compute Platform [machine type](https://cloud.google.com/compute/docs/machine-types) when running using the {ref}`Google Batch <google-batch-executor>` or {ref}`Google Life Sciences <google-lifesciences-executor>` executor, or when using the autopools feature of the {ref}`Azure Batch executor<azurebatch-executor>`.

This directive is optional and if specified overrides the cpus and memory directives:

```nextflow
process foo {
  machineType 'n1-highmem-8'

  script:
  """
  your_command --here
  """
}
```

See also: [cpus](#cpus) and [memory](#memory).

(process-maxsubmitawait)=

### maxSubmitAwait

The `maxSubmitAwait` directive allows you to specify how long a task can remain in submission queue without being executed.
Elapsed this time the task execution will fail.

When used along with `retry` error strategy, it can be useful to re-schedule the task to a difference queue or
resource requirement. For example:

```nextflow
process foo {
  errorStrategy 'retry'
  maxSubmitAwait '10 mins'
  maxRetries 3
  queue "${task.submitAttempt==1 : 'spot-compute' : 'on-demand-compute'}"

  script:
  """
  your_command --here
  """
}
```

In the above example the task is submitted to the `spot-compute` on the first attempt (`task.submitAttempt==1`). If the
task execution does not start in the 10 minutes, a failure is reported and a new submission is attempted using the
queue named `on-demand-compute`.

(process-maxerrors)=

### maxErrors

The `maxErrors` directive allows you to specify the maximum number of times a process can fail when using the `retry` error strategy. By default this directive is disabled, you can set it as shown in the example below:

```nextflow
process retryIfFail {
  errorStrategy 'retry'
  maxErrors 5

  script:
  """
  echo 'do this as that .. '
  """
}
```

:::{note}
This setting considers the **total** errors accumulated for a given process, across all instances. If you want to control the number of times a process **instance** (aka task) can fail, use `maxRetries`.
:::

See also: [errorStrategy](#errorstrategy) and [maxRetries](#maxretries).

(process-maxforks)=

### maxForks

The `maxForks` directive allows you to define the maximum number of process instances that can be executed in parallel. By default this value is equals to the number of CPU cores available minus 1.

If you want to execute a process in a sequential manner, set this directive to one. For example:

```nextflow
process doNotParallelizeIt {
  maxForks 1

  script:
  """
  your_command --here
  """
}
```

(process-maxretries)=

### maxRetries

The `maxRetries` directive allows you to define the maximum number of times a process instance can be re-submitted in case of failure. This value is applied only when using the `retry` error strategy. By default only one retry is allowed, you can increase this value as shown below:

```nextflow
process retryIfFail {
    errorStrategy 'retry'
    maxRetries 3

    script:
    """
    echo 'do this as that .. '
    """
}
```

:::{note}
There is a subtle but important difference between `maxRetries` and the `maxErrors` directive. The latter defines the total number of errors that are allowed during the process execution (the same process can launch different execution instances), while the `maxRetries` defines the maximum number of times the same process execution can be retried in case of an error.
:::

See also: [errorStrategy](#errorstrategy) and [maxErrors](#maxerrors).

(process-memory)=

### memory

The `memory` directive allows you to define how much memory the process is allowed to use. For example:

```nextflow
process big_job {
    memory '2 GB'
    executor 'sge'

    script:
    """
    your_command --here
    """
}
```

The following memory unit suffix can be used when specifying the memory value:

| Unit | Description |
| ---- | ----------- |
| B    | Bytes       |
| KB   | Kilobytes   |
| MB   | Megabytes   |
| GB   | Gigabytes   |
| TB   | Terabytes   |

See {ref}`stdlib-types-memoryunit` for more information.

See also: [cpus](#cpus), [time](#time), [queue](#queue) and {ref}`dynamic-task-resources`.

(process-module)=

### module

[Environment Modules](http://modules.sourceforge.net/) is a package manager that allows you to dynamically configure your execution environment and easily switch between multiple versions of the same software tool.

If it is available in your system you can use it with Nextflow in order to configure the processes execution environment in your pipeline.

In a process definition you can use the `module` directive to load a specific module version to be used in the process execution environment. For example:

```nextflow
process basicExample {
  module 'ncbi-blast/2.2.27'

  script:
  """
  blastp -query <etc..>
  """
}
```

You can repeat the `module` directive for each module you need to load. Alternatively multiple modules can be specified in a single `module` directive by separating all the module names by using a `:` (colon) character as shown below:

```nextflow
process manyModules {
  module 'ncbi-blast/2.2.27:t_coffee/10.0:clustalw/2.1'

  script:
  """
  blastp -query <etc..>
  """
}
```

(process-penv)=

### penv

The `penv` directive allows you to define the parallel environment to be used when submitting a parallel task to the {ref}`SGE <sge-executor>` resource manager. For example:

```nextflow
process big_job {
  cpus 4
  penv 'smp'
  executor 'sge'

  script:
  """
  blastp -query input_sequence -num_threads ${task.cpus}
  """
}
```

This configuration depends on the parallel environment provided by your grid engine installation. Refer to your cluster documentation or contact your admin to learn more about this.

See also: [cpus](#cpus), [memory](#memory), [time](#time)

(process-pod)=

### pod

The `pod` directive allows the definition of pod specific settings, such as environment variables, secrets, and config maps, when using the {ref}`k8s-executor` executor.

For example:

```nextflow
process your_task {
  pod env: 'FOO', value: 'bar'

  script:
  """
  echo $FOO
  """
}
```

The above snippet defines an environment variable named `FOO` whose value is `bar`.

When defined in the Nextflow configuration file, a pod setting can be defined as a map:

```groovy
process {
  pod = [env: 'FOO', value: 'bar']
}
```

Or as a list of maps:

```groovy
process {
  pod = [
    [env: 'FOO', value: 'bar'],
    [secret: 'my-secret/key1', mountPath: '/etc/file.txt']
  ]
}
```

The following options are available:

`affinity: <config>`
: :::{versionadded} 22.01.0-edge
  :::
: Specifies the pod [affinity](https://kubernetes.io/docs/concepts/scheduling-eviction/assign-pod-node/#affinity-and-anti-affinity) with the given configuration.

`annotation: '<name>', value: '<value>'`
: *Can be specified multiple times*
: Defines a pod [annotation](https://kubernetes.io/docs/concepts/overview/working-with-objects/annotations/) with the given name and value.

`automountServiceAccountToken: true | false`
: :::{versionadded} 22.01.0-edge
  :::
: Specifies whether to [automount service account token](https://kubernetes.io/docs/tasks/configure-pod-container/configure-service-account/#opt-out-of-api-credential-automounting) into the pod (default: `true`).

`config: '<configMap>/<key>', mountPath: '</absolute/path>'`
: *Can be specified multiple times*
: Mounts a [ConfigMap](https://kubernetes.io/docs/tasks/configure-pod-container/configure-pod-configmap/) with name and optional key to the given path. If the key is omitted, the path is interpreted as a directory and all entries in the `ConfigMap` are exposed in that path.

`csi: '<config>', mountPath: '</absolute/path>'`
: :::{versionadded} 22.11.0-edge
  :::
: *Can be specified multiple times*
: Mounts a [CSI ephemeral volume](https://kubernetes.io/docs/concepts/storage/ephemeral-volumes/#csi-ephemeral-volumes) with the given configuration to the given path.

`emptyDir: <config>, mountPath: '</absolute/path>'`
: :::{versionadded} 22.11.0-edge
  :::
: *Can be specified multiple times*
: Mounts an [emptyDir](https://kubernetes.io/docs/concepts/storage/volumes/#emptydir) with the given configuration to the given path.

`env: '<name>', config: '<configMap>/<key>'`
: *Can be specified multiple times*
: Defines an environment variable whose value is defined by the given [ConfigMap](https://kubernetes.io/docs/tasks/configure-pod-container/configure-pod-configmap/) and key.

`env: '<name>', fieldPath: '<fieldPath>'`
: :::{versionadded} 21.09.1-edge
  :::
: *Can be specified multiple times*
: Defines an environment variable whose value is defined by the given [field path](https://kubernetes.io/docs/tasks/inject-data-application/environment-variable-expose-pod-information/#use-pod-fields-as-values-for-environment-variables) value.

: For example, the following pod option:

  ```groovy
  pod = [env: 'MY_NODE_NAME', fieldPath: 'spec.nodeName']
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
  pod = [nodeSelector: [disktype: 'ssd', cpu: 'intel']]

  // string
  pod = [nodeSelector: 'disktype=ssd,cpu=intel']
  ```

`priorityClassName: '<name>'`
: :::{versionadded} 22.01.0-edge
  :::
: Specifies the [priority class name](https://kubernetes.io/docs/concepts/scheduling-eviction/pod-priority-preemption/) for pods.

`privileged: true | false`
: :::{versionadded} 22.05.0-edge
  :::
: Specifies whether the pod should run as a *privileged* container (default: `false`).

`runAsUser: '<uid>'`
: Specifies the user ID with which to run the container. Shortcut for the `securityContext` option.

`schedulerName: '<name>'`
: Specifies which [scheduler](https://kubernetes.io/docs/tasks/extend-kubernetes/configure-multiple-schedulers/#specify-schedulers-for-pods) is used to schedule the container.

`secret: '<secret>/<key>', mountPath: '</absolute/path>'`
: *Can be specified multiple times*
: Mounts a [Secret](https://kubernetes.io/docs/concepts/configuration/secret/) with name and optional key to the given path. If the key is omitted, the path is interpreted as a directory and all entries in the `Secret` are exposed in that path.

`securityContext: <config>`
: Specifies the pod [security context](https://kubernetes.io/docs/tasks/configure-pod-container/security-context/) with the given configuration.

`toleration: <config>`
: :::{versionadded} 22.04.0
  :::
: *Can be specified multiple times*
: Specifies the pod [toleration](https://kubernetes.io/docs/concepts/scheduling-eviction/taint-and-toleration/) with the given configuration.

: The configuration should be a map corresponding to a single toleration rule. For example, the following pod options:

  ```groovy
  pod = [
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
: :::{versionadded} 24.02.0-edge
  :::
: Specifies the [TTL mechanism](https://kubernetes.io/docs/concepts/workloads/controllers/job/#ttl-mechanism-for-finished-jobs) for finished jobs in seconds. Applies to both successful and failed jobs.

`volumeClaim: '<name>', mountPath: '</absolute/path>' [, subPath: '<path>', readOnly: true | false]`
: *Can be specified multiple times*
: Mounts a [Persistent volume claim](https://kubernetes.io/docs/concepts/storage/persistent-volumes/) with the given name to the given path.
: The `subPath` option can be used to mount a sub-directory of the volume instead of its root.
: The `readOnly` option can be used to mount the volume as read-only (default: `false`)

(process-publishdir)=

### publishDir

The `publishDir` directive allows you to publish the process output files to a specified folder. For example:

```nextflow
process foo {
    publishDir '/data/chunks'

    output:
    path 'chunk_*'

    script:
    """
    printf 'Hola' | split -b 1 - chunk_
    """
}
```

The above example splits the string `Hola` into file chunks of a single byte. When complete the `chunk_*` output files are published into the `/data/chunks` folder.

:::{note}
Only files that match the declaration in the `output` block are published, not all the outputs of the process.
:::

:::{tip}
The `publishDir` directive can be specified more than once in order to publish output files to different target directories based on different rules.
:::

By default files are published to the target folder creating a *symbolic link* for each process output that links the file produced into the process working directory. This behavior can be modified using the `mode` option, for example:

```nextflow
process foo {
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
Files are copied into the specified directory in an *asynchronous* manner, so they may not be immediately available in the publish directory at the end of the process execution. For this reason, downstream processes should not try to access output files through the publish directory, but through channels.
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
: :::{versionchanged} 24.03.0-edge
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
: When `true` any existing file in the specified folder will be overridden (default: `true` during normal pipeline execution and `false` when pipeline execution is `resumed`).

`path`
: Specifies the directory where files need to be published. **Note**: the syntax `publishDir '/some/dir'` is a shortcut for `publishDir path: '/some/dir'`.

`pattern`
: Specifies a [glob][glob] file pattern that selects which files to publish from the overall set of output files.

`saveAs`
: A closure which, given the name of the file being published, returns the actual file name or a full path where the file is required to be stored. This can be used to rename or change the destination directory of the published files dynamically by using a custom strategy. Return the value `null` from the closure to *not* publish a file. This is useful when the process has multiple output files, but you want to publish only some of them.

`storageClass`
: :::{versionadded} 22.12.0-edge
  :::
: *Experimental: currently only supported for S3.*
: Allow specifying the storage class to be used for the published file.

`tags`
: :::{versionadded} 21.12.0-edge
  :::
: *Experimental: currently only supported for S3.*
: Allow the association of arbitrary tags with the published file e.g. `tags: [FOO: 'Hello world']`.

(process-queue)=

### queue

The `queue` directive allows you to set the `queue` where jobs are scheduled when using a grid based executor in your pipeline. For example:

```nextflow
process grid_job {
    queue 'long'
    executor 'sge'

    script:
    """
    your_command --here
    """
}
```

:::{tip}
Some grid executors support multiple queue names as a comma-separated list:

```nextflow
queue 'short,long,cn-el6'
```

However, this does not generally apply to other executors such as AWS Batch, Azure Batch, Google Batch.
:::

:::{note}
This directive is only used by certain executors. Refer to the {ref}`executor-page` page to see which executors support this directive.
:::

(process-resourcelabels)=

### resourceLabels

:::{versionadded} 22.09.1-edge
:::

The `resourceLabels` directive allows you to specify custom name-value pairs that Nextflow applies to the computing resource used to carry out the process execution. Resource labels can be specified using the syntax shown below:

```nextflow
process my_task {
    resourceLabels region: 'some-region', user: 'some-username'

    script:
    """
    your_command --here
    """
}
```

The limits and the syntax of the corresponding cloud provider should be taken into consideration when using resource labels.

Resource labels are currently supported by the following executors:

- {ref}`awsbatch-executor`
- {ref}`azurebatch-executor`
- {ref}`google-batch-executor`
- {ref}`google-lifesciences-executor`
- {ref}`k8s-executor`

:::{versionadded} 23.09.0-edge
Resource labels are supported for Azure Batch when using automatic pool creation.

Resource labels in Azure are added to pools, rather than jobs, in order to facilitate cost analysis. A new pool will be created for each new set of resource labels, therefore it is recommended to also set `azure.batch.deletePoolsOnCompletion = true` when using process-specific resource labels.
:::

See also: [label](#label)

(process-resourcelimits)=

### resourceLimits

:::{versionadded} 24.04.0
:::

The `resourceLimits` directive allows you to specify environment-specific limits for task resource requests. Resource limits can be specified in a process as follows:

```nextflow
process my_task {
  resourceLimits cpus: 24, memory: 768.GB, time: 72.h

  script:
  """
  your_command --here
  """
}
```

Or in the Nextflow configuration:

```nextflow
process {
    resourceLimits = [ cpus: 24, memory: 768.GB, time: 72.h ]
}
```

Resource limits can be defined for the following directives:

- [cpus](#cpus)
- [disk](#disk)
- [memory](#memory)
- [time](#time)

Resource limits are a useful way to specify environment-specific limits alongside tasks with {ref}`dynamic resources <dynamic-task-resources>`. Normally, if a task requests more resources than can be provisioned (e.g. a task requests 32 cores but the largest node in the cluster has 24), the task will either fail or cause the pipeline to hang forever as it will never be scheduled. If the `resourceLimits` directive is defined with these limits, the task resources will be automatically reduced to comply with these limits before the job is submitted.

(process-scratch)=

### scratch

The `scratch` directive allows you to execute the process in a temporary folder that is local to the execution node.

This is useful when your pipeline is launched by using a grid executor, because it allows you to decrease the NFS overhead by running the pipeline processes in a temporary directory in the local disk of the actual execution node. Only the files declared as output in the process definition will be copied in the pipeline working area.

In its basic form simply specify `true` at the directive value, as shown below:

```nextflow
process simpleTask {
  scratch true

  output:
  path 'data_out'

  script:
  """
  your_command --here
  """
}
```

By doing this, it tries to execute the script in the directory defined by the variable `$TMPDIR` in the execution node. If this variable does not exist, it will create a new temporary directory by using the Linux command `mktemp`.

:::{note}
Cloud-based executors use `scratch = true` by default, since the work directory resides in object storage.
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

(process-directive-shell)=

### shell

The `shell` directive allows you to define a custom shell command for process scripts. By default, script blocks are executed with `/bin/bash -ue`.

```nextflow
process doMoreThings {
    shell '/bin/bash', '-euo', 'pipefail'

    script:
    """
    your_command --here
    """
}
```

The same directive could be specified in your Nextflow configuration as follows:

```groovy
process.shell = ['/bin/bash', '-euo', 'pipefail']
```

(process-spack)=

### spack

The `spack` directive allows for the definition of the process dependencies using the [Spack](https://spack.io) package manager.

Nextflow automatically sets up an environment for the given package names listed by in the `spack` directive. For example:

```nextflow
process foo {
    spack 'bwa@0.7.15'

    script:
    """
    your_command --here
    """
}
```

Multiple packages can be specified separating them with a blank space, e.g. `bwa@0.7.15 fastqc@0.11.5`.

The `spack` directive also allows the specification of a Spack environment file path or the path of an existing environment directory. See the {ref}`spack-page` page for further details.

(process-stageinmode)=

### stageInMode

The `stageInMode` directive defines how input files are staged into the process work directory. The following values are allowed:

`'copy'`
: Input files are staged in the process work directory by creating a copy.

`'link'`
: Input files are staged in the process work directory by creating a hard link for each of them.

`'rellink'`
: Input files are staged in the process work directory by creating a symbolic link with a relative path for each of them.

`'symlink'`
: Input files are staged in the process work directory by creating a symbolic link with an absolute path for each of them (default).

(process-stageoutmode)=

### stageOutMode

The `stageOutMode` directive defines how output files are staged out from the scratch directory to the process work directory. The following values are allowed:

`'copy'`
: Output files are copied from the scratch directory to the work directory.

`'fcp'`
: :::{versionadded} 23.02.0-edge
  :::
: Output files are copied from the scratch directory to the work directory by using the [fcp](https://github.com/Svetlitski/fcp) utility (note: it must be available in your cluster computing nodes).

`'move'`
: Output files are moved from the scratch directory to the work directory.

`'rclone'`
: :::{versionadded} 23.01.0-edge
  :::
: Output files are copied from the scratch directory to the work directory by using the [rclone](https://rclone.org) utility (note: it must be available in your cluster computing nodes).

`'rsync'`
: Output files are copied from the scratch directory to the work directory by using the `rsync` utility.

See also: [scratch](#scratch).

(process-storedir)=

### storeDir

The `storeDir` directive allows you to define a directory that is used as a *permanent* cache for your process results.

In more detail, it affects the process execution in two main ways:

1. The process is executed only if the files declared in the `output` block do not exist in the directory specified by the `storeDir` directive. When the files exist the process execution is skipped and these files are used as the actual process result.
2. Whenever a process is successfully completed the files listed in the `output` block are moved into the directory specified by the `storeDir` directive.

The following example shows how to use the `storeDir` directive to create a directory containing a BLAST database for each species specified by an input parameter:

```nextflow
process formatBlastDatabases {
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
The `storeDir` directive should not be used to publish workflow outputs. Use the [publishDir](#publishdir) directive or the {ref}`workflow output definition <workflow-output-def>` instead.
:::

(process-tag)=

### tag

The `tag` directive allows you to associate each process execution with a custom label, so that it will be easier to identify them in the log file or in the trace execution report. For example:

```nextflow
process foo {
  tag "$code"

  input:
  val code

  script:
  """
  echo $code
  """
}

workflow {
  Channel.of('alpha', 'gamma', 'omega') | foo
}
```

The above snippet will print a log similar to the following one, where process names contain the tag value:

```
[6e/28919b] Submitted process > foo (alpha)
[d2/1c6175] Submitted process > foo (gamma)
[1c/3ef220] Submitted process > foo (omega)
```

See also {ref}`Trace execution report <trace-report>`

(process-time)=

### time

The `time` directive allows you to define how long a process is allowed to run. For example:

```nextflow
process big_job {
    time '1h'

    script:
    """
    your_command --here
    """
}
```

The following time unit suffixes can be used when specifying the duration value:

| Unit                            | Description  |
| ------------------------------- | ------------ |
| `ms`, `milli`, `millis`         | Milliseconds |
| `s`, `sec`, `second`, `seconds` | Seconds      |
| `m`, `min`, `minute`, `minutes` | Minutes      |
| `h`, `hour`, `hours`            | Hours        |
| `d`, `day`, `days`              | Days         |

Multiple units can be used in a single declaration, for example: `'1day 6hours 3minutes 30seconds'`

See {ref}`stdlib-types-duration` for more information.

:::{note}
This directive is only used by certain executors. Refer to the {ref}`executor-page` page to see which executors support this directive.
:::

See also: [cpus](#cpus), [memory](#memory), [queue](#queue) and {ref}`dynamic-task-resources`.

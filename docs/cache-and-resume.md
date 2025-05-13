(cache-resume-page)=

# Caching and resuming

One of the core features of Nextflow is the ability to cache task executions and re-use them in subsequent runs to minimize duplicate work. Resumability is useful both for recovering from errors and for iteratively developing a pipeline. It is similar to [checkpointing](https://en.wikipedia.org/wiki/Application_checkpointing), a common practice used by HPC applications.

You can enable resumability in Nextflow with the `-resume` flag when launching a pipeline with `nextflow run`. In most cases, that is all you need to do and resumability will "just work". This page describes Nextflow's caching behavior in more detail in order to help advanced users understand how the cache works and troubleshoot it when it doesn't work.

## Task cache

All task executions are automatically saved to the task cache, regardless of the `-resume` option (so that you always have the option to resume later). The task cache is a key-value store, where each key-value pair corresponds to a previously-executed task.

The task cache is used in conjunction with the [work directory](#work-directory) to recover cached tasks in a resumed run. It is also used by the {ref}`cli-log` sub-command to query task metadata.

(cache-resume-task-hash)=

### Task hash

The task hash is computed from the following metadata:

- Session ID (see `workflow.sessionId` in the {ref}`stdlib-namespaces-workflow` namespace)
- Task name (see `name` in {ref}`trace-report`)
- Task container image (if applicable)
- Task {ref}`environment modules <process-module>` (if applicable)
- Task {ref}`Conda environment <process-conda>` (if applicable)
- Task {ref}`Spack environment <process-spack>` and {ref}`CPU architecture <process-arch>` (if applicable)
- Task {ref}`inputs <process-input>`
- Task {ref}`script <process-script>`
- Any global variables referenced in the task script
- Any task {ref}`process-ext` properties referenced in the task script
- Any {ref}`bundled scripts <bundling-executables>` used in the task script
- Whether the task is a {ref}`stub run <process-stub>`
- Task attempt

:::{versionchanged} 23.09.2-edge
The {ref}`process-ext` directive was added to the task hash.
:::

Nextflow computes this hash for every task when it is created but before it is executed. If resumability is enabled and there is an entry in the task cache with the same hash, Nextflow tries to recover the previous task execution. A cache hit does not guarantee that the task will be resumed, because it must also recover the task outputs from the [work directory](#work-directory).

Note that files are hashed differently depending on the caching mode. See the {ref}`process-cache` directive for more details.

### Task entry

The task entry is a serialized blob of the task metadata required to resume a task, including the fields used by the {ref}`trace-report` and the task input variables.

### Cache stores

The default cache store uses the `.nextflow/cache` directory, relative to the launch directory (i.e. `workflow.launchDir`), to store the task cache, with a separate subdirectory for each session ID backed by [LevelDB](https://github.com/dain/leveldb).

Due to the limitations of LevelDB, the database for a given session ID can only be accessed by one reader/writer at a time. This means, for example, that you cannot use `nextflow log` to query the task metadata for a pipeline run while it is still running.

:::{versionadded} 23.07.0-edge
:::

The cloud cache is an alternative cache store that uses cloud storage instead of the local cache directory. You can use it by setting the `NXF_CLOUDCACHE_PATH` environment variable to the desired cache path (e.g. `s3://my-bucket/cache`) and providing the necessary credentials.

The cloud cache is particularly useful when launching Nextflow from within the cloud, where the default cache would be lost once the pipeline completes and the VM instance is terminated. Furthermore, because it is backed by cloud storage, it can support multiple readers and writers.

## Work directory

While the [task cache](#task-cache) stores the task metadata for subsequent runs, the work directory stores various files used during a pipeline run.

Each task uses a unique directory based on its hash. When a task is created, Nextflow stages the task input files, script, and other helper files into the task directory. The task writes any output files to this directory during its execution, and Nextflow uses these output files for downstream tasks and/or publishing.

When a previous task is retrieved from the task cache on a resumed run, Nextflow then checks the corresponding task directory in the work directory. If all the required outputs are present and the exit code is valid, then the task is successfully cached; otherwise, the task is re-executed.

For this reason, it is important to preserve both the task cache (`.nextflow/cache`) and work directories in order to resume runs successfully. You can use the {ref}`cli-clean` command to delete specific runs from the cache.

## Troubleshooting

Cache failures occur when a task that was supposed to be cached was re-executed or a task that was supposed to be re-executed was cached. This page provides an overview of common causes for cache failures and strategies to identify them.

Common causes of cache failures include:

- {ref}`Resume not being enabled <cache-failure-resume>`
- {ref}`Non-default cache directives <cache-failure-directives>`
- {ref}`Modified inputs <cache-failure-modified>`
- {ref}`Inconsistent file attributes <cache-failure-inconsistent>`
- {ref}`Race condition on a global variable <cache-global-var-race-condition>`
- {ref}`Non-deterministic process inputs <cache-nondeterministic-inputs>`

The causes of these cache failure and solutions to resolve them are described in detail below.

(cache-failure-resume)=

### Resume not enabled

The `-resume` option is required to resume a pipeline. Ensure `-resume` has been enabled in your run command or your nextflow configuration file.

(cache-failure-directives)=

### Non-default cache directives

The `cache` directive is enabled by default. However, you can disable or modify its behavior for a specific process. For example:

```nextflow
process FOO {
  cache false
  // ...
}
```

Ensure that the cache has not been set to a non-default value. See {ref}`process-cache` for more information about the `cache` directive.

(cache-failure-modified)=

### Modified inputs

Modifying inputs that are used in the task hash will invalidate the cache. Common causes of modified inputs include:

- Changing input files
- Resuming from a different session ID
- Changing the process name
- Changing the task container image or Conda environment
- Changing the task script
- Changing a bundled script used by the task

:::{note}
Changing the value of any directive, except {ref}`process-ext`, will not inactivate the task cache.
:::

A hash for an input file is calculated from the complete file path, the last modified timestamp, and the file size to calculate. If any of these attributes change the task will be re-executed. If a process modifies its input files it cannot be resumed. Processes that modify their own input files are considered to be an anti-pattern and should be avoided.

(cache-failure-inconsistent)=

### Inconsistent file attributes

Some shared file systems, such as NFS, may report inconsistent file timestamps. If you encounter this problem, use the `'lenient'` {ref}`caching mode <process-cache>` to ignore the last modified timestamp and only use the file path.

(cache-global-var-race-condition)=

### Race condition on a global variable

Race conditions can in disrupt caching behavior of your pipeline. For example:

```nextflow
channel.of(1,2,3) | map { v -> X=v; X+=2 } | view { v -> "ch1 = $v" }
channel.of(1,2,3) | map { v -> X=v; X*=2 } | view { v -> "ch2 = $v" }
```

In the above example, `X` is declared in each `map` closure. Without the `def` keyword, or other type qualifier, the variable `X` is global to the entire script. Operators and executed concurrently and, as `X` is global, there is a *race condition* that causes the emitted values to vary depending on the order of the concurrent operations. If these values were passed to a process as inputs the process would execute different tasks during each run due to the race condition.

To resolve this failure type, ensure the variable is not global by using a local variable:

```nextflow
channel.of(1,2,3) | map { v -> def X=v; X+=2 } | view { v -> "ch1 = $v" }
```

Alternatively, remove the variable:

```nextflow
channel.of(1,2,3) | map { v -> v * 2 } | view { v -> "ch2 = $v" }
```

(cache-nondeterministic-inputs)=

### Non-deterministic process inputs

A process that merges inputs from different sources non-deterministically may invalidate the cache. For example:

```nextflow
workflow {
    ch_foo = channel.of( ['1', '1.foo'], ['2', '2.foo'] )
    ch_bar = channel.of( ['2', '2.bar'], ['1', '1.bar'] )
    gather(ch_foo, ch_bar)
}
process gather {
    input:
    tuple val(id), file(foo)
    tuple val(id), file(bar)
    script:
    """
    merge_command $foo $bar
    """
}
```

In the above example, the inputs will be merged without matching. This is the same way method used by the {ref}`operator-merge` operator. When merged, the inputs are incorrect, non-deterministic, and invalidate the cache.

To resolve this failure type, ensure channels are deterministic by joining them before invoking the process:

```nextflow
workflow {
    ch_foo = channel.of( ['1', '1.foo'], ['2', '2.foo'] )
    ch_bar = channel.of( ['2', '2.bar'], ['1', '1.bar'] )
    gather(ch_foo.join(ch_bar))
}
process gather {
    input:
    tuple val(id), file(foo), file(bar)
    script:
    """
    merge_command $foo $bar
    """
}
```

(cache-compare-hashes)=

## Tips

### Resume from a specific run

Nextflow resumes from the previous run by default. If you want to resume from an earlier run, simply specify the session ID for that run with the `-resume` option:

```bash
nextflow run rnaseq-nf -resume 4dc656d2-c410-44c8-bc32-7dd0ea87bebf
```

You can use the {ref}`cli-log` command to view all previous runs as well as the task executions for each run.

### Compare task hashes

By identifying differences between hashes you can detect changes that may be causing cache failures.

To compare the task hashes for a resumed run:

1. Run your pipeline with the `-log` and `-dump-hashes` options:

    ```bash
    nextflow -log run_initial.log run <PIPELINE> -dump-hashes
    ```

2. Run your pipeline with the `-log`, `-dump-hashes`, and `-resume` options:

    ```bash
    nextflow -log run_resumed.log run <PIPELINE> -dump-hashes -resume
    ```

3. Extract the task hash lines from each log:

    ```bash
    cat run_initial.log | grep 'INFO.*TaskProcessor.*cache hash' | cut -d ' ' -f 10- | sort | awk '{ print; print ""; }' > run_initial.tasks.log
    cat run_resumed.log | grep 'INFO.*TaskProcessor.*cache hash' | cut -d ' ' -f 10- | sort | awk '{ print; print ""; }' > run_resumed.tasks.log
    ```

4. Compare the runs:

    ```bash
    diff run_initial.tasks.log run_resumed.tasks.log
    ```

    :::{tip}
    You can also compare the hash lines using a graphical diff viewer.
    :::

:::{versionadded} 23.10.0
:::

Task hashes can also be extracted into a diff using `-dump-hashes json`. The following is an example Bash script to compare two runs and produce a diff:

```bash
nextflow -log run_1.log run $pipeline -dump-hashes json
nextflow -log run_2.log run $pipeline -dump-hashes json -resume
get_hashes() {
    cat $1 \
    | grep 'cache hash:' \
    | cut -d ' ' -f 10- \
    | sort \
    | awk '{ print; print ""; }'
}
get_hashes run_1.log > run_1.tasks.log
get_hashes run_2.log > run_2.tasks.log
diff run_1.tasks.log run_2.tasks.log
```

You can then view the `diff` output or use a graphical diff viewer to compare `run_1.tasks.log` and `run_2.tasks.log`.

:::{versionadded} 25.04.0
Nextflow now has a built-in way to compare two task runs. See the {ref}`data-lineage-page` guide for details.
:::

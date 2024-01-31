(cache-resume-page)=

# Caching and resuming

One of the core features of Nextflow is the ability to cache task executions and re-use them in subsequent runs to minimize duplicate work. Resumability is useful both for recovering from errors and for iteratively developing a pipeline. It is similar to [checkpointing](https://en.wikipedia.org/wiki/Application_checkpointing), a common practice used by HPC applications.

You can enable resumability in Nextflow with the `-resume` flag when launching a pipeline with `nextflow run`. In most cases, that is all you need to do and resumability will "just work". This page describes Nextflow's caching behavior in more detail in order to help advanced users understand how the cache works and troubleshoot it when it doesn't work.

## Task cache

All task executions are automatically saved to the task cache, regardless of the `-resume` option (so that you always have the option to resume later). The task cache is a key-value store, where each key-value pair corresponds to a previously-executed task.

The task cache is used in conjunction with the [work directory](#work-directory) to recover cached tasks in a resumed run. It is also used by the {ref}`cli-log` sub-command to query task metadata.

### Task hash

The task hash is computed from the following metadata:

- Session ID (see `workflow.sessionId` in {ref}`metadata-workflow`)
- Task name (see `name` in {ref}`trace-report`)
- Task container image (if applicable)
- Task {ref}`environment modules <process-module>` (if applicable)
- Task {ref}`Conda environment <process-conda>` (if applicable)
- Task {ref}`Spack environment <process-spack>` and {ref}`CPU architecture <process-arch>` (if applicable)
- Task {ref}`process-ext` directive (if applicable)
- Task {ref}`inputs <process-input>`
- Task {ref}`script <process-script>`
- Any global variables referenced in the task script
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

Cache failures happen when either (1) a task that was supposed to be cached was re-executed, or (2) a task that was supposed to be re-executed was cached.

When this happens, consider the following questions:

- Is resume enabled via `-resume`?
- Is the {ref}`process-cache` directive set to a non-default value?
- Is the task still present in the task cache and work directory?
- Were any of the task inputs changed?

Changing any of the inputs included in the [task hash](#task-hash) will invalidate the cache, for example:

- Resuming from a different session ID
- Changing the process name
- Changing the task container image or Conda environment
- Changing the task script
- Changing an input file or bundled script used by the task

While the following examples would not invalidate the cache:

- Changing the value of a directive (other than {ref}`process-ext`), even if that directive is used in the task script

In many cases, cache failures happen because of a change to the pipeline script or configuration, or because the pipeline itself has some non-deterministic behavior.

Here are some common reasons for cache failures:

### Modified input files

Make sure that your input files have not been changed. Keep in mind that the default caching mode uses the complete file path, the last modified timestamp, and the file size. If any of these attributes change, the task will be re-executed, even if the file content is unchanged.

### Process that modifies its inputs

If a process modifies its own input files, it cannot be resumed for the reasons described in the previous point. As a result, processes that modify their own input files are considered an anti-pattern and should be avoided.

### Inconsistent file attributes

Some shared file systems, such as NFS, may report inconsistent file timestamps, which can invalidate the cache. If you encounter this problem, you can avoid it by using the `'lenient'` {ref}`caching mode <process-cache>`, which ignores the last modified timestamp and uses only the file path and size.

(cache-global-var-race-condition)=

### Race condition on a global variable

While Nextflow tries to make it easy to write safe concurrent code, it is still possible to create race conditions, which can in turn impact the caching behavior of your pipeline.

Consider the following example:

```groovy
Channel.of(1,2,3) | map { it -> X=it; X+=2 } | view { "ch1 = $it" }
Channel.of(1,2,3) | map { it -> X=it; X*=2 } | view { "ch2 = $it" }
```

The problem here is that `X` is declared in each `map` closure without the `def` keyword (or other type qualifier). Using the `def` keyword makes the variable local to the enclosing scope; omitting the `def` keyword makes the variable global to the entire script.

Because `X` is global, and operators are executed concurrently, there is a *race condition* on `X`, which means that the emitted values will vary depending on the particular order of the concurrent operations. If the values were passed as inputs into a process, the process would execute different tasks on each run due to the race condition.

The solution is to not use a global variable where a local variable is enough (or in this simple example, avoid the variable altogether):

```groovy
// local variable
Channel.of(1,2,3) | map { it -> def X=it; X+=2 } | view { "ch1 = $it" }

// no variable
Channel.of(1,2,3) | map { it -> it * 2 } | view { "ch2 = $it" }
```

(cache-nondeterministic-inputs)=

### Non-deterministic process inputs

Sometimes a process needs to merge inputs from different sources. Consider the following example:

```groovy
workflow {
    ch_foo = Channel.of( ['1', '1.foo'], ['2', '2.foo'] )
    ch_bar = Channel.of( ['2', '2.bar'], ['1', '1.bar'] )
    gather(ch_foo, ch_bar)
}

process gather {
    input:
    tuple val(id), file(foo)
    tuple val(id), file(bar)
    """
    merge_command $foo $bar
    """
}
```

It is tempting to assume that the process inputs will be matched by `id` like the {ref}`operator-join` operator. But in reality, they are simply merged like the {ref}`operator-merge` operator. As a result, not only will the process inputs be incorrect, they will also be non-deterministic, thus invalidating the cache.

The solution is to explicitly join the two channels before the process invocation:

```groovy
workflow {
    ch_foo = Channel.of( ['1', '1.foo'], ['2', '2.foo'] )
    ch_bar = Channel.of( ['2', '2.bar'], ['1', '1.bar'] )
    gather(ch_foo.join(ch_bar))
}

process gather {
    input:
    tuple val(id), file(foo), file(bar)
    """
    merge_command $foo $bar
    """
}
```

## Tips

### Resuming from a specific run

Nextflow resumes from the previous run by default. If you want to resume from an earlier run, simply specify the session ID for that run with the `-resume` option:

```bash
nextflow run rnaseq-nf -resume 4dc656d2-c410-44c8-bc32-7dd0ea87bebf
```

You can use the {ref}`cli-log` command to view all previous runs as well as the task executions for each run.

(cache-compare-hashes)=

### Comparing the hashes of two runs

One way to debug a resumed run is to compare the task hashes of each run using the `-dump-hashes` option.

1. Perform an initial run: `nextflow -log run_initial.log run <pipeline> -dump-hashes`
2. Perform a resumed run: `nextflow -log run_resumed.log run <pipeline> -dump-hashes -resume`
3. Extract the task hash lines from each log (search for `cache hash:`)
4. Compare the runs with a diff viewer

While some manual effort is required, the final diff can often reveal the exact change that caused a task to be re-executed.

:::{versionadded} 23.10.0
:::

When using `-dump-hashes json`, the task hashes can be more easily extracted into a diff. Here is an example Bash script to perform two runs and produce a diff:

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

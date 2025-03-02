(cache-resume-page)=

# Caching and resuming

One of the core features of Nextflow is the ability to cache task executions and re-use them in subsequent runs to minimize duplicate work. Resumability is useful both for recovering from errors and for iteratively developing a pipeline. It is similar to [checkpointing](https://en.wikipedia.org/wiki/Application_checkpointing), a common practice used by HPC applications.

You can enable resumability in Nextflow with the `-resume` flag when launching a pipeline with `nextflow run`. In most cases, that is all you need to do and resumability will "just work". This page describes Nextflow's caching behavior in more detail in order to help advanced users understand how the cache works and troubleshoot it when it doesn't work.

## Task cache

All task executions are automatically saved to the task cache, regardless of the `-resume` option (so that you always have the option to resume later). The task cache is a key-value store, where each key-value pair corresponds to a previously-executed task.

The task cache is used in conjunction with the [work directory](#work-directory) to recover cached tasks in a resumed run. It is also used by the {ref}`cli-log` sub-command to query task metadata.

### Task hash

The task hash is computed from the following metadata:

- Session ID (see `workflow.sessionId` in {ref}`stdlib-constants`)
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

## Resume from a specific run

Nextflow resumes from the previous run by default. If you want to resume from an earlier run, simply specify the session ID for that run with the `-resume` option:

```bash
nextflow run rnaseq-nf -resume 4dc656d2-c410-44c8-bc32-7dd0ea87bebf
```

You can use the {ref}`cli-log` command to view all previous runs as well as the task executions for each run.

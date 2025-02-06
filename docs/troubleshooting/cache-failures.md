(cache-failure-page)=

# Cache failures

Cache failures occur when a task that was supposed to be cached was re-executed or a task that was supposed to be re-executed was cached. This page provides an overview of common causes for cache failures and strategies to identify them.

(cache-failure-common)=

## Common causes

Common causes of cache failures include:

- {ref}`Resume not being enabled <cache-failure-resume>`
- {ref}`Non-default cache directives <cache-failure-directives>`
- {ref}`Modified inputs <cache-failure-modified>`
- {ref}`Inconsistent file attributes <cache-failure-inconsistent>`
- {ref}`Race condition on a global variable <cache-failure-race-condition>`
- {ref}`Non-deterministic process inputs <cache-failure-nondeterministic>`

The causes of these cache failure and solutions to resolve them are described in detail below.

(cache-failure-resume)=

### Resume not enabled

The `-resume` option is required to resume a pipeline. Ensure `-resume` has been enabled in your run command or your nextflow configuration file.

(cache-failure-directives)=

### Non-default cache directives

The `cache` directive is enabled by default. However, you can disable or modify it's behavior for a specific process. For example:

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

Some shared file systems, such as NFS, may report inconsistent file timestamps. If you encounter this problem, use the `'lenient'` {ref}`caching mode <process-cache>` to ignore the last modified timestamp and use only the file path.

(cache-failure-race-condition)=

### Race condition on a global variable

Race conditions can in disrupt caching behavior of your pipeline. For example:

```nextflow
Channel.of(1,2,3) | map { v -> X=v; X+=2 } | view { v -> "ch1 = $v" }
Channel.of(1,2,3) | map { v -> X=v; X*=2 } | view { v -> "ch2 = $v" }
```

In the above example, `X` is declared in each `map` closure. Without the `def` keyword, or other type qualifier, the variable `X` is global to the entire script. Operators and executed concurrently and, as `X` is global, there is a *race condition* that causes the emitted values to vary depending on the order of the concurrent operations. If these values were passed to a process as inputs the process would execute different tasks during each run due to the race condition.

To resolve this failure type, ensure the variable is not global by using a local variable:
    
```nextflow
Channel.of(1,2,3) | map { v -> def X=v; X+=2 } | view { v -> "ch1 = $v" }
```

Alternatively, remove the variable:

```nextflow
Channel.of(1,2,3) | map { v -> v * 2 } | view { v -> "ch2 = $v" }
```

(cache-failure-nondeterministic)=

### Non-deterministic process inputs

A process that merges inputs from different sources non-deterministically may invalidate the cache. For example:

```nextflow
workflow {
    ch_foo = Channel.of( ['1', '1.foo'], ['2', '2.foo'] )
    ch_bar = Channel.of( ['2', '2.bar'], ['1', '1.bar'] )
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
    ch_foo = Channel.of( ['1', '1.foo'], ['2', '2.foo'] )
    ch_bar = Channel.of( ['2', '2.bar'], ['1', '1.bar'] )
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

(cache-failure-compare)=

## Compare task hashes

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

Task hashes can also be extracted into a diff using `-dump-hashes json`. See an example Bash script to compare two runs and produce a diff here:

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

(dsl1-page)=

# Migrating from DSL1

In Nextflow version `22.03.0-edge`, DSL2 became the default DSL version. In version `22.12.0-edge`, DSL1 support was removed and the Nextflow documentation was updated to use DSL2 by default. Users who are still using DSL1 should migrate their pipelines to DSL2 to use the latest versions of Nextflow. This page describes the differences between DSL1 and DSL2 and how to migrate to DSL2.

## Processes and workflows

In DSL1, a process definition is also the process invocation. Process inputs and outputs are connected to channels using `from` and `into`. You can see a basic Nextflow script written in DSL1 here:

```nextflow
nextflow.enable.dsl=1

params.str = 'Hello world!'

process splitLetters {
    output:
    file 'chunk_*' into letters

    script:
    """
    printf '${params.str}' | split -b 6 - chunk_
    """
}

process convertToUpper {
    input:
    file x from letters.flatten()

    output:
    stdout result

    script:
    """
    cat $x | tr '[a-z]' '[A-Z]'
    """
}

result.view { it.trim() }
```

To migrate this code to DSL2, you need to move all of your channel logic throughout the script into a `workflow` definition. Additionally, you must call each process explicitly, passing any input channels as arguments (instead of `from ...`) and receiving any output channels as return values (instead of `into ...`).

See {ref}`workflow-page` page to learn how to define a workflow.

You can see the DSL1 Nextflow script from above written in DSL2 here:

```nextflow
params.str = 'Hello world!'

process splitLetters {
    output:
    path 'chunk_*'

    script:
    """
    printf '${params.str}' | split -b 6 - chunk_
    """
}

process convertToUpper {
    input:
    path x

    output:
    stdout

    script:
    """
    cat $x | tr '[a-z]' '[A-Z]'
    """
}

workflow {
    splitLetters | flatten | convertToUpper | view { v -> v.trim() }
}
```

## Channel forking

In DSL1, a channel can be used as an input only once; to use a channel multiple times, the channel must be forked using the `into` operator. In DSL2, channels are automatically forked when connecting two or more consumers. For example:

```nextflow
Channel
    .of('Hello','Hola','Ciao')
    .set { cheers }

cheers
    .map { v -> v.toUpperCase() }
    .view()

cheers
    .map { v -> v.reverse() }
    .view()
```

Similarly, in DSL2, process outputs can be consumed by multiple consumers automatically, which makes workflow scripts much easier to read and write.

## Modules

In DSL1, the entire Nextflow pipeline must be defined in a single file. For example, `main.nf`. This restriction becomes cumbersome as a pipeline grows and hinders the sharing and reuse of pipeline components.

DSL2 introduces the concept of "module scripts" (or "modules" for short), which are Nextflow scripts that can be "included" by other scripts. While modules are not essential to migrating to DSL2, nor are they mandatory in DSL2, modules can help you organize a large pipeline into multiple smaller files and take advantage of modules created by others. See {ref}`module-page` to learn more about modules.

:::{note}
DSL2 scripts cannot exceed 64 KB in size. Split large DSL1 scripts into modules to avoid this limit.
:::

## Deprecations

### Processes

- The `set` process input type is no longer supported, use {ref}`tuple <process-input-tuple>` instead.

- The `set` process output type is no longer supported, use {ref}`tuple <process-out-tuple>` instead.

- The `mode flatten` option for process outputs is no longer available. Use the {ref}`operator-flatten` operator on the corresponding output channel instead.

- Unqualified value and file elements in a tuple declaration are no longer allowed. Use an explicit `val` or `path` qualifier. For example:

  ```nextflow
  process foo {
      input:
      tuple X, 'some-file.sam'

      output:
      tuple X, 'some-file.bam'

      script:
      """
      your_command --in $X some-file.sam > some-file.bam
      """
  }
  ```

  Use:

  ```nextflow
  process foo {
      input:
      tuple val(X), path('some-file.sam')

      output:
      tuple val(X), path('some-file.bam')

      script:
      """
      your_command --in $X some-file.sam > some-file.bam
      """
  }
  ```

### Channels

- Channel method `bind` has been deprecated in DSL2.
- Channel method `<<` has been deprecated in DSL2.
- Channel factory `create` has been deprecated in DSL2.

### Operators

- Operator `choice` has been deprecated in DSL2. Use {ref}`operator-branch` instead.
- Operator `close` has been deprecated in DSL2.
- Operator `countBy` has been deprecated in DSL2.
- Operator `into` has been deprecated in DSL2, as it is no longer needed.
- Operator `fork` has been renamed to {ref}`operator-multimap`.
- Operator `groupBy` has been deprecated in DSL2. Use {ref}`operator-grouptuple` instead.
- Operators `print` and `println` have been deprecated in DSL2. Use {ref}`operator-view` instead.
- Operator `route` has been deprecated in DSL2.
- Operator `separate` has been deprecated in DSL2.
- Operator `spread` has been deprecated in DSL2. Use {ref}`operator-combine` instead.

### DSL2 Preview

An early preview of DSL2 was available in 2020. Note that some of that early DSL2 syntax has since changed.

- The `nextflow.preview.dsl=2` (and `nextflow.enable.dsl=1`) feature flags are no longer needed.

- Anonymous and unwrapped includes are no longer supported. Use an explicit module inclusion instead.

  For example:

  ```nextflow
  include './some/library'
  include bar from './other/library'

  workflow {
      foo()
      bar()
  }
  ```

  Should be replaced with:

  ```nextflow
  include { foo } from './some/library'
  include { bar } from './other/library'

  workflow {
      foo()
      bar()
  }
  ```

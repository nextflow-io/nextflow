(dsl1-page)=

# Migrating from DSL 1

In Nextflow version `22.03.0-edge`, DSL2 became the default DSL version. In version `22.12.0-edge`, DSL1 support was removed, and the Nextflow documentation was updated to use DSL2 by default. Users who are still using DSL1 should migrate their pipelines to DSL2 in order to use the latest versions of Nextflow. This page describes the differences between DSL1 and DSL2, and how to migrate to DSL2.

In Nextflow versions prior to `22.03.0-edge`, you must enable DSL2 explicitly in order to use it. You can either set the feature flag in your pipeline script:

```groovy
nextflow.enable.dsl=2
```

Or set the environment variable where you launch Nextflow:

```bash
export NXF_DEFAULT_DSL=2
```

## Processes and workflows

In DSL1, a process definition is also the process invocation. Process inputs and outputs are connected to channels using `from` and `into`. Here is the {ref}`getstarted-first` example written in DSL1:

```groovy
nextflow.enable.dsl=1

params.str = 'Hello world!'

process splitLetters {
    output:
    file 'chunk_*' into letters

    """
    printf '${params.str}' | split -b 6 - chunk_
    """
}

process convertToUpper {
    input:
    file x from letters.flatten()

    output:
    stdout result

    """
    cat $x | tr '[a-z]' '[A-Z]'
    """
}

result.view { it.trim() }
```

To migrate this code to DSL2, you need to move all of your channel logic throughout the script into a `workflow` definition. Additionally, you must call each process explicitly, passing any input channels as arguments (instead of `from ...`) and receiving any output channels as return values (instead of `into ...`).

Refer to the {ref}`workflow-page` page to learn how to define a workflow. The DSL2 version of the above script is duplicated here for your convenience:

```{literalinclude} snippets/your-first-script.nf
:language: groovy
```

## Channel forking

In DSL1, a channel can be used as an input only once; to use a channel multiple times, the channel must be forked using the `into` operator.

In DSL2, channels are automatically forked when connecting two or more consumers.

For example, this would not work in DSL1 but is not a problem in DSL2:

```groovy
Channel
    .from('Hello','Hola','Ciao')
    .set{ cheers }

cheers
    .map{ it.toUpperCase() }
    .view()

cheers
    .map{ it.reverse() }
    .view()
```

Similarly, process outputs can be consumed by multiple consumers automatically, which makes workflow scripts much easier to read and write.

## Modules

In DSL1, the entire Nextflow pipeline must be defined in a single file (e.g. `main.nf`). This restriction becomes quite cumbersome as a pipeline becomes larger, and it hinders the sharing and reuse of pipeline components.

DSL2 introduces the concept of "module scripts" (or "modules" for short), which are Nextflow scripts that can be "included" by other scripts. While modules are not essential to migrating to DSL2, nor are they mandatory in DSL2 by any means, modules can help you organize a large pipeline into multiple smaller files, and take advantage of modules created by others. Check out the {ref}`module-page` to get started.

:::{note}
With DSL2, the Groovy shell used by Nextflow also imposes a 64KB size limit on pipeline scripts, so if your DSL1 script is very large, you may need to split your script into modules anyway to avoid this limit.
:::

## Deprecations

### Processes

- The `set` process input type is no longer supported, use {ref}`tuple <process-input-tuple>` instead.

- The `set` process output type is no longer supported, use {ref}`tuple <process-out-tuple>` instead.

- The `mode flatten` option for process outputs is no longer available. Use the {ref}`operator-flatten` operator on the corresponding output channel instead.

- Unqualified value and file elements in a tuple declaration are no longer allowed. Use an explicit `val` or `path` qualifier.

  For example:

  ```groovy
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

  ```groovy
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

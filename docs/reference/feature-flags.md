(config-feature-flags)=

# Feature flags

Feature flags introduce experimental or other opt-in features. They must be specified in the pipeline script.

:::warning
Deprecated feature flags may cause pipelines run with newer versions of Nextflow to fail.
:::

`nextflow.enable.configProcessNamesValidation`
: :::{deprecated} 25.10.0
  Nextflow 25.10.0 and later will fail if this flag is present. Remove this flag from your pipeline. Use the {ref}`strict syntax <strict-syntax-page>` instead for more accurate validation of process selectors.
  :::
: When `true`, prints a warning for every `withName:` process selector that doesn't match a process in the pipeline (default: `true`).

`nextflow.enable.dsl`
: :::{deprecated} 25.04.0
  :::
: Specifies the DSL version to use. Options: `1` or `2`.

`nextflow.enable.moduleBinaries`
: When `true`, enables modules with binary scripts. See {ref}`module-binaries` for more information.

`nextflow.enable.strict`
: When `true`, enables "strict" mode. In strict mode, Nextflow:

  - Fails when reading a params file if a dynamic param value references an undefined variable

  - Fails when merging params if the same param is specified from a config file and the command line with different types

  - Fails if an input or output tuple has only one element in a process definition

  - Fails if an output emit name is not a valid identifier in a process definition (i.e., it should match the pattern `/[A-Za-z_][A-Za-z0-9_]*/`)

  - Fails if a the number of elements in a received input tuple does not match the number of elements that were declared

  - Fails if the `storeDir` directive is used with non-file outputs

  - Fails if pipeline params are referenced before they are defined

  - Fails if multiple functions and/or processes with the same name in a module script are defined

  - Sets `failOnDuplicate` to `true` for the `join` operator, regardless of user settings

  - Sets `failOnMismatch` to `true` for the `join` operator (unless `remainder` is `true`), regardless of user settings

  - Sets `failOnError` to `true` for the `publishDir` directive, regardless of user settings

`nextflow.preview.output`
: :::{versionadded} 24.04.0
  :::
: :::{deprecated} 25.10.0
  Nextflow 25.10.0 and later will fail if this flag is present. Remove this flag from your pipeline. Workflow outputs are now enabled by default.
  :::
: When `true`, enables {ref}`workflow outputs <workflow-output-def>`.

`nextflow.preview.recursion`
: *Experimental: may change in a future release.*
: When `true`, enables {ref}`process and workflow recursion <workflow-recursion>`.

`nextflow.preview.topic`
: :::{versionadded} 24.04.0
  :::
: :::{deprecated} 25.04.0
  Nextflow 25.04.0 and later will fail if this flag is present. Remove this flag from your pipeline. Topic channels are now enabled by default.
  :::
: When `true`, enables {ref}`topic channels <channel-topic>`.

`nextflow.preview.types`
: :::{versionadded} 25.10.0
  :::
: When `true`, enables {ref}`typed processes <process-typed-page>`. Must be enabled in every script that uses typed processes. Legacy processes cannot be defined in scripts with this flag enabled.

(config-feature-flags)=

# Feature flags

Feature flags are used to introduce experimental or other opt-in features. They must be specified in the pipeline script.

`nextflow.enable.configProcessNamesValidation`
: :::{deprecated} 25.10.0
  Use the {ref}`strict syntax <strict-syntax-page>` instead. It validates process selectors without producing false warnings.
  :::
: When `true`, prints a warning for every `withName:` process selector that doesn't match a process in the pipeline (default: `true`).

`nextflow.enable.dsl`
: :::{deprecated} 25.04.0
  :::
: Defines the DSL version to use (`1` or `2`).

`nextflow.enable.moduleBinaries`
: When `true`, enables the use of modules with binary scripts. See {ref}`module-binaries` for more information.

`nextflow.enable.strict`
: When `true`, the pipeline is executed in "strict" mode, which introduces the following rules:

  - When reading a params file, Nextflow will fail if a dynamic param value references an undefined variable

  - When merging params from a config file with params from the command line, Nextflow will fail if a param is specified from both sources but with different types

  - When using the `join` operator, the `failOnDuplicate` option is `true` regardless of any user setting

  - When using the `join` operator, the `failOnMismatch` option is `true` (unless `remainder` is also `true`) regardless of any user setting

  - When using the `publishDir` process directive, the `failOnError` option is `true` regardless of any user setting

  - In a process definition, Nextflow will fail if an input or output tuple has only one element

  - In a process definition, Nextflow will fail if an output emit name is not a valid identifier (i.e. it should match the pattern `/[A-Za-z_][A-Za-z0-9_]*/`)

  - During a process execution, Nextflow will fail if a received input tuple does not have the same number of elements as was declared

  - During a process execution, Nextflow will fail if the `storeDir` directive is used with non-file outputs

  - Nextflow will fail if a pipeline param is referenced before it is defined

  - Nextflow will fail if multiple functions and/or processes with the same name are defined in a module script

`nextflow.preview.output`
: :::{versionadded} 24.04.0
  :::
: *Experimental: may change in a future release.*
: When `true`, enables the use of the {ref}`workflow output definition <workflow-output-def>`.

`nextflow.preview.recursion`
: *Experimental: may change in a future release.*
: When `true`, enables {ref}`process and workflow recursion <workflow-recursion>`.

`nextflow.preview.topic`
: :::{versionadded} 24.04.0
  :::
: :::{deprecated} 25.04.0
  This feature flag is no longer required to use topic channels.
  :::
: When `true`, enables {ref}`topic channels <channel-topic>` feature.

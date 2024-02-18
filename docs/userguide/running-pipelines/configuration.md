<!-- FROM basic.md -->

# Configuration basics

Pipeline configuration properties are defined in a file named `nextflow.config` in the pipeline execution directory.

This file can be used to define which executor to use, the process's environment variables, pipeline parameters etc.

A basic configuration file might look like this:

```groovy
process {
  executor = 'sge'
  queue = 'cn-el6'
}
```

Read the {ref}`config-page` section to learn more about the Nextflow configuration file and settings.

<!-- TODO: user-guide content from config.md -->

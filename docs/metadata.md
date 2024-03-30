(metadata-page)=

# Workflow introspection

(metadata-workflow)=

## Runtime metadata

The implicit `workflow` object allows you to access some workflow and runtime metadata in your Nextflow scripts. For example:

```groovy
println "Project : $workflow.projectDir"
println "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"
println "Cmd line: $workflow.commandLine"
println "Manifest's pipeline version: $workflow.manifest.version"
```

:::{tip}
To shortcut access to multiple `workflow` properties, you can use the Groovy [with](<http://docs.groovy-lang.org/latest/html/groovy-jdk/java/lang/Object.html#with(groovy.lang.Closure)>) method.
:::

The following table lists the properties that can be accessed on the `workflow` object:


`workflow.commandLine`
: Command line as entered by the user to launch the workflow execution.

`workflow.commitId`
: Git commit ID of the executed workflow repository.
: When providing a Git tag, branch name, or commit hash using the `-r` CLI option, the associated `workflow.commitId` is also populated.

`workflow.complete`
: *Available only in the `workflow.onComplete` handler*
: Timestamp of workflow when execution is completed.

`workflow.configFiles`
: Configuration files used for the workflow execution.

`workflow.container`
: Docker image used to run workflow tasks. When more than one image is used it returns a map object containing `[process name, image name]` pair entries.

`workflow.containerEngine`
: Returns the name of the container engine (e.g. docker or singularity) or null if no container engine is enabled.

`workflow.duration`
: *Available only in the `workflow.onComplete` handler*
: Time elapsed to complete workflow execution.

`workflow.errorMessage`
: *Available only in the `workflow.onComplete` and `workflow.onError` handlers*
: Error message of the task that caused the workflow execution to fail.

`workflow.errorReport`
: *Available only in the `workflow.onComplete` and `workflow.onError` handlers*
: Detailed error of the task that caused the workflow execution to fail.

`workflow.exitStatus`
: *Available only in the `workflow.onComplete` and `workflow.onError` handlers*
: Exit status of the task that caused the workflow execution to fail.

`workflow.homeDir`
: User system home directory.

`workflow.launchDir`
: Directory where the workflow execution has been launched.

`workflow.manifest`
: Entries of the workflow manifest.

`workflow.profile`
: Used configuration profile.

`workflow.projectDir`
: Directory where the workflow project is stored in the computer.

`workflow.repository`
: Project repository Git remote URL.

`workflow.resume`
: Returns `true` whenever the current instance is resumed from a previous execution.

`workflow.revision`
: Git branch/tag of the executed workflow repository.
: When providing a Git tag or branch name using the `-r` CLI option, the `workflow.revision` is also populated.

`workflow.runName`
: Mnemonic name assigned to this execution instance.

`workflow.scriptFile`
: Project main script file path.

`workflow.scriptId`
: Project main script unique hash ID.

`workflow.scriptName`
: Project main script file name.

`workflow.sessionId`
: Unique identifier (UUID) associated to current execution.

`workflow.start`
: Timestamp of workflow at execution start.

`workflow.stubRun`
: Returns `true` whenever the current instance is a stub-run execution .

`workflow.success`
: *Available only in the `workflow.onComplete` and `workflow.onError` handlers*
: Reports if the execution completed successfully.

`workflow.userName`
: User system account name.

`workflow.workDir`
: Workflow working directory.

(metadata-nextflow)=

## Nextflow metadata

The implicit `nextflow` object allows you to access the metadata information of the Nextflow runtime.

`nextflow.build`
: Nextflow runtime build number.

`nextflow.timestamp`
: Nextflow runtime compile timestamp.

`nextflow.version`
: Nextflow runtime version number.

`nextflow.version.matches()`
: This method allows you to check if the Nextflow runtime satisfies a version requirement for your workflow script. The version requirement string can be prefixed with the usual comparison operators eg `>`, `>=`, `=`, etc. or postfixed with the `+` operator to specify a minimum version requirement. For example:

  ```groovy
  if( !nextflow.version.matches('21.04+') ) {
      println "This workflow requires Nextflow version 21.04 or greater -- You are running version $nextflow.version"
      exit 1
  }
  ```

(metadata-completion-handler)=

## Completion handler

Due to the asynchronous nature of Nextflow the termination of a script does not correspond to the termination of the running workflow. Thus some information, only available on execution completion, needs to be accessed by using an asynchronous handler.

The `onComplete` event handler is invoked by the framework when the workflow execution is completed. It allows one to access the workflow termination status and other useful information. For example:

```groovy
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
```

If you want an e-mail notification on completion, check {ref}`mail-page`.

(metadata-error-handler)=

## Error handler

The `onError` event handler is invoked by Nextflow when a runtime or process error caused the pipeline execution to stop. For example:

```groovy
workflow.onError {
    println "Error: Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
```

:::{note}
Both the `onError` and `onComplete` handlers are invoked when an error condition is encountered. The first is called as soon as the error is raised, while the second is called just before the pipeline execution is about to terminate. When using the `finish` {ref}`process-error-strategy`, there may be a significant gap between the two, depending on the time required to complete any pending job.
:::

## Decoupling metadata

The workflow event handlers can be defined also in the `nextflow.config` file. This is useful to decouple the handling of pipeline events from the main script logic.

When the event handlers are included in a configuration file the only difference is that the `onComplete` and the `onError` closures have to be defined by using the assignment operator as shown below:

```groovy
workflow.onComplete = {
    // any workflow property can be used here
    println "Pipeline complete"
    println "Command line: $workflow.commandLine"
}

workflow.onError = {
    println "Error: something when wrong"
}
```

:::{note}
It is possible to define workflow event handlers both in the pipeline script **and** in the configuration file.
:::

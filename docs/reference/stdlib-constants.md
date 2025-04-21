(stdlib-constants)=

# Constants

The following constants are globally available in a Nextflow script for workflow introspection:

`baseDir: Path`
: :::{deprecated} 20.04.0
  :::
: Alias of `workflow.projectDir`.

`launchDir: Path`
: Alias of `workflow.launchDir`.

`moduleDir: Path`
: The directory where a module script is located (equivalent to `projectDir` if used in the main script).

`nextflow`
: Map of Nextflow runtime information. The following properties are available:

  `build: int`
  : Nextflow runtime build number.

  `timestamp: String`
  : Nextflow runtime compile timestamp.

  `version: VersionNumber`
  : Nextflow runtime version number. See {ref}`VersionNumber <stdlib-types-versionnumber>` for more information.

`params`
: Map of workflow parameters specified in the config file or as command line options.

`projectDir: Path`
: Alias of `workflow.projectDir`.

`secrets: Map<String,String>`
: :::{versionadded} 24.02.0-edge
  :::
: Map of pipeline secrets. See {ref}`secrets-page` for more information.

`workDir: Path`
: Alias of `workflow.workDir`.

`workflow`
: Map of workflow runtime information. The following properties are available:

  `commandLine: String`
  : Command line as entered by the user to launch the workflow execution.

  `commitId: String`
  : Git commit ID of the executed workflow repository.
  : When providing a Git tag, branch name, or commit hash using the `-r` CLI option, the associated `workflow.commitId` is also populated.

  `complete: OffsetDateTime`
  : *Available only in the `workflow.onComplete` handler*
  : Timestamp of workflow when execution is completed.

  `configFiles: List<Path>`
  : Configuration files used for the workflow execution.

  `container: String | Map<String,String>`
  : Docker image used to run workflow tasks, or a map of process names to process containers when multiple images are used.

  `containerEngine: String`
  : Returns the name of the container engine (e.g. docker or singularity) or null if no container engine is enabled.

  `duration: Duration`
  : *Available only in the `workflow.onComplete` handler*
  : Time elapsed to complete workflow execution.

  `errorMessage: String`
  : *Available only in the `workflow.onComplete` and `workflow.onError` handlers*
  : Error message of the task that caused the workflow execution to fail.

  `errorReport: String`
  : *Available only in the `workflow.onComplete` and `workflow.onError` handlers*
  : Detailed error of the task that caused the workflow execution to fail.

  `exitStatus: int`
  : *Available only in the `workflow.onComplete` and `workflow.onError` handlers*
  : Exit status of the task that caused the workflow execution to fail.

  `failOnIgnore: boolean`
  : :::{versionadded} 24.05.0-edge
    :::
  : Whether the `workflow.failOnIgnore` config option was enabled.
  : See also: {ref}`process-error-strategy`

  `fusion`
  : Map of Fusion runtime information. The following properties are available:

  : `enabled: boolean`
    : Whether Fusion is enabled.

  : `version: String`
    : The Fusion version being used.

  `homeDir: Path`
  : User system home directory.

  `launchDir: Path`
  : Directory where the workflow was launched.

  `manifest`
  : Map of properties corresponding to the {ref}`config-manifest` config scope.

  `onComplete( action: Closure )`
  : Define an action to take when the workflow completes (whether successful or not).

  `onError( action: Closure )`
  : Define an action to take if the workflow is terminated due to a runtime error or task failure.

  `outputDir: Path`
  : :::{versionadded} 24.10.0
    :::
  : Workflow output directory.

  `preview: boolean`
  : :::{versionadded} 24.04.0
    :::
  : Whether the current workflow run is a preview run.

  `profile: String`
  : Comma-separated list of active configuration profiles.

  `projectDir: Path`
  : Directory where the workflow project is located.

  `repository: String`
  : Project repository Git remote URL.

  `resume: boolean`
  : Returns `true` whenever the current instance is resumed from a previous execution.

  `revision: String`
  : Git branch/tag of the executed workflow repository.
  : When providing a Git tag or branch name using the `-r` CLI option, the `workflow.revision` is also populated.

  `runName: String`
  : Mnemonic name assigned to this execution instance.

  `scriptFile: Path`
  : Project main script file path.

  `scriptId: String`
  : Project main script unique hash ID.

  `scriptName: String`
  : Project main script file name.

  `sessionId: UUID`
  : Unique identifier (UUID) associated to current execution.

  `start: OffsetDateTime`
  : Timestamp of workflow at execution start.

  `stubRun: boolean`
  : Returns `true` whenever the current instance is a stub-run execution .

  `success: boolean`
  : *Available only in the `workflow.onComplete` and `workflow.onError` handlers*
  : Reports if the execution completed successfully.

  `userName: String`
  : User system account name.

  `wave`
  : Map of Wave runtime information. The following properties are available:

  : `enabled: boolean`
    : Whether Wave is enabled.

  `workDir: Path`
  : The directory where task temporary files are stored.

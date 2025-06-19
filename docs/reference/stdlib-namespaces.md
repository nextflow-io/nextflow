(stdlib-namespaces)=

# Namespaces

This page lists all of the available namespaces in the Nextflow standard library.

(stdlib-namespaces-global)=

## Global namespace

The global namespace contains globally available constants and functions.

**Constants**

`baseDir: Path`
: :::{deprecated} 20.04.0
  :::
: Alias of `workflow.projectDir`.

`launchDir: Path`
: Alias of `workflow.launchDir`.

`moduleDir: Path`
: Directory where a module script is located (equivalent to `projectDir` if used in the main script).

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

**Functions**

`branchCriteria( criteria: Closure ) -> Closure`
: Create a branch criteria to use with the {ref}`operator-branch` operator.

`env( name: String ) -> String`
: :::{versionadded} 24.11.0-edge
  :::
: Get the value of the environment variable with the specified name in the Nextflow launch environment.

`error( message: String = null )`
: Throw a script runtime error with an optional error message.

`exit( exitCode: int = 0, message: String = null )`
: :::{deprecated} 22.06.0-edge
  Use `error()` instead
  :::
: Stop the pipeline execution and return an exit code and optional error message.

`file( filePattern: String, [options] ) -> Path | List<Path>`
: Get a file from a file name or glob pattern. Returns a collection of files if the glob pattern yields zero or multiple files.

: The following options are available:

  `checkIfExists: boolean`
  : When `true`, throws an exception if the specified path does not exist in the file system (default: `false`)

  `followLinks: boolean`
  : When `true`, follows symbolic links when traversing a directory tree, otherwise treats them as files (default: `true`)

  `glob: boolean`
  : When `true`, interprets characters `*`, `?`, `[]` and `{}` as glob wildcards, otherwise handles them as normal characters (default: `true`)

  `hidden: boolean`
  : When `true`, includes hidden files in the resulting paths (default: `false`)

  `maxDepth: int`
  : Maximum number of directory levels to visit (default: *no limit*)

  `type: String`
  : Type of paths returned, can be `'file'`, `'dir'` or `'any'` (default: `'file'`)

: See also: {ref}`channel.fromPath <channel-path>`.

`files( filePattern: String, [options] ) -> List<Path>`
: Get a collection of files from a file name or glob pattern. Supports the same options as `file()`.

`groupKey( key, size: int ) -> GroupKey`
: Create a grouping key to use with the {ref}`operator-grouptuple` operator.

`multiMapCriteria( criteria: Closure ) -> Closure`
: Create a multi-map criteria to use with the {ref}`operator-multiMap` operator.

`print( value )`
: Print a value to standard output.

`printf( format: String, values... )`
: Print a formatted string with the given values to standard output.

`println( value )`
: Print a value to standard output with a newline.

`sendMail( [options] )`
: Send an email. See {ref}`mail-page` for more information.

`sleep( milliseconds: long )`
: Sleep for the given number of milliseconds.

`tuple( collection: List ) -> ArrayTuple`
: Create a tuple object from the given collection.

`tuple( args... ) -> ArrayTuple`
: Create a tuple object from the given arguments.

(stdlib-namespaces-channel)=

## `channel`

The `channel` namespace contains the built-in channel factories. See {ref}`channel-factory` for details.

(stdlib-namespaces-nextflow)=

## `log`

The `log` namepsace contains functions for logging messages to the console.

`error( message: String )`
: Log an error message to the console.
: This function does not terminate the pipeline -- use the global `error()` function instead.

`info( message: String )`
: Log an info message to the console.

`warn( message: String )`
: Log a warning message to the console.

## `nextflow`

The `nextflow` namespace contains information about the current Nextflow runtime.

`build: int`
: Nextflow runtime build number.

`timestamp: String`
: Nextflow runtime compile timestamp.

`version: VersionNumber`
: Nextflow runtime version number. See {ref}`VersionNumber <stdlib-types-versionnumber>` for more information.

(stdlib-namespaces-workflow)=

## `workflow`

The `workflow` namespace contains information about the current workflow run.

**Properties**

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
: Namespace containing information about the current Fusion runtime. The following properties are available:

: `enabled: boolean`
  : Whether Fusion is enabled.

: `version: String`
  : The Fusion version being used.

`homeDir: Path`
: User system home directory.

`launchDir: Path`
: Directory where the workflow was launched.

`manifest`
: Namespace corresponding to the {ref}`config-manifest` config scope.

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
: Namespace containing Wave runtime information. The following properties are available:

: `enabled: boolean`
  : Whether Wave is enabled.

`workDir: Path`
: The directory where task temporary files are stored.

**Functions**

`onComplete( action: Closure )`
: Define an action to take when the workflow completes (whether successful or not).

`onError( action: Closure )`
: Define an action to take if the workflow is terminated due to a runtime error or task failure.

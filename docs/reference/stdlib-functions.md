(stdlib-functions)=

# Functions

The following functions are available in Nextflow scripts:

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

: See also: {ref}`Channel.fromPath <channel-path>`.

`files( filePattern: String, [options] ) -> List<Path>`
: Get a collection of files from a file name or glob pattern. Supports the same options as `file()`.

`groupKey( key, size: int ) -> GroupKey`
: Create a grouping key to use with the {ref}`operator-grouptuple` operator.

`multiMapCriteria( criteria: Closure ) -> Closure`
: Create a multi-map criteria to use with the {ref}`operator-multiMap` operator.

`sendMail( [options] )`
: Send an email. See {ref}`mail-page` for more information.

`sleep( milliseconds: long )`
: Sleep for the given number of milliseconds.

`tuple( collection: List ) -> ArrayTuple`
: Create a tuple object from the given collection.

`tuple( args... ) -> ArrayTuple`
: Create a tuple object from the given arguments.

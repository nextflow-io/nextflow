(stdlib-page)=

# Standard library

This page describes the built-in constants, functions, and types provided by Nextflow.

## Globals

(stdlib-constants)=

### Constants

The following constants are globally available in a Nextflow script:

`baseDir`
: :::{deprecated} 20.04.0
  :::
: Alias of `workflow.projectDir`.

`launchDir`
: Alias of `workflow.launchDir`.

`moduleDir`
: The directory where a module script is located (equivalent to `projectDir` if used in the main script).

`nextflow`
: Map of Nextflow runtime information.

  `nextflow.build`
  : Nextflow runtime build number.

  `nextflow.timestamp`
  : Nextflow runtime compile timestamp.

  `nextflow.version`
  : Nextflow runtime version number.

  `nextflow.version.matches()`
  : Check whether the Nextflow runtime satisfies a version requirement.

  : The version requirement string can be prefixed with the usual comparison operators:
    - `=` or `==`: equal to
    - `<` (`<=`): less than (or equal to)
    - `>` (`>=`): greater than (or equal to)
    - `!=` or `<>`: not equal

    For example:

    ```nextflow
    if( !nextflow.version.matches('>=23.10') ) {
        error "This workflow requires Nextflow version 23.10 or greater -- You are running version $nextflow.version"
    }
    ```

  : Alternatively, the version can be postfixed with `+`, which is similar to `==` but also allows the last version part to be greater. For example, `23.10.1+` is satisfied by `23.10.1` and `23.10.2`, but not `23.11.x` or `23.09.x`. Additionally, `23.10.+` is equivalent to `23.10.0+`. This operator is a useful way to enforce a specific version while allowing for newer patch releases.

`params`
: Map of workflow parameters specified in the config file or as command line options.

`projectDir`
: Alias of `workflow.projectDir`.

`secrets`
: :::{versionadded} 24.02.0-edge
  :::
: Dictionary like object holding workflow secrets. Read the {ref}`secrets-page` page for more information.

`workDir`
: Alias of `workflow.workDir`.

`workflow`
: Map of workflow runtime information.

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

  `workflow.failOnIgnore`
  : :::{versionadded} 24.05.0-edge
    :::
  : Whether the `workflow.failOnIgnore` config option was enabled.
  : See also: {ref}`process-error-strategy`

  `workflow.fusion.enabled`
  : Whether Fusion is enabled.

  `workflow.fusion.version`
  : Fusion version in use.

  `workflow.homeDir`
  : User system home directory.

  `workflow.launchDir`
  : Directory where the workflow was launched.

  `workflow.manifest`
  : Entries of the workflow manifest.

  `workflow.outputDir`
  : :::{versionadded} 24.10.0
    :::
  : Workflow output directory.

  `workflow.preview`
  : :::{versionadded} 24.04.0
    :::
  : Whether the current workflow run is a preview run.

  `workflow.profile`
  : Used configuration profile.

  `workflow.projectDir`
  : Directory where the workflow project is located.

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

  `workflow.wave.enabled`
  : Whether Wave is enabled.

  `workflow.workDir`
  : The directory where task temporary files are stored.

(stdlib-functions)=

### Functions

The following functions are available in Nextflow scripts:

`branchCriteria( closure )`
: Create a branch criteria to use with the {ref}`operator-branch` operator.

`error( message = null )`
: Throw a script runtime error with an optional error message.

`exit( exitCode = 0, message = null )`
: :::{deprecated} 22.06.0-edge
  Use `error()` instead
  :::
: Stop the pipeline execution and return an exit code and optional error message.

`file( filePattern, [options] )`
: Get one or more files from a path or glob pattern. Returns a [Path](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/nio/file/Path.html) or list of Paths if there are multiple files.

: The following options are available:

  `checkIfExists`
  : When `true`, throws an exception if the specified path does not exist in the file system (default: `false`)

  `followLinks`
  : When `true`, follows symbolic links when traversing a directory tree, otherwise treats them as files (default: `true`)

  `glob`
  : When `true`, interprets characters `*`, `?`, `[]` and `{}` as glob wildcards, otherwise handles them as normal characters (default: `true`)

  `hidden`
  : When `true`, includes hidden files in the resulting paths (default: `false`)

  `maxDepth`
  : Maximum number of directory levels to visit (default: *no limit*)

  `type`
  : Type of paths returned, can be `'file'`, `'dir'` or `'any'` (default: `'file'`)

: See also: {ref}`Channel.fromPath <channel-path>`.

`files( filePattern, [options] )`
: Convenience method for `file()` that always returns a list.

`groupKey( key, size )`
: Create a grouping key to use with the {ref}`operator-grouptuple` operator.

`multiMapCriteria( closure )`
: Create a multi-map criteria to use with the {ref}`operator-multiMap` operator.

`sendMail( params )`
: Send an email. See {ref}`mail-page`.

`tuple( collection )`
: Create a tuple object from the given collection.

`tuple( ... args )`
: Create a tuple object from the given arguments.

`workflow.onComplete( closure )`
: Define an action to take when the workflow completes (whether successful or not). Refer to the `workflow` implicit variable to see which additional properties are available in the completion handler.

`workflow.onError( closure )`
: Define an action to take if the workflow is terminated due to a runtime error or task failure. Refer to the `workflow` implicit variable to see which additional properties are available in the error handler.

(stdlib-default-imports)=

## Default imports

The following classes are imported by default in Nextflow scripts:

- `java.io.*`
- `java.lang.*`
- `java.math.BigDecimal`
- `java.math.BigInteger`
- `java.net.*`
- `java.nio.file.Path`
- `java.util.*`
- `groovy.lang.*`
- `groovy.util.*`

## Channel

The `Channel` class provides the channel factory methods. See {ref}`channel-factory` for more information.

(stdlib-types-duration)=

## Duration

A `Duration` represents some duration of time.

You can create a duration by adding a time unit suffix to an integer, e.g. `1.h`. The following suffixes are available:

| Unit                            | Description  |
| ------------------------------- | ------------ |
| `ms`, `milli`, `millis`         | Milliseconds |
| `s`, `sec`, `second`, `seconds` | Seconds      |
| `m`, `min`, `minute`, `minutes` | Minutes      |
| `h`, `hour`, `hours`            | Hours        |
| `d`, `day`, `days`              | Days         |

You can also create a duration with `Duration.of()`:

```nextflow
// integer value (milliseconds)
oneSecond = Duration.of(1000)

// simple string value
oneHour = Duration.of('1h')

// complex string value
complexDuration = Duration.of('1day 6hours 3minutes 30seconds')
```

Durations can be compared like numbers, and they support basic arithmetic operations:

```nextflow
a = 1.h
b = 2.h

assert a < b
assert a + a == b
assert b - a == a
assert a * 2 == b
assert b / 2 == a
```

The following methods are available for a `Duration` object:

`getDays()`, `toDays()`
: Get the duration value in days (rounded down).

`getHours()`, `toHours()`
: Get the duration value in hours (rounded down).

`getMillis()`, `toMillis()`
: Get the duration value in milliseconds.

`getMinutes()`, `toMinutes()`
: Get the duration value in minutes (rounded down).

`getSeconds()`, `toSeconds()`
: Get the duration value in seconds (rounded down).

(stdlib-types-memoryunit)=

## MemoryUnit

A `MemoryUnit` represents a quantity of bytes.

You can create a memory unit by adding a unit suffix to an integer, e.g. `1.GB`. The following suffixes are available:

| Unit | Description |
| ---- | ----------- |
| `B`  | Bytes       |
| `KB` | Kilobytes   |
| `MB` | Megabytes   |
| `GB` | Gigabytes   |
| `TB` | Terabytes   |
| `PB` | Petabytes   |
| `EB` | Exabytes    |
| `ZB` | Zettabytes  |

:::{note}
Technically speaking, a kilobyte is equal to 1000 bytes, whereas 1024 bytes is called a "kibibyte" and abbreviated as "KiB", and so on for the other units. In practice, however, kilobyte is commonly understood to mean 1024 bytes, and Nextflow follows this convention in its implementation as well as this documentation.
:::

You can also create a memory unit with `MemoryUnit.of()`:

```nextflow
// integer value (bytes)
oneKilobyte = MemoryUnit.of(1024)

// string value
oneGigabyte = MemoryUnit.of('1 GB')
```

Memory units can be compared like numbers, and they support basic arithmetic operations:

```nextflow
a = 1.GB
b = 2.GB

assert a < b
assert a + a == b
assert b - a == a
assert a * 2 == b
assert b / 2 == a
```

The following methods are available for a `MemoryUnit` object:

`getBytes()`, `toBytes()`
: Get the memory value in bytes (B).

`getGiga()`, `toGiga()`
: Get the memory value in gigabytes (rounded down), where 1 GB = 1024 MB.

`getKilo()`, `toKilo()`
: Get the memory value in kilobytes (rounded down), where 1 KB = 1024 B.

`getMega()`, `toMega()`
: Get the memory value in megabytes (rounded down), where 1 MB = 1024 KB.

`toUnit( unit )`
: Get the memory value in terms of a given unit (rounded down). The unit can be one of: `'B'`, `'KB'`, `'MB'`, `'GB'`, `'TB'`, `'PB'`, `'EB'`, `'ZB'`.

(stdlib-types-path)=

## Path

The `file()` method returns a [Path](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/nio/file/Path.html), so any method defined for Path can also be used in a Nextflow script.

### Getting attributes

The following methods are useful for getting attributes of a file:

`exists()`
: Returns `true` if the file exists.

`getBaseName()`
: Gets the file name without its extension, e.g. `/some/path/file.tar.gz` -> `file.tar`.

`getExtension()`
: Gets the file extension, e.g. `/some/path/file.txt` -> `txt`.

`getName()`
: Gets the file name, e.g. `/some/path/file.txt` -> `file.txt`.

`getSimpleName()`
: Gets the file name without any extension, e.g. `/some/path/file.tar.gz` -> `file`.

`getParent()`
: Gets the file parent path, e.g. `/some/path/file.txt` -> `/some/path`.

`getScheme()`
: Gets the file URI scheme, e.g. `s3://some-bucket/foo.txt` -> `s3`.

`isDirectory()`
: Returns `true` if the file is a directory.

`isEmpty()`
: Returns `true` if the file is empty or does not exist.

`isFile()`
: Returns `true` if it is a regular file (i.e. not a directory).

`isHidden()`
: Returns `true` if the file is hidden.

`isLink()`
: Returns `true` if the file is a symbolic link.

`lastModified()`
: Returns the file last modified timestamp in Unix time (i.e. milliseconds since January 1, 1970).

`size()`
: Gets the file size in bytes.

`toUriString()`
: Gets the file path along with the protocol scheme:
  ```nextflow
  def ref = file('s3://some-bucket/foo.txt')

  assert ref.toString() == '/some-bucket/foo.txt'
  assert "$ref" == '/some-bucket/foo.txt'
  assert ref.toUriString() == 's3://some-bucket/foo.txt'
  ```

### Reading

The following methods are available for reading files:

`eachByte( closure )`
: Iterates over the file byte by byte, applying the specified {ref}`closure <script-closure>`.

`eachLine( closure )`
: Iterates over the file line by line, applying the specified {ref}`closure <script-closure>`.

`getBytes()`
: Returns the file content as a byte array.

`getText()`
: Returns the file content as a string value.

`newInputStream()`
: Returns an [InputStream](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/io/InputStream.html) object to read a binary file.

`newReader()`
: Returns a [Reader](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/io/Reader.html) object to read a text file.

`readLines()`
: Reads the file line by line and returns the content as a list of strings.

`withInputStream( closure )`
: Opens a file for reading and lets you access it with an [InputStream](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/io/InputStream.html) object.

`withReader( closure )`
: Opens a file for reading and lets you access it with a [Reader](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/io/Reader.html) object.

### Writing

The following methods are available for writing to files:

`append( text )`
: Appends a string value to a file without replacing existing content.

`newOutputStream()`
: Creates an [OutputStream](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/io/OutputStream.html) object that allows you to write binary data to a file.

`newPrintWriter()`
: Creates a [PrintWriter](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/io/PrintWriter.html) object that allows you to write formatted text to a file.

`newWriter()`
: Creates a [Writer](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/io/Writer.html) object that allows you to save text data to a file.

`setBytes( bytes )`
: Writes a byte array to a file. Equivalent to setting the `bytes` property.

`setText( text )`
: Writes a string value to a file. Equivalent to setting the `text` property.

`withOutputStream( closure )`
: Applies the specified closure to an [OutputStream](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/io/OutputStream.html) object, closing it when finished.

`withPrintWriter( closure )`
: Applies the specified closure to a [PrintWriter](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/io/PrintWriter.html) object, closing it when finished.

`withWriter( closure )`
: Applies the specified closure to a [Writer](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/io/Writer.html) object, closing it when finished.

`write( text )`
: Writes a string to a file, replacing any existing content.

### Filesystem operations

The following methods are available for manipulating files and directories in a filesystem:

`copyTo( target )`
: Copies a source file or directory to a target file or directory.

: *When copying a file to another file:* if the target file already exists, it will be replaced.

  ```nextflow
  file('/some/path/my_file.txt').copyTo('/another/path/new_file.txt')
  ```

: *When copying a file to a directory:* the file will be copied into the directory, replacing any file with the same name.

  ```nextflow
  file('/some/path/my_file.txt').copyTo('/another/path')
  ```

: *When copying a directory to another directory:* if the target directory already exists, the source directory will be copied into the target directory, replacing any sub-directory with the same name. If the target path does not exist, it will be created automatically.

  ```nextflow
  file('/any/dir_a').moveTo('/any/dir_b')
  ```

  The result of the above example depends on the existence of the target directory. If the target directory exists, the source is moved into the target directory, resulting in the path `/any/dir_b/dir_a`. If the target directory does not exist, the source is just renamed to the target name, resulting in the path `/any/dir_b`.

: :::{note}
  The `copyTo()` method follows the semantics of the Linux command `cp -r <source> <target>`, with the following caveat: while Linux tools often treat paths ending with a slash (e.g. `/some/path/name/`) as directories, and those not (e.g. `/some/path/name`) as regular files, Nextflow (due to its use of the Java files API) views both of these paths as the same file system object. If the path exists, it is handled according to its actual type (i.e. as a regular file or as a directory). If the path does not exist, it is treated as a regular file, with any missing parent directories created automatically.
  :::

`delete()`
: Deletes the file or directory at the given path, returning `true` if the operation succeeds, and `false` otherwise:

  ```nextflow
  myFile = file('some/file.txt')
  result = myFile.delete()
  println result ? "OK" : "Cannot delete: $myFile"
  ```

  If a directory is not empty, it will not be deleted and `delete()` will return `false`.

`deleteDir()`
: Deletes a directory and all of its contents.

  ```nextflow
  file('any/path').deleteDir()
  ```

`getPermissions()`
: Returns a file's permissions using the [symbolic notation](http://en.wikipedia.org/wiki/File_system_permissions#Symbolic_notation), e.g. `'rw-rw-r--'`.

`list()`
: Returns the first-level elements (files and directories) of a directory as a list of strings.

`listFiles()`
: Returns the first-level elements (files and directories) of a directory as a list of Paths.

`mkdir()`
: Creates a directory at the given path, returning `true` if the directory is created successfully, and `false` otherwise:

  ```nextflow
  myDir = file('any/path')
  result = myDir.mkdir()
  println result ? "OK" : "Cannot create directory: $myDir"
  ```

  If the parent directories do not exist, the directory will not be created and `mkdir()` will return `false`.

`mkdirs()`
: Creates a directory at the given path, including any nonexistent parent directories:

  ```nextflow
  file('any/path').mkdirs()
  ```

`mklink( linkName, [options] )`
: Creates a *filesystem link* to a given path:

  ```nextflow
  myFile = file('/some/path/file.txt')
  myFile.mklink('/user/name/link-to-file.txt')
  ```

  Available options:

  `hard`
  : When `true`, creates a *hard* link, otherwise creates a *soft* (aka *symbolic*) link (default: `false`).

  `overwrite`
  : When `true`, overwrites any existing file with the same name, otherwise throws a [FileAlreadyExistsException](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/nio/file/FileAlreadyExistsException.html) (default: `false`).

`moveTo( target )`
: Moves a source file or directory to a target file or directory. Follows the same semantics as `copyTo()`.

`renameTo( target )`
: Rename a file or directory:

  ```nextflow
  file('my_file.txt').renameTo('new_file_name.txt')
  ```

`setPermissions( permissions )`
: Sets a file's permissions using the [symbolic notation](http://en.wikipedia.org/wiki/File_system_permissions#Symbolic_notation):

  ```nextflow
  myFile.setPermissions('rwxr-xr-x')
  ```

`setPermissions( owner, group, other )`
: Sets a file's permissions using the [numeric notation](http://en.wikipedia.org/wiki/File_system_permissions#Numeric_notation), i.e. as three digits representing the **owner**, **group**, and **other** permissions:

  ```nextflow
  myFile.setPermissions(7,5,5)
  ```

The following methods are available for listing and traversing directories:

`eachDir( closure )`
: Iterates through first-level directories only.

`eachDirMatch( nameFilter, closure )`
: Iterates through directories whose names match the given filter.

`eachDirRecurse( closure )`
: Iterates through directories depth-first (regular files are ignored).

`eachFile( closure )`
: Iterates through first-level files and directories.

`eachFileMatch( nameFilter, closure )`
: Iterates through files and directories whose names match the given filter.

`eachFileRecurse( closure )`
: Iterates through files and directories depth-first.

Refer to the [Groovy documentation](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html) for more details.

### Splitting files

The following methods are available for splitting and counting the records in files:

`countFasta()`
: Counts the number of records in a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file. See the {ref}`operator-splitfasta` operator for available options.

`countFastq()`
: Counts the number of records in a [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file. See the {ref}`operator-splitfastq` operator for available options.

`countJson()`
: Counts the number of records in a JSON file. See the {ref}`operator-splitjson` operator for available options.

`countLines()`
: Counts the number of lines in a text file. See the {ref}`operator-splittext` operator for available options.

`splitCsv()`
: Splits a CSV file into a list of records. See the {ref}`operator-splitcsv` operator for available options.

`splitFasta()`
: Splits a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file into a list of records. See the {ref}`operator-splitfasta` operator for available options.

`splitFastq()`
: Splits a [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file into a list of records. See the {ref}`operator-splitfastq` operator for available options.

`splitJson()`
: Splits a JSON file into a list of records. See the {ref}`operator-splitjson` operator for available options.

`splitText()`
: Splits a text file into a list of lines. See the {ref}`operator-splittext` operator for available options.

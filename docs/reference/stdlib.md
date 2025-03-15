(stdlib-page)=

# Standard library

This page describes the built-in constants, functions, and types provided by Nextflow.

## Globals

(stdlib-constants)=

### Constants

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

(stdlib-functions)=

### Functions

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

## Groovy and Java classes

Any Groovy or Java class can be used in a Nextflow script. The following classes are imported by default:

- `groovy.lang.*`
- `groovy.util.*`
- `java.io.*`
- `java.lang.*`
- `java.math.BigDecimal`
- `java.math.BigInteger`
- `java.net.*`
- `java.util.*`

All other classes must be referenced by their fully-qualified name:

```nextflow
def vals = [1, 2, 3]
println groovy.json.JsonOutput.toJson(vals)
```

The following sections describe the standard types provided by Nextflow.

## Bag\<E\>

A bag is an unordered collection.

The following operations are supported for bags:

`+ : (Bag<E>, Bag<E>) -> Bag<E>`
: Concatenates two bags.

`in, !in : (E, Bag<E>) -> boolean`
: Given a value and a bag, returns `true` if the bag contains the value (or not).

Bags inherit all methods from the [Iterable](#iterablee) trait.

## Channel\<E\>

See {ref}`channel-page` for an overview of channels. See {ref}`channel-factory` and {ref}`operator-page` for the available functions for creating and manipulating channels.

(stdlib-types-duration)=

## Duration

A Duration represents a duration of time with millisecond precision.

A Duration can be created by adding a unit suffix to an integer (e.g. `1.h`), or more explicitly with `Duration.of()`:

```nextflow
// integer with suffix
oneDay = 24.h

// integer value (milliseconds)
oneSecond = Duration.of(1000)

// simple string value
oneHour = Duration.of('1h')

// complex string value
complexDuration = Duration.of('1day 6hours 3minutes 30seconds')
```

The following suffixes are available:

| Unit                            | Description  |
| ------------------------------- | ------------ |
| `ms`, `milli`, `millis`         | Milliseconds |
| `s`, `sec`, `second`, `seconds` | Seconds      |
| `m`, `min`, `minute`, `minutes` | Minutes      |
| `h`, `hour`, `hours`            | Hours        |
| `d`, `day`, `days`              | Days         |

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

The following methods are available for a Duration:

`toDays() -> long`
: Get the duration value in days (rounded down).

`toHours() -> long`
: Get the duration value in hours (rounded down).

`toMillis() -> long`
: Get the duration value in milliseconds.

`toMinutes() -> long`
: Get the duration value in minutes (rounded down).

`toSeconds() -> long`
: Get the duration value in seconds (rounded down).

:::{note}
These methods are also available as `getDays()`, `getHours()`, `getMillis()`, `getMinutes()`, and `getSeconds()`.
:::

## Iterable\<E\>

An iterable is a trait shared by collections that support iteration. The following types are iterables:

- [Bag\<E\>](#bage)
- [List\<E\>](#liste)
- [Set\<E\>](#sete)

Values of these types can be passed to any parameter of type `Iterable`, and they can use all of the methods described below.

The following methods are available for iterables:

`any( condition: (E) -> boolean ) -> boolean`
: Returns `true` if any element in the iterable satisfies the given condition.

`collect( transform: (E) -> R ) -> Iterable<R>`
: Returns a new iterable with each element transformed by the given closure.

`collectMany( transform: (E) -> Iterable<R> ) -> Iterable<R>`
: Transforms each element in the iterable into a collection with the given closure and concatenates the resulting collections into a list.

`contains( value: E ) -> boolean`
: Returns `true` if the iterable contains the given value.

`each( action: (E) -> () )`
: Invoke the given closure for each element in the iterable.

`every( condition: (E) -> boolean ) -> boolean`
: Returns `true` if every element in the iterable satisfies the given condition.

`findAll( condition: (E) -> boolean ) -> Iterable<E>`
: Returns the elements in the iterable that satisfy the given condition.

`groupBy( transform: (E) -> K ) -> Map<K,Iterable<E>>`
: Collect the elements of an iterable into groups based on a matching key. The closure should return the key for a given element.

`inject( accumulator: (E,E) -> E ) -> E`
: Apply the given accumulator to each element in the iterable and return the final accumulated value. The closure should accept two parameters, corresponding to the current accumulated value and the current iterable element, and return the next accumulated value.
: The first element from the iterable is used as the initial accumulated value.

`inject( initialValue: R, accumulator: (R,E) -> R ) -> R`
: Apply the given accumulator to each element in the iterable and return the final accumulated value. The closure should accept two parameters, corresponding to the current accumulated value and the current iterable element, and return the next accumulated value.

`isEmpty() -> boolean`
: Returns `true` if the iterable is empty.

`join( separator: String = '' ) -> String`
: Concatenates the string representation of each element in the iterable, with the given string as the separator between each element.

`max() -> E`
: Returns the maximum element in the iterable.

**`max( comparator: (E) -> R ) -> E`**

`max( comparator: (E,E) -> int ) -> E`
: Returns the maximum element in the iterable according to the given closure.
: The closure should follow the same semantics as the closure parameter of `toSorted()`.

`min() -> E`
: Returns the maximum element in the iterable.

**`min( comparator: (E) -> R ) -> E`**

`min( comparator: (E,E) -> int ) -> E`
: Returns the maximum element in the iterable according to the given closure.
: The closure should follow the same semantics as the closure parameter of `toSorted()`.

`size() -> int`
: Returns the number of elements in the iterable.

`sum() -> E`
: Returns the sum of the elements in the iterable. The elements should support addition (`+`).

`sum( transform: (E) -> R ) -> R`
: Transforms each element in the iterable with the given closure and returns the sum. The values returned by the closure should support addition (`+`).

`toList() -> List<E>`
: Converts the iterable to a list.
: :::{danger}
  Converting an unordered collection to a list can lead to non-deterministic behavior. Consider using `toSorted()` instead to ensure a deterministic ordering. See {ref}`cache-nondeterministic-inputs` for more information.
  :::

`toSet() -> Set<E>`
: Converts the iterable to a set. Duplicate elements are excluded.

`toSorted() -> List<E>`
: Returns a sorted list of the iterable's elements.

`toSorted( comparator: (E) -> R ) -> List<E>`
: Returns the iterable as a list sorted according to the given closure.
: The closure should accept one parameter and transform each element into the value that will be used for comparisons.

`toSorted( comparator: (E,E) -> int ) -> List<E>`
: Returns the iterable as a list sorted according to the given closure.
: The closure should accept two parameters and return a negative integer, zero, or a positive integer to denote whether the first argument is less than, equal to, or greater than the second.

`toUnique() -> Iterable<E>`
: Returns a shallow copy of the iterable with duplicate elements excluded.

**`toUnique( comparator: (E) -> R ) -> Iterable<E>`**

`toUnique( comparator: (E,E) -> int ) -> Iterable<E>`
: Returns a shallow copy of the iterable with duplicate elements excluded.
: The closure should follow the same semantics as the closure parameter of `toSorted()`.

:::{note}
Iterables in Nextflow are backed by the [Java](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/lang/Iterable.html) and [Groovy](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/lang/Iterable.html) standard libraries, which may expose additional methods. Only methods which are recommended for use in Nextflow are documented here.
:::

(stdlib-types-list)=

## List\<E\>

A list is an ordered collection of elements. See {ref}`script-list` for an overview of lists.

The following operations are supported for lists:

`+ : (List<E>, List<E>) -> List<E>`
: Concatenates two lists.

`* : (List<E>, int) -> List<E>`
: Given a list and an integer *n*, repeats the list *n* times.

`[] : (List<E>, int) -> E`
: Given a list and an index, returns the element at the given index in the list, or `null` if the index is out of range.

`in, !in : (E, List<E>) -> boolean`
: Given a value and a list, returns `true` if the list contains the value (or not).

The following methods are available for a list:

`collate( size: int, keepRemainder: boolean = true ) -> List<List<E>>`
: Collates the list into a list of sub-lists of length `size`. If `keepRemainder` is `true`, any remaining elements are included as a partial sub-list, otherwise they are excluded.

: For example:
  ```nextflow
  assert [1, 2, 3, 4, 5, 6, 7].collate(3)        == [[1, 2, 3], [4, 5, 6], [7]]
  assert [1, 2, 3, 4, 5, 6, 7].collate(3, false) == [[1, 2, 3], [4, 5, 6]]
  ```

`collate( size: int, step: int, keepRemainder: boolean = true ) -> List<List<E>>`
: Collates the list into a list of sub-lists of length `size`, stepping through the list `step` elements for each sub-list. If `keepRemainder` is `true`, any remaining elements are included as a partial sub-list, otherwise they are excluded.

: For example:
  ```nextflow
  assert [1, 2, 3, 4].collate(3, 1)        == [[1, 2, 3], [2, 3, 4], [3, 4], [4]]
  assert [1, 2, 3, 4].collate(3, 1, false) == [[1, 2, 3], [2, 3, 4]]
  ```

`find( condition: (E) -> boolean ) -> E`
: Returns the first element in the list that satisfies the given condition.

`first() -> E`
: Returns the first element in the list. Raises an error if the list is empty.

`getIndices() -> List<Integer>`
: Returns the list of integers from 0 to *n - 1*, where *n* is the number of elements in the list.

`head() -> E`
: Equivalent to `first()`.

`indexOf( value: E ) -> int`
: Returns the index of the first occurrence of the given value in the list, or -1 if the list does not contain the value.

`init() -> List<E>`
: Returns a shallow copy of the list with the last element excluded.

`last() -> E`
: Returns the last element in the list. Raises an error if the list is empty.

`reverse() -> List<E>`
: Returns a shallow copy of the list with the elements reversed.

`subList( fromIndex: int, toIndex: int ) -> List<E>`
: Returns the portion of the list between the given `fromIndex` (inclusive) and `toIndex` (exclusive).

`tail() -> List<E>`
: Returns a shallow copy of the list with the first element excluded.

`take( n: int ) -> List<E>`
: Returns the first *n* elements of the list.

`takeWhile( condition: (E) -> boolean ) -> List<E>`
: Returns the longest prefix of the list where each element satisfies the given condition.

`withIndex() -> List<(E,Integer)>`
: Returns a list of 2-tuples corresponding to the value and index of each element in the list.

Lists inherit all methods from the [Iterable](#iterablee) trait.

:::{note}
Lists in Nextflow are backed by the [Java](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/util/List.html) and [Groovy](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/List.html) standard libraries, which may expose additional methods. Only methods which are recommended for use in Nextflow are documented here.
:::

(stdlib-types-map)=

## Map\<K,V\>

A map associates or "maps" keys to values. Each key can map to at most one value -- a map cannot contain duplicate keys. See {ref}`script-map` for an overview of maps.

The following operations are supported for maps:

`+ : (Map<K,V>, Map<K,V>) -> Map<K,V>`
: Concatenates two maps. If a key exists in both maps, the mapping from the right-hand side is used.

`[] : (Map<K,V>, K) -> V`
: Given a map and a key, returns the value for the given key in the map, or `null` if the key is not in the map.

`in, !in : (K, Map<K,V>) -> boolean`
: Given a key and a map, returns `true` if the map contains the key (or not).

The following methods are available for a map:

`any( condition: (K,V) -> boolean ) -> boolean`
: Returns `true` if any key-value pair in the map satisfies the given condition. The closure should accept two parameters corresponding to the key and value of an entry.

`containsKey( key: K ) -> boolean`
: Returns `true` if the map contains a mapping for the given key.

`containsValue( value: V ) -> boolean`
: Returns `true` if the map maps one or more keys to the given value.

`each( action: (K,V) -> () )`
: Invoke the given closure for each key-value pair in the map. The closure should accept two parameters corresponding to the key and value of an entry.

`entrySet() -> Set<(K,V)>`
: Returns a set of the key-value pairs in the map.

`every( condition: (K,V) -> boolean ) -> boolean`
: Returns `true` if every key-value pair in the map satisfies the given condition. The closure should accept two parameters corresponding to the key and value of an entry.

`isEmpty() -> boolean`
: Returns `true` if the map is empty.

`keySet() -> Set<K>`
: Returns a set of the keys in the map.

`size() -> int`
: Returns the number of key-value pairs in the map.

`subMap( keys: Iterable<K> ) -> Map<K,V>`
: Returns a sub-map containing the given keys.

`values() -> Bag<V>`
: Returns a collection of the values in the map.

:::{note}
Maps in Nextflow are backed by the [Java](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/util/Map.html) and [Groovy](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html) standard libraries, which may expose additional methods. Only methods which are recommended for use in Nextflow are documented here.
:::

(stdlib-types-memoryunit)=

## MemoryUnit

A MemoryUnit represents a quantity of bytes.

A MemoryUnit can be created by adding a unit suffix to an integer (e.g. `1.GB`), or more explicitly with `MemoryUnit.of()`:

```nextflow
// integer with suffix
oneMegabyte = 1.MB

// integer value (bytes)
oneKilobyte = MemoryUnit.of(1024)

// string value
oneGigabyte = MemoryUnit.of('1 GB')
```

The following suffixes are available:

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

`toBytes() -> long`
: Get the memory value in bytes (B).

`toGiga() -> long`
: Get the memory value in gigabytes (rounded down), where 1 GB = 1024 MB.

`toKilo() -> long`
: Get the memory value in kilobytes (rounded down), where 1 KB = 1024 B.

`toMega() -> long`
: Get the memory value in megabytes (rounded down), where 1 MB = 1024 KB.

`toUnit( unit: String ) -> long`
: Get the memory value in terms of a given unit (rounded down). The unit can be one of: `'B'`, `'KB'`, `'MB'`, `'GB'`, `'TB'`, `'PB'`, `'EB'`, `'ZB'`.

:::{note}
These methods are also available as `getBytes()`, `getGiga()`, `getKilo()`, and `getMega()`.
:::

(stdlib-types-path)=

## Path

A Path is a handle for hierarchichal paths such as local files and directories, HTTP/FTP URLs, and object storage paths (e.g. Amazon S3).

The `file()` function can be used to get a Path for a given filename or URL:

```nextflow
def hello = file('hello.txt')
println hello.text
```

The `files()` function can be used to get a collection of Paths from a glob pattern:

```nextflow
def inputs = files('*.txt')
inputs.each { input ->
  println "${input.name}: ${input.text}"
}
```

The following sections describe the methods that are available for paths.

:::{note}
Paths in Nextflow are backed by the [Java](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/nio/file/Path.html) and [Groovy](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/nio/file/Path.html) standard libraries, which may expose additional methods. Only methods which are recommended for use in Nextflow are documented here.
:::

### Operations

The following operations are supported for paths:

`/ : (Path, String) -> Path`
: Resolve a relative file name against a directory path.

`<< : (Path, String)`
: Appends text to a file without replacing existing content. Equivalent to `append()`.

### Getting attributes

The following methods are useful for getting attributes of a file:

`exists() -> boolean`
: Returns `true` if the file exists.

`getBaseName() -> String`
: Gets the file name without its extension, e.g. `/some/path/file.tar.gz` -> `file.tar`.

`getExtension() -> String`
: Gets the file extension, e.g. `/some/path/file.txt` -> `txt`.

`getName() -> String`
: Gets the file name, e.g. `/some/path/file.txt` -> `file.txt`.

`getSimpleName() -> String`
: Gets the file name without any extension, e.g. `/some/path/file.tar.gz` -> `file`.

`getParent() -> Path`
: Gets the file parent path, e.g. `/some/path/file.txt` -> `/some/path`.

`getScheme() -> String`
: Gets the file URI scheme, e.g. `s3://some-bucket/foo.txt` -> `s3`.

`isDirectory() -> boolean`
: Returns `true` if the file is a directory.

`isEmpty() -> boolean`
: Returns `true` if the file is empty or does not exist.

`isFile() -> boolean`
: Returns `true` if it is a regular file (i.e. not a directory).

`isHidden() -> boolean`
: Returns `true` if the file is hidden.

`isLink() -> boolean`
: Returns `true` if the file is a symbolic link.

`lastModified() -> long`
: Returns the file last modified timestamp in Unix time (i.e. milliseconds since January 1, 1970).

`size() -> long`
: Gets the file size in bytes.

`toUriString() -> String`
: Gets the file path along with the protocol scheme:
  ```nextflow
  def ref = file('s3://some-bucket/foo.txt')

  assert ref.toString() == '/some-bucket/foo.txt'
  assert "$ref" == '/some-bucket/foo.txt'
  assert ref.toUriString() == 's3://some-bucket/foo.txt'
  ```

### Reading

The following methods are available for reading files:

`eachByte( action: (byte) -> () )`
: Iterates over the file, applying the specified closure to each byte.

`eachLine( action: (String) -> () )`
: Iterates over the file, applying the specified closure to each line.

`getBytes() -> byte[]`
: Returns the file content as a byte array.

`getText() -> String`
: Returns the file content as a string.

`newInputStream() -> InputStream`
: Returns an [InputStream](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/io/InputStream.html) object to read a binary file.

`newReader() -> Reader`
: Returns a [Reader](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/io/Reader.html) object to read a text file.

`readLines() -> List<String>`
: Reads the file line by line and returns the content as a list of strings.

`withInputStream( action: (InputStream) -> () )`
: Opens a file for reading and lets you access it with an [InputStream](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/io/InputStream.html) object.

`withReader( action: (Reader) -> () )`
: Opens a file for reading and lets you access it with a [Reader](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/io/Reader.html) object.

### Writing

The following methods are available for writing to files:

`append( text: String )`
: Appends text to a file without replacing existing content.

`newOutputStream() -> OutputStream`
: Creates an [OutputStream](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/io/OutputStream.html) object that allows you to write binary data to a file.

`newPrintWriter() -> PrintWriter`
: Creates a [PrintWriter](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/io/PrintWriter.html) object that allows you to write formatted text to a file.

`newWriter() -> Writer`
: Creates a [Writer](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/io/Writer.html) object that allows you to save text data to a file.

`setBytes( bytes: byte[] )`
: Writes a byte array to a file. Equivalent to setting the `bytes` property.

`setText( text: String )`
: Writes text to a file. Equivalent to setting the `text` property.

`withOutputStream( action: (OutputStream) -> () )`
: Applies the specified closure to an [OutputStream](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/io/OutputStream.html) object, closing it when finished.

`withPrintWriter( action: (PrintWriter) -> () )`
: Applies the specified closure to a [PrintWriter](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/io/PrintWriter.html) object, closing it when finished.

`withWriter( action: (Writer) -> () )`
: Applies the specified closure to a [Writer](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/io/Writer.html) object, closing it when finished.

`write( text: String )`
: Writes a string to a file, replacing any existing content.

### Filesystem operations

The following methods are available for manipulating files and directories in a filesystem:

`copyTo( target: Path )`
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
  The `copyTo()` function follows the semantics of the Linux command `cp -r <source> <target>`, with the following caveat: while Linux tools often treat paths ending with a slash (e.g. `/some/path/name/`) as directories, and those not (e.g. `/some/path/name`) as regular files, Nextflow (due to its use of the Java files API) views both of these paths as the same file system object. If the path exists, it is handled according to its actual type (i.e. as a regular file or as a directory). If the path does not exist, it is treated as a regular file, with any missing parent directories created automatically.
  :::

`delete() -> boolean`
: Deletes the file or directory at the given path, returning `true` if the operation succeeds, and `false` otherwise:

  ```nextflow
  myFile = file('some/file.txt')
  result = myFile.delete()
  println result ? "OK" : "Cannot delete: $myFile"
  ```

  If a directory is not empty, it will not be deleted and `delete()` will return `false`.

`deleteDir() -> boolean`
: Deletes a directory and all of its contents.

  ```nextflow
  file('any/path').deleteDir()
  ```

`getPermissions() -> String`
: Returns a file's permissions using the [symbolic notation](http://en.wikipedia.org/wiki/File_system_permissions#Symbolic_notation), e.g. `'rw-rw-r--'`.

`list() -> List<String>`
: Returns the first-level elements (files and directories) of a directory as a list of strings.

`listFiles() -> List<Path>`
: Returns the first-level elements (files and directories) of a directory as a list of Paths.

`mkdir() -> boolean`
: Creates a directory at the given path, returning `true` if the directory is created successfully, and `false` otherwise:

  ```nextflow
  myDir = file('any/path')
  result = myDir.mkdir()
  println result ? "OK" : "Cannot create directory: $myDir"
  ```

  If the parent directories do not exist, the directory will not be created and `mkdir()` will return `false`.

`mkdirs() -> boolean`
: Creates a directory at the given path, including any nonexistent parent directories:

  ```nextflow
  file('any/path').mkdirs()
  ```

`mklink( linkName: String, [options] ) -> Path`
: Creates a *filesystem link* to a given path:

  ```nextflow
  myFile = file('/some/path/file.txt')
  myFile.mklink('/user/name/link-to-file.txt')
  ```

  Returns the Path of the link to create.

  Available options:

  `hard: boolean`
  : When `true`, creates a *hard* link, otherwise creates a *soft* (aka *symbolic*) link (default: `false`).

  `overwrite: boolean`
  : When `true`, overwrites any existing file with the same name, otherwise throws a [FileAlreadyExistsException](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/nio/file/FileAlreadyExistsException.html) (default: `false`).

`moveTo( target: Path )`
: Moves a source file or directory to a target file or directory. Follows the same semantics as `copyTo()`.

`renameTo( target: String ) -> boolean`
: Rename a file or directory:

  ```nextflow
  file('my_file.txt').renameTo('new_file_name.txt')
  ```

`setPermissions( permissions: String ) -> boolean`
: Sets a file's permissions using the [symbolic notation](http://en.wikipedia.org/wiki/File_system_permissions#Symbolic_notation):

  ```nextflow
  myFile.setPermissions('rwxr-xr-x')
  ```

`setPermissions( owner: int, group: int, other: int ) -> boolean`
: Sets a file's permissions using the [numeric notation](http://en.wikipedia.org/wiki/File_system_permissions#Numeric_notation), i.e. as three digits representing the **owner**, **group**, and **other** permissions:

  ```nextflow
  myFile.setPermissions(7,5,5)
  ```

The following methods are available for listing and traversing directories:

`eachDir( action: (Path) -> () )`
: Iterates through first-level directories only.

`eachDirMatch( nameFilter: String, action: (Path) -> () )`
: Iterates through directories whose names match the given filter.

`eachDirRecurse( action: (Path) -> () )`
: Iterates through directories depth-first (regular files are ignored).

`eachFile( action: (Path) -> () )`
: Iterates through first-level files and directories.

`eachFileMatch( nameFilter: String, action: (Path) -> () )`
: Iterates through files and directories whose names match the given filter.

`eachFileRecurse( action: (Path) -> () )`
: Iterates through files and directories depth-first.

### Splitting files

The following methods are available for splitting and counting the records in files:

`countFasta() -> long`
: Counts the number of records in a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file. See the {ref}`operator-splitfasta` operator for available options.

`countFastq() -> long`
: Counts the number of records in a [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file. See the {ref}`operator-splitfastq` operator for available options.

`countJson() -> long`
: Counts the number of records in a JSON file. See the {ref}`operator-splitjson` operator for available options.

`countLines() -> long`
: Counts the number of lines in a text file. See the {ref}`operator-splittext` operator for available options.

`splitCsv() -> List<?>`
: Splits a CSV file into a list of records. See the {ref}`operator-splitcsv` operator for available options.

`splitFasta() -> List<?>`
: Splits a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file into a list of records. See the {ref}`operator-splitfasta` operator for available options.

`splitFastq() -> List<?>`
: Splits a [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file into a list of records. See the {ref}`operator-splitfastq` operator for available options.

`splitJson() -> List<?>`
: Splits a JSON file into a list of records. See the {ref}`operator-splitjson` operator for available options.

`splitText() -> List<String>`
: Splits a text file into a list of lines. See the {ref}`operator-splittext` operator for available options.

(stdlib-types-set)=

## Set\<E\>

A set is an unordered collection that cannot contain duplicate elements.

As set literal can be created from a list:

```nextflow
[1, 2, 2, 3].toSet()
// -> [1, 2, 3]
```

The following operations are supported for sets:

`+ : (Set<E>, Iterable<E>) -> Set<E>`
: Given a set and an iterable, returns a new set containing the elements of both collections.

`- : (Set<E>, Iterable<E>) -> Set<E>`
: Given a set and an iterable, returns a shallow copy of the set minus the elements of the iterable.

`in, !in : (E, Set<E>) -> boolean`
: Given a value and a set, returns `true` if the set contains the value (or not).

The following methods are available for a set:

`intersect( right: Iterable<E> ) -> Set<E>`
: Returns the intersection of the set and the given iterable.

Sets inherit all methods from the [Iterable](#iterablee) trait.

:::{note}
Sets in Nextflow are backed by the [Java](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/util/Set.html) and [Groovy](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Set.html) standard libraries, which may expose additional methods. Only methods which are recommended for use in Nextflow are documented here.
:::

(stdlib-types-string)=

## String

A string is an immutable sequence of characters. See {ref}`script-string` for an overview of strings.

The following operations are supported for strings:

`+ : (String, String) -> String`
: Concatenates two strings.

`* : (String, int) -> String`
: Given a string and an integer *n*, repeats the string *n* times.

`[] : (String, int) -> char`
: Given a string and an index, returns the character at the given index in the string.

`in, !in : (String, String) -> boolean`
: Given a substring and a string, returns `true` if the substring occurs anywhere in the string (or not).

`~ : (String) -> Pattern`
: Creates a regular expression from a string.
: See [Pattern](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/util/regex/Pattern.html) in the Java standard library for more information.

`=~ : (String, String) -> Matcher`
: Given a string and a pattern, creates a matcher that is truthy if the pattern occurs anywhere in the string.
: See [Matcher](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/util/regex/Matcher.html) in the Java standard library for more information.

`==~ : (String, String) -> boolean`
: Given a string and a pattern, returns `true` if the string matches the pattern exactly.

The following methods are available for a string:

`endsWith( suffix: String ) -> boolean`
: Returns `true` if the string ends with the given suffix.

`execute() -> Process`
: Execute the string as a command. Returns a [Process](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/lang/Process.html) which provides the exit status and standard input/output/error of the executed command.

`indexOf( str: String ) -> int`
: Returns the index within the string of the first occurrence of the given substring. Returns -1 if the string does not contain the substring.

`indexOf( str: String, fromIndex: int ) -> int`
: Returns the index within the string of the first occurrence of the given substring, starting the search at the given index. Returns -1 if the string does not contain the substring.

`isBlank() -> boolean`
: Returns `true` if the string is empty or contains only whitespace characters.

`isEmpty() -> boolean`
: Returns `true` if the string is empty (i.e. `length()` is 0).

`isFloat() -> boolean`
: Returns `true` if the string can be parsed as a floating-point number.

`isInteger() -> boolean`
: Returns `true` if the string can be parsed as an integer.

`lastIndexOf( str: String ) -> int`
: Returns the index within the string of the last occurrence of the given substring. Returns -1 if the string does not contain the substring.

`lastIndexOf( str: String, fromIndex: int ) -> int`
: Returns the index within the string of the last occurrence of the given substring, searching backwards starting at the given index. Returns -1 if the string does not contain the substring.

`length() -> int`
: Returns the length of the string.

`md5() -> String`
: Returns the MD5 checksum of the string.

`replace( target: String, replacement: String ) -> String`
: Returns a new string in which each occurrence of the target string is replaced with the given replacement string.

`replaceAll( regex: String, replacement: String ) -> String`
: Returns a new string in which each occurrence of the given regular expression is replaced with the given replacement string.

`replaceFirst( regex: String, replacement: String ) -> String`
: Returns a new string in which the first occurrence of the given regular expression is replaced with the given replacement string.

`sha256() -> String`
: Returns the SHA-256 checksum of the string.

`startsWith( prefix: String ) -> boolean`
: Returns `true` if the string ends with the given prefix.

`strip() -> String`
: Returns a copy of the string with all leading and trailing whitespace removed.

`stripIndent() -> String`
: Returns a copy of the string with leading spaces on each line removed.
: The number of spaces to remove is determined by the line with the least number of leading spaces, excluding lines with only whitespace.

`stripLeading() -> String`
: Returns a copy of the string with all leading whitespace removed.

`stripTrailing() -> String`
: Returns a copy of the string with all trailing whitespace removed.

`substring( beginIndex: int ) -> String`
: Returns a substring of this string.

`substring( beginIndex: int, endIndex: int ) -> String`
: Returns a substring of this string.

`toBoolean() -> Boolean`
: Returns `true` if the trimmed string is "true", "y", or "1" (ignoring case).

`toFloat() -> Float`
: Parses the string into a floating-point number.

`toInteger() -> Integer`
: Parses the string into an integer.

`toLowerCase() -> String`
: Returns a copy of this string with all characters converted to lower case.

`toUpperCase() -> String`
: Returns a copy of this string with all characters converted to upper case.

`tokenize( delimiters: String ) -> List<String>`
: Splits the string into a list of substrings using the given delimiters. Each character in the delimiter string is treated as a separate delimiter.

:::{note}
Strings in Nextflow are backed by the [Java](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/lang/String.html) and [Groovy](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/lang/String.html) standard libraries, which may expose additional methods. Only methods which are recommended for use in Nextflow are documented here.
:::

(stdlib-types-versionnumber)=

## VersionNumber

A VersionNumber represents a semantic or calendar version number.

The following methods are available for a VersionNumber:

`getMajor() -> String`
: Get the major version number, i.e. the first version component.

`getMinor() -> String`
: Get the minor version number, i.e. the second version component.

`getPatch() -> String`
: Get the patch version number, i.e. the third version component.

`matches( condition: String ) -> boolean`
: Check whether the version satisfies a version requirement.

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

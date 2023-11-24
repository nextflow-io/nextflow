(channel-page)=

# Channels

Nextflow is based on the dataflow programming model in which processes communicate through channels.

A channel has two major properties:

1. Sending a message is an *asynchronous* (i.e. non-blocking) operation, which means the sender doesn't have to wait for the receiving process.
2. Receiving a message is a *synchronous* (i.e. blocking) operation, which means the receiving process must wait until a message has arrived.

(channel-types)=

## Channel types

In Nextflow there are two kinds of channels: *queue channels* and *value channels*.

(channel-type-queue)=

### Queue channel

A *queue channel* is a non-blocking unidirectional FIFO queue connecting a *producer* process (i.e. outputting a value)
to a consumer process, or an operators.

A queue channel can be created by factory methods ([of](#of), [fromPath](#frompath), etc), operators ({ref}`operator-map`, {ref}`operator-flatmap`, etc), and processes (see {ref}`Process outputs <process-output>`).

(channel-type-value)=

### Value channel

A *value channel* can be bound (i.e. assigned) with one and only one value, and can be consumed any number of times by
a process or an operator.

A value channel can be created with the [value](#value) factory method or by any operator that produces a single value
({ref}`operator-first`, {ref}`operator-collect`, {ref}`operator-reduce`, etc). Additionally, a process will emit value
channels if it is invoked with all value channels, including simple values which are implicitly wrapped in a value channel.

For example:

```groovy
process foo {
  input:
  val x

  output:
  path 'x.txt'

  """
  echo $x > x.txt
  """
}

workflow {
  result = foo(1)
  result.view { "Result: ${it}" }
}
```

In the above example, since the `foo` process is invoked with a simple value instead of a channel, the input is implicitly
wrapped in a value channel, and the output is also emitted as a value channel.

See also: {ref}`process-multiple-input-channels`.

(channel-factory)=

## Channel factories

Channels may be created explicitly using the following channel factory methods.

:::{versionadded} 20.07.0
`channel` was introduced as an alias of `Channel`, allowing factory methods to be specified as `channel.of()` or
`Channel.of()`, and so on.
:::

(channel-empty)=

### empty

The `channel.empty` factory method, by definition, creates a channel that doesn't emit any value.

See also: {ref}`operator-ifempty`.

(channel-from)=

### from

:::{deprecated} 19.09.0-edge
Use [channel.of](#of) or [channel.fromList](#fromlist) instead.
:::

The `channel.from` method allows you to create a channel emitting any sequence of values that are specified as the method argument, for example:

```groovy
ch = channel.from( 1, 3, 5, 7 )
ch.subscribe { println "value: $it" }
```

The first line in this example creates a variable `ch` which holds a channel object. This channel emits the values specified as a parameter in the `from` method. Thus the second line will print the following:

```
value: 1
value: 3
value: 5
value: 7
```

The following example shows how to create a channel from a *range* of numbers or strings:

```groovy
zeroToNine = channel.from( 0..9 )
strings = channel.from( 'A'..'Z' )
```

:::{note}
When the `channel.from` argument is an object implementing the (Java) [Collection](http://docs.oracle.com/javase/7/docs/api/java/util/Collection.html) interface, the resulting channel emits the collection entries as individual items.
:::

Thus the following two declarations produce an identical result even though in the first case the items are specified as multiple arguments while in the second case as a single list object argument:

```groovy
channel.from( 1, 3, 5, 7, 9 )
channel.from( [1, 3, 5, 7, 9] )
```

But when more than one argument is provided, they are always managed as *single* emissions. Thus, the following example creates a channel emitting three entries each of which is a list containing two elements:

```groovy
channel.from( [1, 2], [5,6], [7,9] )
```

(channel-fromlist)=

### fromList

:::{versionadded} 19.10.0
:::

The `channel.fromList` method allows you to create a channel emitting the values provided as a list of elements, for example:

```groovy
channel
    .fromList( ['a', 'b', 'c', 'd'] )
    .view { "value: $it" }
```

Prints:

```
value: a
value: b
value: c
value: d
```

See also: [channel.of](#of) factory method.

(channel-path)=

### fromPath

You can create a channel emitting one or more file paths by using the `channel.fromPath` method and specifying a path
string as an argument. For example:

```groovy
myFileChannel = channel.fromPath( '/data/some/bigfile.txt' )
```

The above line creates a channel and binds it to a [Path](http://docs.oracle.com/javase/7/docs/api/java/nio/file/Path.html)
object for the specified file.

:::{note}
`channel.fromPath` does not check whether the file exists.
:::

Whenever the `channel.fromPath` argument contains a `*` or `?` wildcard character it is interpreted as a [glob][glob] path matcher.
For example:

```groovy
myFileChannel = channel.fromPath( '/data/big/*.txt' )
```

This example creates a channel and emits as many `Path` items as there are files with `txt` extension in the `/data/big` folder.

:::{tip}
Two asterisks, i.e. `**`, works like `*` but crosses directory boundaries. This syntax is generally used for matching complete paths. Curly brackets specify a collection of sub-patterns.
:::

For example:

```groovy
files = channel.fromPath( 'data/**.fa' )
moreFiles = channel.fromPath( 'data/**/*.fa' )
pairFiles = channel.fromPath( 'data/file_{1,2}.fq' )
```

The first line returns a channel emitting the files ending with the suffix `.fa` in the `data` folder *and* recursively in all its sub-folders. While the second one only emits the files which have the same suffix in *any* sub-folder in the `data` path. Finally the last example emits two files: `data/file_1.fq` and `data/file_2.fq`.

:::{note}
As in Linux Bash, the `*` wildcard does not catch hidden files (i.e. files whose name starts with a `.` character).
:::

Multiple paths or glob patterns can be specified using a list:

```groovy
channel.fromPath( ['/some/path/*.fq', '/other/path/*.fastq'] )
```

In order to include hidden files, you need to start your pattern with a period character or specify the `hidden: true` option. For example:

```groovy
expl1 = channel.fromPath( '/path/.*' )
expl2 = channel.fromPath( '/path/.*.fa' )
expl3 = channel.fromPath( '/path/*', hidden: true )
```

The first example returns all hidden files in the specified path. The second one returns all hidden files ending with the `.fa` suffix. Finally the last example returns all files (hidden and non-hidden) in that path.

By default a [glob][glob] pattern only looks for regular file paths that match the specified criteria, i.e. it won't return directory paths.

You can use the `type` option specifying the value `file`, `dir` or `any` in order to define what kind of paths you want. For example:

```groovy
myFileChannel = channel.fromPath( '/path/*b', type: 'dir' )
myFileChannel = channel.fromPath( '/path/a*', type: 'any' )
```

The first example will return all *directory* paths ending with the `b` suffix, while the second will return any file or directory starting with a `a` prefix.

Available options:

`checkIfExists`
: When `true` throws an exception of the specified path do not exist in the file system (default: `false`)

`followLinks`
: When `true` it follows symbolic links during directories tree traversal, otherwise they are managed as files (default: `true`)

`glob`
: When `true` interprets characters `*`, `?`, `[]` and `{}` as glob wildcards, otherwise handles them as normal characters (default: `true`)

`hidden`
: When `true` includes hidden files in the resulting paths (default: `false`)

`maxDepth`
: Maximum number of directory levels to visit (default: no limit)

`relative`
: When `true` returned paths are relative to the top-most common directory (default: `false`)

`type`
: Type of paths returned, either `file`, `dir` or `any` (default: `file`)

(channel-filepairs)=

### fromFilePairs

The `channel.fromFilePairs` method creates a channel emitting the file pairs matching a [glob][glob] pattern provided
by the user. The matching files are emitted as tuples in which the first element is the grouping key of the matching pair and the second element is the list of files (sorted in lexicographical order). For example:

```groovy
channel
    .fromFilePairs('/my/data/SRR*_{1,2}.fastq')
    .view()
```

It will produce an output similar to the following:

```
[SRR493366, [/my/data/SRR493366_1.fastq, /my/data/SRR493366_2.fastq]]
[SRR493367, [/my/data/SRR493367_1.fastq, /my/data/SRR493367_2.fastq]]
[SRR493368, [/my/data/SRR493368_1.fastq, /my/data/SRR493368_2.fastq]]
[SRR493369, [/my/data/SRR493369_1.fastq, /my/data/SRR493369_2.fastq]]
[SRR493370, [/my/data/SRR493370_1.fastq, /my/data/SRR493370_2.fastq]]
[SRR493371, [/my/data/SRR493371_1.fastq, /my/data/SRR493371_2.fastq]]
```

:::{note}
The glob pattern must contain at least one `*` wildcard character.
:::

Multiple glob patterns can be specified using a list:

```groovy
channel.fromFilePairs( ['/some/data/SRR*_{1,2}.fastq', '/other/data/QFF*_{1,2}.fastq'] )
```

Alternatively, it is possible to implement a custom file pair grouping strategy providing a closure which, given the current file as parameter, returns the grouping key. For example:

```groovy
channel
    .fromFilePairs('/some/data/*', size: -1) { file -> file.extension }
    .view { ext, files -> "Files with the extension $ext are $files" }
```

Available options:

`checkIfExists`
: When `true` throws an exception of the specified path do not exist in the file system (default: `false`)

`followLinks`
: When `true` it follows symbolic links during directories tree traversal, otherwise they are managed as files (default: `true`)

`flat`
: When `true` the matching files are produced as sole elements in the emitted tuples (default: `false`).

`hidden`
: When `true` includes hidden files in the resulting paths (default: `false`)

`maxDepth`
: Maximum number of directory levels to visit (default: no limit)

`size`
: Defines the number of files each emitted item is expected to hold (default: 2). Set to `-1` for any.

`type`
: Type of paths returned, either `file`, `dir` or `any` (default: `file`)

(channel-fromsra)=

### fromSRA

:::{versionadded} 19.04.0
:::

The `channel.fromSRA` method queries the [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) database and returns a channel emitting the FASTQ files matching the specified criteria i.e project or accession number(s). For example:

```groovy
channel
    .fromSRA('SRP043510')
    .view()
```

It returns:

```
[SRR1448794, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/004/SRR1448794/SRR1448794.fastq.gz]
[SRR1448795, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/005/SRR1448795/SRR1448795.fastq.gz]
[SRR1448792, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/002/SRR1448792/SRR1448792.fastq.gz]
[SRR1448793, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/003/SRR1448793/SRR1448793.fastq.gz]
[SRR1910483, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/003/SRR1910483/SRR1910483.fastq.gz]
[SRR1910482, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/002/SRR1910482/SRR1910482.fastq.gz]
(remaining omitted)
```

Multiple accession IDs can be specified using a list object:

```groovy
ids = ['ERR908507', 'ERR908506', 'ERR908505']
channel
    .fromSRA(ids)
    .view()
```

```
[ERR908507, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908507/ERR908507_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908507/ERR908507_2.fastq.gz]]
[ERR908506, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908506/ERR908506_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908506/ERR908506_2.fastq.gz]]
[ERR908505, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908505/ERR908505_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908505/ERR908505_2.fastq.gz]]
```

:::{note}
Each read pair is implicitly managed and returned as a list of files.
:::

This method uses the NCBI [ESearch](https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch) API behind the scenes, therefore it allows the use of any query term supported by this API.

To access the ESearch API, you must provide your [NCBI API keys](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities) through one of the following ways:

- The `apiKey` option:
  ```groovy
  channel.fromSRA(ids, apiKey:'0123456789abcdef')
  ```

- The `NCBI_API_KEY` variable in your environment:
  ```bash
  export NCBI_API_KEY=0123456789abcdef
  ```

Available options:

`apiKey`
: NCBI user API key.

`cache`
: Enable/disable the caching API requests (default: `true`).

`max`
: Maximum number of entries that can be retried (default: unlimited) .

`protocol`
: Allow choosing the protocol for the resulting remote URLs. Available choices: `ftp`, `http`, `https` (default: `ftp`).

(channel-of)=

### of

:::{versionadded} 19.10.0
:::

The `channel.of` method allows you to create a channel that emits the arguments provided to it, for example:

```groovy
ch = channel.of( 1, 3, 5, 7 )
ch.view { "value: $it" }
```

The first line in this example creates a variable `ch` which holds a channel object. This channel emits the arguments
supplied to the `of` method. Thus the second line prints the following:

```
value: 1
value: 3
value: 5
value: 7
```

Ranges of values are expanded accordingly:

```groovy
channel
    .of(1..23, 'X', 'Y')
    .view()
```

Prints:

```
1
2
3
4
:
23
X
Y
```

See also: [channel.fromList](#fromlist) factory method.

(channel-topic)=

### topic

:::{versionadded} 23.11.0-edge
:::

:::{note}
This feature requires the `nextflow.preview.topic` feature flag to be enabled.
:::

A *topic* is a channel type introduced as of Nextflow 23.11.0-edge along with {ref}`channel-type-value` and
{ref}`channel-type-queue`.

A *topic channel*, similarly to a *queue channel*, is non-blocking unidirectional FIFO queue, however it connects
multiple *producer* processes with multiple *consumer* processes or operators.

:::{tip}
You can think about it as a channel that is shared across many different process using the same *topic name*.
:::

A process output can be assigned to a topic using the `topic` option on an output, for example:

```groovy
process foo {
  output:
  val('foo'), topic: my_topic
}

process bar {
  output:
  val('bar'), topic: my_topic
}
```

The `channel.topic` method allows referencing the topic channel with the specified name, which can be used as a process
input or operator composition as any other Nextflow channel:

```groovy
channel.topic('my-topic').view()
```

This approach is a convenient way to collect related items from many different sources without explicitly defining
the logic connecting many different queue channels altogether, commonly using the `mix` operator.

:::{warning}
Any process that consumes a channel topic should not send any outputs to that topic, or else the pipeline will hang forever.
:::

See also: {ref}`process-additional-options` for process outputs.

(channel-topic)=

### topic

:::{versionadded} 23.11.0-edge
:::

:::{note}
This feature requires the `nextflow.preview.topic` feature flag to be enabled.
:::

A *topic* is a channel type introduced as of Nextflow 23.11.0-edge along with {ref}`channel-type-value` and
{ref}`channel-type-queue`.

A *topic channel*, similarly to a *queue channel*, is non-blocking unidirectional FIFO queue, however it connects
multiple *producer* processes with multiple *consumer* processes or operators.

:::{tip}
You can think about it as a channel that is shared across many different process using the same *topic name*.
:::

A process output can be assigned to a topic using the `topic` option on an output, for example:

```groovy
process foo {
  output:
  val('foo'), topic: my_topic
}

process bar {
  output:
  val('bar'), topic: my_topic
}
```

The `channel.topic` method allows referencing the topic channel with the specified name, which can be used as a process
input or operator composition as any other Nextflow channel:

```groovy
Channel.topic('my-topic').view()
```

This approach is a convenient way to collect related items from many different sources without explicitly defining
the logic connecting many different queue channels altogether, commonly using the `mix` operator.

:::{warning}
Any process that consumes a channel topic should not send any outputs to that topic, or else the pipeline will hang forever.
:::

See also: {ref}`process-additional-options` for process outputs.

(channel-value)=

### value

The `channel.value` method is used to create a value channel. An optional (not `null`) argument can be specified to bind
the channel to a specific value. For example:

```groovy
expl1 = channel.value()
expl2 = channel.value( 'Hello there' )
expl3 = channel.value( [1,2,3,4,5] )
```

The first line in the example creates an 'empty' variable. The second line creates a channel and binds a string to it.
The third line creates a channel and binds a list object to it that will be emitted as a single value.

(channel-watchpath)=

### watchPath

The `channel.watchPath` method watches a folder for one or more files matching a specified pattern. As soon as there
is a file that meets the specified condition, it is emitted over the channel that is returned by the `watchPath` method.
The condition on files to watch can be specified by using `*` or `?` wildcard characters i.e. by specifying a [glob][glob] path matching criteria.

For example:

```groovy
channel
    .watchPath( '/path/*.fa' )
    .subscribe { println "Fasta file: $it" }
```

By default it watches only for new files created in the specified folder. Optionally, it is possible to provide a second
argument that specifies what event(s) to watch. The supported events are:

- `create`: A new file is created (default)
- `modify`: A file is modified
- `delete`: A file is deleted

You can specify more than one of these events by using a comma separated string as shown below:

```groovy
channel
    .watchPath( '/path/*.fa', 'create,modify' )
    .subscribe { println "File created or modified: $it" }
```

:::{warning}
The `channel.watchPath` factory waits endlessly for files that match the specified pattern and event(s), which means
that it will cause your pipeline to run forever. Consider using the `take` or `until` operator to close the channel when
a certain condition is met (e.g. after receiving 10 files, receiving a file named `DONE`).
:::

See also: [channel.fromPath](#frompath) factory method.

[glob]: http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob

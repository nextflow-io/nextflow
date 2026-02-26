(working-with-files)=

# Working with files

## Retrieving files

Use the `file()` function to obtain a reference to a file by name:

```nextflow
myFile = file('some/path/to/my_file.file')
```

The `file()` function can reference both files and directories.

Use the `files()` function to obtain a list of files. When using the wildcard characters `*`, `?`, `[]` and `{}`, the file name is treated as a [glob](http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob) pattern, returning all files that match the given pattern, or an empty list if no matching files are found:

```nextflow
listOfFiles = files('some/path/*.fa')
```

:::{note}
The `file()` function can also be called with a glob pattern, as long as the pattern is intended to match exactly one file.
:::

A double asterisk (`**`) in a glob pattern works like `*` but also searches through subdirectories:

```nextflow
deeplyNestedFiles = files('some/path/**/*.fa')
```

By default, wildcard characters do not match directories or hidden files. Use the `hidden` option to include hidden files:

```nextflow
listWithHidden = file('some/path/*.fa', hidden: true)
```

Given a file reference, you can use the `resolve()` method or the `/` operator to obtain files relative to that path:

```nextflow
def dir = file('s3://bucket/some/data/path')

dir.resolve('sample.bam')         // correct
dir / 'sample.bam'
file("$dir/sample.bam")           // correct (but verbose)
"$dir/sample.bam"                 // incorrect
```

## Getting file attributes

The `file()` function returns a {ref}`Path <stdlib-types-path>`, which has several methods for retrieving metadata about the file:

```nextflow
def path = file('/some/path/file.txt')

assert path.baseName == 'file'
assert path.extension == 'txt'
assert path.name == 'file.txt'
assert path.parent == '/some/path'
```

See the {ref}`stdlib-types-path` reference for the list of available methods.

## Reading and writing

### Reading and writing an entire file

Reading a file is as easy as using the file's `text` property, which returns the file contents as a string:

```nextflow
print myFile.text
```

Similarly, you can write text to a file by assigning it to the file's `text` property:

```nextflow
myFile.text = 'Hello world!'
```

This approach overwrites any existing file contents, and implicitly creates the file if it doesn't exist.

:::{tip}
The `text` property is shorthand for the `getText()` and `setText()` methods:

```nextflow
println myFile.getText()
myFile.setText('Hello world!')
```
:::

:::{warning}
The above methods read and write the *entire* file contents at once, requiring the entire file to be loaded into memory. Consider using a more memory-efficient approach for large files, such as reading/writing the file line by line.
:::

### Reading a file line by line

You can use the `readLines()` method to read a text file line by line:

```nextflow
file('some/my_file.txt')
    .readLines()
    .each { line ->
        println line
    }
```

The `readLines()` method loads the *entire* file into memory, so it is not ideal for large files.

You can use the `eachLine()` method to read line by line while only loading one line at a time into memory:

```nextflow
count = 0
myFile.eachLine { line ->
    println "line ${count++}: $line"
}
```

The `withReader()` method creates a [Reader](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/io/Reader.html) that you can use to read the file line by line, or even character by character. It is useful when you don't need to read the entire file.

For example, to read only the first line of a file:

```nextflow
myFile.withReader { r ->
    def firstLine = r.readLine()
    println firstLine
}
```

### Writing a file line by line

You can use the `append()` method or left shirt (`<<`) operator to append text to a file without erasing the existing contents:

```nextflow
myFile.append('Add this line\n')
myFile << 'Add a line more\n'
```

For example, the following snippet copies the contents of a source file into a target file, replacing all `U` characters with `X`:

```nextflow
sourceFile.eachLine { line ->
    targetFile << line.replaceAll('U', 'X')
}
```

## Filesystem operations

See the {ref}`stdlib-types-path` reference for the complete list of methods for performing filesystem operations.

### Listing directories

You can use the `listFiles()` method to list the contents of a directory:

```nextflow
children = file('any/path').list()
children.each { file ->
    println file
}
```

:::{versionchanged} 26.04.0
The `listFiles()` method is deprecated -- use `listDirectory()` instead.
:::

You can use the `eachFile()` method to iterate through the contents of a directory:

```nextflow
myDir.eachFile { item ->
    if( item.isFile() ) {
        println "${item.getName()} - size: ${item.size()}"
    }
    else if( item.isDirectory() ) {
        println "${item.getName()} - DIR"
    }
}
```

### Copying files

In general, you should not need to manually copy files, because Nextflow will automatically stage files in and out of the task environment based on the definition of process inputs and outputs. Ideally, any operation which transforms files should be encapsulated in a process, in order to leverage Nextflow's staging capabilities as much as possible.

(remote-files)=

## Remote files

Nextflow works with many types of remote files and objects using the same interface as for local files. The following protocols are supported:

- HTTP(S)/FTP (`http://`, `https://`, `ftp://`)
- Amazon S3 (`s3://`)
- Azure Blob Storage (`az://`)
- Google Cloud Storage (`gs://`)

To reference a remote file, simply specify the URL when opening the file:

```nextflow
pdb = file('http://files.rcsb.org/header/5FID.pdb')
```

It can then be used in the same way as a local file:

```nextflow
println pdb.text
```

:::{note}
Not all operations are supported for all protocols. For example, writing and directory listing is not supported for HTTP(S) and FTP paths.
:::

:::{note}
Additional configuration may be necessary for cloud object storage, such as authenticating with a private bucket. See the documentation for each cloud storage provider for further details.
:::

### Remote file staging

When a process input file resides on a different file system than the work directory, Nextflow copies the file into the work directory using an appropriate Java SDK.

Remote files are staged in a subdirectory of the work directory with the form `stage-<session-id>/<hash>/<filename>`, where `<hash>` is determined by the remote file path. If multiple tasks request the same remote file, the file will be downloaded once and reused by each task. These files can be reused by resumed runs with the same session ID.

:::{note}
Remote file staging can be a bottleneck during large-scale runs, particularly when input files are stored in object storage but need to be staged in a shared filesystem work directory. This bottleneck occurs because Nextflow handles all of these file transfers.

To mitigate this, you can implement a custom process to download the required files, allowing you to stage multiple files efficiently through parallel jobs. Files should be given as a `val` input instead of a `path` input to bypass Nextflow's built-in remote file staging.

Alternatively, use {ref}`fusion-page` with the work directory set to object storage. In this case, tasks can access remote files directly without any prior staging, eliminating the bottleneck.
:::

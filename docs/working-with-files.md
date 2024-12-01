(working-with-files)=

# Working with files

## Opening files

To access and work with files, use the `file()` method, which returns a file system object given a file path string:

```nextflow
myFile = file('some/path/to/my_file.file')
```

The `file()` method can reference both files and directories, depending on what the string path refers to in the file system.

When using the wildcard characters `*`, `?`, `[]` and `{}`, the argument is interpreted as a [glob](http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob) path matcher and the `file()` method returns a list object holding the paths of files whose names match the specified pattern, or an empty list if no match is found:

```nextflow
listOfFiles = file('some/path/*.fa')
```

:::{note}
The `file()` method does not return a list if only one file is matched. Use the `files()` method to always return a list.
:::

:::{note}
A double asterisk (`**`) in a glob pattern works like `*` but also searches through subdirectories.
:::

By default, wildcard characters do not match directories or hidden files. For example, if you want to include hidden files in the result list, enable the `hidden` option:

```nextflow
listWithHidden = file('some/path/*.fa', hidden: true)
```

:::{note}
To compose paths, instead of string interpolation, use the `resolve()` method or the `/` operator:

```nextflow
def dir = file('s3://bucket/some/data/path')
def sample1 = dir.resolve('sample.bam')         // correct
def sample2 = dir / 'sample.bam'
def sample3 = file("$dir/sample.bam")           // correct (but verbose)
def sample4 = "$dir/sample.bam"                 // incorrect
```
:::

## Getting file attributes

The `file()` method returns a [Path](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/nio/file/Path.html), which has several methods for retrieving metadata about the file:

```nextflow
def path = file('/some/path/file.txt')

assert path.baseName == 'file'
assert path.extension == 'txt'
assert path.name == 'file.txt'
assert path.parent == '/some/path'
```

:::{tip}
When calling an object method, any method that looks like `get*()` can also be accessed as a field. For example, `path.getName()` is equivalent to `path.name`, `path.getBaseName()` is equivalent to `path.baseName`, and so on.
:::

See the {ref}`stdlib-types-path` reference for the list of available methods.

## Reading and writing

### Reading and writing an entire file

Given a file variable, created with the `file()` method as shown previously, reading a file is as easy as getting the file's `text` property, which returns the file content as a string:

```nextflow
print myFile.text
```

Similarly, you can save a string to a file by assigning it to the file's `text` property:

```nextflow
myFile.text = 'Hello world!'
```

Binary data can managed in the same way, just using the file property `bytes` instead of `text`. Thus, the following example reads the file and returns its content as a byte array:

```nextflow
binaryContent = myFile.bytes
```

Or you can save a byte array to a file:

```nextflow
myFile.bytes = binaryContent
```

:::{note}
The above assignment overwrites any existing file contents, and implicitly creates the file if it doesn't exist.
:::

:::{warning}
The above methods read and write the **entire** file contents at once, in a single variable or buffer. For this reason, when dealing with large files it is recommended that you use a more memory efficient approach, such as reading/writing a file line by line or using a fixed size buffer.
:::

### Appending to a file

In order to append a string value to a file without erasing existing content, you can use the `append()` method:

```nextflow
myFile.append('Add this line\n')
```

Or use the left shift operator, a more idiomatic way to append text content to a file:

```nextflow
myFile << 'Add a line more\n'
```

### Reading a file line by line

In order to read a text file line by line you can use the method `readLines()` provided by the file object, which returns the file content as a list of strings:

```nextflow
myFile = file('some/my_file.txt')
allLines = myFile.readLines()
for( line : allLines ) {
    println line
}
```

This can also be written in a more idiomatic syntax:

```nextflow
file('some/my_file.txt')
    .readLines()
    .each { println it }
```

:::{warning}
The method `readLines()` reads the **entire** file at once and returns a list containing all the lines. For this reason, do not use it to read big files.
:::

To process a big file, use the method `eachLine()`, which reads only a single line at a time into memory:

```nextflow
count = 0
myFile.eachLine { str ->
    println "line ${count++}: $str"
}
```

### Advanced file reading

The classes `Reader` and `InputStream` provide fine-grained control for reading text and binary files, respectively.

The method `newReader()` creates a [Reader](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/io/Reader.html) object for the given file that allows you to read the content as single characters, lines or arrays of characters:

```nextflow
myReader = myFile.newReader()
String line
while( line = myReader.readLine() ) {
    println line
}
myReader.close()
```

The method `withReader()` works similarly, but automatically calls the `close()` method for you when you have finished processing the file. So, the previous example can be written more simply as:

```nextflow
myFile.withReader {
    String line
    while( line = it.readLine() ) {
        println line
    }
}
```

The methods `newInputStream()` and `withInputStream()` work similarly. The main difference is that they create an [InputStream](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/io/InputStream.html) object useful for writing binary data.

See the {ref}`stdlib-types-path` reference for the list of available methods.

### Advanced file writing

The `Writer` and `OutputStream` classes provide fine-grained control for writing text and binary files, respectively, including low-level operations for single characters or bytes, and support for big files.

For example, given two file objects `sourceFile` and `targetFile`, the following code copies the first file's content into the second file, replacing all `U` characters with `X`:

```nextflow
sourceFile.withReader { source ->
    targetFile.withWriter { target ->
        String line
        while( line=source.readLine() ) {
            target << line.replaceAll('U','X')
        }
    }
}
```

See the {ref}`stdlib-types-path` reference for the list of available methods.

## Filesystem operations

Methods for performing filesystem operations such as copying, deleting, and directory listing are documented in the {ref}`stdlib-types-path` reference.

### Listing directories

The simplest way to list a directory is to use `list()` or `listFiles()`, which return a collection of first-level elements (files and directories) of a directory:

```nextflow
for( def file : file('any/path').list() ) {
    println file
}
```

Additionally, the `eachFile()` method allows you to iterate through the first-level elements only (just like `listFiles()`). As with other `each*()` methods, `eachFile()` takes a closure as a parameter:

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

## Remote files

Nextflow can work with many kinds of remote files and objects using the same interface as for local files. The following protocols are supported:

- HTTP(S) / FTP (`http://`, `https://`, `ftp://`)
- Amazon S3 (`s3://`)
- Azure Blob Storage (`az://`)
- Google Cloud Storage (`gs://`)

To reference a remote file, simple specify the URL when opening the file:

```nextflow
pdb = file('http://files.rcsb.org/header/5FID.pdb')
```

You can then access it as a local file as described previously:

```nextflow
println pdb.text
```

:::{note}
Not all operations are supported for all protocols. In particular, writing and directory listing are not supported for HTTP(S) and FTP paths.
:::

:::{note}
Additional configuration may be required to work with cloud object storage (e.g. to authenticate with a private bucket). Refer to the respective page for each cloud storage provider for more information.
:::

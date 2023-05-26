(script-page)=

# Nextflow scripting

The Nextflow scripting language is an extension of the Groovy programming language. Groovy is a powerful programming language for the Java virtual machine. The Nextflow syntax has been specialized to ease the writing of computational pipelines in a declarative manner.

Nextflow can execute any piece of Groovy code or use any library for the JVM platform.

For a detailed description of the Groovy programming language, reference these links:

- [Groovy User Guide](http://groovy-lang.org/documentation.html)
- [Groovy Cheat sheet](http://www.cheat-sheets.org/saved-copy/rc015-groovy_online.pdf)
- [Groovy in Action](http://www.manning.com/koenig2/)

Below you can find a crash course in the most important language constructs used in the Nextflow scripting language.

:::{warning}
Nextflow uses UTF-8 as the default character encoding for source files. Make sure to use UTF-8 encoding when editing Nextflow scripts with your preferred text editor.
:::

## Language basics

### Hello world

To print something is as easy as using one of the `print` or `println` methods.

```groovy
println "Hello, World!"
```

The only difference between the two is that the `println` method implicitly appends a newline character to the printed string.

### Variables

To define a variable, simply assign a value to it:

```groovy
x = 1
println x

x = new java.util.Date()
println x

x = -3.1499392
println x

x = false
println x

x = "Hi"
println x
```

### Lists

A List object can be defined by placing the list items in square brackets:

```groovy
myList = [1776, -1, 33, 99, 0, 928734928763]
```

You can access a given item in the list with square-bracket notation (indexes start at 0):

```groovy
println myList[0]
```

In order to get the length of the list use the `size` method:

```groovy
println myList.size()
```

Learn more about lists:

- [Groovy Lists tutorial](http://groovy-lang.org/groovy-dev-kit.html#Collections-Lists)
- [Groovy List SDK](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/List.html)
- [Java List SDK](http://docs.oracle.com/javase/7/docs/api/java/util/List.html)

### Maps

Maps are used to store *associative arrays* (also known as *dictionaries*). They are unordered collections of heterogeneous, named data:

```groovy
scores = [ "Brett":100, "Pete":"Did not finish", "Andrew":86.87934 ]
```

Note that each of the values stored in the map can be of a different type. `Brett` is an integer, `Pete` is a string, and `Andrew` is a floating-point number.

We can access the values in a map in two main ways:

```groovy
println scores["Pete"]
println scores.Pete
```

To add data to or modify a map, the syntax is similar to adding values to list:

```groovy
scores["Pete"] = 3
scores["Cedric"] = 120
```

You can also use the `+` operator to add two maps together:

```groovy
new_scores = scores + ["Pete": 3, "Cedric": 120]
```

When adding two maps, the first map is copied and then appended with the keys from the second map. Any conflicting keys are overwritten by the second map.

:::{tip}
Appending an "update" map is a safer way to modify maps in Nextflow, specifically when passing maps through channels. This way, any references to the original map elsewhere in the pipeline won't be modified.
:::

Learn more about maps:

- [Groovy Maps tutorial](http://groovy-lang.org/groovy-dev-kit.html#Collections-Maps)
- [Groovy Map SDK](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html)
- [Java Map SDK](http://docs.oracle.com/javase/7/docs/api/java/util/Map.html)

(script-multiple-assignment)=

### Multiple assignment

An array or a list object can used to assign to multiple variables at once:

```groovy
(a, b, c) = [10, 20, 'foo']
assert a == 10 && b == 20 && c == 'foo'
```

The three variables on the left of the assignment operator are initialized by the corresponding item in the list.

Read more about [Multiple assignment](http://www.groovy-lang.org/semantics.html#_multiple_assignment) in the Groovy documentation.

### Conditional Execution

One of the most important features of any programming language is the ability to execute different code under different conditions. The simplest way to do this is to use the `if` construct:

```groovy
x = Math.random()
if( x < 0.5 ) {
    println "You lost."
}
else {
    println "You won!"
}
```

### Strings

Strings can be defined by enclosing text in single or double quotes (`'` or `"` characters):

```groovy
println "he said 'cheese' once"
println 'he said "cheese!" again'
```

Strings can be concatenated with `+`:

```groovy
a = "world"
print "hello " + a + "\n"
```

(string-interpolation)=

### String interpolation

There is an important difference between single-quoted and double-quoted strings: Double-quoted strings support variable interpolations, while single-quoted strings do not.

In practice, double-quoted strings can contain the value of an arbitrary variable by prefixing its name with the `$` character, or the value of any expression by using the `${expression}` syntax, similar to Bash/shell scripts:

```groovy
foxtype = 'quick'
foxcolor = ['b', 'r', 'o', 'w', 'n']
println "The $foxtype ${foxcolor.join()} fox"

x = 'Hello'
println '$x + $y'
```

This code prints:

```
The quick brown fox
$x + $y
```

### Multi-line strings

A block of text that span multiple lines can be defined by delimiting it with triple single or double quotes:

```groovy
text = """
    hello there James
    how are you today?
    """
```

:::{note}
Like before, multi-line strings inside double quotes support variable interpolation, while single-quoted multi-line strings do not.
:::

As in Bash/shell scripts, terminating a line in a multi-line string with a `\` character prevents a newline character from separating that line from the one that follows:

```groovy
myLongCmdline = """
    blastp \
    -in $input_query \
    -out $output_file \
    -db $blast_database \
    -html
    """

result = myLongCmdline.execute().text
```

In the preceding example, `blastp` and its `-in`, `-out`, `-db` and `-html` switches and their arguments are effectively a single line.

(implicit-variables)=

## Implicit variables

### Script implicit variables

The following variables are implicitly defined in the script global execution scope:

`baseDir`
: :::{deprecated} 20.04.0
  Use `projectDir` instead
  :::
: The directory where the main workflow script is located.

`launchDir`
: :::{versionadded} 20.04.0
  :::
: The directory where the workflow is run.

`moduleDir`
: :::{versionadded} 20.04.0
  :::
: The directory where a module script is located for DSL2 modules or the same as `projectDir` for a non-module script.

`nextflow`
: Dictionary like object representing nextflow runtime information (see {ref}`metadata-nextflow`).

`params`
: Dictionary like object holding workflow parameters specifying in the config file or as command line options.

`projectDir`
: :::{versionadded} 20.04.0
  :::
: The directory where the main script is located.

`workDir`
: The directory where tasks temporary files are created.

`workflow`
: Dictionary like object representing workflow runtime information (see {ref}`metadata-workflow`).

### Configuration implicit variables

The following variables are implicitly defined in the Nextflow configuration file:

`baseDir`
: :::{deprecated} 20.04.0
  Use `projectDir` instead
  :::
: The directory where the main workflow script is located.

`launchDir`
: :::{versionadded} 20.04.0
  :::
: The directory where the workflow is run.

`projectDir`
: :::{versionadded} 20.04.0
  :::
: The directory where the main script is located.

### Process implicit variables

The following variables are implicitly defined in the `task` object of each process:

`attempt`
: The current task attempt

`hash`
: *Available only in `exec:` blocks*
: The task unique hash ID

`index`
: The task index (corresponds to `task_id` in the execution trace)

`name`
: *Available only in `exec:` blocks*
: The current task name

`process`
: The current process name

`workDir`
: *Available only in `exec:` blocks*
: The task unique directory

The `task` object also contains the values of all process directives for the given task, which allows you to access these settings at runtime. For examples:

```groovy
process foo {
  script:
  """
  some_tool --cpus $task.cpus --mem $task.memory
  """
}
```

In the above snippet the `task.cpus` holds the value for the {ref}`cpus directive<process-cpus>` and the `task.memory` the current value for {ref}`memory directive<process-memory>` depending on the actual setting given in the workflow configuration file.

See {ref}`Process directives <process-directives>` for details.

(script-closure)=

## Closures

Briefly, a closure is a block of code that can be passed as an argument to a function. Thus, you can define a chunk of code and then pass it around as if it were a string or an integer.

More formally, you can create functions that are defined as *first-class objects*.

```groovy
square = { it * it }
```

The curly brackets around the expression `it * it` tells the script interpreter to treat this expression as code. The `it` identifier is an implicit variable that represents the value that is passed to the function when it is invoked.

Once compiled the function object is assigned to the variable `square` as any other variable assignments shown previously. Now we can do something like this:

```groovy
println square(9)
```

and get the value 81.

This is not very interesting until we find that we can pass the function `square` as an argument to other functions or methods. Some built-in functions take a function like this as an argument. One example is the `collect` method on lists:

```groovy
[ 1, 2, 3, 4 ].collect(square)
```

This expression says: Create an array with the values 1, 2, 3 and 4, then call its `collect` method, passing in the closure we defined above. The `collect` method runs through each item in the array, calls the closure on the item, then puts the result in a new array, resulting in:

```groovy
[ 1, 4, 9, 16 ]
```

For more methods that you can call with closures as arguments, see the [Groovy GDK documentation](http://docs.groovy-lang.org/latest/html/groovy-jdk/).

By default, closures take a single parameter called `it`, but you can also create closures with multiple, custom-named parameters. For example, the method `Map.each()` can take a closure with two arguments, to which it binds the `key` and the associated `value` for each key-value pair in the `Map`. Here, we use the obvious variable names `key` and `value` in our closure:

```groovy
printMapClosure = { key, value ->
    println "$key = $value"
}

[ "Yue" : "Wu", "Mark" : "Williams", "Sudha" : "Kumari" ].each(printMapClosure)
```

Prints:

```
Yue = Wu
Mark = Williams
Sudha = Kumari
```

A closure has two other important features. First, it can access variables in the scope where it is defined, so that it can interact with them.

Second, a closure can be defined in an anonymous manner, meaning that it is not given a name, and is defined in the place where it needs to be used.

As an example showing both these features, see the following code fragment:

```groovy
myMap = ["China": 1 , "India" : 2, "USA" : 3]

result = 0
myMap.keySet().each( { result+= myMap[it] } )

println result
```

Learn more about closures in the [Groovy documentation](http://groovy-lang.org/closures.html)

(script-regexp)=

## Regular expressions

Regular expressions are the Swiss Army knife of text processing. They provide the programmer with the ability to match and extract patterns from strings.

Regular expressions are available via the `~/pattern/` syntax and the `=~` and `==~` operators.

Use `=~` to check whether a given pattern occurs anywhere in a string:

```groovy
assert 'foo' =~ /foo/       // return TRUE
assert 'foobar' =~ /foo/    // return TRUE
```

Use `==~` to check whether a string matches a given regular expression pattern exactly.

```groovy
assert 'foo' ==~ /foo/       // return TRUE
assert 'foobar' ==~ /foo/    // return FALSE
```

It is worth noting that the `~` operator creates a Java `Pattern` object from the given string, while the `=~` operator creates a Java `Matcher` object.

```groovy
x = ~/abc/
println x.class
// prints java.util.regex.Pattern

y = 'some string' =~ /abc/
println y.class
// prints java.util.regex.Matcher
```

Regular expression support is imported from Java. Java's regular expression language and API is documented in the [Pattern Java documentation](http://download.oracle.com/javase/7/docs/api/java/util/regex/Pattern.html).

You may also be interested in this post: [Groovy: Don't Fear the RegExp](https://web.archive.org/web/20170621185113/http://www.naleid.com/blog/2008/05/19/dont-fear-the-regexp).

### String replacement

To replace pattern occurrences in a given string, use the `replaceFirst` and `replaceAll` methods:

```groovy
x = "colour".replaceFirst(/ou/, "o")
println x
// prints: color

y = "cheesecheese".replaceAll(/cheese/, "nice")
println y
// prints: nicenice
```

### Capturing groups

You can match a pattern that includes groups. First create a matcher object with the `=~` operator. Then, you can index the matcher object to find the matches: `matcher[0]` returns a list representing the first match of the regular expression in the string. The first list element is the string that matches the entire regular expression, and the remaining elements are the strings that match each group.

Here's how it works:

```groovy
programVersion = '2.7.3-beta'
m = programVersion =~ /(\d+)\.(\d+)\.(\d+)-?(.+)/

assert m[0] == ['2.7.3-beta', '2', '7', '3', 'beta']
assert m[0][1] == '2'
assert m[0][2] == '7'
assert m[0][3] == '3'
assert m[0][4] == 'beta'
```

Applying some syntactic sugar, you can do the same in just one line of code:

```groovy
programVersion = '2.7.3-beta'
(full, major, minor, patch, flavor) = (programVersion =~ /(\d+)\.(\d+)\.(\d+)-?(.+)/)[0]

println full    // 2.7.3-beta
println major   // 2
println minor   // 7
println patch   // 3
println flavor  // beta
```

### Removing part of a string

You can remove part of a `String` value using a regular expression pattern. The first match found is replaced with an empty String:

```groovy
// define the regexp pattern
wordStartsWithGr = ~/(?i)\s+Gr\w+/

// apply and verify the result
('Hello Groovy world!' - wordStartsWithGr) == 'Hello world!'
('Hi Grails users' - wordStartsWithGr) == 'Hi users'
```

Remove the first 5-character word from a string:

```groovy
assert ('Remove first match of 5 letter word' - ~/\b\w{5}\b/) == 'Remove match of 5 letter word'
```

Remove the first number with its trailing whitespace from a string:

```groovy
assert ('Line contains 20 characters' - ~/\d+\s+/) == 'Line contains characters'
```

(script-file-io)=

## Files and I/O

### Opening files

To access and work with files, use the `file` method, which returns a file system object given a file path string:

```groovy
myFile = file('some/path/to/my_file.file')
```

The `file` method can reference both files and directories, depending on what the string path refers to in the file system.

When using the wildcard characters `*`, `?`, `[]` and `{}`, the argument is interpreted as a [glob][glob] path matcher and the `file` method returns a list object holding the paths of files whose names match the specified pattern, or an empty list if no match is found:

```groovy
listOfFiles = file('some/path/*.fa')
```

:::{note}
Two asterisks (`**`) in a glob pattern works like `*` but also searches through subdirectories.
:::

By default, wildcard characters do not match directories or hidden files. For example, if you want to include hidden files in the result list, add the optional parameter `hidden`:

```groovy
listWithHidden = file('some/path/*.fa', hidden: true)
```

Here are `file`'s available options:

| Name          | Description                                                                                                                                |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| glob          | When `true` interprets characters `*`, `?`, `[]` and `{}` as glob wildcards, otherwise handles them as normal characters (default: `true`) |
| type          | Type of paths returned, either `file`, `dir` or `any` (default: `file`)                                                                    |
| hidden        | When `true` includes hidden files in the resulting paths (default: `false`)                                                                |
| maxDepth      | Maximum number of directory levels to visit (default: `no limit`)                                                                          |
| followLinks   | When `true` follows symbolic links during directory tree traversal, otherwise treats them as files (default: `true`)                       |
| checkIfExists | When `true` throws an exception of the specified path do not exist in the file system (default: `false`)                                   |

:::{note}
Nextflow also provides a `files()` method, which is identical to `file()` except that it always returns a list, whereas `file()` only returns a list if it matches multiple files.
:::

:::{tip}
If you are a Java geek, you might be interested to know that the `file` method returns a [Path](http://docs.oracle.com/javase/8/docs/api/java/nio/file/Path.html) object, which allows you to use the same methods you would use in a Java program.
:::

:::{warning}
When a file reference is created using a URL path and it's converted to a string, the protocol schema is not included in the resulting string.

You can check this behavior with the code below:

```groovy
def ref = file('s3://some-bucket/foo.txt')
assert ref.toString() == '/some-bucket/foo.txt'
assert "$ref" == '/some-bucket/foo.txt'
```

To include the schema, the `toUriString()` method should be used instead:

```groovy
assert ref.toUriString() == 's3://some-bucket/foo.txt'
```

Also, instead of composing paths through string interpolation, the `resolve()` method or the `/` operator should be used instead::

```groovy
def dir = file('s3://bucket/some/data/path')
def sample = "$dir/sample.bam"                // don't do this
def sample1 = dir.resolve('sample.bam')
def sample2 = dir / 'sample.bam'              // the `/` operator can be used to compose paths
```
:::

See also: {ref}`Channel.fromPath <channel-path>`.

### Basic read/write

Given a file variable, declared using the `file` method as shown in the previous example, reading a file is as easy as getting the value of the file's `text` property, which returns the file content as a string value:

```groovy
print myFile.text
```

Similarly, you can save a string value to a file by simply assigning it to the file's `text` property:

```groovy
myFile.text = 'Hello world!'
```

:::{note}
The above assignment overwrites any existing file contents, and implicitly creates the file if it doesn't exist.
:::

In order to append a string value to a file without erasing existing content, you can use the `append` method:

```groovy
myFile.append('Add this line\n')
```

Or use the left shift operator, a more idiomatic way to append text content to a file:

```groovy
myFile << 'Add a line more\n'
```

Binary data can managed in the same way, just using the file property `bytes` instead of `text`. Thus, the following example reads the file and returns its content as a byte array:

```groovy
binaryContent = myFile.bytes
```

Or you can save a byte array data buffer to a file, by simply writing:

```groovy
myFile.bytes = binaryBuffer
```

:::{warning}
The above methods read and write the **entire** file contents at once, in a single variable or buffer. For this reason, when dealing with large files it is recommended that you use a more memory efficient approach, such as reading/writing a file line by line or using a fixed size buffer.
:::

### Read a file line by line

In order to read a text file line by line you can use the method `readLines()` provided by the file object, which returns the file content as a list of strings:

```groovy
myFile = file('some/my_file.txt')
allLines = myFile.readLines()
for( line : allLines ) {
    println line
}
```

This can also be written in a more idiomatic syntax:

```groovy
file('some/my_file.txt')
    .readLines()
    .each { println it }
```

:::{warning}
The method `readLines()` reads the **entire** file at once and returns a list containing all the lines. For this reason, do not use it to read big files.
:::

To process a big file, use the method `eachLine`, which reads only a single line at a time into memory:

```groovy
count = 0
myFile.eachLine { str ->
    println "line ${count++}: $str"
}
```

### Advanced file reading operations

The classes `Reader` and `InputStream` provide fine control for reading text and binary files, respectively.\_

The method `newReader` creates a [Reader](http://docs.oracle.com/javase/7/docs/api/java/io/Reader.html) object for the given file that allows you to read the content as single characters, lines or arrays of characters:

```groovy
myReader = myFile.newReader()
String line
while( line = myReader.readLine() ) {
    println line
}
myReader.close()
```

The method `withReader` works similarly, but automatically calls the `close` method for you when you have finished processing the file. So, the previous example can be written more simply as:

```groovy
myFile.withReader {
    String line
    while( line = it.readLine() ) {
        println line
    }
}
```

The methods `newInputStream` and `withInputStream` work similarly. The main difference is that they create an [InputStream](http://docs.oracle.com/javase/7/docs/api/java/io/InputStream.html) object useful for writing binary data.

Here are the most important methods for reading from files:

| Name            | Description                                                                                                                                     |
| --------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| getText         | Returns the file content as a string value                                                                                                      |
| getBytes        | Returns the file content as byte array                                                                                                          |
| readLines       | Reads the file line by line and returns the content as a list of strings                                                                        |
| eachLine        | Iterates over the file line by line, applying the specified {ref}`closure <script-closure>`                                                     |
| eachByte        | Iterates over the file byte by byte, applying the specified {ref}`closure <script-closure>`                                                     |
| withReader      | Opens a file for reading and lets you access it with a [Reader](http://docs.oracle.com/javase/7/docs/api/java/io/Reader.html) object            |
| withInputStream | Opens a file for reading and lets you access it with an [InputStream](http://docs.oracle.com/javase/7/docs/api/java/io/InputStream.html) object |
| newReader       | Returns a [Reader](http://docs.oracle.com/javase/7/docs/api/java/io/Reader.html) object to read a text file                                     |
| newInputStream  | Returns an [InputStream](http://docs.oracle.com/javase/7/docs/api/java/io/InputStream.html) object to read a binary file                        |

Read the Java documentation for [Reader](http://docs.oracle.com/javase/7/docs/api/java/io/Reader.html) and [InputStream](http://docs.oracle.com/javase/7/docs/api/java/io/InputStream.html) classes to learn more about methods available for reading data from files.

### Advanced file writing operations

The `Writer` and `OutputStream` classes provide fine control for writing text and binary files, respectively, including low-level operations for single characters or bytes, and support for big files.

For example, given two file objects `sourceFile` and `targetFile`, the following code copies the first file's content into the second file, replacing all `U` characters with `X`:

```groovy
sourceFile.withReader { source ->
    targetFile.withWriter { target ->
        String line
        while( line=source.readLine() ) {
            target << line.replaceAll('U','X')
        }
    }
}
```

Here are the most important methods for writing to files:

| Name             | Description                                                                                                                                             |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| setText          | Writes a string value to a file                                                                                                                         |
| setBytes         | Writes a byte array to a file                                                                                                                           |
| write            | Writes a string to a file, replacing any existing content                                                                                               |
| append           | Appends a string value to a file without replacing existing content                                                                                     |
| newWriter        | Creates a [Writer](http://docs.oracle.com/javase/7/docs/api/java/io/Writer.html) object that allows you to save text data to a file                     |
| newPrintWriter   | Creates a [PrintWriter](http://docs.oracle.com/javase/7/docs/api/java/io/PrintWriter.html) object that allows you to write formatted text to a file     |
| newOutputStream  | Creates an [OutputStream](http://docs.oracle.com/javase/7/docs/api/java/io/OutputStream.html) object that allows you to write binary data to a file     |
| withWriter       | Applies the specified closure to a [Writer](http://docs.oracle.com/javase/7/docs/api/java/io/Writer.html) object, closing it when finished              |
| withPrintWriter  | Applies the specified closure to a [PrintWriter](http://docs.oracle.com/javase/7/docs/api/java/io/PrintWriter.html) object, closing it when finished    |
| withOutputStream | Applies the specified closure to an [OutputStream](http://docs.oracle.com/javase/7/docs/api/java/io/OutputStream.html) object, closing it when finished |

Read the Java documentation for the [Writer](http://docs.oracle.com/javase/7/docs/api/java/io/Writer.html), [PrintWriter](http://docs.oracle.com/javase/7/docs/api/java/io/PrintWriter.html) and [OutputStream](http://docs.oracle.com/javase/7/docs/api/java/io/OutputStream.html) classes to learn more about methods available for writing data to files.

### List directory content

Let's assume that you need to walk through a directory of your choice. You can define the `myDir` variable that points to it:

```groovy
myDir = file('any/path')
```

The simplest way to get a directory list is by using the methods `list` or `listFiles`, which return a collection of first-level elements (files and directories) of a directory:

```groovy
allFiles = myDir.list()
for( def file : allFiles ) {
    println file
}
```

:::{note}
The only difference between `list` and `listFiles` is that the former returns a list of strings, and the latter returns a list of file objects that allow you to access file metadata (size, last modified time, etc).
:::

The `eachFile` method allows you to iterate through the first-level elements only (just like `listFiles`). As with other `each-` methods, `eachFiles` takes a closure as a parameter:

```groovy
myDir.eachFile { item ->
    if( item.isFile() ) {
        println "${item.getName()} - size: ${item.size()}"
    }
    else if( item.isDirectory() ) {
        println "${item.getName()} - DIR"
    }
}
```

Several variants of the above method are available. See the table below for a complete list.

| Name            | Description                                                                                                                                                                                                    |
| --------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| eachFile        | Iterates through first-level elements (files and directories). [Read more](<http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachFile(groovy.io.FileType,%20groovy.lang.Closure)>)         |
| eachDir         | Iterates through first-level directories only. [Read more](<http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachDir(groovy.lang.Closure)>)                                                |
| eachFileMatch   | Iterates through files and dirs whose names match the given filter. [Read more](<http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachFileMatch(java.lang.Object,%20groovy.lang.Closure)>) |
| eachDirMatch    | Iterates through directories whose names match the given filter. [Read more](<http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachDirMatch(java.lang.Object,%20groovy.lang.Closure)>)     |
| eachFileRecurse | Iterates through directory elements depth-first. [Read more](<http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachFileRecurse(groovy.lang.Closure)>)                                      |
| eachDirRecurse  | Iterates through directories depth-first (regular files are ignored). [Read more](<http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachDirRecurse(groovy.lang.Closure)>)                  |

See also: Channel {ref}`channel-path` method.

### Create directories

Given a file variable representing a nonexistent directory, like the following:

```groovy
myDir = file('any/path')
```

the `mkdir` method creates a directory at the given path, returning `true` if the directory is created successfully, and `false` otherwise:

```groovy
result = myDir.mkdir()
println result ? "OK" : "Cannot create directory: $myDir"
```

:::{note}
If the parent directories do not exist, the above method will fail and return `false`.
:::

The `mkdirs` method creates the directory named by the file object, including any nonexistent parent directories:

```groovy
myDir.mkdirs()
```

### Create links

Given a file, the `mklink` method creates a *file system link* for that file using the path specified as a parameter:

```groovy
myFile = file('/some/path/file.txt')
myFile.mklink('/user/name/link-to-file.txt')
```

Table of optional parameters:

| Name      | Description                                                                                                                                                                                                             |
| --------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| hard      | When `true` creates a *hard* link, otherwise creates a *soft* (aka *symbolic*) link. (default: `false`)                                                                                                                 |
| overwrite | When `true` overwrites any existing file with the same name, otherwise throws a [FileAlreadyExistsException](http://docs.oracle.com/javase/8/docs/api/java/nio/file/FileAlreadyExistsException.html) (default: `false`) |

### Copy files

The `copyTo` method copies a file into a new file or into a directory, or copies a directory to a new directory:

```groovy
myFile.copyTo('new_name.txt')
```

:::{note}
If the target file already exists, it will be replaced by the new one. Note also that, if the target is a directory, the source file will be copied into that directory, maintaining the file's original name.
:::

When the source file is a directory, all its content is copied to the target directory:

```groovy
myDir = file('/some/path')
myDir.copyTo('/some/new/path')
```

If the target path does not exist, it will be created automatically.

:::{note}
The `copyTo` method mimics the semantics of the Linux command `cp -r <source> <target>`, with the following caveat: while Linux tools often treat paths ending with a slash (e.g. `/some/path/name/`) as directories, and those not (e.g. `/some/path/name`) as regular files, Nextflow (due to its use of the Java files API) views both these paths as the same file system object. If the path exists, it is handled according to its actual type (i.e. as a regular file or as a directory). If the path does not exist, it is treated as a regular file, with any missing parent directories created automatically.
:::

### Move files

You can move a file by using the `moveTo` method:

```groovy
myFile = file('/some/path/file.txt')
myFile.moveTo('/another/path/new_file.txt')
```

:::{note}
When a file with the same name as the target already exists, it will be replaced by the source. Note also that, when the target is a directory, the file will be moved to (or within) that directory, maintaining the file's original name.
:::

When the source is a directory, all the directory content is moved to the target directory:

```groovy
myDir = file('/any/dir_a')
myDir.moveTo('/any/dir_b')
```

Please note that the result of the above example depends on the existence of the target directory. If the target directory exists, the source is moved into the target directory, resulting in the path:

```
/any/dir_b/dir_a
```

If the target directory does not exist, the source is just renamed to the target name, resulting in the path:

```
/any/dir_b
```

:::{note}
The `moveTo` method mimics the semantics of the Linux command `mv <source> <target>`, with the same caveat as that given above for `copyTo`.
:::

### Rename files

You can rename a file or directory by simply using the `renameTo` method:

```groovy
myFile = file('my_file.txt')
myFile.renameTo('new_file_name.txt')
```

### Delete files

The `delete` method deletes the file or directory at the given path, returning `true` if the operation succeeds, and `false` otherwise:

```groovy
myFile = file('some/file.txt')
result = myFile.delete()
println result ? "OK" : "Cannot delete: $myFile"
```

:::{note}
This method deletes a directory **only** if it does not contain any files or sub-directories. To delete a directory and **all** its contents (i.e. removing all the files and sub-directories it may contain), use the method `deleteDir`.
:::

### Check file attributes

The following methods can be used on a file variable created by using the `file` method:

| Name          | Description                                                                          |
| ------------- | ------------------------------------------------------------------------------------ |
| getName       | Gets the file name e.g. `/some/path/file.txt` -> `file.txt`                          |
| getBaseName   | Gets the file name without its extension e.g. `/some/path/file.tar.gz` -> `file.tar` |
| getSimpleName | Gets the file name without any extension e.g. `/some/path/file.tar.gz` -> `file`     |
| getExtension  | Gets the file extension e.g. `/some/path/file.txt` -> `txt`                          |
| getParent     | Gets the file parent path e.g. `/some/path/file.txt` -> `/some/path`                 |
| size          | Gets the file size in bytes                                                          |
| exists        | Returns `true` if the file exists, or `false` otherwise                              |
| isEmpty       | Returns `true` if the file is zero length or does not exist, `false` otherwise       |
| isFile        | Returns `true` if it is a regular file e.g. not a directory                          |
| isDirectory   | Returns `true` if the file is a directory                                            |
| isHidden      | Returns `true` if the file is hidden                                                 |
| lastModified  | Returns the file last modified timestamp i.e. a long as Linux epoch time             |

For example, the following line prints a file name and size:

```groovy
println "File ${myFile.getName()} size: ${myFile.size()}"
```

:::{tip}
The invocation of any method name starting with the `get` prefix can be shortcut by omitting the `get` prefix and `()` parentheses. Therefore, writing `myFile.getName()` is exactly the same as `myFile.name` and `myFile.getBaseName()` is the same as `myFile.baseName` and so on.
:::

### Get and modify file permissions

Given a file variable representing a file (or directory), the `getPermissions` method returns a 9-character string representing the file's permissions using the [Linux symbolic notation](http://en.wikipedia.org/wiki/File_system_permissions#Symbolic_notation), e.g. `rw-rw-r--`:

```groovy
permissions = myFile.getPermissions()
```

Similarly, the `setPermissions` method sets the file's permissions using the same notation:

```groovy
myFile.setPermissions('rwxr-xr-x')
```

A second version of the `setPermissions` method sets a file's permissions given three digits representing, respectively, the `owner`, `group` and `other` permissions:

```groovy
myFile.setPermissions(7,5,5)
```

Learn more about [File permissions numeric notation](http://en.wikipedia.org/wiki/File_system_permissions#Numeric_notation).

### HTTP/FTP files

Nextflow provides transparent integration of HTTP/S and FTP protocols for handling remote resources as local file system objects. Simply specify the resource URL as the argument of the `file` object:

```groovy
pdb = file('http://files.rcsb.org/header/5FID.pdb')
```

Then, you can access it as a local file as described in the previous sections:

```groovy
println pdb.text
```

The above one-liner prints the content of the remote PDB file. Previous sections provide code examples showing how to stream or copy the content of files.

:::{note}
Write and list operations are not supported for HTTP/S and FTP files.
:::

### Counting records

#### countLines

The `countLines` method counts the lines in a text files.

```groovy
def sample = file('/data/sample.txt')
println sample.countLines()
```

Files whose name ends with the `.gz` suffix are expected to be GZIP compressed and automatically uncompressed.

#### countFasta

The `countFasta` method counts the number of records in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) formatted file.

```groovy
def sample = file('/data/sample.fasta')
println sample.countFasta()
```

Files whose name ends with the `.gz` suffix are expected to be GZIP compressed and automatically uncompressed.

#### countFastq

The `countFastq` method counts the number of records in a [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) formatted file.

```groovy
def sample = file('/data/sample.fastq')
println sample.countFastq()
```

Files whose name ends with the `.gz` suffix are expected to be GZIP compressed and automatically uncompressed.

[glob]: http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob

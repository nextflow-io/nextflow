(script-page)=

# Scripts

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

:::{warning}
Nextflow scripts have a maximum size of 64 KiB. To avoid this limit for large pipelines, consider moving pipeline components into separate files and including them as modules.
:::

## Groovy basics

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
- [Groovy List API](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/List.html)
- [Java List API](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/util/List.html)

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
Copying a map with the `+` operator is a safer way to modify maps in Nextflow, specifically when passing maps through channels. This way, a new instance of the map will be created, and any references to the original map won't be affected.
:::

Learn more about maps:

- [Groovy Maps tutorial](http://groovy-lang.org/groovy-dev-kit.html#Collections-Maps)
- [Groovy Map API](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html)
- [Java Map API](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/util/Map.html)

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

:::{warning}
When using backslashes to continue a multi-line command, make sure to not put any spaces after the backslash, otherwise it will be interpreted by the Groovy lexer as an escaped space instead of a backslash, which will make your script incorrect. It will also print this warning:

```
unknown recognition error type: groovyjarjarantlr4.v4.runtime.LexerNoViableAltException
```
:::

(script-regexp)=

### Regular expressions

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

Regular expression support is imported from Java. Java's regular expression language and API is documented in the [Pattern](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/util/regex/Pattern.html) class.

You may also be interested in this post: [Groovy: Don't Fear the RegExp](https://web.archive.org/web/20170621185113/http://www.naleid.com/blog/2008/05/19/dont-fear-the-regexp).

#### String replacement

To replace pattern occurrences in a given string, use the `replaceFirst` and `replaceAll` methods:

```groovy
x = "colour".replaceFirst(/ou/, "o")
println x
// prints: color

y = "cheesecheese".replaceAll(/cheese/, "nice")
println y
// prints: nicenice
```

#### Capturing groups

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

#### Removing part of a string

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

### Functions

Functions can be defined using the following syntax:

```groovy
def <function name> ( arg1, arg, .. ) {
    <function body>
}
```

For example:

```groovy
def foo() {
    'Hello world'
}

def bar(alpha, omega) {
    alpha + omega
}
```

The above snippet defines two simple functions, that can be invoked in the workflow script as `foo()`, which returns `'Hello world'`, and `bar(10, 20)`, which returns the sum of two parameters (`30` in this case).

Functions implicitly return the result of the last statement. Additionally, the `return` keyword can be used to explicitly exit from a function and return the specified value. For example:

```groovy
def fib( x ) {
    if( x <= 1 )
        return x

    fib(x-1) + fib(x-2)
}
```

(script-closure)=

### Closures

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

Closures can also access variables outside of their scope, and they can be used anonymously, that is without assigning them to a variable. Here is an example that demonstrates both of these things:

```groovy
myMap = ["China": 1, "India": 2, "USA": 3]

result = 0
myMap.keySet().each { result += myMap[it] }

println result
```

A closure can also declare local variables that exist only for the lifetime of the closure:

```groovy
result = 0
myMap.keySet().each {
  def count = myMap[it]
  result += count
}
```

:::{warning}
Local variables should be declared using a qualifier such as `def` or a type name, otherwise they will be interpreted as global variables, which could lead to a {ref}`race condition <cache-global-var-race-condition>`.
:::

Learn more about closures in the [Groovy documentation](http://groovy-lang.org/closures.html)

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

(implicit-functions)=

## Implicit functions

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

`file( filePattern, options = [:] )`
: Get one or more files from a path or glob pattern. Returns a [Path](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/nio/file/Path.html) or list of Paths if there are multiple files. See [Files and I/O](#files-and-io).

`files( filePattern, options = [:] )`
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

(implicit-classes)=

## Implicit classes

The following classes are imported by default in Nextflow scripts:

- `java.lang.*`
- `java.util.*`
- `java.io.*`
- `java.net.*`
- `groovy.lang.*`
- `groovy.util.*`
- `java.math.BigInteger`
- `java.math.BigDecimal`
- `java.nio.file.Path`

Additionally, Nextflow imports several new classes which are described below.

### Channel

The `Channel` class provides the channel factory methods. See {ref}`channel-factory` for more information.

(implicit-classes-duration)=

### Duration

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

```groovy
// integer value (milliseconds)
oneSecond = Duration.of(1000)

// simple string value
oneHour = Duration.of('1h')

// complex string value
complexDuration = Duration.of('1day 6hours 3minutes 30seconds')
```

Durations can be compared like numbers, and they support basic arithmetic operations:

```groovy
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

(implicit-classes-memoryunit)=

### MemoryUnit

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

```groovy
// integer value (bytes)
oneKilobyte = MemoryUnit.of(1024)

// string value
oneGigabyte = MemoryUnit.of('1 GB')
```

Memory units can be compared like numbers, and they support basic arithmetic operations:

```groovy
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

### ValueObject

`ValueObject` is an AST transformation for classes and enums, which simply combines [AutoClone](http://docs.groovy-lang.org/latest/html/gapi/groovy/transform/AutoClone.html) and [Immutable](https://docs.groovy-lang.org/latest/html/gapi/groovy/transform/Immutable.html). It is useful for defining custom "record" types.

(script-file-io)=

## Files and I/O

### Opening files

To access and work with files, use the `file()` method, which returns a file system object given a file path string:

```groovy
myFile = file('some/path/to/my_file.file')
```

The `file()` method can reference both files and directories, depending on what the string path refers to in the file system.

When using the wildcard characters `*`, `?`, `[]` and `{}`, the argument is interpreted as a [glob](http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob) path matcher and the `file()` method returns a list object holding the paths of files whose names match the specified pattern, or an empty list if no match is found:

```groovy
listOfFiles = file('some/path/*.fa')
```

:::{note}
The `file()` method does not return a list if only one file is matched. Use the `files()` method to always return a list.
:::

:::{note}
A double asterisk (`**`) in a glob pattern works like `*` but also searches through subdirectories.
:::

By default, wildcard characters do not match directories or hidden files. For example, if you want to include hidden files in the result list, enable the `hidden` option:

```groovy
listWithHidden = file('some/path/*.fa', hidden: true)
```

:::{note}
To compose paths, instead of string interpolation, use the `resolve()` method or the `/` operator:

```groovy
def dir = file('s3://bucket/some/data/path')
def sample1 = dir.resolve('sample.bam')         // correct
def sample2 = dir / 'sample.bam'
def sample3 = file("$dir/sample.bam")           // correct (but verbose)
def sample4 = "$dir/sample.bam"                 // incorrect
```
:::

The following options are available:

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

See also: {ref}`Channel.fromPath <channel-path>`.

### Getting file attributes

The `file()` method returns a [Path](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/nio/file/Path.html), so any method defined for Path can also be used in a Nextflow script.

Additionally, the following methods are also defined for Paths in Nextflow:

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
  ```groovy
  def ref = file('s3://some-bucket/foo.txt')

  assert ref.toString() == '/some-bucket/foo.txt'
  assert "$ref" == '/some-bucket/foo.txt'
  assert ref.toUriString() == 's3://some-bucket/foo.txt'
  ```

:::{tip}
In Groovy, any method that looks like `get*()` can also be accessed as a field. For example, `myFile.getName()` is equivalent to `myFile.name`, `myFile.getBaseName()` is equivalent to `myFile.baseName`, and so on.
:::

### Reading and writing

#### Reading and writing an entire file

Given a file variable, created with the `file()` method as shown previously, reading a file is as easy as getting the file's `text` property, which returns the file content as a string:

```groovy
print myFile.text
```

Similarly, you can save a string to a file by assigning it to the file's `text` property:

```groovy
myFile.text = 'Hello world!'
```

Binary data can managed in the same way, just using the file property `bytes` instead of `text`. Thus, the following example reads the file and returns its content as a byte array:

```groovy
binaryContent = myFile.bytes
```

Or you can save a byte array to a file:

```groovy
myFile.bytes = binaryContent
```

:::{note}
The above assignment overwrites any existing file contents, and implicitly creates the file if it doesn't exist.
:::

:::{warning}
The above methods read and write the **entire** file contents at once, in a single variable or buffer. For this reason, when dealing with large files it is recommended that you use a more memory efficient approach, such as reading/writing a file line by line or using a fixed size buffer.
:::

#### Appending to a file

In order to append a string value to a file without erasing existing content, you can use the `append()` method:

```groovy
myFile.append('Add this line\n')
```

Or use the left shift operator, a more idiomatic way to append text content to a file:

```groovy
myFile << 'Add a line more\n'
```

#### Reading a file line by line

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

To process a big file, use the method `eachLine()`, which reads only a single line at a time into memory:

```groovy
count = 0
myFile.eachLine { str ->
    println "line ${count++}: $str"
}
```

#### Advanced file reading

The classes `Reader` and `InputStream` provide fine-grained control for reading text and binary files, respectively.

The method `newReader()` creates a [Reader](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/io/Reader.html) object for the given file that allows you to read the content as single characters, lines or arrays of characters:

```groovy
myReader = myFile.newReader()
String line
while( line = myReader.readLine() ) {
    println line
}
myReader.close()
```

The method `withReader()` works similarly, but automatically calls the `close()` method for you when you have finished processing the file. So, the previous example can be written more simply as:

```groovy
myFile.withReader {
    String line
    while( line = it.readLine() ) {
        println line
    }
}
```

The methods `newInputStream()` and `withInputStream()` work similarly. The main difference is that they create an [InputStream](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/io/InputStream.html) object useful for writing binary data.

The following methods are useful for reading files:

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

#### Advanced file writing

The `Writer` and `OutputStream` classes provide fine-grained control for writing text and binary files, respectively, including low-level operations for single characters or bytes, and support for big files.

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

  ```groovy
  file('/some/path/my_file.txt').copyTo('/another/path/new_file.txt')
  ```

: *When copying a file to a directory:* the file will be copied into the directory, replacing any file with the same name.

  ```groovy
  file('/some/path/my_file.txt').copyTo('/another/path')
  ```

: *When copying a directory to another directory:* if the target directory already exists, the source directory will be copied into the target directory, replacing any sub-directory with the same name. If the target path does not exist, it will be created automatically.

  ```groovy
  file('/any/dir_a').moveTo('/any/dir_b')
  ```

  The result of the above example depends on the existence of the target directory. If the target directory exists, the source is moved into the target directory, resulting in the path `/any/dir_b/dir_a`. If the target directory does not exist, the source is just renamed to the target name, resulting in the path `/any/dir_b`.

: :::{note}
  The `copyTo()` method follows the semantics of the Linux command `cp -r <source> <target>`, with the following caveat: while Linux tools often treat paths ending with a slash (e.g. `/some/path/name/`) as directories, and those not (e.g. `/some/path/name`) as regular files, Nextflow (due to its use of the Java files API) views both of these paths as the same file system object. If the path exists, it is handled according to its actual type (i.e. as a regular file or as a directory). If the path does not exist, it is treated as a regular file, with any missing parent directories created automatically.
  :::

`delete()`
: Deletes the file or directory at the given path, returning `true` if the operation succeeds, and `false` otherwise:

  ```groovy
  myFile = file('some/file.txt')
  result = myFile.delete()
  println result ? "OK" : "Cannot delete: $myFile"
  ```

  If a directory is not empty, it will not be deleted and `delete()` will return `false`.

`deleteDir()`
: Deletes a directory and all of its contents.

  ```groovy
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

  ```groovy
  myDir = file('any/path')
  result = myDir.mkdir()
  println result ? "OK" : "Cannot create directory: $myDir"
  ```

  If the parent directories do not exist, the directory will not be created and `mkdir()` will return `false`.

`mkdirs()`
: Creates a directory at the given path, including any nonexistent parent directories:

  ```groovy
  file('any/path').mkdirs()
  ```

`mklink( linkName, options = [:] )`
: Creates a *filesystem link* to a given path:

  ```groovy
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

  ```groovy
  file('my_file.txt').renameTo('new_file_name.txt')
  ```

`setPermissions( permissions )`
: Sets a file's permissions using the [symbolic notation](http://en.wikipedia.org/wiki/File_system_permissions#Symbolic_notation):

  ```groovy
  myFile.setPermissions('rwxr-xr-x')
  ```

`setPermissions( owner, group, other )`
: Sets a file's permissions using the [numeric notation](http://en.wikipedia.org/wiki/File_system_permissions#Numeric_notation), i.e. as three digits representing the **owner**, **group**, and **other** permissions:

  ```groovy
  myFile.setPermissions(7,5,5)
  ```

#### Listing directories

The simplest way to list a directory is to use `list()` or `listFiles()`, which return a collection of first-level elements (files and directories) of a directory:

```groovy
for( def file : file('any/path').list() ) {
    println file
}
```

Additionally, the `eachFile()` method allows you to iterate through the first-level elements only (just like `listFiles()`). As with other `each*()` methods, `eachFile()` takes a closure as a parameter:

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

The following methods are available for listing and traversing directories:

`eachDir( closure )`
: Iterates through first-level directories only. [Read more](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachDir(groovy.lang.Closure))

`eachDirMatch( nameFilter, closure )`
: Iterates through directories whose names match the given filter. [Read more](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachDirMatch(java.lang.Object,%20groovy.lang.Closure))

`eachDirRecurse( closure )`
: Iterates through directories depth-first (regular files are ignored). [Read more](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachDirRecurse(groovy.lang.Closure))

`eachFile( closure )`
: Iterates through first-level files and directories. [Read more](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachFile(groovy.lang.Closure))

`eachFileMatch( nameFilter, closure )`
: Iterates through files and directories whose names match the given filter. [Read more](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachFileMatch(java.lang.Object,%20groovy.lang.Closure))

`eachFileRecurse( closure )`
: Iterates through files and directories depth-first. [Read more](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachFileRecurse(groovy.lang.Closure))

See also: {ref}`Channel.fromPath <channel-path>`.

### Fetching HTTP/FTP files

Nextflow integrates seamlessly with the HTTP(S) and FTP protocols for handling remote resources the same as local files. Simply specify the resource URL when opening the file:

```groovy
pdb = file('http://files.rcsb.org/header/5FID.pdb')
```

Then, you can access it as a local file as described previously:

```groovy
println pdb.text
```

The above one-liner prints the content of the remote PDB file. Previous sections provide code examples showing how to stream or copy the content of files.

:::{note}
Write and list operations are not supported for HTTP(S) and FTP files.
:::

### Splitting and counting records

The following methods are defined for Paths for splitting and counting records:

`countFasta()`
: Counts the number of records in a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file. See the {ref}`operator-splitfasta` operator for available options.

`countFastq()`
: Counts the number of records in a [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file. See the {ref}`operator-splitfastq` operator for available options.

`countJson()`
: Counts the number of records in a JSON file. See the {ref}`operator-splitjson` operator for available options.

`countLines()`
: Counts the number of lines in a text file. See the {ref}`operator-splittext` operator for available options.

`splitFasta()`
: Splits a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file into a list of records. See the {ref}`operator-splitfasta` operator for available options.

`splitFastq()`
: Splits a [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file into a list of records. See the {ref}`operator-splitfastq` operator for available options.

`splitJson()`
: Splits a JSON file into a list of records. See the {ref}`operator-splitjson` operator for available options.

`splitLines()`
: Splits a text file into a list of lines. See the {ref}`operator-splittext` operator for available options.

(script-page)=

# Scripts

Nextflow is a workflow language that runs on the Java virtual machine (JVM). Nextflow's syntax is very similar to [Groovy](https://groovy-lang.org/), a scripting language for the JVM. However, Nextflow is specialized for writing computational pipelines in a declarative manner. See {ref}`syntax-page` for a full description of the Nextflow language.

Nextflow scripts can also make full use of the Java and Groovy standard libraries. See {ref}`stdlib-page` for more information.

:::{warning}
Nextflow uses UTF-8 as the default character encoding for source files. Make sure to use UTF-8 encoding when editing Nextflow scripts with your preferred text editor.
:::

:::{warning}
Nextflow scripts have a maximum size of 64 KiB. To avoid this limit for large pipelines, consider moving pipeline components into separate files and including them as modules.
:::

## Hello world

You can use the `println` function to print to the console:

```nextflow
println 'Hello, World!'
```


## Variables

Variables are declared using the `def` keyword:

```nextflow
def num = 1
println num

def date = new java.util.Date()
println date

def x = -3.1499392
println x

def flag = false
println flag

def str = "Hi"
println str
```

:::{warning}
Variables can also be declared without `def` in some cases. However, this practice is discouraged outside of simple code snippets because it can lead to a {ref}`race condition <cache-global-var-race-condition>`.
:::

## Lists

Lists are defined using square brackets:

```nextflow
def myList = [1776, -1, 33, 99, 0, 928734928763]
```

You can access a given item in the list with square-bracket notation (indexes start at 0):

```nextflow
println myList[0]
```

In order to get the length of the list use the `size` method:

```nextflow
println myList.size()
```

Refer to the [Java](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/util/List.html) and [Groovy](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/List.html) standard libraries for the set of available list operations.

## Maps

Maps are used to store *associative arrays* (also known as *dictionaries*). They are unordered collections of heterogeneous, named data:

```nextflow
def scores = ["Brett": 100, "Pete": "Did not finish", "Andrew": 86.87934]
```

Note that each of the values stored in the map can be of a different type. `Brett` is an integer, `Pete` is a string, and `Andrew` is a floating-point number.

We can access the values in a map in two main ways:

```nextflow
println scores["Pete"]
println scores.Pete
```

To add data to or modify a map, the syntax is similar to adding values to list:

```nextflow
scores["Pete"] = 3
scores["Cedric"] = 120
```

You can also use the `+` operator to add two maps together:

```nextflow
def new_scores = scores + ["Pete": 3, "Cedric": 120]
```

When adding two maps, the first map is copied and then appended with the keys from the second map. Any conflicting keys are overwritten by the second map.

:::{tip}
Copying a map with the `+` operator is a safer way to modify maps in Nextflow, specifically when passing maps through channels. This way, a new instance of the map will be created, and any references to the original map won't be affected.
:::

See the [Java](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/util/Map.html) and [Groovy](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html) standard libraries for the set of available map operations.

## Conditional execution

One of the most important features of any programming language is the ability to execute different code under different conditions. The simplest way to do this is to use the `if` construct:

```nextflow
def x = Math.random()
if( x < 0.5 ) {
    println "You lost."
}
else {
    println "You won!"
}
```

## Strings

Strings can be defined by enclosing text in single or double quotes (`'` or `"` characters):

```nextflow
println "he said 'cheese' once"
println 'he said "cheese!" again'
```

Strings can be concatenated with `+`:

```nextflow
def a = "world"
print "hello " + a + "\n"
```

(string-interpolation)=

### String interpolation

There is an important difference between single-quoted and double-quoted strings: Double-quoted strings support variable interpolations, while single-quoted strings do not.

In practice, double-quoted strings can contain the value of an arbitrary variable by prefixing its name with the `$` character, or the value of any expression by using the `${expression}` syntax, similar to Bash/shell scripts:

```nextflow
def foxtype = 'quick'
def foxcolor = ['b', 'r', 'o', 'w', 'n']
println "The $foxtype ${foxcolor.join()} fox"

def x = 'Hello'
println '$x + $y'
```

This code prints:

```
The quick brown fox
$x + $y
```

### Multi-line strings

A block of text that span multiple lines can be defined by delimiting it with triple single or double quotes:

```nextflow
def text = """
    hello there James
    how are you today?
    """
```

:::{note}
Like before, multi-line strings inside double quotes support variable interpolation, while single-quoted multi-line strings do not.
:::

As in Bash/shell scripts, terminating a line in a multi-line string with a `\` character prevents a newline character from separating that line from the one that follows:

```nextflow
def myLongCmdline = """
    blastp \
    -in $input_query \
    -out $output_file \
    -db $blast_database \
    -html
    """

def result = myLongCmdline.execute().text
```

In the preceding example, `blastp` and its `-in`, `-out`, `-db` and `-html` switches and their arguments are effectively a single line.

:::{warning}
Do not put any spaces after the backslash when using backslashes to continue a multi-line command. Spaces after the backslash will be interpreted as an escaped space and will make your script incorrect. It will also print this warning:

```
unknown recognition error type: groovyjarjarantlr4.v4.runtime.LexerNoViableAltException
```
:::

(script-regexp)=

## Regular expressions

Regular expressions are the Swiss Army knife of text processing. They provide the programmer with the ability to match and extract patterns from strings.

Regular expressions are available via the `~/pattern/` syntax and the `=~` and `==~` operators.

Use `=~` to check whether a given pattern occurs anywhere in a string:

```nextflow
assert 'foo' =~ /foo/       // return TRUE
assert 'foobar' =~ /foo/    // return TRUE
```

Use `==~` to check whether a string matches a given regular expression pattern exactly.

```nextflow
assert 'foo' ==~ /foo/       // return TRUE
assert 'foobar' ==~ /foo/    // return FALSE
```

The `~` operator creates a [Pattern](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/util/regex/Pattern.html) from the given string, while the `=~` operator creates a [Matcher](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/util/regex/Matcher.html):

```nextflow
x = ~/abc/
println x.class
// prints java.util.regex.Pattern

y = 'some string' =~ /abc/
println y.class
// prints java.util.regex.Matcher
```

See the linked Java documentation for the available operations for these classes.

### String replacement

To replace pattern occurrences in a given string, use the `replaceFirst` and `replaceAll` methods:

```nextflow
def x = "colour".replaceFirst(/ou/, "o")
println x
// prints: color

def y = "cheesecheese".replaceAll(/cheese/, "nice")
println y
// prints: nicenice
```

To remove part of a string, simply replace it with a blank string:

```nextflow
def z = 'Hello World!'.replaceFirst(/(?i)\s+Wo\w+/, '')
println z
// prints: Hello!
```

### Capturing groups

You can match a pattern that includes groups. First create a matcher object with the `=~` operator. Then, you can index the matcher object to find the matches: `matcher[0]` returns a list representing the first match of the regular expression in the string. The first list element is the string that matches the entire regular expression, and the remaining elements are the strings that match each group.

Here's how it works:

```nextflow
def programVersion = '2.7.3-beta'
def m = programVersion =~ /(\d+)\.(\d+)\.(\d+)-?(.+)/

assert m[0] == ['2.7.3-beta', '2', '7', '3', 'beta']
assert m[0][1] == '2'
assert m[0][2] == '7'
assert m[0][3] == '3'
assert m[0][4] == 'beta'
```

Applying some syntactic sugar, you can do the same in just one line of code:

```nextflow
def programVersion = '2.7.3-beta'
def (full, major, minor, patch, flavor) = (programVersion =~ /(\d+)\.(\d+)\.(\d+)-?(.+)/)[0]

println full    // 2.7.3-beta
println major   // 2
println minor   // 7
println patch   // 3
println flavor  // beta
```

(script-closure)=

## Closures

A closure is a function that can be used like a regular value. Typically, closures are passed as arguments to *higher-order functions* to express computations in a declarative manner.

For example:

```nextflow
def square = { v -> v * v }
```

The above example defines a closure, which takes one parameter named `v` and returns the "square" of `v` (`v * v`). The closure is assigned to the variable `square`.

`square` can now be called like a function:

```nextflow
println square(9)
```

The above example prints `81`.

The main use case for a closure is as an argument to a higher-order function:

```nextflow
[ 1, 2, 3, 4 ].collect(square)
```

The `collect` method of a list applies a mapping function to each value in the list and produces a new list. The above example produces:

```nextflow
[ 1, 4, 9, 16 ]
```

The example can be expressed more concisely as:

```nextflow
[ 1, 2, 3, 4 ].collect { v -> v * v }
```

Another example is the `each` method of a map, which takes a closure with two arguments corresponding to the key and value of each map entry:

```nextflow
[ "Yue" : "Wu", "Mark" : "Williams", "Sudha" : "Kumari" ].each { key, value ->
    println "$key = $value"
}
```

Prints:

```
Yue = Wu
Mark = Williams
Sudha = Kumari
```

Closures can access variables outside of their scope:

```nextflow
def counts = ["China": 1, "India": 2, "USA": 3]

def result = 0
counts.keySet().each { v ->
    result += counts[v]
}

println result
```

A closure can also declare local variables that exist only for the lifetime of each closure invocation:

```nextflow
def result = 0
myMap.keySet().each { v ->
    def count = myMap[v]
    result += count
}
```

While the `each` method is a convenient way to iterate through a collection and build up some result, a more idiomatic way to do this is to use the `inject` method:

```nextflow
def result = counts.values().inject { sum, v -> sum + v }
```

This way, the closure is fully "self-contained" because it doesn't access or mutate any variables outside of its scope.

:::{note}
When a closure takes a single parameter, the parameter can be omitted, in which case the implicit `it` parameter will be used:

```nextflow
[1, 2, 3].each { println it }
```
:::

## Script definitions

So far, we have been focusing on the basic building blocks of Nextflow code, like variables, lists, strings, and closures.

In practice, however, Nextflow scripts are composed of *workflows*, *processes*, and *functions* (collectively known as *definitions*), and can *include*  definitions from other scripts.

To transition a code snippet into a proper workflow script, simply wrap it in a `workflow` block:

```nextflow
workflow {
    println 'Hello!'
}
```

This block is called the *entry workflow*. It serves as the entrypoint when the script is executed. A script can only have one entry workflow. Whenever a script contains only simple statements like `println 'Hello!'`, Nextflow simply treats it as an entry workflow.

You can also break up code into functions, for example:

```nextflow
def sayHello() {
    println 'Hello!'
}

def add(a, b) {
    a + b
}

workflow {
    sayHello()
    println "2 + 2 = ${add(2, 2)}!"
}
```

See {ref}`workflow-page`, {ref}`process-page`, and {ref}`module-page`  for more information about how to use these features in your Nextflow scripts.

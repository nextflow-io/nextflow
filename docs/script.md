(script-page)=

# Scripts

Nextflow is a workflow language that runs on the Java virtual machine (JVM). Nextflow's syntax is very similar to [Groovy](https://groovy-lang.org/), a scripting language for the JVM, but Nextflow is specialized for writing computational pipelines in a declarative manner.

Nextflow scripts can also make full use of the Java and Groovy standard libraries; see the {ref}`stdlib-page` page for more information.

:::{warning}
Nextflow uses UTF-8 as the default character encoding for source files. Make sure to use UTF-8 encoding when editing Nextflow scripts with your preferred text editor.
:::

:::{warning}
Nextflow scripts have a maximum size of 64 KiB. To avoid this limit for large pipelines, consider moving pipeline components into separate files and including them as modules.
:::

## Hello world

To print something is as easy as using one of the `print` or `println` methods.

```groovy
println "Hello, World!"
```

The only difference between the two is that the `println` method implicitly appends a newline character to the printed string.

## Variables

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

## Lists

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

## Maps

Maps are used to store *associative arrays* (also known as *dictionaries*). They are unordered collections of heterogeneous, named data:

```groovy
scores = ["Brett": 100, "Pete": "Did not finish", "Andrew": 86.87934]
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

## Multiple assignment

An array or a list object can used to assign to multiple variables at once:

```groovy
(a, b, c) = [10, 20, 'foo']
assert a == 10 && b == 20 && c == 'foo'
```

The three variables on the left of the assignment operator are initialized by the corresponding item in the list.

Read more about [Multiple assignment](http://www.groovy-lang.org/semantics.html#_multiple_assignment) in the Groovy documentation.

## Conditional execution

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

## Strings

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

## String interpolation

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

## Multi-line strings

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
When using backslashes to continue a multi-line command, make sure to not put any spaces after the backslash, otherwise it will be interpreted as an escaped space instead of a backslash, which will make your script incorrect. It will also print this warning:

```
unknown recognition error type: groovyjarjarantlr4.v4.runtime.LexerNoViableAltException
```
:::

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

Regular expression support is imported from Java. Java's regular expression language and API is documented in the [Pattern](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/util/regex/Pattern.html) class.

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

To remove part of a string, simply replace it with a blank string:

```groovy
z = 'Hello World!'.replaceFirst(/(?i)\s+Wo\w+/, '')
println z
// prints: Hello!
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

(script-functions)=

## Functions

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

## Closures

A closure is a function that can be used like a regular value. Typically, closures are passed as arguments to *higher-order functions* to express computations in a declarative manner.

For example:

```groovy
square = { v -> v * v }
```

The above example defines a closure, which takes one parameter named `v` and returns the "square" of `v` (`v * v`), and assigns the closure to the variable `square`.

Now we can call `square` like a function:

```groovy
println square(9)
```

which prints `81`.

The main use case for a closure, however, is as an argument to a higher-order function:

```groovy
[ 1, 2, 3, 4 ].collect(square)
```

The `collect` method of a list applies a mapping function to each value in the list and produces a new list. The above example produces:

```groovy
[ 1, 4, 9, 16 ]
```

The example can be expressed more concisely as:

```groovy
[ 1, 2, 3, 4 ].collect { v -> v * v }
```

Another example is the `each` method of a map, which takes a closure with two arguments corresponding to the key and value of each map entry:

```groovy
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

```groovy
counts = ["China": 1, "India": 2, "USA": 3]

result = 0
counts.keySet().each { v ->
    result += counts[v]
}

println result
```

A closure can also declare local variables that exist only for the lifetime of each closure invocation:

```groovy
result = 0
myMap.keySet().each { v ->
    def count = myMap[v]
    result += count
}
```

:::{warning}
Local variables should be declared using `def`, otherwise they will be interpreted as global variables, which could lead to a {ref}`race condition <cache-global-var-race-condition>`.
:::

While the `each` method is a convenient way to iterate through a collection and build up some result, a more idiomatic way to do this is to use the `inject` method:

```groovy
result = counts.values().inject { sum, v -> sum + v }
```

This way, the closure is fully "self-contained" because it doesn't access or mutate any variables outside of its scope.

Learn more about closures in the [Groovy documentation](http://groovy-lang.org/closures.html)

## Syntax sugar

Nextflow provides several forms of "syntax sugar", or shorthands that can make your code easier to read.

Some programming languages require every statement to be terminated by a semi-colon. In Nextflow, semi-colons are optional, but they can still be used to write multiple statements on the same line:

```groovy
println 'Hello!' ; println 'Hello again!'
```

When calling a function, the parentheses around the function arguments are optional:

```groovy
// full syntax
printf('Hello %s!\n', 'World')

// shorthand
printf 'Hello %s!\n', 'World'
```

It is especially useful when calling a function with a closure parameter:

```groovy
// full syntax
[1, 2, 3].each({ v -> println v })

// shorthand
[1, 2, 3].each { v -> println v }
```

If the last argument is a closure, the closure can be written outside of the parentheses:

```groovy
// full syntax
[1, 2, 3].inject('result:', { acc, v -> acc + ' ' + v })

// shorthand
[1, 2, 3].inject('result:') { acc, v -> acc + ' ' + v }
```

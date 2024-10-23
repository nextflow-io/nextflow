(script-page)=

# Scripts

Nextflow is a domain-specific language (DSL) based on Groovy, a general-purpose programming language for the Java virtual machine. Nextflow extends the Groovy syntax with features that ease the writing of computational pipelines in a declarative manner.

For more background on Groovy, refer to these resources:

- [Groovy User Guide](http://groovy-lang.org/documentation.html)
- [Groovy Cheat sheet](http://www.cheat-sheets.org/saved-copy/rc015-groovy_online.pdf)
- [Groovy in Action](http://www.manning.com/koenig2/)

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
When using backslashes to continue a multi-line command, make sure to not put any spaces after the backslash, otherwise it will be interpreted by the Groovy lexer as an escaped space instead of a backslash, which will make your script incorrect. It will also print this warning:

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

## Syntax sugar

Groovy provides several forms of "syntax sugar", or shorthands that can make your code easier to read.

Some programming languages require every statement to be terminated by a semi-colon. In Groovy, semi-colons are optional, but they can still be used to write multiple statements on the same line:

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
[1, 2, 3].each({ println it })

// shorthand
[1, 2, 3].each { println it }
```

If the last argument is a closure, the closure can be written outside of the parentheses:

```groovy
// full syntax
[1, 2, 3].inject('result:', { accum, v -> accum + ' ' + v })

// shorthand
[1, 2, 3].inject('result:') { accum, v -> accum + ' ' + v }
```

:::{note}
In some cases, you might not be able to omit the parentheses because it would be syntactically ambiguous. You can use the `groovysh` REPL console to play around with Groovy and figure out what works.
:::

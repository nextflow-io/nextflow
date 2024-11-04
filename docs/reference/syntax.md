(syntax-page)=

# Syntax

This page provides a comprehensive description of the Nextflow language.

## Comments

A line comment starts with `//` and includes the rest of the line.

```nextflow
println 'Hello world!' // line comment
```

A block comment starts with `/*` and includes all subsequent characters up to the first `*/`.

```nextflow
/*
 * block comment
 */
println 'Hello again!'
```

## Top-level declarations

A Nextflow script may contain the following top-level declarations:

- Shebang
- Feature flags
- Includes
- Parameter definitions
- Workflow definitions
- Process definitions
- Function definitions
- Enum types
- Output block

Script declarations are in turn composed of statements and expressions.

If there are no top-level declarations, a script may contain one or more [statements](#statements), in which case the entire script is treated as an entry workflow. For example:

```nextflow
println 'Hello world!'
```

Is equivalent to:

```nextflow
workflow {
    println 'Hello world!'
}
```

:::{warning}
Statements and top-level declarations can not be mixed at the same level. If your script has top-level declarations, all statements must be contained within top-level declarations such as the entry workflow.
:::

### Shebang

The first line of a script can be a [shebang](https://en.wikipedia.org/wiki/Shebang_(Unix)):

```sh
#!/usr/bin/env nextflow
```

### Feature flag

A feature flag declaration is an assignment. The target should be a valid {ref}`feature flag <config-feature-flags>` and the source should be a literal (i.e. number, string, boolean):

```nextflow
nextflow.preview.topic = true
```

### Include

An include declaration consists of an *include source* and one or more *include clauses*:

```nextflow
include { foo as bar } from './some/module'
```

The include source should be a string literal and should refer to either a local path (e.g. `./module.nf`) or a plugin (e.g. `plugin/nf-hello`). Each include clause should specify a name, and may also specify an *alias*. In the above example, `foo` is included under the alias `bar`.

Include clauses can be separated by semi-colons or newlines:

```nextflow
// semi-colons
include { foo ; bar as baz } from './some/module'

// newlines
include {
    foo
    bar as baz
} from './some/module'
```

Include clauses can also be specified as separate includes:

```nextflow
include { foo } from './some/module'
include { bar as baz } from './some/module'
```

The following definitions can be included:

- Functions
- Processes
- Named workflows

### Parameter

A parameter declaration is an assignment. The target should be a pipeline parameter and the source should be an expression:

```nextflow
params.message = 'Hello world!'
```

Parameters supplied via command line options, params files, and config files take precedence over parameter definitions in a script.

(syntax-workflow)=

### Workflow

A workflow can be a *named workflow* or an *entry workflow*.

A *named workflow* consists of a name and a body, and may consist of a *take*, *main*, *emit*, and *publish* section:

```nextflow
workflow greet {
    take:
    greetings

    main:
    messages = greetings.map { v -> "$v world!" }

    emit:
    messages
}
```

- The take, emit, and publish sections are optional. The `main:` section label can be omitted if they are not specified.

- The take section consists of one or more parameters.

- The main section consists of one or more [statements](#statements).

- The emit section consists of one or more *emit statements*. An emit statement can be a [variable name](#variable), an [assignment](#assignment), or an [expression statement](#expression-statement). If an emit statement is an expression statement, it must be the only emit.

- The publish section can be specified but is intended to be used in the entry workflow (see below).


An *entry workflow* has no name and may consist of a *main* and *publish* section:

```nextflow
workflow {
    main:
    greetings = Channel.of('Bonjour', 'Ciao', 'Hello', 'Hola')
    messages = greetings.map { v -> "$v world!" }
    greetings.view { it -> '$it world!' }

    publish:
    messages >> 'messages'
}
```

- Only one entry workflow may be defined in a script.

- The `main:` section label can be omitted if the publish section is not specified.

- The publish section consists of one or more *publish statements*. A publish statement is a [right-shift expression](#binary-expressions), where the left-hand side is an expression that refers to a value in the workflow body, and the right-hand side is an expression that returns a string.

In order for a script to be executable, it must either define an entry workflow or use the implicit workflow syntax described [above](#top-level-declarations).

Entry workflow definitions are ignored when a script is included as a module. This way, the same script can be included as a module or executed as a pipeline.

(syntax-process)=

### Process

A process consists of a name and a body. The process body consists of one or more [statements](#statements). A minimal process definition must return a string:

```nextflow
process sayHello {
    """
    echo 'Hello world!'
    """
}
```

A process may define additional sections for *directives*, *inputs*, *outputs*, *script*, *shell*, *exec*, and *stub*:

```nextflow
process greet {
    // directives
    errorStrategy 'retry'
    tag { "${greeting}/${name}" }

    input: 
    val greeting
    val name

    output:
    stdout

    script: // or shell: or exec:
    """
    echo '${greeting}, ${name}!'
    """

    stub:
    """
    # do nothing
    """
}
```

- A process must define a script, shell, or exec section (see below). All other sections are optional. Directives do not have an explicit section label, but must be defined first.

- The `script:` section label can be omitted only when there are no other sections in the body.

- Sections must be defined in the order shown above, with the exception of the output section, which can also be specified after the script and stub.

Each section may contain one or more statements. For directives, inputs, and outputs, these statements must be [function calls](#function-call). See {ref}`process-reference` for the set of available input qualifiers, output qualifiers, and directives.

The script section can be substituted with a shell or exec section:

```nextflow
process greetShell {
    input: 
    val greeting

    shell:
    '''
    echo '!{greeting}, ${USER}!'
    '''
}

process greetExec {
    input: 
    val greeting
    val name

    exec:
    message = "${greeting}, ${name}!"

    output:
    val message
}
```

The script, shell, and stub sections must return a string in the same manner as a [function](#function).

See {ref}`process-page` for more information on the semantics of each process section.

(syntax-function)=

### Function

A function consists of a name, parameter list, and a body:

```nextflow
def greet(greeting, name) {
    println "${greeting}, ${name}!"
}
```

The function body consists of one or more [statements](#statements). The last statement is implicitly treated as a return statement if it is an [expression statement](#expression-statement) that returns a value.

The [return statement](#return) can be used to explicitly return from a function:

```nextflow
// return with no value
def greet(greeting, name) {
    if( !greeting || !name )
        return
    println "${greeting}, ${name}!"
}

// return a value
def fib(x) {
    if( x <= 1 )
        return x
    fib(x - 1) + fib(x - 2)
}
```

### Enum type

An enum type declaration consists of a name and a body. The body consists of a comma-separated list of identifiers:

```nextflow
enum Day {
    MONDAY,
    TUESDAY,
    WEDNESDAY,
    THURSDAY,
    FRIDAY,
    SATURDAY,
    SUNDAY
}
```

Enum values in the above example can be accessed as `Day.MONDAY`, `Day.TUESDAY`, and so on.

:::{note}
Enum types cannot be included across modules at this time.
:::

### Output block

The output block consists of one or more *target blocks*. A target block consists of a *target name* and one or more *target directives* for configuring the corresponding publish target:

```nextflow
output {
    'fastq' {
        path 'samples'
        index {
            path 'index.csv'
        }
    }
}
```

Only one output block may be defined in a script. See {ref}`workflow-output-def` for the set of available target directives.

## Statements

Statements can be separated by either a newline or a semi-colon:

```nextflow
// newline
println 'Hello!'
println 'Hello again!'

// semi-colon
println 'Hello!' ; println 'Hello again!'
```

### Variable declaration

Variables can be declared with the `def` keyword:

```nextflow
def x = 42
```

Multiple variables can be declared in a single statement if the initializer is a [list literal](#list) with the same number of elements and declared variables:

```nextflow
def (x, y) = [ 1, 2 ]
```

Each variable has a *scope*, which is the region of code in which the variable can be used.

Variables declared in a function, as well as the parameters of that function, exist for the duration of that function call. The same applies to closures.

Workflow inputs exist for the entire workflow body. Variables declared in the main section exist for the main, emit, and publish sections. Named outputs are not considered variable declarations and therefore do not have any scope.

Process input variables exist for the entire process body. Variables declared in the process script, shell, exec, and stub sections exist only in their respective section, with one exception -- variables declared without the `def` keyword also exist in the output section.

Variables declared in an if or else branch exist only within that branch:

```nextflow
if( true )
    def x = 'foo'
println x           // error: `x` is undefined

// solution: declare `x` outside of if branch
def x
if( true )
    x = 'foo'
println x
```

A variable cannot be declared with the same name as another variable in the same scope or an enclosing scope:

```nextflow
def clash(x) {
    def x           // error: `x` is already declared
    if( true )
        def x       // error: `x` is already declared
}
```

### Assignment

An assignment statement consists of a *target* expression and a *source* expression separated by an equals sign:

```nextflow
v = 42
list[0] = 'first'
map.key = 'value'
```

The target expression must be a [variable](#variable), [index](#binary-expressions), or [property](#binary-expressions) expression. The source expression can be any expression.

Multiple variables can be assigned in a single statement as long as the source expression is a [list literal](#list) with the same number of elements and assigned variables:

```nextflow
(x, y) = [ 1, 2 ]
```

### Expression statement

Any [expression](#expressions) can be a statement.

In general, the only expressions that can have any effect as expression statements are function calls that have side effects (e.g. `println`) or an implicit return statement in a [function](#function) or [closure](#closure).

### assert

An assert statement consists of the `assert` keyword followed by a boolean expression, with an optional error message separated by a colon:

```nextflow
assert 2 + 2 == 4 : 'The math broke!'
```

If the condition is false, an error will be raised with the given error message.

### if/else

An if/else statement consists of an *if branch* and an optional *else branch*. Each branch consists of a boolean expression in parentheses, followed by either a single statement or a *block statement* (one or more statements in curly braces). For example:

```nextflow
def x = Math.random()
if( x < 0.5 ) {
    println 'You lost.'
}
else {
    println 'You won!'
}
```

If the condition is true, the if branch will be executed, otherwise the else branch will be executed.

If/else statements can be chained any number of times by making the else branch another if/else statement:

```nextflow
def grade = 89
if( grade >= 90 )
    println 'You get an A!'
else if( grade >= 80 )
    println 'You get a B!'
else if( grade >= 70 )
    println 'You get a C!'
else if( grade >= 60 )
    println 'You get a D!'
else
    println 'You failed.'
```

A more verbose way to write the same code is:

```nextflow
def grade = 89
if( grade >= 90 ) {
    println 'You get an A!'
}
else {
    if( grade >= 80 ) {
        println 'You get a B!'
    }
    else {
        if( grade >= 70 ) {
            println 'You get a C!'
        }
        else {
            if( grade >= 60 ) {
                println 'You get a D!'
            }
            else {
                println 'You failed.'
            }
        }
    }
}
```

### return

A return statement consists of the `return` keyword with an optional expression:

```nextflow
def add(a, b) {
    return a + b
}

def sayHello(name) {
    if( !name )
        return
    println "Hello, ${name}!"
}
```

Return statements can only be used in functions and closures. In the case of a nested closure, the return statement will return from the nearest enclosing closure.

If a function or closure has multiple return statements (including implicit returns), all of the return statements should either return a value or return nothing. If a function or closure does return a value, it should do so for every conditional branch.

```nextflow
def isEven1(n) {
    if( n % 2 == 1 )
        return          // error: return value is required here
    return true
}

def isEven2(n) {
    if( n % 2 == 0 )
        return true
                        // error: return value is required here
}
```

:::{note}
If the last statement is not a return or expression statement (implicit return), it is equivalent to appending an empty return.
:::

### throw

A throw statement consists of the `throw` keyword followed by an expression that returns an error type:

```nextflow
throw new Exception('something failed!')
```

:::{note}
In general, the appropriate way to raise an error is to use the {ref}`error <stdlib-functions>` function:
```nextflow
error 'something failed!'
```
:::

### try/catch

A try/catch statement consists of a *try block* followed by any number of *catch clauses*:

```nextflow
def text = null
try {
    text = file('foo.txt').text
}
catch( IOException e ) {
    log.warn "Could not load foo.txt"
}
```

The try block will be executed, and if an error is raised and matches the expected error type of a catch clause, the code in that catch clause will be executed. If no catch clause is matched, the error will be raised to the next enclosing try/catch statement, or to the Nextflow runtime.

## Expressions

An expression represents a value. A *literal* value is an expression whose value is known at compile-time, such as a number, string, or boolean. All other expressions must be evaluated at run-time.

Every expression has a *type*, which may be resolved at compile-time or run-time.

### Variable

A variable expression is a reference to a variable or other named value:

```nextflow
def x = 42

x
// -> 42
```

### Number

A number literal can be an integer or floating-point number, and can be positive or negative. Integers can specified in binary with `0b`, octal with `0`, or hexadecimal with `0x`. Floating-point numbers can use scientific notation with the `e` or `E` prefix. Underscores can be used as thousands separators to make long numbers more readable.

```nextflow
// integer
42
-1
0b1001  // -> 9
031     // -> 25
0xabcd  // -> 43981

// real
3.14
-0.1
1.59e7  // -> 15_900_000
1.59e-7 // -> 0.000000159
```

### Boolean

A boolean literal can be `true` or `false`:

```nextflow
assert true != false
assert !true == false
assert true == !false
```

### Null

The null literal is specified as `null`. It can be used to represent an "empty" value:

```nextflow
def x = null
x = 42
```

:::{note}
Using a null value in certain expressions (e.g. the object of a property expression or method call) will lead to a "null reference" error. It is best to avoid the use of `null` where possible.
:::

### String

A string literal consists of arbitrary text enclosed by single or double quotes:

```nextflow
println "I said 'hello'"
println 'I said "hello" again!'
```

A triple-quoted string can span multiple lines:

```nextflow
println '''
    Hello,
    How are you today?
    '''

println """
    We don't have to escape quotes anymore!
    Even "double" quotes!
    """
```

A *slashy string* is enclosed by slashes instead of quotes:

```nextflow
/no escape!/
```

Slashy strings can also span multiple lines:

```nextflow
/
Patterns in the code,
Symbols dance to match and find,
Logic unconfined.
/
```

:::{note}
A slashy string cannot be empty because it would become a line comment.
:::

### Dynamic string

Double-quoted strings can be interpolated using the `${}` placeholder with an expression:

```nextflow
def names = ['Thing 1', 'Thing 2']
println "Hello, ${names.join(' and ')}!"
// -> Hello, Thing 1 and Thing 2!
```

If the expression is a name or simple property expression (one or more identifiers separated by dots), the curly braces can be omitted:

```nextflow
def name = [first: '<FIRST_NAME>', last: '<LAST_NAME>']
println "Hello, ${name.first} ${name.last}!"
// -> Hello, <FIRST_NAME> <LAST_NAME>!
```

Multi-line double-quoted strings can also be interpolated:

```nextflow
"""
blastp \
    -in $input \
    -out $output \
    -db $blast_db \
    -html
"""
```

Single-quoted strings are not interpolated:

```nextflow
println 'Hello, ${names.join(" and ")}!'
// -> Hello, ${names.join(" and ")}!
```

### List

A list literal consists of a comma-separated list of zero or more expressions, enclosed in square brackets:

```nextflow
[1, 2, 3]
```

### Map

A map literal consists of a comma-separated list of one or more *map entries*, enclosed in square brackets. Each map entry consists of a *key expression* and *value expression* separated by a colon:

```nextflow
[foo: 1, bar: 2, baz: 3]
```

An empty map is specified with a single colon to distinguish it from an empty list:

```nextflow
[:]
```

Both the key and value can be any expression. Identifier keys are treated as string literals (i.e. the quotes can be omitted). A variable can be used as a key by enclosing it in parentheses:

```nextflow
def x = 'foo'
[(x): 1]
// -> ['foo': 1]
```

### Closure

A closure, also known as an anonymous function, consists of a parameter list followed by zero or more statements, enclosed in curly braces:

```nextflow
{ a, b -> a + b }
```

The above closure takes two arguments and returns their sum.

The closure body is identical to that of a [function](#function). Statements should be separated by newlines or semi-colons, and the last statement is implicitly treated as a [return statement](#return):

```nextflow
{ v ->
    println 'Hello!'
    println "We're in a closure!"
    println 'Goodbye...'
    v * v
}
```

Closures can access variables outside of their scope:

```nextflow
def factor = 2
println [1, 2, 3].collect { v -> factor * v }
// -> [2, 4, 6]
```

Closures can declare local variables that exist only for the lifetime of each closure invocation:

```nextflow
def result = 0
[1, 2, 3].each { v ->
    def squared = v * v
    result += squared
}

println result
// -> 14
```

See {ref}`standard library <stdlib-page>` and {ref}`operator <operator-page>` for more examples of how closures are used in practice.

### Index expression

An index expression consists of a *left expression* and a *right expression*, with the right expression enclosed in square brackets:

```nextflow
myList[0]
```

### Property expression

A property expression consists of an *object expression* and a *property*, separated by a dot:

```nextflow
file.text
```

The property must be an identifier or string literal.

### Function call

A function call consists of a name and argument list:

```nextflow
printf('Hello %s!\n', 'World')
```

A *method call* consists of an *object expression* and a function call separated by a dot:

```nextflow
myList.size()
```

The argument list may contain any number of *positional arguments* and *named arguments*:

```nextflow
file('hello.txt', checkIfExists: true)
```

The named arguments are collected into a map and provided as the first positional argument to the function. The above function call can be rewritten as:

```nextflow
file([checkIfExists: true], 'hello.txt')
```

The argument name must be an identifier or string literal.

The parentheses can be omitted when the function call is also an [expression statement](#expression-statement) and there is at least one argument:

```nextflow
// positional args
printf 'Hello %s!\n', 'World'

// positional and named args
file 'hello.txt', checkIfExists: true
```

If the last argument is a closure, it can be specified outside of the parentheses:

```nextflow
// closure arg with additional args
[1, 2, 3].inject('result:') { acc, v -> acc + ' ' + v }

// single closure arg
[1, 2, 3].each() { v -> println v }

// single closure arg without parentheses
[1, 2, 3].each { v -> println v }
```

### Constructor call

A constructor call consists of the `new` keyword followed by a *type name* and an argument list enclosed in parentheses:

```nextflow
new java.util.Date()
```

If the type is implicitly available in the script, the *fully-qualified type name* can be elided to the *simple type name*:

```nextflow
new Date()
```

See {ref}`stdlib-default-imports` for the set of types which are implicitly available in Nextflow scripts.

### Unary expressions

A unary expression consists of a *unary operator* followed by an expression:

```nextflow
!(2 + 2 == 4)
```

The following unary operators are available:

- `~`: bitwise NOT
- `!`: logical NOT
- `+`: unary plus
- `-`: unary minus

### Binary expressions

A binary expression consists of a *left expression* and a *right expression* separated by a *binary operator*:

```nextflow
2 + 2
```

The following binary operators are available:

- `**`: power (i.e. exponentiation)
- `*`: multiplication
- `/`: division
- `%`: remainder (i.e. modulo)
- `+`: addition
- `-`: subtraction
- `<<`: left shift
- `>>`: right shift
- `>>>`: unsigned right shift
- `..`: inclusive range
- `..<`: right-exclusive range
- `as`: type cast
- `instanceof`: type relation
- `!instanceof`: negated type relation
- `<`: less than
- `>`: greater than
- `<=`: less than or equals
- `>=`: greater than or equals
- `in`: membership
- `!in`: negated membership
- `==`: equals
- `!=`: negated equals
- `<=>`: spaceship (i.e. three-way comparison)
- `=~`: regex find
- `==~`: regex match
- `&`: bitwise AND
- `^`: bitwise XOR (exclusive or)
- `|`: bitwise OR
- `&&`: logical AND
- `||`: logical OR
- `?:` elvis (i.e. short ternary)

### Ternary expression

A ternary expression consists of a *test expression*, a *true expression*, and a *false expression*, separated by a question mark and a colon:

```nextflow
println x % 2 == 0 ? 'x is even!' : 'x is odd!'
```

If the test expression is true, the true expression is evaluated, otherwise the false expression is evaluated.

### Parentheses

Any expression can be enclosed in parentheses:

```nextflow
1 + 2 * 3
// -> 1 + 6 -> 7

(1 + 2) * 3
// -> 3 * 3 -> 9
```

### Precedence

Compound expressions are evaluated in the following order:

- parentheses
- property expressions
- function calls
- index expressions
- `~`,  `!`
- `**`
- `+`, `-` (unary)
- `*`, `/`, `%`
- `+`, `-` (binary)
- `<<`, `>>>`, `>>`, `..`, `..<`
- `as`
- `instanceof`, `!instanceof`
- `<`, `>`, `<=`, `>=`, `in`, `!in`
- `==`, `!=`, `<=>`
- `=~`, `==~`
- `&`
- `^`
- `|`
- `&&`
- `||`
- `?:` (ternary)
- `?:` (elvis)

## Deprecations

The following legacy features were excluded from this page because they are deprecated:

- The `addParams` and `params` clauses of include declarations. See {ref}`module-params` for more information.
- The `when:` section of a process definition. See {ref}`process-when` for more information.
- The implicit `it` closure parameter. See {ref}`script-closure` for more information.

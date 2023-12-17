(process-page)=

# Processes

In Nextflow, a **process** is the basic processing primitive to execute a user script.

The process definition starts with the keyword `process`, followed by process name and finally the process body delimited by curly braces. The process body must contain a string which represents the command or, more generally, a script that is executed by it. A basic process looks like the following example:

```groovy
process sayHello {
    """
    echo 'Hello world!' > file
    """
}
```

A process may contain any of the following definition blocks: directives, inputs, outputs, when clause, and the process script. The syntax is defined as follows:

```
process < name > {

  [ directives ]

  input:
    < process inputs >

  output:
    < process outputs >

  when:
    < condition >

  [script|shell|exec]:
    < user script to be executed >

}
```

(process-script)=

## Script

The `script` block defines, as a string expression, the script that is executed by the process.

A process may contain only one script, and if the `script` guard is not explicitly declared, the script must be the final statement in the process block.

The script string is executed as a [Bash](<http://en.wikipedia.org/wiki/Bash_(Unix_shell)>) script in the host environment. It can be any command or script that you would normally execute on the command line or in a Bash script. Naturally, the script may only use commands that are available in the host environment.

The script block can be a simple string or a multi-line string. The latter approach makes it easier to write scripts with multiple commands spanning multiple lines. For example:

```groovy
process doMoreThings {
  """
  blastp -db $db -query query.fa -outfmt 6 > blast_result
  cat blast_result | head -n 10 | cut -f 2 > top_hits
  blastdbcmd -db $db -entry_batch top_hits > sequences
  """
}
```

As explained in the script tutorial section, strings can be defined using single-quotes or double-quotes, and multi-line strings are defined by three single-quote or three double-quote characters.

There is a subtle but important difference between them. Like in Bash, strings delimited by a `"` character support variable substitutions, while strings delimited by `'` do not.

In the above code fragment, the `$db` variable is replaced by the actual value defined elsewhere in the pipeline script.

:::{warning}
Since Nextflow uses the same Bash syntax for variable substitutions in strings, you must manage them carefully depending on whether you want to evaluate a *Nextflow* variable or a *Bash* variable.
:::

When you need to access a system environment variable in your script, you have two options.

If you don't need to access any Nextflow variables, you can define your script block with single-quotes:

```groovy
process printPath {
  '''
  echo The path is: $PATH
  '''
}
```

Otherwise, you can define your script with double-quotes and escape the system environment variables by prefixing them with a back-slash `\` character, as shown in the following example:

```groovy
process doOtherThings {
  """
  blastp -db \$DB -query query.fa -outfmt 6 > blast_result
  cat blast_result | head -n $MAX | cut -f 2 > top_hits
  blastdbcmd -db \$DB -entry_batch top_hits > sequences
  """
}
```

In this example, `$MAX` is a Nextflow variable that must be defined elsewhere in the pipeline script. Nextflow replaces it with the actual value before executing the script. Meanwhile, `$DB` is a Bash variable that must exist in the execution environment, and Bash will replace it with the actual value during execution.

:::{tip}
Alternatively, you can use the {ref}`process-shell` block definition, which allows a script to contain both Bash and Nextflow variables without having to escape the first.
:::

### Scripts *à la carte*

The process script is interpreted by Nextflow as a Bash script by default, but you are not limited to Bash.

You can use your favourite scripting language (Perl, Python, R, etc), or even mix them in the same pipeline.

A pipeline may be composed of processes that execute very different tasks. With Nextflow, you can choose the scripting language that best fits the task performed by a given process. For example, for some processes R might be more useful than Perl, whereas for others you may need to use Python because it provides better access to a library or an API, etc.

To use a language other than Bash, simply start your process script with the corresponding [shebang](<http://en.wikipedia.org/wiki/Shebang_(Unix)>). For example:

```groovy
process perlTask {
    """
    #!/usr/bin/perl

    print 'Hi there!' . '\n';
    """
}

process pythonTask {
    """
    #!/usr/bin/python

    x = 'Hello'
    y = 'world!'
    print "%s - %s" % (x,y)
    """
}

workflow {
    perlTask()
    pythonTask()
}
```

:::{tip}
Since the actual location of the interpreter binary file can differ across platforms, it is wise to use the `env` command followed by the interpreter name, e.g. `#!/usr/bin/env perl`, instead of the absolute path, in order to make your script more portable.
:::

### Conditional scripts

So far, our `script` block has always been a simple string expression, but in reality, the `script` block is just Groovy code that returns a string. This means that you can write arbitrary Groovy code to determine the script to execute, as long as the final statement is a string (remember that the `return` keyword is optional in Groovy).

For example, you can use flow control statements (`if`, `switch`, etc) to execute a different script based on the process inputs. The only difference here is that you must explicitly declare the `script` guard, whereas before it was not required. Here is an example:

```groovy
mode = 'tcoffee'

process align {
    input:
    path sequences

    script:
    if( mode == 'tcoffee' )
        """
        t_coffee -in $sequences > out_file
        """

    else if( mode == 'mafft' )
        """
        mafft --anysymbol --parttree --quiet $sequences > out_file
        """

    else if( mode == 'clustalo' )
        """
        clustalo -i $sequences -o out_file
        """

    else
        error "Invalid alignment mode: ${mode}"
}
```

In the above example, the process will execute one of the script fragments depending on the value of the `mode` parameter. By default it will execute the `tcoffee` command, but changing the `mode` variable will cause a different branch to be executed.

(process-template)=

### Template

Process scripts can be externalised to **template** files, which can be reused across different processes and tested independently from the overall pipeline execution.

A template is simply a shell script file that Nextflow is able to execute by using the `template` function as shown below:

```groovy
process templateExample {
    input:
    val STR

    script:
    template 'my_script.sh'
}

workflow {
    Channel.of('this', 'that') | templateExample
}
```

By default, Nextflow looks for the `my_script.sh` template file in the `templates` directory located alongside the Nextflow script and/or the module script in which the process is defined. Any other location can be specified by using an absolute template path.

The template script may contain any code that can be executed by the underlying environment. For example:

```bash
#!/bin/bash
echo "process started at `date`"
echo $STR
echo "process completed"
```

:::{tip}
The dollar character (`$`) is interpreted as a Nextflow variable when the script is run as a Nextflow template, whereas it is evaluated as a Bash variable when run as a Bash script. This can be very useful for testing your script independently from Nextflow execution. You only need to provide a Bash environment variable for each of the Nextflow variables that are referenced in your script. For example, it would be possible to execute the above script with the following command in the terminal: `STR='foo' bash templates/my_script.sh`
:::

:::{tip}
As a best practice, the template script should not contain any `\$` escaped variables, because these variables will not be evaluated properly when the script is executed directly.
:::

(process-shell)=

### Shell

The `shell` block is a string expression that defines the script that is executed by the process. It is an alternative to the {ref}`process-script` definition with one important difference: it uses the exclamation mark `!` character, instead of the usual dollar `$` character, to denote Nextflow variables.

This way, it is possible to use both Nextflow and Bash variables in the same script without having to escape the latter, which makes process scripts easier to read and maintain. For example:

```groovy
process myTask {
    input:
    val str

    shell:
    '''
    echo "User $USER says !{str}"
    '''
}

workflow {
    Channel.of('Hello', 'Hola', 'Bonjour') | myTask
}
```

In the above example, `$USER` is treated as a Bash variable, while `!{str}` is treated as a Nextflow variable.

:::{note}
- Shell script definitions require the use of single-quote `'` delimited strings. When using double-quote `"` delimited strings, dollar variables are interpreted as Nextflow variables as usual. See {ref}`string-interpolation`.
- Variables prefixed with `!` must always be enclosed in curly brackets, i.e. `!{str}` is a valid variable whereas `!str` is ignored.
- Shell scripts support the use of the {ref}`process-template` mechanism. The same rules are applied to the variables defined in the script template.
:::

(process-native)=

### Native execution

Nextflow processes can also execute native Groovy code as the task itself, using the `exec` block. Whereas the `script` block defines a script to be executed, the `exec` block defines Groovy code to be executed directly.

For example:

```groovy
process simpleSum {
    input:
    val x

    exec:
    println "Hello Mr. $x"
}

workflow {
    Channel.of('a', 'b', 'c') | simpleSum
}
```

will display:

```
Hello Mr. b
Hello Mr. a
Hello Mr. c
```

(process-stub)=

## Stub

:::{versionadded} 20.11.0-edge
:::

You can define a command *stub*, which replaces the actual process command when the `-stub-run` or `-stub` command-line option is enabled:

```groovy
process INDEX {
  input:
    path transcriptome

  output:
    path 'index'

  script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """

  stub:
    """
    mkdir index
    touch index/seq.bin
    touch index/info.json
    touch index/refseq.bin
    """
}
```

The `stub` block can be defined before or after the `script` block. When the pipeline is executed with the `-stub-run` option and a process's `stub` is not defined, the `script` block is executed.

This feature makes it easier to quickly prototype the workflow logic without using the real commands. The developer can use it to provide a dummy script that mimics the execution of the real one in a quicker manner. In other words, it is a way to perform a dry-run.

(process-input)=

## Inputs

The `input` block allows you to define the input channels of a process, similar to function arguments. A process may have at most one input block, and it must contain at least one input.

The input block follows the syntax shown below:

```
input:
  <input qualifier> <input name>
```

An input definition consists of a *qualifier* and a *name*. The input qualifier defines the type of data to be received. This information is used by Nextflow to apply the semantic rules associated with each qualifier, and handle it properly depending on the target execution platform (grid, cloud, etc).

When a process is invoked in a workflow block, it must be provided a channel for each channel in the process input block, similar to calling a function with specific arguments. The examples provided in the following sections demonstrate how a process is invoked with input channels.

The following input qualifiers are available:

- `val`: Access the input value by name in the process script.
- `file`: (DEPRECATED) Handle the input value as a file, staging it properly in the execution context.
- `path`: Handle the input value as a path, staging the file properly in the execution context.
- `env`: Use the input value to set an environment variable in the process script.
- `stdin`: Forward the input value to the process `stdin` special file.
- `tuple`: Handle a group of input values having any of the above qualifiers.
- `each`: Execute the process for each element in the input collection.

### Input type `val`

The `val` qualifier accepts any data type. It can be accessed in the process script by using the specified input name, as shown in the following example:

```groovy
process basicExample {
  input:
  val x

  "echo process job $x"
}

workflow {
  def num = Channel.of(1,2,3)
  basicExample(num)
}
```

In the above example, the process is executed three times: once for each value emitted by the `num` channel. The resulting output is similar to the one shown below:

```
process job 3
process job 1
process job 2
```

:::{note}
While channels do emit items in the order that they are received, *processes* do not necessarily *process* items in the order that they are received. In the above example, the value `3` was processed before the others.
:::

:::{note}
When the process declares exactly one input, the pipe `|` operator can be used to provide inputs to the process, instead of passing it as a parameter. Both methods have identical semantics:

```groovy
process basicExample {
  input:
  val x

  "echo process job $x"
}

workflow {
  Channel.of(1,2,3) | basicExample
}
```
:::

### Input type `file`

:::{note}
The `file` qualifier was the standard way to handle input files prior to Nextflow 19.10.0. In later versions of Nextflow, the `path` qualifier should be preferred over `file`.
:::

The `file` qualifier is identical to `path`, with one important difference. When a `file` input receives a value that is not a file, it automatically converts the value to a string and saves it to a temporary file. This behavior is useful in some cases, but tends to be confusing in general. The `path` qualifier instead interprets string values as the path location of the input file and automatically converts to a file object.

(process-input-path)=

### Input type `path`

The `path` qualifier allows you to provide input files to the process execution context. Nextflow will stage the files into the process execution directory, and they can be accessed in the script by using the specified input name. For example:

```groovy
process blastThemAll {
  input:
  path query_file

  "blastp -query ${query_file} -db nr"
}

workflow {
  def proteins = Channel.fromPath( '/some/path/*.fa' )
  blastThemAll(proteins)
}
```

In the above example, all the files ending with the suffix `.fa` are sent over the channel `proteins`. These files are received by the process, which executes a BLAST query on each of them.

It's worth noting that in the above example, the name of the file in the file-system is not used. You can access the file without even knowing its name, because you can reference it in the process script by the input name.

There may be cases where your task needs to use a file whose name is fixed, i.e. it does not have to change along with the actual provided file. In this case, you can specify a fixed name with the `name` attribute in the input file parameter definition, as shown in the following example:

```groovy
input:
path query_file, name: 'query.fa'
```

or, using a shorter syntax:

```groovy
input:
path 'query.fa'
```

The previous example can be re-written as shown below:

```groovy
process blastThemAll {
  input:
  path 'query.fa'

  "blastp -query query.fa -db nr"
}

workflow {
  def proteins = Channel.fromPath( '/some/path/*.fa' )
  blastThemAll(proteins)
}
```

In this example, each file received by the process is staged with the name `query.fa` in a different execution context (i.e. the folder where a task is executed).

:::{tip}
This feature allows you to execute the process command multiple times without worrying about the file names changing. In other words, Nextflow helps you write pipeline tasks that are self-contained and decoupled from the execution environment. As a best practice, you should avoid referencing files in your process script other than those defined in your input block.
:::

Channel factories like `Channel.fromPath` produce file objects, but a `path` input can also accept a string literal path. The string value should be an absolute path, i.e. it must be prefixed with a `/` character or a supported URI protocol (`file://`, `http://`, `s3://`, etc), and it cannot contain special characters (`\n`, etc).

```groovy
process foo {
  input:
  path x

  """
  your_command --in $x
  """
}

workflow {
  foo('/some/data/file.txt')
}
```

Available options:

`arity`
: :::{versionadded} 23.09.0-edge
  :::
: Specify the number of expected files. Can be a number or a range:

  ```groovy
  input:
      path('one.txt', arity: '1')         // exactly one file is expected
      path('pair_*.txt', arity: '2')      // exactly two files are expected
      path('many_*.txt', arity: '1..*')   // one or more files are expected
  ```

  When a task is created, Nextflow will check whether the received files for each path input match the declared arity, and fail if they do not.

`stageAs`
: Specify how the file should be named in the task work directory:

  ```groovy
  process foo {
    input:
    path x, stageAs: 'data.txt'

    """
    your_command --in data.txt
    """
  }

  workflow {
    foo('/some/data/file.txt')
  }
  ```

  Can be a name or a pattern as described in the [Multiple input files](#multiple-input-files) section.

:::{note}
Process `path` inputs have nearly the same interface as described in {ref}`script-file-io`, with one difference which is relevant when files are staged into a subdirectory. Given the following input:

```groovy
path x, stageAs: 'my-dir/*'
```

In this case, `x.name` returns the file name with the parent directory (e.g. `my-dir/file.txt`), whereas normally it would return the file name (e.g. `file.txt`). You can use `x.fileName.name` to get the file name.
:::

### Multiple input files

A `path` input can also accept a collection of files instead of a single value. In this case, the input variable will be a Groovy list, and you can use it as such.

When the input has a fixed file name and a collection of files is received by the process, the file name will be appended with a numerical suffix representing its ordinal position in the list. For example:

```groovy
process blastThemAll {
    input:
    path 'seq'

    "echo seq*"
}

workflow {
    def fasta = Channel.fromPath( "/some/path/*.fa" ).buffer(size: 3)
    blastThemAll(fasta)
}
```

will output:

```
seq1 seq2 seq3
seq1 seq2 seq3
...
```

The target input file name may contain the `*` and `?` wildcards, which can be used to control the name of staged files. The following table shows how the wildcards are replaced depending on the cardinality of the received input collection.

| Cardinality | Name pattern | Staged file names                                                                                       |
| ----------- | ------------ | ------------------------------------------------------------------------------------------------------- |
| any         | `*`          | named as the source file                                                                                |
| 1           | `file*.ext`  | `file.ext`                                                                                              |
| 1           | `file?.ext`  | `file1.ext`                                                                                             |
| 1           | `file??.ext` | `file01.ext`                                                                                            |
| many        | `file*.ext`  | `file1.ext`, `file2.ext`, `file3.ext`, ..                                                               |
| many        | `file?.ext`  | `file1.ext`, `file2.ext`, `file3.ext`, ..                                                               |
| many        | `file??.ext` | `file01.ext`, `file02.ext`, `file03.ext`, ..                                                            |
| many        | `dir/*`      | named as the source file, created in `dir` subdirectory                                                 |
| many        | `dir??/*`    | named as the source file, created in a progressively indexed subdirectory e.g. `dir01/`, `dir02/`, etc. |
| many        | `dir*/*`     | (as above)                                                                                              |

The following example shows how a wildcard can be used in the input file definition:

```groovy
process blastThemAll {
    input:
    path 'seq?.fa'

    "cat seq1.fa seq2.fa seq3.fa"
}

workflow {
    def fasta = Channel.fromPath( "/some/path/*.fa" ).buffer(size: 3)
    blastThemAll(fasta)
}
```

:::{note}
Rewriting input file names according to a named pattern is an extra feature and not at all required. The normal file input syntax introduced in the {ref}`process-input-path` section is valid for collections of multiple files as well. To handle multiple input files while preserving the original file names, use a variable identifier or the `*` wildcard.
:::

### Dynamic input file names

When the input file name is specified by using the `name` option or a string literal, you can also use other input values as variables in the file name string. For example:

```groovy
process simpleCount {
  input:
  val x
  path "${x}.fa"

  """
  cat ${x}.fa | grep '>'
  """
}
```

In the above example, the input file name is determined by the current value of the `x` input value.

This approach allows input files to be staged in the task directory with a name that is coherent with the current execution context.

:::{tip}
In most cases, you won't need to use dynamic file names, because each task is executed in its own directory, and input files are automatically staged into this directory by Nextflow. This behavior guarantees that input files with the same name won't overwrite each other.

An example of when you may have to deal with that is when you have many input files in a task, and some of these files may have the same filename. In this case, a solution would be to use the `stageAs` option.
:::

### Input type `env`

The `env` qualifier allows you to define an environment variable in the process execution context based on the input value. For example:

```groovy
process printEnv {
    input:
    env HELLO

    '''
    echo $HELLO world!
    '''
}

workflow {
    Channel.of('hello', 'hola', 'bonjour', 'ciao') | printEnv
}
```

```
hello world!
ciao world!
bonjour world!
hola world!
```

### Input type `stdin`

The `stdin` qualifier allows you to forward the input value to the [standard input](http://en.wikipedia.org/wiki/Standard_streams#Standard_input_.28stdin.29) of the process script. For example:

```groovy
process printAll {
  input:
  stdin str

  """
  cat -
  """
}

workflow {
  Channel.of('hello', 'hola', 'bonjour', 'ciao')
    | map { it + '\n' }
    | printAll
}
```

will output:

```
hola
bonjour
ciao
hello
```

(process-input-set)=

### Input type `set`

:::{deprecated} 19.08.1-edge
Use `tuple` instead.
:::

(process-input-tuple)=

### Input type `tuple`

The `tuple` qualifier allows you to group multiple values into a single input definition. It can be useful when a channel emits tuples of values that need to be handled separately. Each element in the tuple is associated with a corresponding element in the `tuple` definition. For example:

```groovy
process tupleExample {
    input:
    tuple val(x), path('latin.txt')

    """
    echo "Processing $x"
    cat - latin.txt > copy
    """
}

workflow {
  Channel.of( [1, 'alpha'], [2, 'beta'], [3, 'delta'] ) | tupleExample
}
```

In the above example, the `tuple` input consists of the value `x` and the file `latin.txt`.

A `tuple` definition may contain any of the following qualifiers, as previously described: `val`, `env`, `path` and `stdin`. Files specified with the `path` qualifier are treated exactly the same as standalone `path` inputs.

### Input repeaters (`each`)

The `each` qualifier allows you to repeat the execution of a process for each item in a collection, each time a new value is received. For example:

```groovy
process alignSequences {
  input:
  path seq
  each mode

  """
  t_coffee -in $seq -mode $mode > result
  """
}

workflow {
  sequences = Channel.fromPath('*.fa')
  methods = ['regular', 'espresso', 'psicoffee']

  alignSequences(sequences, methods)
}
```

In the above example, each time a file of sequences is emitted from the `sequences` channel, the process executes *three* tasks, each running a T-coffee alignment with a different value for the `mode` parameter. This behavior is useful when you need to repeat the same task over a given set of parameters.

Input repeaters can be applied to files as well. For example:

```groovy
process alignSequences {
  input:
  path seq
  each mode
  each path(lib)

  """
  t_coffee -in $seq -mode $mode -lib $lib > result
  """
}

workflow {
  sequences = Channel.fromPath('*.fa')
  methods = ['regular', 'espresso']
  libraries = [ file('PQ001.lib'), file('PQ002.lib'), file('PQ003.lib') ]

  alignSequences(sequences, methods, libraries)
}
```

In the above example, each sequence input file emitted by the `sequences` channel triggers six alignment tasks, three with the `regular` method against each library file, and three with the `espresso` method.

:::{note}
When multiple repeaters are defined, the process is executed for each *combination* of them.
:::

:::{note}
Input repeaters currently do not support tuples. However, you can emulate an input repeater on a channel of tuples by using the {ref}`operator-combine` or {ref}`operator-cross` operator with other input channels to produce all of the desired input combinations.
:::

(process-multiple-input-channels)=

### Multiple input channels

A key feature of processes is the ability to handle inputs from multiple channels.

When two or more channels are declared as process inputs, the process waits until there is a complete input configuration, i.e. until it receives a value from each input channel. When this condition is satisfied, the process consumes a value from each channel and launches a new task, repeating this logic until one or more channels are empty.

As a result, channel values are consumed sequentially and any empty channel will cause the process to wait, even if the other channels have values.

For example:

```groovy
process foo {
  input:
  val x
  val y

  script:
  """
  echo $x and $y
  """
}

workflow {
  x = Channel.of(1, 2)
  y = Channel.of('a', 'b', 'c')
  foo(x, y)
}
```

The process `foo` is executed two times because the `x` channel emits only two values, therefore the `c` element is discarded. It outputs:

```
1 and a
2 and b
```

A different semantic is applied when using a {ref}`value channel <channel-type-value>`. This kind of channel is created by the {ref}`Channel.value <channel-value>` factory method or implicitly when a process is invoked with an argument that is not a channel. By definition, a value channel is bound to a single value and it can be read an unlimited number of times without consuming its content. Therefore, when mixing a value channel with one or more (queue) channels, it does not affect the process termination because the underlying value is applied repeatedly.

To better understand this behavior, compare the previous example with the following one:

```groovy
process bar {
  input:
  val x
  val y

  script:
  """
  echo $x and $y
  """
}

workflow {
  x = Channel.value(1)
  y = Channel.of('a', 'b', 'c')
  foo(x, y)
}
```

The above example executes the `bar` process three times because `x` is a value channel, therefore its value can be read as many times as needed. The process termination is determined by the contents of `y`. It outputs:

```
1 and a
1 and b
1 and c
```

:::{note}
In general, multiple input channels should be used to process *combinations* of different inputs, using the `each` qualifier or value channels. Having multiple queue channels as inputs is equivalent to using the {ref}`operator-merge` operator, which is not recommended as it may lead to {ref}`non-deterministic process inputs <cache-nondeterministic-inputs>`.
:::

See also: {ref}`channel-types`.

(process-output)=

## Outputs

The `output` block allows you to define the output channels of a process, similar to function outputs. A process may have at most one output block, and it must contain at least one output.

The output block follows the syntax shown below:

```
output:
  <output qualifier> <output name> [, <option>: <option value>]
```

An output definition consists of a *qualifier* and a *name*. Some optional attributes can also be specified.

When a process is invoked, each process output is returned as a channel. The examples provided in the following sections demonstrate how to access the output channels of a process.

The following output qualifiers are available:

- `val`: Emit the variable with the specified name.
- `file`: (DEPRECATED) Emit a file produced by the process with the specified name.
- `path`: Emit a file produced by the process with the specified name.
- `env`: Emit the variable defined in the process environment with the specified name.
- `stdout`: Emit the `stdout` of the executed process.
- `tuple`: Emit multiple values.

### Output type `val`

The `val` qualifier allows you to output any Nextflow variable defined in the process. A common use case is to output a variable that was defined in the `input` block, as shown in the following example:

```groovy
process foo {
  input:
  each x

  output:
  val x

  """
  echo $x > file
  """
}

workflow {
  methods = ['prot', 'dna', 'rna']

  receiver = foo(methods)
  receiver.view { "Received: $it" }
}
```

The output value can be a value literal, an input variable, any other Nextflow variable in the process scope, or a value expression. For example:

```groovy
process foo {
  input:
  path infile

  output:
  val x
  val 'BB11'
  val "${infile.baseName}.out"

  script:
  x = infile.name
  """
  cat $x > file
  """
}

workflow {
  ch_dummy = Channel.fromPath('*').first()
  (ch_var, ch_str, ch_exp) = foo(ch_dummy)

  ch_var.view { "ch_var: $it" }
  ch_str.view { "ch_str: $it" }
  ch_exp.view { "ch_exp: $it" }
}
```

### Output type `file`

:::{note}
The `file` qualifier was the standard way to handle input files prior to Nextflow 19.10.0. In later versions of Nextflow, the `path` qualifier should be preferred over `file`.
:::

The `file` qualifier is similar to `path`, but with some differences. The `file` qualifier interprets `:` as a path separator, therefore `file 'foo:bar'` captures two files named `foo` and `bar`, whereas `path 'foo:bar'` captures a single file named `foo:bar`. Additionally, `file` does not support all of the extra options provided by `path`.

### Output type `path`

The `path` qualifier allows you to output one or more files produced by the process. For example:

```groovy
process randomNum {
  output:
  path 'result.txt'

  '''
  echo $RANDOM > result.txt
  '''
}

workflow {
  numbers = randomNum()
  numbers.view { "Received: ${it.text}" }
}
```

In the above example, the `randomNum` process creates a file named `result.txt` which contains a random number. Since a `path` output with the same name is declared, that file is emitted by the corresponding output channel. A downstream process with a compatible input channel will be able to receive it.

Available options:

`arity`
: :::{versionadded} 23.09.0-edge
  :::
: Specify the number of expected files. Can be a number or a range:

  ```groovy
  output:
      path('one.txt', arity: '1')         // exactly one file is expected
      path('pair_*.txt', arity: '2')      // exactly two files are expected
      path('many_*.txt', arity: '1..*')   // one or more files are expected
  ```

  When a task completes, Nextflow will check whether the produced files for each path output match the declared arity,
  and fail if they do not. If the arity is `1`, a sole file object will be emitted. Otherwise, a list will always be emitted,
  even if only one file is produced.

`followLinks`
: When `true` target files are return in place of any matching symlink (default: `true`)

`glob`
: When `true` the specified name is interpreted as a glob pattern (default: `true`)

`hidden`
: When `true` hidden files are included in the matching output files (default: `false`)

`includeInputs`
: When `true` any input files matching an output file glob pattern are included.

`maxDepth`
: Maximum number of directory levels to visit (default: no limit)

`type`
: Type of paths returned, either `file`, `dir` or `any` (default: `any`, or `file` if the specified file name pattern contains a double star (`**`))

### Multiple output files

When an output file name contains a `*` or `?` wildcard character, it is interpreted as a [glob][glob] path matcher. This allows you to capture multiple files into a list and emit the list as a single value. For example:

```groovy
process splitLetters {
    output:
    path 'chunk_*'

    '''
    printf 'Hola' | split -b 1 - chunk_
    '''
}

workflow {
    splitLetters
        | flatten
        | view { "File: ${it.name} => ${it.text}" }
}
```

It prints:

```
File: chunk_aa => H
File: chunk_ab => o
File: chunk_ac => l
File: chunk_ad => a
```

By default, all the files matching the specified glob pattern are emitted as a single list. However, as the above example demonstrates, the {ref}`operator-flatten` operator can be used to transform the list of files into a channel that emits each file individually.

Some caveats on glob pattern behavior:

- Input files are not included (unless `includeInputs` is `true`)
- Directories are included, unless the `**` pattern is used to recurse through directories

:::{warning}
Although the input files matching a glob output declaration are not included in the resulting output channel, these files may still be transferred from the task scratch directory to the original task work directory. Therefore, to avoid unnecessary file copies, avoid using loose wildcards when defining output files, e.g. `path '*'`. Instead, use a prefix or a suffix to restrict the set of matching files to only the expected ones, e.g. `path 'prefix_*.sorted.bam'`.
:::

Read more about glob syntax at the following link [What is a glob?][what is a glob?]

### Dynamic output file names

When an output file name needs to be expressed dynamically, it is possible to define it using a dynamic string which references variables in the `input` block or in the script global context. For example:

```groovy
process align {
  input:
  val species
  path seq

  output:
  path "${species}.aln"

  """
  t_coffee -in $seq > ${species}.aln
  """
}
```

In the above example, each process execution produces an alignment file whose name depends on the actual value of the `species` input.

:::{tip}
The management of output files in Nextflow is often misunderstood.

With other tools it is generally necessary to organize the output files into some kind of directory structure or to guarantee a unique file name scheme, so that result files don't overwrite each other and so they can be referenced unequivocally by downstream tasks.

With Nextflow, in most cases, you don't need to manage the naming of output files, because each task is executed in its own unique directory, so files produced by different tasks can't overwrite each other. Also, metadata can be associated with outputs by using the {ref}`tuple output <process-out-tuple>` qualifier, instead of including them in the output file name.

One example in which you'd need to manage the naming of output files is when you use the `publishDir` directive to have output files also in a specific path of your choice. If two tasks have the same filename for their output and you want them to be in the same path specified by `publishDir`, the last task to finish will overwrite the output of the task that finished before. You can dynamically change that by adding the `saveAs` option to your `publishDir` directive.

To sum up, the use of output files with static names over dynamic ones is preferable whenever possible, because it will result in simpler and more portable code.
:::

(process-env)=

### Output type `env`

The `env` qualifier allows you to output a variable defined in the process execution environment:

```groovy
process myTask {
    output:
    env FOO

    script:
    '''
    FOO=$(ls -la)
    '''
}

workflow {
    myTask | view { "directory contents: $it" }
}
```

(process-stdout)=

### Output type `stdout`

The `stdout` qualifier allows you to output the `stdout` of the executed process:

```groovy
process sayHello {
    output:
    stdout

    """
    echo Hello world!
    """
}

workflow {
    sayHello | view { "I say... $it" }
}
```

(process-set)=

### Output type `set`

:::{deprecated} 19.08.1-edge
Use `tuple` instead.
:::

(process-out-tuple)=

### Output type `tuple`

The `tuple` qualifier allows you to output multiple values in a single channel. It is useful when you need to associate outputs with metadata, for example:

```groovy
process blast {
  input:
    val species
    path query

  output:
    tuple val(species), path('result')

  script:
    """
    blast -db nr -query $query > result
    """
}

workflow {
  ch_species = Channel.from('human', 'cow', 'horse')
  ch_query = Channel.fromPath('*.fa')

  blast(ch_species, ch_query)
}
```

In the above example, a `blast` task is executed for each pair of `species` and `query` that are received. Each task produces a new tuple containing the value for `species` and the file `result`.

A `tuple` definition may contain any of the following qualifiers, as previously described: `val`, `path`, `env` and `stdout`. Files specified with the `path` qualifier are treated exactly the same as standalone `path` inputs.

:::{note}
While parentheses for input and output qualifiers are generally optional, they are required when specifying elements in an input/output tuple.

Here's an example with a single path output (parentheses optional):

```groovy
process foo {
    output:
    path 'result.txt', hidden: true

    '''
    echo 'another new line' >> result.txt
    '''
}
```

And here's an example with a tuple output (parentheses required):

```groovy
process foo {
    output:
    tuple path('last_result.txt'), path('result.txt', hidden: true)

    '''
    echo 'another new line' >> result.txt
    echo 'another new line' > last_result.txt
    '''
}
```
:::

(process-additional-options)=

### Additional options

The following options are available for all process outputs:

`emit: <name>`

: Defines the name of the output channel, which can be used to access the channel by name from the process output:

  ```groovy
  process FOO {
      output:
      path 'hello.txt', emit: hello
      path 'bye.txt', emit: bye

      """
      echo "hello" > hello.txt
      echo "bye" > bye.txt
      """
  }

  workflow {
      FOO()
      FOO.out.hello.view()
  }
  ```

  See {ref}`workflow-process-invocation` for more details.

`optional: true | false`

: Normally, if a specified output is not produced by the task, the task will fail. Setting `optional: true` will cause the task to not fail, and instead emit nothing to the given output channel.

  ```groovy
  output:
  path("output.txt"), optional: true
  ```

  In this example, the process is normally expected to produce an `output.txt` file, but in the cases where the file is missing, the task will not fail. The output channel will only contain values for those tasks that produced `output.txt`.

: :::{note}
  While this option can be used with any process output, it cannot be applied to individual elements of a [tuple](#output-type-tuple) output. The entire tuple must be optional or not optional.
  :::

`topic: <name>`

: :::{versionadded} 23.11.0-edge
  :::

: *Experimental: may change in a future release.*

: Defines the {ref}`channel topic <channel-topic>` to which the output will be sent.

## When

The `when` block allows you to define a condition that must be satisfied in order to execute the process. The condition can be any expression that returns a boolean value.

It can be useful to enable/disable the process execution depending on the state of various inputs and parameters. For example:

```groovy
process find {
  input:
  path proteins
  val dbtype

  when:
  proteins.name =~ /^BB11.*/ && dbtype == 'nr'

  script:
  """
  blastp -query $proteins -db nr
  """
}
```

:::{tip}
As a best practice, it is better to define such control flow logic in the workflow block, i.e. with an `if` statement or with channel operators, to make the process more portable.
:::

(process-directives)=

## Directives

Directives are optional settings that affect the execution of the current process.

They must be entered at the top of the process body, before any other declaration blocks (`input`, `output`, etc), and have the following syntax:

```groovy
// directive with simple value
name value

// directive with list value
name arg1, arg2, arg3

// directive with map value
name key1: val1, key2: val2

// directive with value and options
name arg, opt1: val1, opt2: val2
```

By default, directives are evaluated when the process is defined. However, if the value is a dynamic string or closure, it will be evaluated separately for each task, which allows task-specific variables like `task` and `val` inputs to be used.

Some directives are generally available to all processes, while others depend on the `executor` currently defined.

(process-accelerator)=

### accelerator

:::{versionadded} 19.09.0-edge
:::

The `accelerator` directive allows you to request hardware accelerators (e.g. GPUs) for the task execution. For example:

```groovy
process foo {
    accelerator 4, type: 'nvidia-tesla-k80'

    script:
    """
    your_gpu_enabled --command --line
    """
}
```

The above examples will request 4 GPUs of type `nvidia-tesla-k80`.

:::{note}
This directive is only used by certain executors. Refer to the {ref}`executor-page` page to see which executors support this directive.
:::

:::{note}
The accelerator `type` option depends on the target execution platform. Refer to the platform-specific documentation for details on the available accelerators:

- [Google Cloud](https://cloud.google.com/compute/docs/gpus/)
- [Kubernetes](https://kubernetes.io/docs/tasks/manage-gpus/scheduling-gpus/#clusters-containing-different-types-of-gpus)

The accelerator `type` option is not supported for AWS Batch. You can control the accelerator type indirectly through the allowed instance types in your Compute Environment. See the [AWS Batch FAQs](https://aws.amazon.com/batch/faqs/?#GPU_Scheduling_) for more information.
:::

(process-afterscript)=

### afterScript

The `afterScript` directive allows you to execute a custom (Bash) snippet immediately *after* the main process has run. This may be useful to clean up your staging area.

:::{note}
When combined with the {ref}`container directive <process-container>`, the `afterScript` will be executed outside the specified container. In other words, the `afterScript` is always executed in the host environment.
:::

(process-arch)=

### arch

The `arch` directive allows you to define the CPU architecture to build the software in use by the process' task. For example:

```groovy
process cpu_task {
    spack 'blast-plus@2.13.0'
    arch 'linux/x86_64', target: 'cascadelake'

    """
    blastp -query input_sequence -num_threads ${task.cpus}
    """
}
```

The example above declares that the CPU generic architecture is `linux/x86_64` (X86 64 bit), and more specifically that the microarchitecture is `cascadelake` (a specific generation of Intel CPUs).

This directive is currently used by the following Nextflow functionalities:

- by the [spack](#spack) directive, to build microarchitecture-optimised applications;
- by the {ref}`wave-page` service, to build containers for one of the generic families of CPU architectures (see below);
- by the `spack` strategy within {ref}`wave-page`, to optimise the container builds for specific CPU microarchitectures.

Allowed values for the `arch` directive are as follows, grouped by equivalent family (choices available for the sake of compatibility):
- X86 64 bit: `linux/x86_64`, `x86_64`, `linux/amd64`, `amd64`
- ARM 64 bit: `linux/aarch64`, `aarch64`, `linux/arm64`, `arm64`, `linux/arm64/v8`
- ARM 64 bit, older generation: `linux/arm64/v7`

Examples of values for the architecture `target` option are `cascadelake`, `icelake`, `zen2` and `zen3`. See the Spack documentation for the full and up-to-date [list of meaningful targets](https://spack.readthedocs.io/en/latest/basic_usage.html#support-for-specific-microarchitectures).

(process-beforescript)=

### beforeScript

The `beforeScript` directive allows you to execute a custom (Bash) snippet *before* the main process script is run. This may be useful to initialise the underlying cluster environment or for other custom initialisation.

For example:

```groovy
process foo {
  beforeScript 'source /cluster/bin/setup'

  """
  echo bar
  """
}
```

:::{note}
When combined with the {ref}`container directive <process-container>`, the `beforeScript` will be executed outside the specified container. In other words, the `beforeScript` is always executed in the host environment.
:::

(process-cache)=

### cache

The `cache` directive allows you to store the process results to a local cache. When the cache is enabled *and* the pipeline is launched with the {ref}`resume <getstarted-resume>` option, any task executions that are already cached will be re-used. See the {ref}`cache-resume-page` page for more information about how the cache works.

The cache is enabled by default, but you can disable it for a specific process by setting the `cache` directive to `false`. For example:

```groovy
process noCacheThis {
  cache false

  script:
  <your command string here>
}
```

The following options are available:

`false`
: Disable caching.

`true` (default)
: Enable caching. Input file metadata (name, size, last updated timestamp) are included in the cache keys.

`'deep'`
: Enable caching. Input file content is included in the cache keys.

`'lenient'`
: Enable caching. Minimal input file metadata (name and size only) are included in the cache keys.
: This strategy provides a workaround for incorrect caching invalidation observed on shared file systems due to inconsistent file timestamps.

(process-clusteroptions)=

### clusterOptions

The `clusterOptions` directive allows the usage of any native configuration option accepted by your cluster submit command. You can use it to request non-standard resources or use settings that are specific to your cluster and not supported out of the box by Nextflow.

:::{note}
This directive is only used by grid executors. Refer to the {ref}`executor-page` page to see which executors support this directive.
:::

:::{warning}
While you can use the `clusterOptions` directive to specify options that are supported as process directives (`queue`, `memory`, `time`, etc), you should not use both at the same time, as it will cause undefined behavior. Most HPC schedulers will either fail or simply ignore one or the other.
:::

(process-conda)=

### conda

The `conda` directive allows for the definition of the process dependencies using the [Conda](https://conda.io) package manager.

Nextflow automatically sets up an environment for the given package names listed by in the `conda` directive. For example:

```groovy
process foo {
  conda 'bwa=0.7.15'

  '''
  your_command --here
  '''
}
```

Multiple packages can be specified separating them with a blank space e.g. `bwa=0.7.15 fastqc=0.11.5`. The name of the channel from where a specific package needs to be downloaded can be specified using the usual Conda notation i.e. prefixing the package with the channel name as shown here `bioconda::bwa=0.7.15`.

The `conda` directive also allows the specification of a Conda environment file path or the path of an existing environment directory. See the {ref}`conda-page` page for further details.

(process-container)=

### container

The `container` directive allows you to execute the process script in a [Docker](http://docker.io) container.

It requires the Docker daemon to be running in machine where the pipeline is executed, i.e. the local machine when using the *local* executor or the cluster nodes when the pipeline is deployed through a *grid* executor.

For example:

```groovy
process runThisInDocker {
  container 'dockerbox:tag'

  """
  <your holy script here>
  """
}
```

Simply replace in the above script `dockerbox:tag` with the name of the Docker image you want to use.

:::{tip}
Containers are a very useful way to execute your scripts in a reproducible self-contained environment or to run your pipeline in the cloud.
:::

:::{note}
This directive is ignored for processes that are {ref}`executed natively <process-native>`.
:::

(process-containeroptions)=

### containerOptions

The `containerOptions` directive allows you to specify any container execution option supported by the underlying container engine (ie. Docker, Singularity, etc). This can be useful to provide container settings only for a specific process e.g. mount a custom path:

```groovy
process runThisWithDocker {
    container 'busybox:latest'
    containerOptions '--volume /data/db:/db'

    output:
    path 'output.txt'

    '''
    your_command --data /db > output.txt
    '''
}
```

:::{warning}
This feature is not supported by the {ref}`k8s-executor` and {ref}`google-lifesciences-executor` executors.
:::

(process-cpus)=

### cpus

The `cpus` directive allows you to define the number of (logical) CPU required by the process' task. For example:

```groovy
process big_job {
  cpus 8
  executor 'sge'

  """
  blastp -query input_sequence -num_threads ${task.cpus}
  """
}
```

This directive is required for tasks that execute multi-process or multi-threaded commands/tools and it is meant to reserve enough CPUs when a pipeline task is executed through a cluster resource manager.

See also: [penv](#penv), [memory](#memory), [time](#time), [queue](#queue), [maxForks](#maxforks)

(process-debug)=

### debug

By default the `stdout` produced by the commands executed in all processes is ignored. By setting the `debug` directive to `true`, you can forward the process `stdout` to the current top running process `stdout` file, showing it in the shell terminal.

For example:

```groovy
process sayHello {
  debug true

  script:
  "echo Hello"
}
```

```
Hello
```

Without specifying `debug true`, you won't see the `Hello` string printed out when executing the above example.

(process-disk)=

### disk

The `disk` directive allows you to define how much local disk storage the process is allowed to use. For example:

```groovy
process big_job {
    disk '2 GB'
    executor 'cirrus'

    """
    your task script here
    """
}
```

The following memory unit suffix can be used when specifying the disk value:

| Unit | Description |
| ---- | ----------- |
| B    | Bytes       |
| KB   | Kilobytes   |
| MB   | Megabytes   |
| GB   | Gigabytes   |
| TB   | Terabytes   |

See {ref}`implicit-classes-memoryunit` for more information.

:::{note}
This directive is only used by certain executors. Refer to the {ref}`executor-page` page to see which executors support this directive.
:::

See also: [cpus](#cpus), [memory](#memory) [time](#time), [queue](#queue) and [Dynamic computing resources](#dynamic-computing-resources).

(process-echo)=

### echo

:::{deprecated} 22.04.0
Use `debug` instead
:::

(process-error-strategy)=

### errorStrategy

The `errorStrategy` directive allows you to define how an error condition is managed by the process. By default when an error status is returned by the executed script, the process stops immediately. This in turn forces the entire pipeline to terminate.

The following error strategies are available:

`terminate` (default)
: Terminate the execution as soon as an error condition is reported. Pending jobs are killed.

`finish`
: Initiate an orderly pipeline shutdown when an error condition is raised, waiting for the completion of any submitted jobs.

`ignore`
: Ignore process execution errors.

`retry`
: Re-submit any process that returns an error condition.

When setting the `errorStrategy` directive to `ignore` the process doesn't stop on an error condition, it just reports a message notifying you of the error event.

For example:

```groovy
process ignoreAnyError {
  errorStrategy 'ignore'

  script:
  <your command string here>
}
```

:::{note}
By definition, a command script fails when it ends with a non-zero exit status.
:::

The `retry` error strategy allows you to re-submit for execution a process returning an error condition. For example:

```groovy
process retryIfFail {
  errorStrategy 'retry'

  script:
  <your command string here>
}
```

The number of times a failing process is re-executed is defined by the [maxRetries](#maxretries) and [maxErrors](#maxerrors) directives.

:::{tip}
More complex strategies depending on the task exit status or other parametric values can be defined using a dynamic `errorStrategy`. See the [Dynamic directives](#dynamic-directives) section for details.
:::

See also: [maxErrors](#maxerrors), [maxRetries](#maxretries) and [Dynamic computing resources](#dynamic-computing-resources).

(process-executor)=

### executor

The `executor` defines the underlying system where processes are executed. By default a process uses the executor defined globally in the `nextflow.config` file.

The `executor` directive allows you to configure what executor has to be used by the process, overriding the default configuration.

The following executors are available:

| Name                  | Executor                                                                                    |
| --------------------- | ------------------------------------------------------------------------------------------- |
| `awsbatch`            | [AWS Batch](https://aws.amazon.com/batch/) service                                          |
| `azurebatch`          | [Azure Batch](https://azure.microsoft.com/en-us/services/batch/) service                    |
| `condor`              | [HTCondor](https://research.cs.wisc.edu/htcondor/) job scheduler                            |
| `google-lifesciences` | [Google Genomics Pipelines](https://cloud.google.com/life-sciences) service                 |
| `ignite`              | [Apache Ignite](https://ignite.apache.org/) cluster                                         |
| `k8s`                 | [Kubernetes](https://kubernetes.io/) cluster                                                |
| `local`               | The computer where `Nextflow` is launched                                                   |
| `lsf`                 | [Platform LSF](http://en.wikipedia.org/wiki/Platform_LSF) job scheduler                     |
| `moab`                | [Moab](http://www.adaptivecomputing.com/moab-hpc-basic-edition/) job scheduler              |
| `nqsii`               | [NQSII](https://www.rz.uni-kiel.de/en/our-portfolio/hiperf/nec-linux-cluster) job scheduler |
| `oge`                 | Alias for the `sge` executor                                                                |
| `pbs`                 | [PBS/Torque](http://en.wikipedia.org/wiki/Portable_Batch_System) job scheduler              |
| `pbspro`              | [PBS Pro](https://www.pbsworks.com/) job scheduler                                          |
| `sge`                 | Sun Grid Engine / [Open Grid Engine](http://gridscheduler.sourceforge.net/)                 |
| `slurm`               | [SLURM](https://en.wikipedia.org/wiki/Slurm_Workload_Manager) workload manager              |
| `tes`                 | [GA4GH TES](https://github.com/ga4gh/task-execution-schemas) service                        |
| `uge`                 | Alias for the `sge` executor                                                                |

The following example shows how to set the process's executor:

```groovy
process doSomething {
  executor 'sge'

  script:
  <your script here>
}
```

:::{note}
Each executor supports additional directives and `executor` configuration options. Refer to the {ref}`executor-page` page to see what each executor supports.
:::

(process-ext)=

### ext

The `ext` is a special directive used as *namespace* for user custom process directives. This can be useful for advanced configuration options. For example:

```groovy
process mapping {
  container "biocontainers/star:${task.ext.version}"

  input:
  path genome
  tuple val(sampleId), path(reads)

  """
  STAR --genomeDir $genome --readFilesIn $reads ${task.ext.args ?: ''}
  """
}
```

In the above example, the process container version is controlled by `ext.version`, and the script supports additional command line arguments through `ext.args`.

The `ext` directive can be set in the process definition:

```groovy
process mapping {
  ext version: '2.5.3', args: '--foo --bar'
}
```

Or in the Nextflow configuration:

```groovy
process.ext.version = '2.5.3'
process.ext.args = '--foo --bar'
```

(process-fair)=

### fair

:::{versionadded} 22.12.0-edge
:::

The `fair` directive, when enabled, guarantees that process outputs will be emitted in the order in which they were received. For example:

```groovy
process foo {
    fair true

    input:
    val x
    output:
    tuple val(task.index), val(x)

    script:
    """
    sleep \$((RANDOM % 3))
    """
}

workflow {
    channel.of('A','B','C','D') | foo | view
}
```

The above example produces:

```
[1, A]
[2, B]
[3, C]
[4, D]
```

(process-label)=

### label

The `label` directive allows the annotation of processes with mnemonic identifier of your choice. For example:

```groovy
process bigTask {
  label 'big_mem'

  '''
  <task script>
  '''
}
```

The same label can be applied to more than a process and multiple labels can be applied to the same process using the `label` directive more than one time.

:::{note}
A label must consist of alphanumeric characters or `_`, must start with an alphabetic character and must end with an alphanumeric character.
:::

Labels are useful to organise workflow processes in separate groups which can be referenced in the configuration file to select and configure subset of processes having similar computing requirements. See the {ref}`config-process-selectors` documentation for details.

See also: [resourceLabels](#resourcelabels)

(process-machinetype)=

### machineType

:::{versionadded} 19.07.0
:::

The `machineType` can be used to specify a predefined Google Compute Platform [machine type](https://cloud.google.com/compute/docs/machine-types) when running using the {ref}`Google Life Sciences <google-lifesciences-executor>` executor.

This directive is optional and if specified overrides the cpus and memory directives:

```groovy
process foo {
  machineType 'n1-highmem-8'

  """
  <your script here>
  """
}
```

See also: [cpus](#cpus) and [memory](#memory).

(process-maxsubmitawait)=

### maxSubmitAwait (experimental)

The `maxSubmitAwait` directives allows you to specify how long a task can remain in submission queue without being executed.
Elapsed this time the task execution will fail.

When used along with `retry` error strategy, it can be useful to re-schedule the task to a difference queue or
resource requirement. For example:

```groovy
process foo {
  errorStrategy 'retry'
  maxSubmitAwait '10 mins'
  maxRetries 3
  queue "${task.submitAttempt==1 : 'spot-compute' : 'on-demand-compute'}"
  script:
  '''
  your_job --here
  '''
}
```

In the above example the task is submitted to the `spot-compute` on the first attempt (`task.submitAttempt==1`). If the
task execution does not start in the 10 minutes, a failure is reported and a new submission is attempted using the
queue named `on-demand-compute`. 

(process-maxerrors)=

### maxErrors

The `maxErrors` directive allows you to specify the maximum number of times a process can fail when using the `retry` error strategy. By default this directive is disabled, you can set it as shown in the example below:

```groovy
process retryIfFail {
  errorStrategy 'retry'
  maxErrors 5

  """
  echo 'do this as that .. '
  """
}
```

:::{note}
This setting considers the **total** errors accumulated for a given process, across all instances. If you want to control the number of times a process **instance** (aka task) can fail, use `maxRetries`.
:::

See also: [errorStrategy](#errorstrategy) and [maxRetries](#maxretries).

(process-maxforks)=

### maxForks

The `maxForks` directive allows you to define the maximum number of process instances that can be executed in parallel. By default this value is equals to the number of CPU cores available minus 1.

If you want to execute a process in a sequential manner, set this directive to one. For example:

```groovy
process doNotParallelizeIt {
  maxForks 1

  '''
  <your script here>
  '''
}
```

(process-maxretries)=

### maxRetries

The `maxRetries` directive allows you to define the maximum number of times a process instance can be re-submitted in case of failure. This value is applied only when using the `retry` error strategy. By default only one retry is allowed, you can increase this value as shown below:

```groovy
process retryIfFail {
    errorStrategy 'retry'
    maxRetries 3

    """
    echo 'do this as that .. '
    """
}
```

:::{note}
There is a subtle but important difference between `maxRetries` and the `maxErrors` directive. The latter defines the total number of errors that are allowed during the process execution (the same process can launch different execution instances), while the `maxRetries` defines the maximum number of times the same process execution can be retried in case of an error.
:::

See also: [errorStrategy](#errorstrategy) and [maxErrors](#maxerrors).

(process-memory)=

### memory

The `memory` directive allows you to define how much memory the process is allowed to use. For example:

```groovy
process big_job {
    memory '2 GB'
    executor 'sge'

    """
    your task script here
    """
}
```

The following memory unit suffix can be used when specifying the memory value:

| Unit | Description |
| ---- | ----------- |
| B    | Bytes       |
| KB   | Kilobytes   |
| MB   | Megabytes   |
| GB   | Gigabytes   |
| TB   | Terabytes   |

See {ref}`implicit-classes-memoryunit` for more information.

See also: [cpus](#cpus), [time](#time), [queue](#queue) and [Dynamic computing resources](#dynamic-computing-resources).

(process-module)=

### module

[Environment Modules](http://modules.sourceforge.net/) is a package manager that allows you to dynamically configure your execution environment and easily switch between multiple versions of the same software tool.

If it is available in your system you can use it with Nextflow in order to configure the processes execution environment in your pipeline.

In a process definition you can use the `module` directive to load a specific module version to be used in the process execution environment. For example:

```groovy
process basicExample {
  module 'ncbi-blast/2.2.27'

  """
  blastp -query <etc..>
  """
}
```

You can repeat the `module` directive for each module you need to load. Alternatively multiple modules can be specified in a single `module` directive by separating all the module names by using a `:` (colon) character as shown below:

```groovy
 process manyModules {

   module 'ncbi-blast/2.2.27:t_coffee/10.0:clustalw/2.1'

   """
   blastp -query <etc..>
   """
}
```

(process-penv)=

### penv

The `penv` directive allows you to define the parallel environment to be used when submitting a parallel task to the {ref}`SGE <sge-executor>` resource manager. For example:

```groovy
process big_job {
  cpus 4
  penv 'smp'
  executor 'sge'

  """
  blastp -query input_sequence -num_threads ${task.cpus}
  """
}
```

This configuration depends on the parallel environment provided by your grid engine installation. Refer to your cluster documentation or contact your admin to learn more about this.

See also: [cpus](#cpus), [memory](#memory), [time](#time)

(process-pod)=

### pod

The `pod` directive allows the definition of pod specific settings, such as environment variables, secrets, and config maps, when using the {ref}`k8s-executor` executor.

For example:

```groovy
process your_task {
  pod env: 'FOO', value: 'bar'

  '''
  echo $FOO
  '''
}
```

The above snippet defines an environment variable named `FOO` whose value is `bar`.

When defined in the Nextflow configuration file, a pod setting can be defined as a map:

```groovy
process {
  pod = [env: 'FOO', value: 'bar']
}
```

Or as a list of maps:

```groovy
process {
  pod = [
    [env: 'FOO', value: 'bar'],
    [secret: 'my-secret/key1', mountPath: '/etc/file.txt']
  ]
}
```

The following options are available:

`affinity: <config>`
: :::{versionadded} 22.01.0-edge
  :::
: Specifies the pod [affinity](https://kubernetes.io/docs/concepts/scheduling-eviction/assign-pod-node/#affinity-and-anti-affinity) with the given configuration.

`annotation: '<name>', value: '<value>'`
: *Can be specified multiple times*
: Defines a pod [annotation](https://kubernetes.io/docs/concepts/overview/working-with-objects/annotations/) with the given name and value.

`automountServiceAccountToken: true | false`
: :::{versionadded} 22.01.0-edge
  :::
: Specifies whether to [automount service account token](https://kubernetes.io/docs/tasks/configure-pod-container/configure-service-account/#opt-out-of-api-credential-automounting) into the pod (default: `true`).

`config: '<configMap>/<key>', mountPath: '</absolute/path>'`
: *Can be specified multiple times*
: Mounts a [ConfigMap](https://kubernetes.io/docs/tasks/configure-pod-container/configure-pod-configmap/) with name and optional key to the given path. If the key is omitted, the path is interpreted as a directory and all entries in the `ConfigMap` are exposed in that path.

`csi: '<name>', mountPath: '</absolute/path>'`
: :::{versionadded} 22.11.0-edge
  :::
: *Can be specified multiple times*
: Mounts a [CSI ephemeral volume](https://kubernetes.io/docs/concepts/storage/ephemeral-volumes/#csi-ephemeral-volumes) by name to the given path.

`emptyDir: <config>, mountPath: '</absolute/path>'`
: :::{versionadded} 22.11.0-edge
  :::
: *Can be specified multiple times*
: Mounts an [emptyDir](https://kubernetes.io/docs/concepts/storage/volumes/#emptydir) with the given configuration to the given path.

`env: '<name>', config: '<configMap>/<key>'`
: *Can be specified multiple times*
: Defines an environment variable whose value is defined by the given [ConfigMap](https://kubernetes.io/docs/tasks/configure-pod-container/configure-pod-configmap/) and key.

`env: '<name>', fieldPath: '<fieldPath>'`
: :::{versionadded} 21.09.1-edge
  :::
: *Can be specified multiple times*
: Defines an environment variable whose value is defined by the given [field path](https://kubernetes.io/docs/tasks/inject-data-application/environment-variable-expose-pod-information/#use-pod-fields-as-values-for-environment-variables) value.

: For example, the following pod option:

  ```groovy
  pod = [env: 'MY_NODE_NAME', fieldPath: 'spec.nodeName']
  ```

  Maps to the following pod spec:

  ```yaml
  env:
    - name: MY_NODE_NAME
      valueFrom:
        fieldRef:
          fieldPath: spec.nodeName
  ```

`env: '<name>', secret: '<secret>/<key>'`
: *Can be specified multiple times*
: Defines an environment variable whose value is defined by the given [Secret](https://kubernetes.io/docs/concepts/configuration/secret/) and key.

`env: '<name>', value: '<value>'`
: *Can be specified multiple times*
: Defines an environment variable with the given name and value.

`hostPath: '/host/absolute/path', mountPath: '</pod/absolute/path>'`
: :::{versionadded} 23.10.0
  :::
: *Can be specified multiple times*
: Allows creating [hostPath](https://kubernetes.io/docs/concepts/storage/volumes/#hostpath) volume and access it with the specified `mountPath` in the pod.

`imagePullPolicy: 'IfNotPresent' | 'Always' | 'Never'`
: Specifies the [image pull policy](https://kubernetes.io/docs/concepts/containers/images/#image-pull-policy) used by the pod to pull the container image.

`imagePullSecret: '<name>'`
: Specifies the [image pull secret](https://kubernetes.io/docs/concepts/containers/images/#specifying-imagepullsecrets-on-a-pod) used to access a private container image registry.

`label: '<name>', value: '<value>'`
: *Can be specified multiple times*
: Defines a pod [label](https://kubernetes.io/docs/concepts/overview/working-with-objects/labels/) with the given name and value.

`nodeSelector: <config>`
: Specifies the [node selector](https://kubernetes.io/docs/concepts/scheduling-eviction/assign-pod-node/#nodeselector) with the given configuration.

: The configuration can be a map or a string:

  ```groovy
  // map
  pod = [nodeSelector: [disktype: 'ssd', cpu: 'intel']]

  // string
  pod = [nodeSelector: 'disktype=ssd,cpu=intel']
  ```

`priorityClassName: '<name>'`
: :::{versionadded} 22.01.0-edge
  :::
: Specifies the [priority class name](https://kubernetes.io/docs/concepts/scheduling-eviction/pod-priority-preemption/) for pods.

`privileged: true | false`
: :::{versionadded} 22.05.0-edge
  :::
: Specifies whether the pod should run as a *privileged* container (default: `false`).

`runAsUser: '<uid>'`
: Specifies the user ID with which to run the container. Shortcut for the `securityContext` option.

`schedulerName: '<name>'`
: Specifies which [scheduler](https://kubernetes.io/docs/tasks/extend-kubernetes/configure-multiple-schedulers/#specify-schedulers-for-pods) is used to schedule the container. 

`secret: '<secret>/<key>', mountPath: '</absolute/path>'`
: *Can be specified multiple times*
: Mounts a [Secret](https://kubernetes.io/docs/concepts/configuration/secret/) with name and optional key to the given path. If the key is omitted, the path is interpreted as a directory and all entries in the `Secret` are exposed in that path.

`securityContext: <config>`
: Specifies the pod [security context](https://kubernetes.io/docs/tasks/configure-pod-container/security-context/) with the given configuration.

`toleration: <config>`
: :::{versionadded} 22.04.0
  :::
: *Can be specified multiple times*
: Specifies the pod [toleration](https://kubernetes.io/docs/concepts/scheduling-eviction/taint-and-toleration/) with the given configuration.

: The configuration should be a map corresponding to a single toleration rule. For example, the following pod options:

  ```groovy
  pod = [
      [toleration: [key: 'key1', operator: 'Equal', value: 'value1', effect: 'NoSchedule']],
      [toleration: [key: 'key1', operator: 'Exists', effect: 'NoSchedule']],
  ]
  ```

  Maps to the following pod spec:

  ```yaml
  tolerations:
    - key: "key1"
      operator: "Equal"
      value: "value1"
      effect: "NoSchedule"
    - key: "key1"
      operator: "Exists"
      effect: "NoSchedule"
  ```

`volumeClaim: '<name>', mountPath: '</absolute/path>' [, subPath: '<path>', readOnly: true | false]`
: *Can be specified multiple times*
: Mounts a [Persistent volume claim](https://kubernetes.io/docs/concepts/storage/persistent-volumes/) with the given name to the given path.
: The `subPath` option can be used to mount a sub-directory of the volume instead of its root.
: The `readOnly` option can be used to mount the volume as read-only (default: `false`)

(process-publishdir)=

### publishDir

The `publishDir` directive allows you to publish the process output files to a specified folder. For example:

```groovy
process foo {
    publishDir '/data/chunks'

    output:
    path 'chunk_*'

    '''
    printf 'Hola' | split -b 1 - chunk_
    '''
}
```

The above example splits the string `Hola` into file chunks of a single byte. When complete the `chunk_*` output files are published into the `/data/chunks` folder.

:::{note}
Only files that match the declaration in the `output` block are published, not all the outputs of the process.
:::

:::{tip}
The `publishDir` directive can be specified more than once in order to publish output files to different target directories based on different rules.
:::

By default files are published to the target folder creating a *symbolic link* for each process output that links the file produced into the process working directory. This behavior can be modified using the `mode` option, for example:

```groovy
process foo {
    publishDir '/data/chunks', mode: 'copy', overwrite: false

    output:
    path 'chunk_*'

    '''
    printf 'Hola' | split -b 1 - chunk_
    '''
}
```

:::{warning}
Files are copied into the specified directory in an *asynchronous* manner, so they may not be immediately available in the publish directory at the end of the process execution. For this reason, downstream processes should not try to access output files through the publish directory, but through channels.
:::

Available options:

`contentType`
: :::{versionadded} 22.10.0
  :::
: *Experimental: currently only supported for S3.*
: Allow specifying the media content type of the published file a.k.a. [MIME type](https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/MIME_Types). If set to `true`, the content type is inferred from the file extension (default: `false`).

`enabled`
: Enable or disable the publish rule depending on the boolean value specified (default: `true`).

`failOnError`
: When `true` abort the execution if some file can't be published to the specified target directory or bucket for any cause (default: `false`)

`mode`
: The file publishing method. Can be one of the following values:

  - `'copy'`: Copies the output files into the publish directory.
  - `'copyNoFollow'`: Copies the output files into the publish directory without following symlinks ie. copies the links themselves.
  - `'link'`: Creates a hard link in the publish directory for each output file.
  - `'move'`: Moves the output files into the publish directory. **Note**: this is only supposed to be used for a *terminal* process i.e. a process whose output is not consumed by any other downstream process.
  - `'rellink'`: Creates a relative symbolic link in the publish directory for each output file.
  - `'symlink'`: Creates an absolute symbolic link in the publish directory for each output file (default).

`overwrite`
: When `true` any existing file in the specified folder will be overridden (default: `true` during normal pipeline execution and `false` when pipeline execution is `resumed`).

`path`
: Specifies the directory where files need to be published. **Note**: the syntax `publishDir '/some/dir'` is a shortcut for `publishDir path: '/some/dir'`.

`pattern`
: Specifies a [glob][glob] file pattern that selects which files to publish from the overall set of output files.

`saveAs`
: A closure which, given the name of the file being published, returns the actual file name or a full path where the file is required to be stored. This can be used to rename or change the destination directory of the published files dynamically by using a custom strategy. Return the value `null` from the closure to *not* publish a file. This is useful when the process has multiple output files, but you want to publish only some of them.

`storageClass`
: :::{versionadded} 22.12.0-edge
  :::
: *Experimental: currently only supported for S3.*
: Allow specifying the storage class to be used for the published file.

`tags`
: :::{versionadded} 21.12.0-edge
  :::
: *Experimental: currently only supported for S3.*
: Allow the association of arbitrary tags with the published file e.g. `tags: [FOO: 'Hello world']`.

(process-queue)=

### queue

The `queue` directive allows you to set the `queue` where jobs are scheduled when using a grid based executor in your pipeline. For example:

```groovy
process grid_job {
    queue 'long'
    executor 'sge'

    """
    your task script here
    """
}
```

Multiple queues can be specified by separating their names with a comma for example:

```groovy
process grid_job {
    queue 'short,long,cn-el6'
    executor 'sge'

    """
    your task script here
    """
}
```

:::{note}
This directive is only used by certain executors. Refer to the {ref}`executor-page` page to see which executors support this directive.
:::

(process-resourcelabels)=

### resourceLabels

:::{versionadded} 22.09.1-edge
:::

The `resourceLabels` directive allows you to specify custom name-value pairs that Nextflow applies to the computing resource used to carry out the process execution. Resource labels can be specified using the syntax shown below:

```groovy
process my_task {
    resourceLabels region: 'some-region', user: 'some-username'

    '''
    <task script>
    '''
}
```

The limits and the syntax of the corresponding cloud provider should be taken into consideration when using resource labels.

Resource labels are currently supported by the following executors:

- {ref}`awsbatch-executor`
- {ref}`azurebatch-executor`
- {ref}`google-batch-executor`
- {ref}`google-lifesciences-executor`
- {ref}`k8s-executor`

:::{versionadded} 23.09.0-edge
Resource labels are supported for Azure Batch when using automatic pool creation.

Resource labels in Azure are added to pools, rather than jobs, in order to facilitate cost analysis. A new pool will be created for each new set of resource labels, therefore it is recommended to also set `azure.batch.deletePoolsOnCompletion = true` when using process-specific resource labels.
:::

See also: [label](#label)

(process-scratch)=

### scratch

The `scratch` directive allows you to execute the process in a temporary folder that is local to the execution node.

This is useful when your pipeline is launched by using a grid executor, because it allows you to decrease the NFS overhead by running the pipeline processes in a temporary directory in the local disk of the actual execution node. Only the files declared as output in the process definition will be copied in the pipeline working area.

In its basic form simply specify `true` at the directive value, as shown below:

```groovy
process simpleTask {
  scratch true

  output:
  path 'data_out'

  '''
  <task script>
  '''
}
```

By doing this, it tries to execute the script in the directory defined by the variable `$TMPDIR` in the execution node. If this variable does not exist, it will create a new temporary directory by using the Linux command `mktemp`.

:::{note}
Cloud-based executors use `scratch = true` by default, since the work directory resides in object storage.
:::

The following values are supported:

`false`
: Do not use a scratch directory.

`true`
: Create a scratch directory in the directory defined by the `$TMPDIR` environment variable, or `$(mktemp /tmp)` if `$TMPDIR` is not set.

`'$YOUR_VAR'`
: Create a scratch directory in the directory defined by the given environment variable, or `$(mktemp /tmp)` if that variable is not set. The value must use single quotes, otherwise the environment variable will be evaluated in the pipeline script context.

`'/my/tmp/path'`
: Create a scratch directory in the specified directory.

`'ram-disk'`
: Create a scratch directory in the RAM disk `/dev/shm/`.

(process-directive-shell)=

### shell

The `shell` directive allows you to define a custom shell command for process scripts. By default, script blocks are executed with `/bin/bash -ue`.

```groovy
process doMoreThings {
    shell '/bin/bash', '-euo', 'pipefail'

    '''
    your_command_here
    '''
}
```

The same directive could be specified in your Nextflow configuration as follows:

```groovy
process.shell = ['/bin/bash', '-euo', 'pipefail']
```

(process-spack)=

### spack

The `spack` directive allows for the definition of the process dependencies using the [Spack](https://spack.io) package manager.

Nextflow automatically sets up an environment for the given package names listed by in the `spack` directive. For example:

```groovy
process foo {
    spack 'bwa@0.7.15'

    '''
    your_command --here
    '''
}
```

Multiple packages can be specified separating them with a blank space, e.g. `bwa@0.7.15 fastqc@0.11.5`.

The `spack` directive also allows the specification of a Spack environment file path or the path of an existing environment directory. See the {ref}`spack-page` page for further details.

(process-stageinmode)=

### stageInMode

The `stageInMode` directive defines how input files are staged into the process work directory. The following values are allowed:

`'copy'`
: Input files are staged in the process work directory by creating a copy.

`'link'`
: Input files are staged in the process work directory by creating a hard link for each of them.

`'rellink'`
: Input files are staged in the process work directory by creating a symbolic link with a relative path for each of them.

`'symlink'`
: Input files are staged in the process work directory by creating a symbolic link with an absolute path for each of them (default).

(process-stageoutmode)=

### stageOutMode

The `stageOutMode` directive defines how output files are staged out from the scratch directory to the process work directory. The following values are allowed:

`'copy'`
: Output files are copied from the scratch directory to the work directory.

`'fcp'`
: :::{versionadded} 23.02.0-edge
  :::
: Output files are copied from the scratch directory to the work directory by using the [fcp](https://github.com/Svetlitski/fcp) utility (note: it must be available in your cluster computing nodes).

`'move'`
: Output files are moved from the scratch directory to the work directory.

`'rclone'`
: :::{versionadded} 23.01.0-edge
  :::
: Output files are copied from the scratch directory to the work directory by using the [rclone](https://rclone.org) utility (note: it must be available in your cluster computing nodes).

`'rsync'`
: Output files are copied from the scratch directory to the work directory by using the `rsync` utility.

See also: [scratch](#scratch).

(process-storedir)=

### storeDir

The `storeDir` directive allows you to define a directory that is used as a *permanent* cache for your process results.

In more detail, it affects the process execution in two main ways:

1. The process is executed only if the files declared in the `output` block do not exist in the directory specified by the `storeDir` directive. When the files exist the process execution is skipped and these files are used as the actual process result.
2. Whenever a process is successfully completed the files listed in the `output` block are moved into the directory specified by the `storeDir` directive.

The following example shows how to use the `storeDir` directive to create a directory containing a BLAST database for each species specified by an input parameter:

```groovy
process formatBlastDatabases {
  storeDir '/db/genomes'

  input:
  path species

  output:
  path "${dbName}.*"

  script:
  dbName = species.baseName
  """
  makeblastdb -dbtype nucl -in ${species} -out ${dbName}
  """
}
```

:::{warning}
The `storeDir` directive is meant for long-term process caching and should not be used to publish output files or organize outputs into a semantic directory structure. In those cases, use the [publishDir](#publishdir) directive instead.
:::

:::{note}
The use of AWS S3 paths is supported, however it requires the installation of the [AWS CLI](https://aws.amazon.com/cli/) (i.e. `aws`) in the target compute node.
:::

(process-tag)=

### tag

The `tag` directive allows you to associate each process execution with a custom label, so that it will be easier to identify them in the log file or in the trace execution report. For example:

```groovy
process foo {
  tag "$code"

  input:
  val code

  """
  echo $code
  """
}

workflow {
  Channel.of('alpha', 'gamma', 'omega') | foo
}
```

The above snippet will print a log similar to the following one, where process names contain the tag value:

```
[6e/28919b] Submitted process > foo (alpha)
[d2/1c6175] Submitted process > foo (gamma)
[1c/3ef220] Submitted process > foo (omega)
```

See also {ref}`Trace execution report <trace-report>`

(process-time)=

### time

The `time` directive allows you to define how long a process is allowed to run. For example:

```groovy
process big_job {
    time '1h'

    """
    your task script here
    """
}
```

The following time unit suffixes can be used when specifying the duration value:

| Unit                            | Description  |
| ------------------------------- | ------------ |
| `ms`, `milli`, `millis`         | Milliseconds |
| `s`, `sec`, `second`, `seconds` | Seconds      |
| `m`, `min`, `minute`, `minutes` | Minutes      |
| `h`, `hour`, `hours`            | Hours        |
| `d`, `day`, `days`              | Days         |

Multiple units can be used in a single declaration, for example: `'1day 6hours 3minutes 30seconds'`

See {ref}`implicit-classes-duration` for more information.

:::{note}
This directive is only used by certain executors. Refer to the {ref}`executor-page` page to see which executors support this directive.
:::

See also: [cpus](#cpus), [memory](#memory), [queue](#queue) and [Dynamic computing resources](#dynamic-computing-resources).

### Dynamic directives

A directive can be assigned *dynamically*, during the process execution, so that its actual value can be evaluated based on the process inputs.

In order to be defined in a dynamic manner, the directive's value needs to be expressed using a {ref}`closure <script-closure>`, as in the following example:

```groovy
process foo {
  executor 'sge'
  queue { entries > 100 ? 'long' : 'short' }

  input:
  tuple val(entries), path('data.txt')

  script:
  """
  < your job here >
  """
}
```

In the above example, the [queue](#queue) directive is evaluated dynamically, depending on the input value `entries`. When it is larger than 100, jobs will be submitted to the `long` queue, otherwise the `short` queue will be used.

All directives can be assigned a dynamic value except the following:

- [executor](#executor)
- [label](#label)
- [maxForks](#maxforks)

:::{tip}
Assigning a string value with one or more variables is always resolved in a dynamic manner, and therefore is equivalent to the above syntax. For example, the above directive can also be written as:

```groovy
queue "${ entries > 100 ? 'long' : 'short' }"
```

Note, however, that the latter syntax can be used both for a directive's main argument (as in the above example) and for a directive's optional named attributes, whereas the closure syntax is only resolved dynamically for a directive's main argument.
:::

:::{tip}
You can retrieve the current value of a dynamic directive in the process script by using the implicit variable `task`, which holds the directive values defined in the current task. For example:

```groovy
process foo {
  queue { entries > 100 ? 'long' : 'short' }

  input:
  tuple val(entries), path('data.txt')

  script:
  """
  echo Current queue: ${task.queue}
  """
}
```
:::

### Dynamic computing resources

It's a very common scenario that different instances of the same process may have very different needs in terms of computing resources. In such situations requesting, for example, an amount of memory too low will cause some tasks to fail. Instead, using a higher limit that fits all the tasks in your execution could significantly decrease the execution priority of your jobs.

The [Dynamic directives](#dynamic-directives) evaluation feature can be used to modify the amount of computing resources requested in case of a process failure and try to re-execute it using a higher limit. For example:

```groovy
process foo {
    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    script:
    <your job here>
}
```

In the above example the [memory](#memory) and execution [time](#time) limits are defined dynamically. The first time the process is executed the `task.attempt` is set to `1`, thus it will request a two GB of memory and one hour of maximum execution time.

If the task execution fail reporting an exit status in the range between 137 and 140, the task is re-submitted (otherwise terminates immediately). This time the value of `task.attempt` is `2`, thus increasing the amount of the memory to four GB and the time to 2 hours, and so on.

The directive [maxRetries](#maxretries) set the maximum number of time the same task can be re-executed.

### Dynamic Retry with backoff

There are cases in which the required execution resources may be temporary unavailable e.g. network congestion. In these cases immediately re-executing the task will likely result in the identical error. A retry with an exponential backoff delay can better recover these error conditions:

```groovy
process foo {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  maxRetries 5

  script:
  '''
  your_command --here
  '''
}
```

[glob]: http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
[what is a glob?]: http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob

(process-page)=

# Processes

In Nextflow, a **process** is a specialized function for executing scripts in a scalable and portable manner.

Here is an example process definition:

```nextflow
process sayHello {
    output:
    path 'hello.txt'

    script:
    """
    echo 'Hello world!' > hello.txt
    """
}
```

See {ref}`syntax-process` for a full description of the process syntax.

(process-script)=

## Script

The `script` block defines, as a string expression, the script that is executed by the process.

A process may contain only one script, and if the `script` guard is not explicitly declared, the script must be the final statement in the process block.

The script string is executed as a [Bash](<http://en.wikipedia.org/wiki/Bash_(Unix_shell)>) script in the host environment. It can be any command or script that you would normally execute on the command line or in a Bash script. Naturally, the script may only use commands that are available in the host environment.

The script block can be a simple string or a multi-line string. The latter approach makes it easier to write scripts with multiple commands spanning multiple lines. For example:

```nextflow
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

```nextflow
process printPath {
  '''
  echo The path is: $PATH
  '''
}
```

Otherwise, you can define your script with double-quotes and escape the system environment variables by prefixing them with a back-slash `\` character, as shown in the following example:

```nextflow
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

```nextflow
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

The `script` block is like a function that returns a string. This means that you can write arbitrary code to determine the script, as long as the final statement is a string.

If-else statements based on task inputs can be used to produce a different script. For example:

```nextflow
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

In the above example, the process will execute one of several scripts depending on the value of the `mode` parameter. By default it will execute the `tcoffee` command.

(process-template)=

### Template

Process scripts can be externalised to **template** files, which can be reused across different processes and tested independently from the overall pipeline execution.

A template is simply a shell script file that Nextflow is able to execute by using the `template` function as shown below:

```nextflow
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

```nextflow
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

The `exec` block executes the given code without launching a job.

For example:

```nextflow
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

A native process is very similar to a {ref}`function <syntax-function>`. However, it provides additional capabilities such as parallelism, caching, and progress logging.

(process-stub)=

## Stub

:::{versionadded} 20.11.0-edge
:::

You can define a command *stub*, which replaces the actual process command when the `-stub-run` or `-stub` command-line option is enabled:

```nextflow
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
- `path`: Handle the input value as a path, staging the file properly in the execution context.
- `env`: Use the input value to set an environment variable in the process script.
- `stdin`: Forward the input value to the process `stdin` special file.
- `tuple`: Handle a group of input values having any of the above qualifiers.
- `each`: Execute the process for each element in the input collection.

See {ref}`process reference <process-reference-inputs>` for the full list of input methods and options.

### Input variables (`val`)

The `val` qualifier accepts any data type. It can be accessed in the process script by using the specified input name, as shown in the following example:

```nextflow
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

```nextflow
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

(process-input-path)=

### Input files (`path`)

The `path` qualifier allows you to provide input files to the process execution context. Nextflow will stage the files into the process execution directory, and they can be accessed in the script by using the specified input name. For example:

```nextflow
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

```nextflow
input:
path query_file, name: 'query.fa'
```

or, using a shorter syntax:

```nextflow
input:
path 'query.fa'
```

The previous example can be re-written as shown below:

```nextflow
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

```nextflow
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

:::{versionadded} 23.09.0-edge
:::

By default, `path` inputs will accept any number of files and stage them accordingly. The `arity` option can be used to enforce the expected number of files, either as a number or a range.

For example:

```nextflow
input:
    path('one.txt', arity: '1')         // exactly one file is expected
    path('pair_*.txt', arity: '2')      // exactly two files are expected
    path('many_*.txt', arity: '1..*')   // one or more files are expected
```

When a task is executed, Nextflow will check whether the received files for each path input match the declared arity, and fail if they do not.

:::{note}
Process `path` inputs have nearly the same interface as described in {ref}`stdlib-types-path`, with one difference which is relevant when files are staged into a subdirectory. Given the following input:

```nextflow
path x, name: 'my-dir/*'
```

In this case, `x.name` returns the file name with the parent directory (e.g. `my-dir/file.txt`), whereas normally it would return the file name (e.g. `file.txt`). You can use `x.fileName.name` to get the file name.
:::

### Multiple input files

A `path` input can also accept a collection of files instead of a single value. In this case, the input variable will be a list.

When the input has a fixed file name and a collection of files is received by the process, the file name will be appended with a numerical suffix representing its ordinal position in the list. For example:

```nextflow
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

```nextflow
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

```nextflow
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
In most cases, you won't need to use dynamic file names, because each task is executed in its own directory, and input files are automatically staged into this directory by Nextflow. This behavior guarantees that input files with the same name won't overwrite each other. The above example is useful specifically when there are potential file name conflicts within a single task.
:::

### Input environment variables (`env`)

The `env` qualifier allows you to define an environment variable in the process execution context based on the input value. For example:

```nextflow
process printEnv {
    input:
    env 'HELLO'

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

### Standard input (`stdin`)

The `stdin` qualifier allows you to forward the input value to the [standard input](http://en.wikipedia.org/wiki/Standard_streams#Standard_input_.28stdin.29) of the process script. For example:

```nextflow
process printAll {
  input:
  stdin

  """
  cat -
  """
}

workflow {
  Channel.of('hello', 'hola', 'bonjour', 'ciao')
    | map { v -> v + '\n' }
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

(process-input-tuple)=

### Input tuples (`tuple`)

The `tuple` qualifier allows you to group multiple values into a single input definition. It can be useful when a channel emits tuples of values that need to be handled separately. Each element in the tuple is associated with a corresponding element in the `tuple` definition. For example:

```nextflow
process tupleExample {
    input:
    tuple val(x), path('input.txt')

    """
    echo "Processing $x"
    cat input.txt > copy
    """
}

workflow {
  Channel.of( [1, 'alpha.txt'], [2, 'beta.txt'], [3, 'delta.txt'] ) | tupleExample
}
```

In the above example, the `tuple` input consists of the value `x` and the file `input.txt`.

A `tuple` definition may contain any of the following qualifiers, as previously described: `val`, `env`, `path` and `stdin`. Files specified with the `path` qualifier are treated exactly the same as standalone `path` inputs.

### Input repeaters (`each`)

The `each` qualifier allows you to repeat the execution of a process for each item in a collection, each time a new value is received. For example:

```nextflow
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

```nextflow
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

```nextflow
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

```nextflow
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
- `path`: Emit a file produced by the process with the specified name.
- `env`: Emit the variable defined in the process environment with the specified name.
- `stdout`: Emit the `stdout` of the executed process.
- `tuple`: Emit multiple values.
- `eval`: Emit the result of a script or command evaluated in the task execution context.

Refer to the {ref}`process reference <process-reference-outputs>` for the full list of available output methods and options.

### Output variables (`val`)

The `val` qualifier allows you to output any Nextflow variable defined in the process. A common use case is to output a variable that was defined in the `input` block, as shown in the following example:

```nextflow
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
  receiver.view { method -> "Received: $method" }
}
```

The output value can be a value literal, an input variable, any other Nextflow variable in the process scope, or a value expression. For example:

```nextflow
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

  ch_var.view { var -> "ch_var: $var" }
  ch_str.view { str -> "ch_str: $str" }
  ch_exp.view { exp -> "ch_exp: $exp" }
}
```

### Output files (`path`)

The `path` qualifier allows you to output one or more files produced by the process. For example:

```nextflow
process randomNum {
  output:
  path 'result.txt'

  '''
  echo $RANDOM > result.txt
  '''
}

workflow {
  numbers = randomNum()
  numbers.view { file -> "Received: ${file.text}" }
}
```

In the above example, the `randomNum` process creates a file named `result.txt` which contains a random number. Since a `path` output with the same name is declared, that file is emitted by the corresponding output channel. A downstream process with a compatible input channel will be able to receive it.

Refer to the {ref}`process reference <process-reference-outputs>` for the list of available options for `path` outputs.

:::{versionadded} 23.09.0-edge
:::

By default, `path` outputs will accept any number of matching files from the task directory. The `arity` option can be used to enforce the expected number of files, either as a number or a range.

For example:

```nextflow
output:
path('one.txt', arity: '1')         // exactly one file is expected
path('pair_*.txt', arity: '2')      // exactly two files are expected
path('many_*.txt', arity: '1..*')   // one or more files are expected
```

When a task completes, Nextflow will check whether the produced files for each path output match the declared arity, and fail if they do not.

### Multiple output files

When an output file name contains a `*` or `?` wildcard character, it is interpreted as a [glob][glob] path matcher. This allows you to capture multiple files into a list and emit the list as a single value. For example:

```nextflow
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
        | view { chunk -> "File: ${chunk.name} => ${chunk.text}" }
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

Read more about glob syntax at the following link [What is a glob?][glob]

### Dynamic output file names

When an output file name needs to be expressed dynamically, it is possible to define it using a dynamic string which references variables in the `input` block or in the script global context. For example:

```nextflow
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

### Output environment variables (`env`)

The `env` qualifier allows you to output a variable defined in the process execution environment:

```{literalinclude} snippets/process-out-env.nf
:language: nextflow
```

(process-stdout)=

### Standard output (`stdout`)

The `stdout` qualifier allows you to output the `stdout` of the executed process:

```{literalinclude} snippets/process-stdout.nf
:language: nextflow
```

(process-out-eval)=

### Eval output (`eval`)

:::{versionadded} 24.02.0-edge
:::

The `eval` qualifier allows you to capture the standard output of an arbitrary command evaluated the task shell interpreter context:

```{literalinclude} snippets/process-out-eval.nf
:language: nextflow
```

Only one-line Bash commands are supported. You can use a semi-colon `;` to specify multiple Bash commands on a single line, and many interpreters can execute arbitrary code on the command line, e.g. `python -c 'print("Hello world!")'`.

If the command fails, the task will also fail. In Bash, you can append `|| true` to a command to suppress any command failure.

(process-out-tuple)=

### Output tuples (`tuple`)

The `tuple` qualifier allows you to output multiple values in a single channel. It is useful when you need to associate outputs with metadata, for example:

```nextflow
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
  ch_species = Channel.of('human', 'cow', 'horse')
  ch_query = Channel.fromPath('*.fa')

  blast(ch_species, ch_query)
}
```

In the above example, a `blast` task is executed for each pair of `species` and `query` that are received. Each task produces a new tuple containing the value for `species` and the file `result`.

A `tuple` definition may contain any of the following qualifiers, as previously described: `val`, `path`, `env` and `stdout`. Files specified with the `path` qualifier are treated exactly the same as standalone `path` inputs.

:::{note}
While parentheses for input and output qualifiers are generally optional, they are required when specifying elements in an input/output tuple.

Here's an example with a single path output (parentheses optional):

```nextflow
process foo {
    output:
    path 'result.txt', hidden: true

    '''
    echo 'another new line' >> result.txt
    '''
}
```

And here's an example with a tuple output (parentheses required):

```nextflow
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

(process-naming-outputs)=

### Naming outputs

The `emit` option can be used on a process output to define a name for the corresponding output channel, which can be used to access the channel by name from the process output. For example:

```nextflow
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

### Optional outputs

Normally, if a specified output is not produced by the task, the task will fail. Setting `optional: true` will cause the task to not fail, and instead emit nothing to the given output channel.

```nextflow
output:
path("output.txt"), optional: true
```

In this example, the process is normally expected to produce an `output.txt` file, but in this case, if the file is missing, the task will not fail. The output channel will only contain values for those tasks that produced `output.txt`.

:::{note}
While this option can be used with any process output, it cannot be applied to individual elements of a [tuple](#output-tuples-tuple) output. The entire tuple must be optional or not optional.
:::

(process-when)=

## When

:::{deprecated} 24.10.0
Use conditional logic (e.g. `if` statement, {ref}`operator-filter` operator) in the calling workflow instead.
:::

The `when` block allows you to define a condition that must be satisfied in order to execute the process. The condition can be any expression that returns a boolean value.

It can be useful to enable/disable the process execution depending on the state of various inputs and parameters. For example:

```nextflow
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

(process-directives)=

## Directives

Directives are optional settings that affect the execution of the current process.

By default, directives are evaluated when the process is defined. However, if the value is a dynamic string or closure, it will be evaluated separately for each task, which allows task-specific variables like `task` and `val` inputs to be used.

Some directives are only supported by specific executors. Refer to the {ref}`executor-page` page for more information about each executor.

Refer to the {ref}`process reference <process-reference-directives>` for the full list of process directives. If you are new to Nextflow, here are some commonly-used operators to learn first:

General:
- {ref}`process-error-strategy`: strategy for handling task failures
- {ref}`process-executor`: the {ref}`executor <executor-page>` with which to execute tasks
- {ref}`process-tag`: a semantic name used to differentiate between task executions of the same process

Resource requirements:
- {ref}`process-cpus`: the number of CPUs to request for each task
- {ref}`process-memory`: the amount of memory to request for each task
- {ref}`process-time`: the amount of walltime to request for each task

Software dependencies:
- {ref}`process-conda`: list of conda packages to provision for tasks
- {ref}`process-container`: container image to use for tasks

### Using task directive values

The `task` object also contains the values of all process directives for the given task, which allows you to access these settings at runtime. For examples:

```nextflow
process foo {
  script:
  """
  some_tool --cpus $task.cpus --mem $task.memory
  """
}
```

In the above snippet, `task.cpus` and `task.memory` hold the values for the {ref}`cpus directive<process-cpus>` and {ref}`memory directive<process-memory>` directives, respectively, which were resolved for this task based on the process configuration.

(dynamic-directives)=

### Dynamic directives

A directive can be assigned *dynamically*, during the process execution, so that its actual value can be evaluated based on the process inputs.

To be defined dynamically, the directive's value needs to be expressed using a {ref}`closure <script-closure>`. For example:

```nextflow
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

In the above example, the {ref}`process-queue` directive is evaluated dynamically, depending on the input value `entries`. When it is larger than 100, jobs will be submitted to the `long` queue, otherwise the `short` queue will be used.

All directives can be assigned a dynamic value except the following:

- {ref}`process-executor`
- {ref}`process-label`
- {ref}`process-maxforks`

:::{tip}
Assigning a string value with one or more variables is always resolved in a dynamic manner, and therefore is equivalent to the above syntax. For example, the above directive can also be written as:

```nextflow
queue "${ entries > 100 ? 'long' : 'short' }"
```

Note, however, that the latter syntax can be used both for a directive's main argument (as in the above example) and for a directive's optional named attributes, whereas the closure syntax is only resolved dynamically for a directive's main argument.
:::

(dynamic-task-resources)=

### Dynamic task resources

It's a very common scenario that different instances of the same process may have very different needs in terms of computing resources. In such situations requesting, for example, an amount of memory too low will cause some tasks to fail. Instead, using a higher limit that fits all the tasks in your execution could significantly decrease the execution priority of your jobs.

The [Dynamic directives](#dynamic-directives) evaluation feature can be used to modify the amount of computing resources requested in case of a process failure and try to re-execute it using a higher limit. For example:

```nextflow
process foo {
    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    script:
    <your job here>
}
```

In the above example the {ref}`process-memory` and execution {ref}`process-time` limits are defined dynamically. The first time the process is executed the `task.attempt` is set to `1`, thus it will request a two GB of memory and one hour of maximum execution time.

If the task execution fail reporting an exit status in the range between 137 and 140, the task is re-submitted (otherwise terminates immediately). This time the value of `task.attempt` is `2`, thus increasing the amount of the memory to four GB and the time to 2 hours, and so on.

The directive {ref}`process-maxretries` set the maximum number of time the same task can be re-executed.

### Dynamic task resources with previous execution trace
:::{versionadded} 24.10.0
:::

Task resource requests can be updated relative to the {ref}`trace record <trace-report>` metrics of the previous task attempt. The metrics can be accessed through the `task.previousTrace` variable. For example:

```nextflow
process foo {
    memory { task.attempt > 1 ? task.previousTrace.memory * 2 : (1.GB) }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    script:
    <your job here>
}
```
In the above example, the {ref}`process-memory` is set according to previous trace record metrics. In the first attempt, when no trace metrics are available, it is set to one GB. In the subsequent attempts, it doubles the previously allocated memory. See {ref}`trace-report` for more information about trace records.


### Dynamic retry with backoff

There are cases in which the required execution resources may be temporary unavailable e.g. network congestion. In these cases immediately re-executing the task will likely result in the identical error. A retry with an exponential backoff delay can better recover these error conditions:

```nextflow
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

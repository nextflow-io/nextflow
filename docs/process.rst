.. _process-page:

*********
Processes
*********

In Nextflow a `process` is the basic processing `primitive` to execute a user script.

The process definition starts with keyword the ``process``, followed by process name and finally the process `body`
delimited by curly brackets. The process body must contain a string which represents the command or, more generally,
a script that is executed by it. A basic process looks like the following example::

  process sayHello {

      """
      echo 'Hello world!' > file
      """

  }


A process may contain five definition blocks, respectively: directives,
inputs, outputs, when clause and finally the process script. The syntax is defined as follows:

::

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


.. _process-script:

Script
======

The `script` block is a string statement that defines the command that is executed by the process to carry out its task.

A process contains one and only one script block, and it must be the last statement when the process contains
input and output declarations.

The entered string is executed as a `Bash <http://en.wikipedia.org/wiki/Bash_(Unix_shell)>`_ script in the
`host` system. It can be any command, script or combination of them, that you would normally use in terminal shell
or in a common Bash script.

The only limitation to the commands that can be used in the script statement is given by the availability of those
programs in the target execution system.


The script block can be a simple string or multi-line string. The latter simplifies the writing of non trivial scripts
composed by multiple commands spanning over multiple lines. For example::

    process doMoreThings {

      """
      blastp -db $db -query query.fa -outfmt 6 > blast_result
      cat blast_result | head -n 10 | cut -f 2 > top_hits
      blastdbcmd -db $db -entry_batch top_hits > sequences
      """

    }

As explained in the script tutorial section, strings can be defined by using a single-quote
or a double-quote, and multi-line strings are defined by three single-quote or three double-quote characters.

There is a subtle but important difference between them. Like in Bash, strings delimited by a ``"`` character support
variable substitutions, while strings delimited by ``'`` do not.

In the above code fragment the ``$db`` variable is replaced by the actual value defined somewhere in the
pipeline script.

.. warning:: Since Nextflow uses the same Bash syntax for variable substitutions in strings, you need to manage them
  carefully depending on if you want to evaluate a variable in the Nextflow context - or - in the Bash environment execution.

When you need to access a system environment variable  in your script you have two options. The first choice is as
easy as defining your script block by using a single-quote string. For example::

    process printPath {

       '''
       echo The path is: $PATH
       '''

    }

The drawback of this solution is that you will not able to access variables defined in the pipeline script context,
in your script block.

To fix this, define your script by using a double-quote string and `escape` the system environment variables by
prefixing them with a back-slash ``\`` character, as shown in the following example::


    process doOtherThings {

      """
      blastp -db \$DB -query query.fa -outfmt 6 > blast_result
      cat blast_result | head -n $MAX | cut -f 2 > top_hits
      blastdbcmd -db \$DB -entry_batch top_hits > sequences
      """

    }

In this example the ``$MAX`` variable has to be defined somewhere before, in the pipeline script.
`Nextflow` replaces it with the actual value before executing the script. Instead, the ``$DB`` variable
must exist in the script execution environment and the Bash interpreter will replace it with the actual value.

.. tip::
  Alternatively you can use the :ref:`process-shell` block definition which allows a script to contain both
  Bash and Nextflow variables without having to escape the first.

Scripts `à la carte`
--------------------

The process script is interpreted by Nextflow as a Bash script by default, but you are not limited to it.

You can use your favourite scripting language (e.g. Perl, Python, Ruby, R, etc), or even mix them in the same pipeline.

A pipeline may be composed by processes that execute very different tasks. Using `Nextflow` you can choose the scripting
language that better fits the task carried out by a specified process. For example for some processes `R` could be
more useful than `Perl`, in other you may need to use `Python` because it provides better access to a library or an API, etc.

To use a scripting other than Bash, simply start your process script with the corresponding
`shebang <http://en.wikipedia.org/wiki/Shebang_(Unix)>`_ declaration. For example::

    process perlStuff {

        """
        #!/usr/bin/perl

        print 'Hi there!' . '\n';
        """

    }

    process pyStuff {

        """
        #!/usr/bin/python

        x = 'Hello'
        y = 'world!'
        print "%s - %s" % (x,y)
        """

    }


.. tip:: Since the actual location of the interpreter binary file can change across platforms, to make your scripts
   more portable it is wise to use the ``env`` shell command followed by the interpreter's name, instead of the absolute
   path of it. Thus, the `shebang` declaration for a Perl script, for example,
   would look like: ``#!/usr/bin/env perl`` instead of the one in the above pipeline fragment.


Conditional scripts
-------------------

Complex process scripts may need to evaluate conditions on the input parameters or use traditional flow control
statements (i.e. ``if``, ``switch``, etc) in order to execute specific script commands, depending on the current
inputs configuration.

Process scripts can contain conditional statements by simply prefixing the script block with the keyword ``script:``.
By doing that the interpreter will evaluate all the following statements as a code block that must return the
script string to be executed. It's much easier to use than to explain, for example::


    seq_to_align = ...
    mode = 'tcoffee'

    process align {
        input:
        file seq_to_aln from sequences

        script:
        if( mode == 'tcoffee' )
            """
            t_coffee -in $seq_to_aln > out_file
            """

        else if( mode == 'mafft' )
            """
            mafft --anysymbol --parttree --quiet $seq_to_aln > out_file
            """

        else if( mode == 'clustalo' )
            """
            clustalo -i $seq_to_aln -o out_file
            """

        else
            error "Invalid alignment mode: ${mode}"

    }


In the above example the process will execute the script fragment depending on the value of the ``mode`` parameter.
By default it will execute the ``tcoffee`` command, changing the ``mode`` variable to ``mafft`` or ``clustalo`` value,
the other branches will be executed.

.. _process-template:

Template
--------

Process script can be externalised by using *template* files which can be reused across different processes and tested
independently by the overall pipeline execution.

A template is simply a shell script file that Nextflow is able to execute by using the ``template`` function
as shown below::

    process template_example {

        input:
        val STR from 'this', 'that'

        script:
        template 'my_script.sh'

    }


Nextflow looks for the ``my_script.sh`` template file in the directory ``templates`` that must exist in the same folder
where the Nextflow script file is located (any other location can be provided by using an absolute template path).

The template script can contain any piece of code that can be executed by the underlying system. For example::

  #!/bin/bash
  echo "process started at `date`"
  echo $STR
  :
  echo "process completed"



.. tip::
  Note that the dollar character (``$``) is interpreted as a Nextflow variable placeholder, when the script is run as a
  Nextflow template, while it is evaluated as a Bash variable when it is run alone. This can be very useful to test
  your script autonomously, i.e. independently from Nextflow execution. You only need to provide a Bash environment
  variable for each the Nextflow variable existing in your script. For example, it would be possible to execute the above
  script entering the following command in the shell terminal: ``STR='foo' bash templates/my_script.sh``


.. _process-shell:

Shell
-----

.. warning:: This is an incubating feature. It may change in future Nextflow releases.


The ``shell`` block is a string statement that defines the *shell* command executed by the process to carry out its task.
It is an alternative to the :ref:`process-script` definition with an important difference, it uses
the exclamation mark ``!`` character as the variable placeholder for Nextflow variables in place of the usual dollar character.

In this way it is possible to use both Nextflow and Bash variables in the same piece of code without having to escape
the latter and making process scripts more readable and easy to maintain. For example::

    process myTask {

        input:
        val str from 'Hello', 'Hola', 'Bonjour'

        shell:
        '''
        echo User $USER says !{str}
        '''

    }



In the above trivial example the ``$USER`` variable is managed by the Bash interpreter, while ``!{str}`` is handled
as a process input variable managed by Nextflow.

.. note::

    - Shell script definition requires the use of single-quote ``'`` delimited strings. When using double-quote ``"``
      delimited strings, dollar variables are interpreted as Nextflow variables as usual. See :ref:`string-interpolation`.

    - Exclamation mark prefixed variables always need to be enclosed in curly brackets i.e. ``!{str}`` is a valid 
      variable while ``!str`` is ignored.

    - Shell script supports the use of the file :ref:`process-template` mechanism. The same rules are applied to the variables
      defined in the script template.

.. _process-native:

Native execution
----------------

Nextflow processes can execute native code other than system scripts as shown in the previous paragraphs.

This means that instead of specifying the process command to be executed as a string script, you can
define it by providing one or more language statements, as you would do in the rest of the pipeline script.
Simply starting the script definition block with the ``exec:`` keyword, for example::

    x = Channel.from( 'a', 'b', 'c')

    process simpleSum {
        input:
        val x

        exec:
        println "Hello Mr. $x"
    }

Will display::

    Hello Mr. b
    Hello Mr. a
    Hello Mr. c



.. _process-input:

Inputs
======

Nextflow processes are isolated from each other but can communicate between themselves sending values through channels.

The `input` block defines which channels the process is expecting to receive inputs data from. You can only define one
input block at a time and it must contain one or more inputs declarations.

The input block follows the syntax shown below::

    input:
      <input qualifier> <input name> [from <source channel>] [attributes]


An input definition starts with an input `qualifier` and the input `name`, followed by the keyword ``from`` and
the actual channel over which inputs are received. Finally some input optional attributes can be specified.

.. note:: When the input name is the same as the channel name, the ``from`` part of the declaration can be omitted.

The input qualifier declares the `type` of data to be received. This information is used by Nextflow to apply the
semantic rules associated to each qualifier and handle it properly depending on the target execution platform
(grid, cloud, etc).

The qualifiers available are the ones listed in the following table:

=========== =============
Qualifier   Semantic
=========== =============
val         Lets you access the received input value by its name in the process script.
env         Lets you use the received value to set an environment variable named
            as the specified input name.
file        Lets you handle the received value as a file, staging it properly in the execution context.
stdin       Lets you forward the received value to the process `stdin` special file.
set         Lets you handle a group of input values having one of the above qualifiers.
each        Lets you execute the process for each entry in the input collection.
=========== =============


Input of generic values
-----------------------

The ``val`` qualifier allows you to receive data of any type as input. It can be accessed in the process script
by using the specified input name, as shown in the following example::

    num = Channel.from( 1, 2, 3 )

    process basicExample {
      input:
      val x from num

      "echo process job $x"

    }


In the above example the process is executed three times, each time a value is received from the channel ``num``
and used to process the script. Thus, it results in an output similar to the one shown below::

    process job 3
    process job 1
    process job 2

.. note:: The `channel` guarantees that items are delivered in the same order as they have been sent - but -
  since the process is executed in a parallel manner, there is no guarantee that they are processed in the
  same order as they are received. In fact, in the above example, value ``3`` is processed before the others.


When the ``val`` has the same name as the channel from where the data is received, the ``from`` part can be omitted.
Thus the above example can be written as shown below::

    num = Channel.from( 1, 2, 3 )

    process basicExample {
      input:
      val num

      "echo process job $num"

    }


Input of files
--------------

The ``file`` qualifier allows the handling of file values in the process execution context. This means that
Nextflow will stage it in the process execution directory, and it can be access in the script by using the name
specified in the input declaration. For example::

    proteins = Channel.fromPath( '/some/path/*.fa' )

    process blastThemAll {
      input:
      file query_file from proteins

      "blastp -query ${query_file} -db nr"

    }

In the above example all the files ending with the suffix ``.fa`` are sent over the channel ``proteins``.
Then, these files are received by the process which will execute a `BLAST` query on each of them.

When the file input name is the same as the channel name, the ``from`` part of the input declaration can be omitted.
Thus, the above example could be written as shown below::

    proteins = Channel.fromPath( '/some/path/*.fa' )

    process blastThemAll {
      input:
      file proteins

      "blastp -query $proteins -db nr"

    }


It's worth noting that in the above examples, the name of the file in the file-system is not touched, you can
access the file even without knowing its name because you can reference it in the process script by using the
variable whose name is specified in the input file parameter declaration.

There may be cases where your task needs to use a file whose name is fixed, it does not have to change along
with the actual provided file. In this case you can specify its name by specifying the ``name`` attribute in the
input file parameter declaration, as shown in the following example::

    input:
        file query_file name 'query.fa' from proteins


Or alternatively using a shorter syntax::

    input:
        file 'query.fa' from proteins


Using this, the previous example can be re-written as shown below::

    proteins = Channel.fromPath( '/some/path/*.fa' )

    process blastThemAll {
      input:
      file 'query.fa' from proteins

      "blastp -query query.fa -db nr"

    }


What happens in this example is that each file, that the process receives, is staged with the name ``query.fa``
in a different execution context (i.e. the folder where the job is executed) and an independent process
execution is launched.

.. tip:: This allows you to execute the process command various time without worrying the files names changing.
  In other words, `Nextflow` helps you write pipeline tasks that are self-contained and decoupled by the execution
  environment. This is also the reason why you should avoid whenever possible to use absolute or relative paths
  referencing files in your pipeline processes.


.. TODO describe that file can handle channels containing any data type not only file


Multiple input files
--------------------

A process can declare as input file a channel that emits a collection of values, instead of a simple value.

In this case, the script variable defined by the input file parameter will hold a list of files. You can
use it as shown before, referring to all the files in the list, or by accessing a specific entry using the
usual square brackets notation.

When a target file name is defined in the input parameter and a collection of files is received by the process,
the file name will be appended by a numerical suffix representing its ordinal position in the list. For example::

    fasta = Channel.fromPath( "/some/path/*.fa" ).buffer(size:3)

    process blastThemAll {
        input:
        file 'seq' from fasta

        "echo seq*"

    }

Will output::

    seq1 seq2 seq3
    seq1 seq2 seq3
    ...

The target input file name can contain the ``*`` and ``?`` wildcards, that can be used
to control the name of staged files. The following table shows how the wildcards are
replaced depending on the cardinality of the received input collection.

============ ============== ==================================================
Cardinality   Name pattern     Staged file names
============ ============== ==================================================
 any         ``*``           named as the source file
 1           ``file*.ext``   ``file.ext``
 1           ``file?.ext``   ``file1.ext``
 1           ``file??.ext``  ``file01.ext``
 many        ``file*.ext``   ``file1.ext``, ``file2.ext``, ``file3.ext``, ..
 many        ``file?.ext``   ``file1.ext``, ``file2.ext``, ``file3.ext``, ..
 many        ``file??.ext``  ``file01.ext``, ``file02.ext``, ``file03.ext``, ..
 many        ``dir/*``       named as the source file, created in ``dir`` subdirectory
 many        ``dir??/*``     named as the source file, created in a progressively indexed subdirectory e.g. ``dir01/``, ``dir02/``, etc.
 many        ``dir*/*``      (as above)
============ ============== ==================================================

The following fragment shows how a wildcard can be used in the input file declaration::


    fasta = Channel.fromPath( "/some/path/*.fa" ).buffer(size:3)

    process blastThemAll {
        input:
        file 'seq?.fa' from fasta

        "cat seq1.fa seq2.fa seq3.fa"

    }


.. note:: Rewriting input file names according to a named pattern is an extra feature and not at all obligatory.
  The normal file input constructs introduced in the `Input of files`_ section are valid for collections of
  multiple files as well. To handle multiple input files preserving the original file names, use the ``*`` wildcard as
  name pattern or a variable identifier.

Dynamic input file names
------------------------

When the input file name is specified by using the ``name`` file clause or the short `string` notation, you
are allowed to use other input values as variables in the file name string. For example::


  process simpleCount {
    input:
    val x from species
    file "${x}.fa" from genomes

    """
    cat ${x}.fa | grep '>'
    """
  }


In the above example, the input file name is set by using the current value of the ``x`` input value.

This allows the input files to be staged in the script working directory with a name that is coherent
with the current execution context.

.. tip:: In most cases, you won't need to use dynamic file names, because each process is executed in its 
  own private temporary directory, and input files are automatically staged to this directory by Nextflow. 
  This guarantees that input files with the same name won't overwrite each other.



Input of type 'stdin'
---------------------

The ``stdin`` input qualifier allows you the forwarding of the value received from a channel to the
`standard input <http://en.wikipedia.org/wiki/Standard_streams#Standard_input_.28stdin.29>`_
of the command executed by the process. For example::

    str = Channel.from('hello', 'hola', 'bonjour', 'ciao').map { it+'\n' }

    process printAll {
       input:
       stdin str

       """
       cat -
       """

    }

It will output::

    hola
    bonjour
    ciao
    hello




Input of type 'env'
-------------------

The ``env`` qualifier allows you to define an environment variable in the process execution context based
on the value received from the channel. For example::

    str = Channel.from('hello', 'hola', 'bonjour', 'ciao')

    process printEnv {

        input:
        env HELLO from str

        '''
        echo $HELLO world!
        '''

    }

::

    hello world!
    ciao world!
    bonjour world!
    hola world!



Input of type 'set'
-------------------

The ``set`` qualifier allows you to group multiple parameters in a single parameter definition. It can be useful
when a process receives, in input, tuples of values that need to be handled separately. Each element in the tuple
is associated to a corresponding element with the ``set`` definition. For example::

     tuple = Channel.from( [1, 'alpha'], [2, 'beta'], [3, 'delta'] )

     process setExample {
         input:
         set val(x), file('latin.txt')  from tuple

         """
         echo Processing $x
         cat - latin.txt > copy
         """

     }


In the above example the ``set`` parameter is used to define the value ``x`` and the file ``latin.txt``,
which will receive a value from the same channel.

In the ``set`` declaration items can be defined by using the following qualifiers: ``val``, ``env``, ``file`` and ``stdin``.

A shorter notation can be used by applying the following substitution rules:

============== =======
long            short
============== =======
val(x)          x
file(x)         (not supported)
file('name')    'name'
file(x:'name')  x:'name'
stdin           '-'
env(x)          (not supported)
============== =======

Thus the previous example could be rewritten as follows::

      tuple = Channel.from( [1, 'alpha'], [2, 'beta'], [3, 'delta'] )

      process setExample {
          input:
          set x, 'latin.txt' from tuple

          """
          echo Processing $x
          cat - latin.txt > copy
          """

      }

File names can be defined in *dynamic* manner as explained in the `Dynamic input file names`_ section.


Input repeaters
---------------

The ``each`` qualifier allows you to repeat the execution of a process for each item in a collection,
every time a new data is received. For example::

  sequences = Channel.fromPath('*.fa')
  methods = ['regular', 'expresso', 'psicoffee']

  process alignSequences {
    input:
    file seq from sequences
    each mode from methods

    """
    t_coffee -in $seq -mode $mode > result
    """
  }


In the above example every time a file of sequences is received as input by the process,
it executes *three* tasks running a T-coffee alignment with a different value for the ``mode`` parameter.
This is useful when you need to `repeat` the same task for a given set of parameters.

Since version 0.25+ input repeaters can be applied to files as well. For example::

    sequences = Channel.fromPath('*.fa')
    methods = ['regular', 'expresso']
    libraries = [ file('PQ001.lib'), file('PQ002.lib'), file('PQ003.lib') ]

    process alignSequences {
      input:
      file seq from sequences
      each mode from methods
      each file(lib) from libraries

      """
      t_coffee -in $seq -mode $mode -lib $lib > result
      """
    }


.. note:: When multiple repeaters are declared, the process is executed for each *combination* of them.

In the latter example for any sequence input file emitted by the ``sequences`` channel are executed 6 alignments,
3 using the ``regular`` method against each library files, and other 3 by using the ``expresso`` method always
against the same library files.


.. hint:: If you need to repeat the execution of a process over n-tuple of elements instead a simple values or files,
  create a channel combining the input values as needed to trigger the process execution multiple times.
  In this regard, see the :ref:`operator-combine`, :ref:`operator-cross` and :ref:`operator-phase` operators.

.. _process-understand-how-multiple-input-channels-work:

Understand how multiple input channels work
-------------------------------------------

A key feature of processes is the ability to handle inputs from multiple channels.

When two or more channels are declared as process inputs, the process stops until
there's a complete input configuration ie. it receives an input value from all the channels declared
as input.

When this condition is verified, it consumes the input values coming from the respective channels,
and spawns a task execution, then repeat the same logic until one or more channels have no more content.

This means channel values are consumed serially one after another and the first empty channel
cause the process execution to stop even if there are other values in other channels.

For example::

  process foo {
    echo true
    input:
    val x from Channel.from(1,2)
    val y from Channel.from('a','b','c')
    script:
     """
     echo $x and $y
     """
  }


The process ``foo`` is executed two times because the first input channel only provides two values and therefore
the ``c`` element is discarded. It prints::

    1 and a
    2 and b


.. warning:: A different semantic is a applied when using *Value channel* a.k.a. *Singleton channel*.

This kind of channel is created by the :ref:`Channel.value <channel-value>` factory method or implicitly
when a process input specifies a simple value in the ``from`` clause.

By definition, a *Value channel* is bound to a single value and it can be read unlimited times without
consuming its content.

These properties make that when mixing a *value channel* with one or more (queue) channels,
it does not affect the process termination which only depends by the other channels and its
content is applied repeatedly.

To better understand this behavior compare the previous example with the following one::

  process bar {
    echo true
    input:
    val x from Channel.value(1)
    val y from Channel.from('a','b','c')
    script:
     """
     echo $x and $y
     """
  }

The above snippet executes the ``bar`` process three times because the first input is a *value channel*, therefore
its content can be read as many times as needed. The process termination is determined by the content of the second
channel. It prints::


  1 and a
  1 and b
  1 and c

See also: :ref:`channel-types`.

Outputs
=======

The `output` declaration block allows to define the channels used by the process to send out the results produced.

It can be defined at most one output block and it can contain one or more outputs declarations.
The output block follows the syntax shown below::

    output:
      <output qualifier> <output name> [into <target channel>[,channel,..]] [attribute [,..]]

Output definitions start by an output `qualifier` and the output `name`, followed by the keyword ``into`` and
one or more channels over which outputs are sent. Finally some optional attributes can be specified.

.. note:: When the output name is the same as the channel name, the ``into`` part of the declaration can be omitted.


.. TODO the channel is implicitly create if does not exist

The qualifiers that can be used in the output declaration block are the ones listed in the following table:

=========== =============
Qualifier   Semantic
=========== =============
val         Sends variable's with the name specified over the output channel.
file        Sends a file produced by the process with the name specified over the output channel.
stdout      Sends the executed process `stdout` over the output channel.
set         Lets to send multiple values over the same output channel.
=========== =============


Output values
-------------

The ``val`` qualifier allows to output a `value` defined in the script context. In a common usage scenario,
this is a value which has been defined in the `input` declaration block, as shown in the following example::

   methods = ['prot','dna', 'rna']

   process foo {
     input:
     val x from methods

     output:
     val x into receiver

     """
     echo $x > file
     """

   }

   receiver.println { "Received: $it" }


Valid output values are value literals, input values identifiers, variables accessible in the process scope and
value expressions. For example::

    process foo {
      input:
      file fasta from 'dummy'

      output:
      val x into var_channel
      val 'BB11' into str_channel
      val "${fasta.baseName}.out" into exp_channel

      script:
      x = fasta.name
      """
      cat $x > file
      """
    }




Output files
------------

The ``file`` qualifier allows to output one or more files, produced by the process, over the specified channel.
For example::


    process randomNum {

       output:
       file 'result.txt' into numbers

       '''
       echo $RANDOM > result.txt
       '''

    }

    numbers.subscribe { println "Received: " + it.text }


In the above example the process, when executed, creates a file named ``result.txt`` containing a random number.
Since a file parameter using the same name is declared between the outputs, when the task is completed that
file is sent over the ``numbers`` channel. A downstream `process` declaring the same channel as `input` will
be able to receive it.

.. note:: If the channel specified as output has not been previously declared in the pipeline script, it
  will implicitly created by the output declaration itself.


.. TODO explain Path object

Multiple output files
---------------------

When an output file name contains a ``*`` or ``?`` wildcard character it is interpreted as a `glob`_ path matcher.
This allows to *capture* multiple files into a list object and output them as a sole emission. For example::

    process splitLetters {

        output:
        file 'chunk_*' into letters

        '''
        printf 'Hola' | split -b 1 - chunk_
        '''
    }

    letters
        .flatMap()
        .subscribe { println "File: ${it.name} => ${it.text}" }

It prints::

    File: chunk_aa => H
    File: chunk_ab => o
    File: chunk_ac => l
    File: chunk_ad => a

.. note:: In the above example the operator :ref:`operator-flatmap` is used to transform the list of files emitted by
  the ``letters`` channel into a channel that emits each file object independently.

Some caveats on glob pattern behavior:

* Input files are not included in the list of possible matches.
* Glob pattern matches against both files and directories path.
* When a two stars pattern ``**`` is used to recourse across directories, only file paths are matched
  i.e. directories are not included in the result list.


.. tip::
    By default all the files matching the specified glob pattern are emitted by the channel as a sole (list) item.
    It is also possible to emit each file as a sole item by adding the ``mode flatten`` attribute in the output file
    declaration.

By using the `mode` attribute the previous example can be re-written as show below::

    process splitLetters {

        output:
        file 'chunk_*' into letters mode flatten

        '''
        printf 'Hola' | split -b 1 - chunk_
        '''
    }

    letters .subscribe { println "File: ${it.name} => ${it.text}" }



Read more about glob syntax at the following link `What is a glob?`_

.. _glob: http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
.. _What is a glob?: http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob

.. _process-dynoutname:

Dynamic output file names
-------------------------

When an output file name needs to be expressed dynamically, it is possible to define it using a dynamic evaluated
string which references values defined in the input declaration block or in the script global context.
For example::


  process align {
    input:
    val x from species
    file seq from sequences

    output:
    file "${x}.aln" into genomes

    """
    t_coffee -in $seq > ${x}.aln
    """
  }

In the above example, each time the process is executed an alignment file is produced whose name depends
on the actual value of the ``x`` input.

.. tip:: The management of output files is a very common misunderstanding when using Nextflow. 
  With other tools, it is generally necessary to organize the outputs files into some kind of directory 
  structure or to guarantee a unique file name scheme, so that result files won't overwrite each other 
  and that they can be referenced univocally by downstream tasks.

  With Nextflow, in most cases, you don't need to take care of naming output files, because each task is executed 
  in its own unique temporary directory, so files produced by different tasks can never override each other.
  Also meta-data can be associated with outputs by using the :ref:`set output <process-set>` qualifier, instead of
  including them in the output file name.

  To sum up, the use of output files with static names over dynamic ones is preferable whenever possible, 
  because it will result in a simpler and more portable code.


.. _process-stdout:

Output 'stdout' special file
----------------------------

The ``stdout`` qualifier allows to `capture` the `stdout` output of the executed process and send it over
the channel specified in the output parameter declaration. For example::

    process echoSomething {
        output:
        stdout channel

        "echo Hello world!"

    }

    channel.subscribe { print "I say..  $it" }



.. _process-set:

Output 'set' of values
----------------------

The ``set`` qualifier allows to send multiple values into a single channel. This feature is useful
when you need to `group together` the results of multiple executions of the same process, as shown in the following
example::

    query = Channel.fromPath '*.fa'
    species = Channel.from 'human', 'cow', 'horse'

    process blast {

    input:
        val species
        file query

    output:
        set val(species), file('result') into blastOuts


    "blast -db nr -query $query" > result

    }


In the above example a `BLAST` task is executed for each pair of ``species`` and ``query`` that are received.
When the task completes a new tuple containing the value for ``species`` and the file ``result`` is sent to the ``blastOuts`` channel.


A `set` declaration can contain any combination of the following qualifiers, previously described: ``val``, ``file`` and ``stdout``.

.. tip:: Variable identifiers are interpreted as `values` while strings literals are interpreted as `files` by default,
  thus the above output `set` can be rewritten using a short notation as shown below.

::

    output:
        set species, 'result' into blastOuts



File names can be defined in a dynamic manner as explained in the :ref:`process-dynoutname` section.

Optional Output
---------------

In most cases a process is expected to generate output that is added to the output channel.  However, there are situations where it is valid for a process to `not` generate output. In these cases ``optional true`` may be added to the output declaration, which tells Nextflow not to fail the process if the declared output is not created.

::

    output:
        file("output.txt") optional true into outChannel

In this example, the process is normally expected to generate an ``output.txt`` file, but in the cases where the file is legitimately missing, the process does not fail. ``outChannel`` is only populated by those processes that do generate ``output.txt``. 


When
====

The ``when`` declaration allows you to define a condition that must be verified in order to execute the process.
This can be any expression that evaluates a boolean value.

It is useful to enable/disable the process execution depending the state of various inputs and parameters. For example::


    process find {
      input:
      file proteins
      val type from dbtype

      when:
      proteins.name =~ /^BB11.*/ && type == 'nr'

      script:
      """
      blastp -query $proteins -db nr
      """

    }


.. _process-directives:

Directives
==========

Using the `directive` declarations block you can provide optional settings that will affect the execution of the current
process.

They must be entered at the top of the process `body`, before any other declaration blocks (i.e. ``input``, ``output``, etc) 
and have the following syntax::

    name value [, value2 [,..]]

Some directives are generally available to all processes, some others depends on the `executor` currently defined.

The directives are:

* `afterScript`_
* `beforeScript`_
* `cache`_
* `cpus`_
* `conda`_
* `container`_
* `containerOptions`_
* `clusterOptions`_
* `disk`_
* `echo`_
* `errorStrategy`_
* `executor`_
* `ext`_
* `label`_
* `maxErrors`_
* `maxForks`_
* `maxRetries`_
* `memory`_
* `module`_
* `penv`_
* `pod`_
* `publishDir`_
* `queue`_
* `scratch`_
* `stageInMode`_
* `stageOutMode`_
* `storeDir`_
* `tag`_
* `time`_
* `validExitStatus`_

afterScript
-----------

The ``afterScript`` directive allows you to execute a custom (Bash) snippet immediately *after* the main process has run.
This may be useful to clean up your staging area.

beforeScript
------------

The ``beforeScript`` directive allows you to execute a custom (Bash) snippet *before* the main process script is run.
This may be useful to initialise the underlying cluster environment or for other custom initialisation.

For example::

    process foo {

      beforeScript 'source /cluster/bin/setup'

      """
      echo bar
      """

    }


cache
-----

The ``cache`` directive allows you to store the process results to a local cache. When the cache is enabled *and*
the pipeline is launched with the :ref:`resume <getstart-resume>` option, any following attempt to execute the process,
along with the same inputs, will cause the process execution to be skipped, producing the stored data as
the actual results.

The caching feature generates a unique `key` by indexing the process script and inputs. This key is used
identify univocally the outputs produced by the process execution.


The cache is enabled by default, you can disable it for a specific process by setting the ``cache``
directive to ``false``. For example:: 

  process noCacheThis {
    cache false

    script:
    <your command string here>
  }

The ``cache`` directive possible values are shown in the following table:

===================== =================
Value                 Description
===================== =================
``false``             Disable cache feature.
``true`` (default)    Enable caching. Cache keys are created indexing input files meta-data information (name, size and last update timestamp attributes).
``'deep'``            Enable caching. Cache keys are created indexing input files content.
``'lenient'``         Enable caching. Cache keys are created indexing input files path and size attributes (this policy provides a workaround for incorrect caching invalidation observed on shared file systems due to inconsistent files timestamps; requires version 0.32.x or later).
===================== =================


.. _process-conda:

conda
-----

The ``conda`` directive allows for the definition of the process dependencies using the `Conda <https://conda.io>`_
package manager.

Nextflow automatically sets up an environment for the given package names listed by in the ``conda`` directive.
For example::

  process foo {
    conda 'bwa=0.7.15'

    '''
    your_command --here
    '''
  }


Multiple packages can be specified separating them with a blank space eg. ``bwa=0.7.15 fastqc=0.11.5``.
The name of the channel from where a specific package needs to be downloaded can be specified using the usual
Conda notation i.e. prefixing the package with the channel name as shown here ``bioconda::bwa=0.7.15``.

The ``conda`` directory also allows the specification of a Conda environment file
path or the path of an existing environment directory. See the :ref:`conda-page` page for further details.

.. _process-container:

container
---------

The ``container`` directive allows you to execute the process script in a `Docker <http://docker.io>`_ container.

It requires the Docker daemon to be running in machine where the pipeline is executed, i.e. the local machine when using the
*local* executor or the cluster nodes when the pipeline is deployed through a *grid* executor.

For example::


    process runThisInDocker {

      container 'dockerbox:tag'

      """
      <your holy script here>
      """

    }


Simply replace in the above script ``dockerbox:tag`` with the Docker image name you want to use.

.. tip:: This can be very useful to execute your scripts into a replicable self-contained environment or to deploy your pipeline in the cloud.

.. note:: This directive is ignore for processes :ref:`executed natively <process-native>`.


.. _process-containerOptions:

containerOptions
----------------

The ``containerOptions`` directive allows you to specify any container execution option supported by the underlying
container engine (ie. Docker, Singularity, etc). This can be useful to provide container settings
only for a specific process e.g. mount a custom path::


  process runThisWithDocker {

      container 'busybox:latest'
      containerOptions '--volume /data/db:/db'

      output: file 'output.txt'

      '''
      your_command --data /db > output.txt
      '''
  }


.. warning:: This feature is not supported by :ref:`awsbatch-executor` and :ref:`k8s-executor` executors.

.. _process-cpus:

cpus
----

The ``cpus`` directive allows you to define the number of (logical) CPU required by the process' task.
For example::

    process big_job {

      cpus 8
      executor 'sge'

      """
      blastp -query input_sequence -num_threads ${task.cpus}
      """
    }


This directive is required for tasks that execute multi-process or multi-threaded commands/tools and it is meant
to reserve enough CPUs when a pipeline task is executed through a cluster resource manager.

See also: `penv`_, `memory`_, `time`_, `queue`_, `maxForks`_

.. _process-clusterOptions:

clusterOptions
--------------

The ``clusterOptions`` directive allows the usage of any `native` configuration option accepted by your cluster submit command.
You can use it to request non-standard resources or use settings that are specific to your cluster and not supported
out of the box by Nextflow.

.. note:: This directive is taken in account only when using a grid based executor:
  :ref:`sge-executor`, :ref:`lsf-executor`, :ref:`slurm-executor`, :ref:`pbs-executor` and
  :ref:`condor-executor` executors.

.. _process-disk:

disk
----

The ``disk`` directive allows you to define how much local disk storage the process is allowed to use. For example::

    process big_job {

        disk '2 GB'
        executor 'cirrus'

        """
        your task script here
        """
    }

The following memory unit suffix can be used when specifying the disk value:

======= =============
Unit    Description
======= =============
B       Bytes
KB      Kilobytes
MB      Megabytes
GB      Gigabytes
TB      Terabytes
======= =============

.. note:: This directive currently is taken in account only by the :ref:`ignite-executor`
  and the :ref:`condor-executor` executors.

See also: `cpus`_, `memory`_ `time`_, `queue`_ and `Dynamic computing resources`_.

.. _process-echo:

echo
----

By default the `stdout` produced by the commands executed in all processes is ignored.
Setting the ``echo`` directive to ``true`` you can forward the process `stdout` to the current top
running process `stdout` file, showing it in the shell terminal.

For example::

    process sayHello {
      echo true

      script:
      "echo Hello"
    }

::

    Hello

Without specifying ``echo true`` you won't see the ``Hello`` string printed out when executing the above example.


.. _process-page-error-strategy:

errorStrategy
-------------

The ``errorStrategy`` directive allows you to define how an error condition is managed by the process. By default when
an error status is returned by the executed script, the process stops immediately. This in turn forces the entire pipeline
to terminate.

Table of available error strategies:

============== ==================
Name            Executor
============== ==================
``terminate``   Terminates the execution as soon as an error condition is reported. Pending jobs are killed (default)
``finish``      Initiates an orderly pipeline shutdown when an error condition is raised, waiting the completion of any submitted job.
``ignore``      Ignores processes execution errors.
``retry``       Re-submit for execution a process returning an error condition.
============== ==================


When setting the ``errorStrategy`` directive to ``ignore`` the process doesn't stop on an error condition,
it just reports a message notifying you of the error event.

For example::

    process ignoreAnyError {
       errorStrategy 'ignore'

       script:
       <your command string here>
    }

.. tip:: By definition a command script fails when it ends with a non-zero exit status. To change this behavior
  see `validExitStatus`_.

The ``retry`` `error strategy`, allows you to re-submit for execution a process
returning an error condition. For example::

    process retryIfFail {
       errorStrategy 'retry'

       script:
       <your command string here>
    }


The number of times a failing process is re-executed is defined by the `maxRetries`_ and `maxErrors`_ directives.

.. _process-executor:

executor
--------

The `executor` defines the underlying system where processes are executed. By default a process uses the executor
defined globally in the ``nextflow.config`` file.

The ``executor`` directive allows you to configure what executor has to be used by the process, overriding the default
configuration. The following values can be used:

============== ==================
Name            Executor
============== ==================
``local``      The process is executed in the computer where `Nextflow` is launched.
``sge``        The process is executed using the Sun Grid Engine / `Open Grid Engine <http://gridscheduler.sourceforge.net/>`_.
``uge``        The process is executed using the `Univa Grid Engine <https://en.wikipedia.org/wiki/Univa_Grid_Engine/>`_ job scheduler.
``lsf``        The process is executed using the `Platform LSF <http://en.wikipedia.org/wiki/Platform_LSF>`_ job scheduler.
``slurm``      The process is executed using the SLURM job scheduler.
``pbs``        The process is executed using the `PBS/Torque <http://en.wikipedia.org/wiki/Portable_Batch_System>`_ job scheduler.
``condor``     The process is executed using the `HTCondor <https://research.cs.wisc.edu/htcondor/>`_ job scheduler.
``nqsii``      The process is executed using the `NQSII <https://www.rz.uni-kiel.de/en/our-portfolio/hiperf/nec-linux-cluster>`_ job scheduler.
``ignite``     The process is executed using the `Apache Ignite <https://ignite.apache.org/>`_ cluster.
``k8s``        The process is executed using the `Kubernetes <https://kubernetes.io/>`_ cluster.
============== ==================

The following example shows how to set the process's executor::


   process doSomething {

      executor 'sge'

      script:
      <your script here>

   }


.. note:: Each executor provides its own set of configuration options that can set be in the `directive` declarations block.
   See :ref:`executor-page` section to read about specific executor directives.

.. _process-ext:

ext
---

The ``ext`` is a special directive used as *namespace* for user custom process directives. This can be useful for
advanced configuration options. For example::

    process mapping {
      container "biocontainers/star:${task.ext.version}"

      input:
      file genome from genome_file
      set sampleId, file(reads) from reads_ch

      """
      STAR --genomeDir $genome --readFilesIn $reads
      """
    }

In the above example, the process uses a container whose version is controlled by the ``ext.version`` property.
This can be defined in the ``nextflow.config`` file as shown below::

    process.ext.version = '2.5.3'



.. _process-maxErrors:

maxErrors
---------

The ``maxErrors`` directive allows you to specify the maximum number of times a process can fail when using the ``retry`` `error strategy`.
By default this directive is disabled, you can set it as shown in the example below::

    process retryIfFail {
      errorStrategy 'retry'
      maxErrors 5

      """
      echo 'do this as that .. '
      """
    }
    
.. note:: This setting considers the **total** errors accumulated for a given process, across all instances. If you want
  to control the number of times a process **instance** (aka task) can fail, use ``maxRetries``.

See also: `errorStrategy`_ and `maxRetries`_.

.. _process-maxForks:

maxForks
--------

The ``maxForks`` directive allows you to define the maximum number of process instances that can be executed in parallel.
By default this value is equals to the number of CPU cores available minus 1.

If you want to execute a process in a sequential manner, set this directive to one. For example::

    process doNotParallelizeIt {

       maxForks 1

       '''
       <your script here>
       '''

    }

.. _process-maxRetries:

maxRetries
----------

The ``maxRetries`` directive allows you to define the maximum number of times a process instance can be
re-submitted in case of failure. This value is applied only when using the ``retry`` `error strategy`. By default
only one retry is allowed, you can increase this value as shown below::

    process retryIfFail {
        errorStrategy 'retry'
        maxRetries 3

        """
        echo 'do this as that .. '
        """
    }


.. note:: There is a subtle but important difference between ``maxRetries`` and the ``maxErrors`` directive.
    The latter defines the total number of errors that are allowed during the process execution (the same process can
    launch different execution instances), while the ``maxRetries`` defines the maximum number of times the same process
    execution can be retried in case of an error.

See also: `errorStrategy`_ and `maxErrors`_.


.. _process-memory:

memory
------

The ``memory`` directive allows you to define how much memory the process is allowed to use. For example::

    process big_job {

        memory '2 GB'
        executor 'sge'

        """
        your task script here
        """
    }


The following memory unit suffix can be used when specifying the memory value:

======= =============
Unit    Description
======= =============
B       Bytes
KB      Kilobytes
MB      Megabytes
GB      Gigabytes
TB      Terabytes
======= =============

.. This setting is equivalent to set the ``qsub -l virtual_free=<mem>`` command line option.

See also: `cpus`_, `time`_, `queue`_ and `Dynamic computing resources`_.


.. _process-module:

module
------

`Environment Modules <http://modules.sourceforge.net/>`_ is a package manager that allows you to dynamically configure
your execution environment and easily switch between multiple versions of the same software tool.

If it is available in your system you can use it with Nextflow in order to configure the processes execution
environment in your pipeline.

In a process definition you can use the ``module`` directive to load a specific module version to be used in the
process execution environment. For example::

  process basicExample {

    module 'ncbi-blast/2.2.27'

    """
    blastp -query <etc..>
    """
  }

You can repeat the ``module`` directive for each module you need to load. Alternatively multiple modules
can be specified in a single ``module`` directive by separating all the module names by using a ``:``
(colon) character as shown below::

   process manyModules {

     module 'ncbi-blast/2.2.27:t_coffee/10.0:clustalw/2.1'

     """
     blastp -query <etc..>
     """
  }


.. _process-penv:

penv
----

The ``penv`` directive  allows you to define the `parallel environment` to be used when submitting a parallel task to the
:ref:`SGE <sge-executor>` resource manager. For example::

    process big_job {

      cpus 4
      penv 'smp'
      executor 'sge'

      """
      blastp -query input_sequence -num_threads ${task.cpus}
      """
    }

This configuration depends on the parallel environment provided by your grid engine installation. Refer to your
cluster documentation or contact your admin to lean more about this.

.. note:: This setting is available when using the :ref:`sge-executor` executor.

See also: `cpus`_, `memory`_, `time`_

.. _process-pod:

pod
---

The ``pod`` directive allows the definition of pods specific settings, such as environment variables, secrets
and config maps when using the :ref:`k8s-executor` executor.

For example::

  process your_task {
    pod env: 'FOO', value: 'bar'

    '''
    echo $FOO
    '''
  }

The above snippet defines an environment variable named ``FOO`` which value is ``bar``.

The ``pod`` directive allows the definition of the following options:

================================================= =================================================
``label: <K>, value: <V>``                        Defines a pod label with key ``K`` and value ``V``.
``env: <E>, value: <V>``                          Defines an environment variable with name ``E`` and whose value is given by the ``V`` string.
``env: <E>, config: <C/K>``                       Defines an environment variable with name ``E`` and whose value is given by the entry associated to the key with name ``K`` in the `ConfigMap <https://kubernetes.io/docs/tasks/configure-pod-container/configure-pod-configmap/>`_ with name ``C``.
``env: <E>, secret: <S/K>``                       Defines an environment variable with name ``E`` and whose value is given by the entry associated to the key with name ``K`` in the `Secret <https://kubernetes.io/docs/concepts/configuration/secret/>`_ with name ``S``.
``config: <C/K>, mountPath: </absolute/path>``    The content of the `ConfigMap <https://kubernetes.io/docs/tasks/configure-pod-container/configure-pod-configmap/>`_ with name ``C`` with key ``K`` is made available to the path ``/absolute/path``. When the key component is omitted the path is interpreted as a directory and all the `ConfigMap` entries are exposed in that path.
``secret: <S/K>, mountPath: </absolute/path>``    The content of the `Secret <https://kubernetes.io/docs/concepts/configuration/secret/>`_ with name ``S`` with key ``K`` is made available to the path ``/absolute/path``. When the key component is omitted the path is interpreted as a directory and all the `Secret` entries are exposed in that path.
``volumeClaim: <V>, mountPath: </absolute/path>`` Mounts a `Persistent volume claim <https://kubernetes.io/docs/concepts/storage/persistent-volumes/>`_ with name ``V`` to the specified path location. Use the optional `subPath` parameter to mount a directory inside the referenced volume instead of its root.
``imagePullPolicy: <V>``                          Specifies the strategy to be used to pull the container image e.g. ``imagePullPolicy: 'Always'``.
``imagePullSecret: <V>``                          Specifies the secret name to access a private container image registry. See `Kubernetes documentation <https://kubernetes.io/docs/concepts/containers/images/#specifying-imagepullsecrets-on-a-pod>`_ for details.
``runAsUser: <UID>``                              Specifies the user ID to be used to run the container.
================================================= =================================================

When defined in the Nextflow configuration file, a pod setting can be defined using the canonical
associative array syntax. For example::

  process {
    pod = [env: 'FOO', value: 'bar']
  }

When more than one setting needs to be provides they must be enclosed in a list definition as shown below::

  process {
    pod = [ [env: 'FOO', value: 'bar'], [secret: 'my-secret/key1', mountPath: '/etc/file.txt'] ]
  }


.. _process-publishDir:

publishDir
----------

The ``publishDir`` directive allows you to publish the process output files to a specified folder. For example::


    process foo {

        publishDir '/data/chunks'

        output:
        file 'chunk_*' into letters

        '''
        printf 'Hola' | split -b 1 - chunk_
        '''
    }


The above example splits the string ``Hola`` into file chunks of a single byte. When complete the ``chunk_*`` output files
are published into the ``/data/chunks`` folder.

.. tip:: The ``publishDir`` directive can be specified more than one time in to publish the output files
  to different target directories. This feature requires version 0.29.0 or higher.

By default files are published to the target folder creating a *symbolic link* for each process output that links
the file produced into the process working directory. This behavior can be modified using the ``mode`` parameter.

Table of optional parameters that can be used with the ``publishDir`` directive:

=============== =================
Name            Description
=============== =================
mode            The file publishing method. See the following table for possible values.
overwrite       When ``true`` any existing file in the specified folder will be overridden (default: ``true`` during normal
                pipeline execution and ``false`` when pipeline execution is `resumed`).
pattern         Specifies a `glob`_ file pattern that selects which files to publish from the overall set of output files.
path            Specifies the directory where files need to be published. **Note**: the syntax ``publishDir '/some/dir'`` is a shortcut for ``publishDir path: '/some/dir'``.
saveAs          A closure which, given the name of the file being published, returns the actual file name or a full path where the file is required to be stored.
                This can be used to rename or change the destination directory of the published files dynamically by using
                a custom strategy.
                Return the value ``null`` from the closure to *not* publish a file.
                This is useful when the process has multiple output files, but you want to publish only some of them.
=============== =================

Table of publish modes:

=============== =================
 Mode           Description
=============== =================
symlink         Creates an absolute `symbolic link` in the published directory for each process output file (default).
rellink         Creates a relative `symbolic link` in the published directory for each process output file.
link            Creates a `hard link` in the published directory for each process output file.
copy            Copies the output files into the published directory.
copyNoFollow    Copies the output files into the published directory without following symlinks ie. copies the links themselves. 
move            Moves the output files into the published directory. **Note**: this is only supposed to be used for a `terminating` process i.e. a process whose output is not consumed by any other downstream process.
=============== =================

.. note:: The `mode` value needs to be specified as a string literal i.e. enclosed by quote characters. Multiple parameters
  need to be separated by a colon character. For example:

::

    process foo {

        publishDir '/data/chunks', mode: 'copy', overwrite: false

        output:
        file 'chunk_*' into letters

        '''
        printf 'Hola' | split -b 1 - chunk_
        '''
    }


.. warning:: Files are copied into the specified directory in an *asynchronous* manner, thus they may not be immediately
  available in the published directory at the end of the process execution. For this reason files published by a process
  must not be accessed by other downstream processes.


.. _process-queue:

queue
-----

The ``queue`` directory allows you to set the `queue` where jobs are scheduled when using a grid based executor
in your pipeline. For example::

    process grid_job {

        queue 'long'
        executor 'sge'

        """
        your task script here
        """
    }


Multiple queues can be specified by separating their names with a comma for example::

    process grid_job {

        queue 'short,long,cn-el6'
        executor 'sge'

        """
        your task script here
        """
    }


.. note:: This directive is taken in account only by the following executors: :ref:`sge-executor`, :ref:`lsf-executor`,
  :ref:`slurm-executor` and :ref:`pbs-executor` executors.


.. _process-label:

label
-----

The ``label`` directive allows the annotation of processes with mnemonic identifier of your choice.
For example::

  process bigTask {

    label 'big_mem'

    '''
    <task script>
    '''
  }


The same label can be applied to more than a process and multiple labels can be applied to the same
process using the ``label`` directive more than one time.

.. note:: A label must consist of alphanumeric characters or ``_``, must start with an alphabetic character
  and must end with an alphanumeric character.

Labels are useful to organise workflow processes in separate groups which can be referenced
in the configuration file to select and configure subset of processes having similar computing requirements.

See the :ref:`config-process-selectors` documentation for details.


.. _process-scratch:

scratch
-------

The ``scratch`` directive allows you to execute the process in a temporary folder that is local to the execution node.

This is useful when your pipeline is launched by using a `grid` executor, because it permits to decrease the NFS
overhead by running the pipeline processes in a temporary directory in the local disk of the actual execution node.
Only the files declared as output in the process definition will be copied in the pipeline working area.

In its basic form simply specify ``true`` at the directive value, as shown below::

  process simpleTask {

    scratch true

    output:
    file 'data_out'

    '''
    <task script>
    '''
  }


By doing this, it tries to execute the script in the directory defined by the variable ``$TMPDIR`` in the execution node.
If this variable does not exist, it will create a new temporary directory by using the Linux command ``mktemp``.

A custom environment variable, other than ``$TMPDIR``, can be specified by simply using it as the scratch value, for
example::

  scratch '$MY_GRID_TMP'

Note, it must be wrapped by single quotation characters, otherwise the variable will be evaluated in the
pipeline script context.

You can also provide a specific folder path as scratch value, for example::

  scratch '/tmp/my/path'

By doing this, a new temporary directory will be created in the specified path each time a process is executed.

Finally, when the ``ram-disk`` string is provided as ``scratch`` value, the process will be execute in the node
RAM virtual disk.

Summary of allowed values:

=========== ==================
scratch     Description
=========== ==================
false       Do not use the scratch folder.
true        Creates a scratch folder in the directory defined by the ``$TMPDIR`` variable; fallback to ``mktemp /tmp`` if that variable do not exists.
$YOUR_VAR   Creates a scratch folder in the directory defined by the ``$YOUR_VAR`` environment variable; fallback to ``mktemp /tmp`` if that variable do not exists.
/my/tmp     Creates a scratch folder in the specified directory.
ram-disk    Creates a scratch folder in the RAM disk ``/dev/shm/`` (experimental).
=========== ==================

.. _process-storeDir:

storeDir
--------

The ``storeDir`` directive allows you to define a directory that is used as `permanent` cache for your process results.

In more detail, it affects the process execution in two main ways:

#. The process is executed only if the files declared in the `output` clause do not exist in the directory specified by
   the ``storeDir`` directive. When the files exist the process execution is skipped and these files are used as
   the actual process result.

#. Whenever a process is successfully completed the files listed in the `output` declaration block are copied into the directory
   specified by the ``storeDir`` directive.

The following example shows how to use the ``storeDir`` directive to create a directory containing a BLAST database
for each species specified by an input parameter::

  genomes = Channel.fromPath(params.genomes)

  process formatBlastDatabases {

    storeDir '/db/genomes'

    input:
    file species from genomes

    output:
    file "${dbName}.*" into blastDb

    script:
    dbName = species.baseName
    """
    makeblastdb -dbtype nucl -in ${species} -out ${dbName}
    """

  }


.. warning:: The ``storeDir`` directive is meant for long term process caching and should not be used to
    output the files produced by a process to a specific folder or organise result data in `semantic` directory structure.
    In these cases you may use the `publishDir`_ directive instead.

.. note:: The use of AWS S3 path is supported however it requires the installation of the `AWS CLI tool <https://aws.amazon.com/cli/>`_
  (ie. ``aws``) in the target computing node.

.. _process-stageInMode:

stageInMode
-----------

The ``stageInMode`` directive defines how input files are staged-in to the process work directory. The following values
are allowed:

======= ==================
Value   Description
======= ==================
copy    Input files are staged in the process work directory by creating a copy.
link    Input files are staged in the process work directory by creating an (hard) link for each of them.
symlink Input files are staged in the process work directory by creating a symbolic link with an absolute path for each of them (default).
rellink Input files are staged in the process work directory by creating a symbolic link with a relative path for each of them.
======= ==================


.. _process-stageOutMode:

stageOutMode
------------

The ``stageOutMode`` directive defines how output files are staged-out from the scratch directory to the process work
directory. The following values are allowed:

======= ==================
Value   Description
======= ==================
copy    Output files are copied from the scratch directory to the work directory.
move    Output files are moved from the scratch directory to the work directory.
rsync   Output files are copied from the scratch directory to the work directory by using the ``rsync`` utility.
======= ==================

See also: `scratch`_.


.. _process-tag:

tag
---

The ``tag`` directive allows you to associate each process executions with a custom label, so that it will be easier
to identify them in the log file or in the trace execution report. For example::

    process foo {
      tag "$code"

      input:
      val code from 'alpha', 'gamma', 'omega'

      """
      echo $code
      """
    }

The above snippet will print a log similar to the following one, where process names contain the tag value::

    [6e/28919b] Submitted process > foo (alpha)
    [d2/1c6175] Submitted process > foo (gamma)
    [1c/3ef220] Submitted process > foo (omega)


See also :ref:`Trace execution report <trace-report>`


.. _process-time:

time
----

The ``time`` directive allows you to define how long a process is allowed to run. For example::

    process big_job {

        time '1h'

        """
        your task script here
        """
    }



The following time unit suffix can be used when specifying the duration value:

======= =============
Unit    Description
======= =============
s       Seconds
m       Minutes
h       Hours
d       Days
======= =============

.. note:: This directive is taken in account only when using one of the following grid based executors:
  :ref:`sge-executor`, :ref:`lsf-executor`, :ref:`slurm-executor`, :ref:`pbs-executor`,
  :ref:`condor-executor` and :ref:`awsbatch-executor` executors.

See also: `cpus`_, `memory`_, `queue`_ and `Dynamic computing resources`_.


.. _process-validExitStatus:

validExitStatus
---------------

A process is terminated when the executed command returns an error exit status. By default any error status
other than ``0`` is interpreted as an error condition.

The ``validExitStatus`` directive allows you to fine control which error status will represent a successful command execution.
You can specify a single value or multiple values as shown in the following example::


    process returnOk {
        validExitStatus 0,1,2

         script:
         """
         echo Hello
         exit 1
         """
    }


In the above example, although the command script ends with a ``1`` exit status, the process
will not return an error condition because the value ``1`` is declared as a `valid` status in
the ``validExitStatus`` directive.



Dynamic directives
------------------

A directive can be assigned *dynamically*, during the process execution, so that its actual value can be evaluated
depending on the value of one, or more, process' input values.

In order to be defined in a dynamic manner the directive's value needs to be expressed by using a :ref:`closure <script-closure>`
statement, as in the following example::

    process foo {

      executor 'sge'
      queue { entries > 100 ? 'long' : 'short' }

      input:
      set entries, file(x) from data

      script:
      """
      < your job here >
      """
    }

In the above example the `queue`_ directive is evaluated dynamically, depending on the input value ``entries``. When it is
bigger than 100, jobs will be submitted to the queue ``long``, otherwise the ``short`` one will be used.

All directives can be assigned to a dynamic value except the following:

* `executor`_
* `maxForks`_


.. note:: You can retrieve the current value of a dynamic directive in the process script by using the implicit variable ``task``
  which holds the directive values defined in the current process instance.

For example::


   process foo {

      queue { entries > 100 ? 'long' : 'short' }

      input:
      set entries, file(x) from data

      script:
      """
      echo Current queue: ${task.queue}
      """
    }


Dynamic computing resources
---------------------------

It's a very common scenario that different instances of the same process may have very different needs in terms of computing resources. 
In such situations requesting, for example, an amount of memory too low will cause some tasks to fail. 
Instead, using a higher limit that fits all the tasks in your execution could significantly decrease the execution priority of your jobs.

The `Dynamic directives`_ evaluation feature can be used to modify the amount of computing resources requested in case
of a process failure and try to re-execute it using a higher limit. For example::


    process foo {

        memory { 2.GB * task.attempt }
        time { 1.hour * task.attempt }

        errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }
        maxRetries 3

        script:
        <your job here>

    }


In the above example the `memory`_ and execution `time`_ limits are defined dynamically. The first time the process
is executed the ``task.attempt`` is set to ``1``, thus it will request a two GB of memory and one hour of maximum execution
time.

If the task execution fail reporting an exit status equals ``140``, the task is re-submitted (otherwise terminates immediately).
This time the value of ``task.attempt`` is ``2``, thus increasing the amount of the memory to four GB and the time to 2 hours, and so on.

The directive `maxRetries`_ set the maximum number of time the same task can be re-executed.

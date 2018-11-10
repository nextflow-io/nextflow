.. _example-page:

********
Examples
********

Basic pipeline
--------------

This example shows a pipeline that is made of two processes. The first process receives a
`FASTA formatted <http://en.wikipedia.org/wiki/FASTA_format>`_ file and splits it into file chunks whose names start with
the prefix ``seq_``.

The process that follows, receives these files and it simply `reverses` their content by using the ``rev`` command line tool.

.. raw:: html

    <style type='text/css'>
    .gist {font-size:13px;line-height:18px;margin-bottom:20px;width:100%}
    .gist pre {font-family:Menlo,Monaco,'Bitstream Vera Sans Mono','Courier New',monospace !important}
    .gist-meta {font-family:Helvetica,Arial,sans-serif;font-size:13px !important}
    .gist-meta a {color:#26a !important;text-decoration:none}
    .gist-meta a:hover{color:#0e4071 !important}
    .gist .file-data { background-color: white; }
    </style>
    <script src="https://gist.github.com/pditommaso/92fea4525cd66c286904.js"></script>


In more detail:

* line 1: The script starts with a `shebang <http://en.wikipedia.org/wiki/Shebang_(Unix)>`_ declaration. This allows you
  to launch your pipeline, as any other Bash script

* line 3: Declares a pipeline parameter named ``params.in`` that is initialized with the value ``$HOME/sample.fa``.This value
  can be overridden when launching the pipeline, by simply adding the option ``--in <value>`` to the script command line

* line 5: Defines a variable ``sequences`` holding a reference for the file whose name is specified by the ``params.in``
  parameter

* line 6: Defines a variable ``SPLIT`` whose value is ``gcsplit`` when the script is executed on a Mac OSX or ``csplit``
  when it runs on Linux. This is the name of the tool that is used to split the file.

* lines 8-20: The process that splits the provided file.

* line 10: Opens the `input` declaration block. The lines following this clause are interpreted as input definitions.

* line 11: Defines the process input file. This file is received from the variable ``sequences`` and will be named ``input.fa``.

* line 13: Opens the `output` declaration block. Lines following this clause are interpreted as output definitions.

* line 14: Defines that the process outputs files whose names match the pattern ``seq_*``. These files are sent over the
  channel ``records``.

* lines 16-18: The actual script executed by the process to split the provided file.

* lines 22-33: Defines the second process, that receives the splits produced by the previous process and reverses their
  content.

* line 24: Opens the `input` declaration block. Lines following this clause are interpreted as input definitions.

* line 25: Defines the process input file. This file is received through the channel ``records``.

* line 27: Opens the `output` declaration block. Lines following this clause are interpreted as output definitions.

* line 28: The `standard output` of the executed script is declared as the process output. This output is sent over the
  channel ``result``.

* lines 30-32: The actual script executed by the process to `reverse` the content of the received files.

* line 35: Prints a `result` each time a new item is received on the ``result`` channel.


.. tip:: The above example can manage only a single file at a time. If you want to execute it for two (or more) different files
   you will need to launch it several times.

   It is possible to modify it in such a way that it can handle any number of input files, as shown below.


In order to make the above script able to handle any number of files simply replace `line 3` with the following line::

  sequences = Channel.fromPath(params.in)


By doing this the ``sequences`` variable is assigned to the channel created by the :ref:`channel-path` method. This
channel emits all the files that match the pattern specified by the parameter ``params.in``.

Given that you saved the script to a file named ``example.nf`` and you have a list of `FASTA` files in a folder
named ``dataset/``, you can execute it by entering this command::

  nextflow example.nf --in 'dataset/*.fa'


.. warning:: Make sure you enclose the ``dataset/*.fa`` parameter value in single-quotation characters,
  otherwise the Bash environment will expand the ``*`` symbol to the actual file names and the example won't work.

More examples
-------------

You can find at `this link <https://github.com/nextflow-io/examples>`_ a collection of examples introducing Nextflow
scripting.

Check also `Awesome Nextflow <https://github.com/nextflow-io/awesome-nextflow/>`_ for a list
of pipelines developed by the Nextflow community.
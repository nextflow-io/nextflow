.. _channel-page:

********
Channels
********

Nextflow is based on the Dataflow programming model in which processes communicate through channels.

A channel has two major properties:

#. Sending a message is an `asynchronous` operation which completes immediately,
   without having to wait for the receiving process.

#. Receiving data is a blocking operation which stops the receiving process until the message has arrived.


.. _channel-types:

Channel types
=============

In Nextflow there are two kinds of channels: `queue channels` and `value channels`.


.. _channel-type-queue:

Queue channel
-------------

A `queue channel` is a non-blocking unidirectional FIFO queue which connects two processes,
channel factories, or operators.

A queue channel can be created by factory methods (`of`_, `fromPath`_, etc),
operators (:ref:`operator-map`, :ref:`operator-flatmap`, etc), and
processes (see :ref:`Process outputs <process-output>`).


.. _channel-type-value:

Value channel
-------------

A `value channel` a.k.a. *singleton channel* is bound to a single value and can be read any
number of times without being consumed.

A value channel can be created with the `value`_ factory method or by any operator that produces
a single value (:ref:`operator-first`, :ref:`operator-collect`, :ref:`operator-reduce`, etc). Additionally,
a process will emit value channels if it is invoked with all value channels, including
simple values which are implicitly wrapped in a value channel.

For example::

    process foo {
      input:
      val x

      output:
      path 'x.txt'

      """
      echo $x > x.txt
      """
    }

    workflow {
      result = foo(1)
      result.view { "Result: ${it}" }
    }

In the above example, since the ``foo`` process is invoked with a simple value instead of a channel,
the input is implicitly converted to a value channel, and the output is also emitted as a value channel.

See also: :ref:`process-multiple-input-channels`.


.. _channel-factory:

Channel factory
===============

Channels may be created explicitly using the following channel factory methods.

.. note::
  As of version 20.07.0, ``channel`` has been introduced as an alias of ``Channel``, therefore factory
  methods can be specified either as ``channel.of()`` or ``Channel.of()``, and so on.


.. _channel-empty:

empty
-----

The ``empty`` factory method, by definition, creates a channel that doesn't emit any value.

See also: :ref:`operator-ifempty` operator.


.. _channel-from:

from
----

.. warning::
    This method is deprecated. Use `of`_ or `fromList`_ instead.

The ``from`` method allows you to create a channel emitting any sequence of values that are specified as the method argument,
for example::

    ch = Channel.from( 1, 3, 5, 7 )
    ch.subscribe { println "value: $it" }

The first line in this example creates a variable ``ch`` which holds a channel object. This channel emits the values
specified as a parameter in the ``from`` method. Thus the second line will print the following::

    value: 1
    value: 3
    value: 5
    value: 7

The following example shows how to create a channel from a `range` of numbers or strings::

    zeroToNine = Channel.from( 0..9 )
    strings = Channel.from( 'A'..'Z' )

.. note::
  When the ``from`` argument is an object implementing the (Java)
  `Collection <http://docs.oracle.com/javase/7/docs/api/java/util/Collection.html>`_ interface, the resulting channel
  emits the collection entries as individual items.

Thus the following two declarations produce an identical result even though in the first case the items are specified
as multiple arguments while in the second case as a single list object argument::

    Channel.from( 1, 3, 5, 7, 9 )
    Channel.from( [1, 3, 5, 7, 9] )

But when more than one argument is provided, they are always managed as `single` emissions. Thus, the following example
creates a channel emitting three entries each of which is a list containing two elements::

    Channel.from( [1, 2], [5,6], [7,9] )


.. _channel-fromlist:

fromList
--------

.. note::
  This feature requires Nextflow version 19.10.0 or later.

The ``fromList`` method allows you to create a channel emitting the values provided as a list of elements,
for example::

    Channel
        .fromList( ['a', 'b', 'c', 'd'] )
        .view { "value: $it" }

Prints::

    value: a
    value: b
    value: c
    value: d

See also: `of`_ factory method.


.. _channel-path:

fromPath
--------

You can create a channel emitting one or more file paths by using the ``fromPath`` method and specifying a path string
as an argument. For example::

    myFileChannel = Channel.fromPath( '/data/some/bigfile.txt' )

The above line creates a channel and binds it to a `Path <http://docs.oracle.com/javase/7/docs/api/java/nio/file/Path.html>`_
object for the specified file.

.. note::
    ``fromPath`` does not check whether the file exists.

Whenever the ``fromPath`` argument contains a ``*`` or ``?`` wildcard character it is interpreted as a `glob`_ path matcher.
For example::

    myFileChannel = Channel.fromPath( '/data/big/*.txt' )

This example creates a channel and emits as many ``Path`` items as there are files with ``txt`` extension in the ``/data/big`` folder.

.. tip::
  Two asterisks, i.e. ``**``, works like ``*`` but crosses directory boundaries.
  This syntax is generally used for matching complete paths. Curly brackets specify a collection of sub-patterns.

For example::

    files = Channel.fromPath( 'data/**.fa' )
    moreFiles = Channel.fromPath( 'data/**/*.fa' )
    pairFiles = Channel.fromPath( 'data/file_{1,2}.fq' )

The first line returns a channel emitting the files ending with the suffix ``.fa`` in the ``data`` folder `and` recursively
in all its sub-folders. While the second one only emits the files which have the same suffix in `any` sub-folder in the ``data`` path.
Finally the last example emits two files: ``data/file_1.fq`` and ``data/file_2.fq``.

.. note::
    As in Linux Bash, the ``*`` wildcard does not catch hidden files (i.e. files whose name starts with a ``.`` character).

In order to include hidden files, you need to start your pattern with a period character or specify the ``hidden: true`` option. For example::

    expl1 = Channel.fromPath( '/path/.*' )
    expl2 = Channel.fromPath( '/path/.*.fa' )
    expl3 = Channel.fromPath( '/path/*', hidden: true )

The first example returns all hidden files in the specified path. The second one returns all hidden files
ending with the ``.fa`` suffix. Finally the last example returns all files (hidden and non-hidden) in that path.

By default a `glob`_ pattern only looks for `regular file` paths that match the specified criteria, i.e.
it won't return directory paths.

You may use the parameter ``type`` specifying the value ``file``, ``dir`` or ``any`` in order to define what kind of paths
you want. For example::

    myFileChannel = Channel.fromPath( '/path/*b', type: 'dir' )
    myFileChannel = Channel.fromPath( '/path/a*', type: 'any' )

The first example will return all `directory` paths ending with the ``b`` suffix, while the second will return any file
and directory starting with a ``a`` prefix.

=============== ===================
Name            Description
=============== ===================
glob            When ``true`` interprets characters ``*``, ``?``, ``[]`` and ``{}`` as glob wildcards, otherwise handles them as normal characters (default: ``true``)
type            Type of paths returned, either ``file``, ``dir`` or ``any`` (default: ``file``)
hidden          When ``true`` includes hidden files in the resulting paths (default: ``false``)
maxDepth        Maximum number of directory levels to visit (default: `no limit`)
followLinks     When ``true`` it follows symbolic links during directories tree traversal, otherwise they are managed as files (default: ``true``)
relative        When ``true`` returned paths are relative to the top-most common directory (default: ``false``)
checkIfExists   When ``true`` throws an exception of the specified path do not exist in the file system (default: ``false``)
=============== ===================

.. note::
  Multiple paths or glob patterns can be specified using a list::

      Channel.fromPath( ['/some/path/*.fq', '/other/path/*.fastq'] )


.. _channel-filepairs:

fromFilePairs
-------------

The ``fromFilePairs`` method creates a channel emitting the file pairs matching a `glob`_ pattern provided by the user.
The matching files are emitted as tuples in which the first element is the grouping key of the matching
pair and the second element is the list of files (sorted in lexicographical order). For example::

    Channel
        .fromFilePairs('/my/data/SRR*_{1,2}.fastq')
        .view()

It will produce an output similar to the following::

    [SRR493366, [/my/data/SRR493366_1.fastq, /my/data/SRR493366_2.fastq]]
    [SRR493367, [/my/data/SRR493367_1.fastq, /my/data/SRR493367_2.fastq]]
    [SRR493368, [/my/data/SRR493368_1.fastq, /my/data/SRR493368_2.fastq]]
    [SRR493369, [/my/data/SRR493369_1.fastq, /my/data/SRR493369_2.fastq]]
    [SRR493370, [/my/data/SRR493370_1.fastq, /my/data/SRR493370_2.fastq]]
    [SRR493371, [/my/data/SRR493371_1.fastq, /my/data/SRR493371_2.fastq]]

.. note::
    The glob pattern must contain at least one ``*`` wildcard character.

Alternatively it is possible to implement a custom file pair grouping strategy providing a closure which,
given the current file as parameter, returns the grouping key.
For example::

    Channel
        .fromFilePairs('/some/data/*', size: -1) { file -> file.extension }
        .view { ext, files -> "Files with the extension $ext are $files" }

Table of optional parameters available:

=============== ===================
Name            Description
=============== ===================
type            Type of paths returned, either ``file``, ``dir`` or ``any`` (default: ``file``)
hidden          When ``true`` includes hidden files in the resulting paths (default: ``false``)
maxDepth        Maximum number of directory levels to visit (default: `no limit`)
followLinks     When ``true`` it follows symbolic links during directories tree traversal, otherwise they are managed as files (default: ``true``)
size            Defines the number of files each emitted item is expected to hold (default: 2). Set to ``-1`` for any.
flat            When ``true`` the matching files are produced as sole elements in the emitted tuples (default: ``false``).
checkIfExists   When ``true`` throws an exception of the specified path do not exist in the file system (default: ``false``)
=============== ===================

.. note::
  Multiple glob patterns can be specified using a list::

      Channel.fromFilePairs( ['/some/data/SRR*_{1,2}.fastq', '/other/data/QFF*_{1,2}.fastq'] )


.. _channel-fromsra:

fromSRA
-------

.. note:: This feature requires Nextflow version 19.04.0 or later.

The ``fromSRA`` method queries the `NCBI SRA <https://www.ncbi.nlm.nih.gov/sra>`_ database and returns a channel emitting
the FASTQ files matching the specified criteria i.e project or accession number(s). For example::

    Channel
        .fromSRA('SRP043510')
        .view()

It returns::

    [SRR1448794, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/004/SRR1448794/SRR1448794.fastq.gz]
    [SRR1448795, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/005/SRR1448795/SRR1448795.fastq.gz]
    [SRR1448792, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/002/SRR1448792/SRR1448792.fastq.gz]
    [SRR1448793, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/003/SRR1448793/SRR1448793.fastq.gz]
    [SRR1910483, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/003/SRR1910483/SRR1910483.fastq.gz]
    [SRR1910482, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/002/SRR1910482/SRR1910482.fastq.gz]
    (remaining omitted)

Multiple accession IDs can be specified using a list object::

    ids = ['ERR908507', 'ERR908506', 'ERR908505']
    Channel
        .fromSRA(ids)
        .view()

::

    [ERR908507, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908507/ERR908507_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908507/ERR908507_2.fastq.gz]]
    [ERR908506, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908506/ERR908506_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908506/ERR908506_2.fastq.gz]]
    [ERR908505, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908505/ERR908505_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908505/ERR908505_2.fastq.gz]]

.. note::
  Each read pair is implicitly managed and returned as a list of files.

.. tip::
  This method uses the NCBI `ESearch <https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch>`_
  API behind the scenes, therefore it allows the use of any query term supported by this API.

Table of optional parameters available:

=============== ===================
Name            Description
=============== ===================
apiKey          NCBI user API key.
cache           Enable/disable the caching API requests (default: ``true``).
max             Maximum number of entries that can be retried (default: unlimited) .
protocol        Allow choosing the protocol for the resulting remote URLs. Available choices: ``ftp``, ``http``, ``https`` (default: ``ftp``).
=============== ===================

To access the NCBI search service the `NCBI API keys <https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities>`_
should be provided either:

* Using the ``apiKey`` optional parameter e.g. ``Channel.fromSRA(ids, apiKey:'0123456789abcdef')``.
* Exporting the ``NCBI_API_KEY`` variable in your environment e.g. ``export NCBI_API_KEY=0123456789abcdef``.


.. _channel-of:

of
--

.. note::
  This feature requires Nextflow version 19.10.0 of later.

The ``of`` method allows you to create a channel that emits the arguments provided to it,
for example::

    ch = Channel.of( 1, 3, 5, 7 )
    ch.view { "value: $it" }

The first line in this example creates a variable ``ch`` which holds a channel object. This channel emits the arguments
supplied to the ``of`` method. Thus the second line prints the following::

    value: 1
    value: 3
    value: 5
    value: 7

Ranges of values are expanded accordingly::

    Channel
        .of(1..23, 'X', 'Y')
        .view()

Prints::

    1
    2
    3
    4
    :
    23
    X
    Y

See also: `fromList`_ factory method.


.. _channel-value:

value
-----

The ``value`` method is used to create a value channel. An optional (not ``null``) argument
can be specified to bind the channel to a specific value. For example::

    expl1 = Channel.value()
    expl2 = Channel.value( 'Hello there' )
    expl3 = Channel.value( [1,2,3,4,5] )

The first line in the example creates an 'empty' variable. The second line creates a channel and binds a string to it.
The third line creates a channel and binds a list object to it that will be emitted as a single value.


.. _channel-watchpath:

watchPath
---------

The ``watchPath`` method watches a folder for one or more files matching a specified pattern. As soon as
there is a file that meets the specified condition, it is emitted over the channel that is returned by the ``watchPath``
method. The condition on files to watch can be specified by using ``*`` or ``?`` wildcard characters i.e. by specifying
a `glob`_ path matching criteria.

For example::

    Channel
        .watchPath( '/path/*.fa' )
        .subscribe { println "Fasta file: $it" }

By default it watches only for new files created in the specified folder. Optionally, it is possible to provide a
second argument that specifies what event(s) to watch. The supported events are:

=========== ================
Name        Description
=========== ================
``create``  A new file is created (default)
``modify``  A file is modified
``delete``  A file is deleted
=========== ================

You can specify more than one of these events by using a comma separated string as shown below::

    Channel
        .watchPath( '/path/*.fa', 'create,modify' )
        .subscribe { println "File created or modified: $it" }

.. warning::
    The ``watchPath`` factory waits endlessly for files that match the specified pattern and event(s),
    which means that it will cause your pipeline to run forever. Consider using the ``until`` operator
    to close the channel when a certain condition is met (e.g. receiving a file named ``DONE``).

See also: `fromPath`_ factory method.


.. _glob: http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob

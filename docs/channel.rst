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

Nextflow distinguish two different kinds of channels: `queue channels` and `value channels`.

Queue channel
-------------

A `queue channel` is a non-blocking unidirectional FIFO queue which connects two processes or operators.

A queue channel is usually created using a factory method such as a `from`_, `fromPath`_, etc.
or chaining it with a channel operator such as :ref:`operator-map`, :ref:`operator-flatmap`, etc.

Queue channels are also created by process output declarations using the ``into`` clause.

.. note:: The definition implies that the same queue channel cannot be used more than one time as process
 output and more then one time as process input.

In you need to connect a process output channel to more then one process or operator use the
:ref:`operator-into` operator to create two (or more) copies of the same channel and use each
of them to connect a separate process.


Value channel
-------------

A `value channel` a.k.a. *singleton channel* by definition is bound to a single value and it can be read
unlimited times without consuming its content.

.. tip:: For this reason a value channel can be used as input by more than one process.

A value channel is created using the `value`_ factory method or by operators returning
a single value, such us :ref:`operator-first`, :ref:`operator-last`, :ref:`operator-collect`,
:ref:`operator-count`, :ref:`operator-min`, :ref:`operator-max`, :ref:`operator-reduce`, :ref:`operator-sum`.


.. note:: A value channel is implicitly created by a process when an input specifies a simple value
  in the ``from`` clause.
  Moreover, a value channel is also implicitly created as output for a process whose
  inputs are only value channels.

For example::

    process foo {
      input:
      val x from 1
      output:
      file 'x.txt' into result

      """
      echo $x > x.txt
      """
    }

The process in the above snippet declare a single input which implicitly is a value channel.
Therefore also the ``result`` output is a value channel that can be read by more than one process.

See also: :ref:`process-understand-how-multiple-input-channels-work`.

.. _channel-factory:

Channel factory
===============

Channels may be created implicitly by the process output(s) declaration or explicitly using the following channel
factory methods.

The available factory methods are:

* `create`_
* `empty`_
* `from`_
* `fromPath`_
* `fromFilePairs`_
* `value`_
* `watchPath`_

.. _channel-create:

create
------

Creates a new `channel` by using the ``create`` method, as shown below::

    channelObj = Channel.create()


.. _channel-from:

from
----

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



.. note:: Note that when the ``from`` argument is an object implementing the (Java)
  `Collection <http://docs.oracle.com/javase/7/docs/api/java/util/Collection.html>`_ interface, the resulting channel
  emits the collection entries as individual emissions.

Thus the following two declarations produce an identical result even tough in the first case the items are specified
as multiple arguments while in the second case as a single list object argument::

    Channel.from( 1, 3, 5, 7, 9 )
    Channel.from( [1, 3, 5, 7, 9] )


But when more than one argument is provided, they are always managed as `single` emissions. Thus, the following example
creates a channel emitting three entries each of which is a list containing two elements::

    Channel.from( [1, 2], [5,6], [7,9] )



.. _channel-value:

value
-----

The `value` factory method is used to create a *value* channel. An optional not ``null`` argument
can be specified to bind the channel to a specific value. For example::


    expl1 = Channel.value()
    expl2 = Channel.value( 'Hello there' )
    expl3 = Channel.value( [1,2,3,4,5] )


The first line in the example creates an 'empty' variable. The second line creates a channel and binds a string to it.
Finally the last one creates a channel and binds a list object to it that will be emitted as a sole emission.

.. _channel-path:

fromPath
--------

You can create a channel emitting one or more file paths by using the ``fromPath`` method and specifying a path string
as an argument. For example::

    myFileChannel = Channel.fromPath( '/data/some/bigfile.txt' )

The above line creates a channel and binds to it a `Path <http://docs.oracle.com/javase/7/docs/api/java/nio/file/Path.html>`_
item referring the specified file.

.. note:: It does not check the file existence.

Whenever the ``fromPath`` argument contains a ``*`` or ``?`` wildcard character it is interpreted as a `glob`_ path matcher.
For example::

    myFileChannel = Channel.fromPath( '/data/big/*.txt' )


This example creates a channel and emits as many ``Path`` items as there are files with ``txt`` extension in the ``/data/big`` folder.

.. tip:: Two asterisks, i.e. ``**``, works like ``*`` but crosses directory boundaries.
  This syntax is generally used for matching complete paths. Curly brackets specify a collection of sub-patterns.

For example::

    files = Channel.fromPath( 'data/**.fa' )
    moreFiles = Channel.fromPath( 'data/**/*.fa' )
    pairFiles = Channel.fromPath( 'data/file_{1,2}.fq' )

The first line returns a channel emitting the files ending with the suffix ``.fa`` in the ``data`` folder `and` recursively
in all its sub-folders. While the second one only emits the files which have the same suffix in `any` sub-folder in the ``data`` path.
Finally the last example emits two files: ``data/file_1.fq`` and ``data/file_2.fq``.

.. note:: As in Linux Bash the ``*`` wildcard does not match against hidden files (i.e. files whose name start with a ``.`` character).

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

.. note:: More than one path or glob pattern can be specified using a list as argument::

      Channel.fromPath( ['/some/path/*.fq', '/other/path/*.fastq'] )

  (requires version 0.31.x or later)

.. _channel-filepairs:

fromFilePairs
-------------

The ``fromFilePairs`` method creates a channel emitting the file pairs matching a `glob`_ pattern provided by the user.
The matching files are emitted as tuples in which the first element is the grouping key of the matching
pair and the second element is the list of files (sorted in lexicographical order). For example::

    Channel
        .fromFilePairs('/my/data/SRR*_{1,2}.fastq')
        .println()

It will produce an output similar to the following::

    [SRR493366, [/my/data/SRR493366_1.fastq, /my/data/SRR493366_2.fastq]]
    [SRR493367, [/my/data/SRR493367_1.fastq, /my/data/SRR493367_2.fastq]]
    [SRR493368, [/my/data/SRR493368_1.fastq, /my/data/SRR493368_2.fastq]]
    [SRR493369, [/my/data/SRR493369_1.fastq, /my/data/SRR493369_2.fastq]]
    [SRR493370, [/my/data/SRR493370_1.fastq, /my/data/SRR493370_2.fastq]]
    [SRR493371, [/my/data/SRR493371_1.fastq, /my/data/SRR493371_2.fastq]]


.. note::
    The glob pattern must contain at least a star wildcard character.

Alternatively it is possible to implement a custom file pair grouping strategy providing a closure which,
given the current file as parameter, returns the grouping key.
For example::

    Channel
        .fromFilePairs('/some/data/*', size: -1) { file -> file.extension }
        .println { ext, files -> "Files with the extension $ext are $files" }


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

.. note:: More than one glob pattern can be specified using a list as argument::

      Channel.fromFilePairs( ['/some/data/SRR*_{1,2}.fastq', '/other/data/QFF*_{1,2}.fastq'] )

  (requires version 0.31.x or later)


.. _channel-watch:

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

You can specified more than one of these events by using a comma separated string as shown below::

     Channel
        .watchPath( '/path/*.fa', 'create,modify' )
        .subscribe { println "File created or modified: $it" }


.. warning:: The ``watchPath`` factory waits endlessly for files that match the specified pattern and event(s).
  Thus, whenever you use it in your script, the resulting pipeline will never finish.

See also: `fromPath`_ factory method.


.. _channel-empty:

empty
-----

The ``empty`` factory method, by definition, creates a channel that doesn't emit any value.

See also: :ref:`operator-ifempty` and :ref:`operator-close` operators.


Binding values
==============

Since in `Nextflow` channels are implemented using `dataflow` variables or queues. Thus sending a message
is equivalent to `bind` a value to object representing the communication channel.

bind( )
-------

Channel objects provide a `bind( )` method which is the basic operation to send a message over the channel.
For example::

    myChannel = Channel.create()
    myChannel.bind( 'Hello world' )


operator <<
-----------

The operator ``<<`` is just a syntax sugar for the `bind( )` method. Thus, the following example produce
an identical result as the previous one::

    myChannel = Channel.create()
    myChannel << 'Hello world'



Observing events
================


.. _channel-subscribe:

subscribe( )
------------

The ``subscribe( )`` method permits to execute a user define function each time a new value is emitted by the source channel.

The emitted value is passed implicitly to the specified function. For example::

    // define a channel emitting three values
    source = Channel.from ( 'alpha', 'beta', 'delta' )

    // subscribe a function to the channel printing the emitted values
    source.subscribe {  println "Got: $it"  }

::

    Got: alpha
    Got: beta
    Got: delta


.. note:: Formally the user defined function is a ``Closure`` as defined by the Groovy programming language on which
  the `Nextflow` scripts are based on.

If needed the closure parameter can be defined explicitly, using a name other than ``it`` and, optionally,
specifying the expected value type, as showed in the following example::

    Channel
        .from( 'alpha', 'beta', 'lambda' )
        .subscribe { String str ->
            println "Got: ${str}; len: ${str.size()}"
         }

::

    Got: alpha; len: 5
    Got: beta; len: 4
    Got: lambda; len: 6

Read :ref:`script-closure` paragraph to learn more about `closure` feature.


onNext, onComplete, and onError
-------------------------------

The ``subscribe()`` method may accept one or more of the following event handlers:

* ``onNext``: registers a function that is invoked whenever the channel emits a value.
  This is the same as using the ``subscribe( )`` with a `plain` closure as describe in the examples above.

* ``onComplete``: registers a function that is invoked after the `last` value is emitted by the channel.

* ``onError``: registers a function that it is invoked when an exception is raised while handling the
  ``onNext`` event. It will not make further calls to ``onNext`` or ``onComplete``.
  The ``onError`` method takes as its parameter the ``Throwable`` that caused the error.


For example::

    Channel
        .from( 1, 2, 3 )
        .subscribe onNext: { println it }, onComplete: { println 'Done.' }

::

    1
    2
    3
    Done.


.. Special messages
.. STOP
.. VOID



.. _glob: http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob

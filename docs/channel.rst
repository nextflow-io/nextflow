.. _channel-page:

**********
Channels
**********

Nextflow is based on the Dataflow programming model in which processes communicate through channels.

A `channel` is a non-blocking unidirectional FIFO queue which connects two processes.

The channel has the property that sending a message is an `asynchronous` operation which completes immediately,
without any interaction with or waiting for the receiving process. While receiving data is blocking
operation which stops the receiving process until the message has arrived.



.. _channel-factory:

Channel factory
===============

Channels may be created implicitly by the process output(s) declaration or explicitly using the following channel
factory methods.

The available factory methods are:

* `create`_
* `from`_
* `just`_
* `fromPath`_
* `watchPath`_


.. _channel-create:

create
---------

Creates a new `channel` by using the method ``create`` method, as showed below::

    channelObj = Channel.create()


.. _channel-from:

from
-------

The method allows you to create a channel emitting any sequence of values that are specified as the method argument,
For example::

    ch = Channel.from( 1, 3, 5, 7 )
    ch.subscribe { println "value: $it" }

The first line in this example creates variable ``ch`` which holds a channel object. This channel emits the values
specified as parameter in the ``from`` method. Thus the second line will print the following::

    value: 1
    value: 3
    value: 5
    value: 7


The following example shows how create channel from a `range` of numbers or strings::

    zeroToNine = Channel.from( 0..9 )
    strings = Channel.from( 'A'..'Z' )



.. note:: Note that when the ``from`` argument is an object implement the (Java)
  `Collection <http://docs.oracle.com/javase/8/docs/api/java/util/Collection.html>`_ interface, the resulting channel
  emits the collections entries as individual emissions.

Thus, the following two declarations produce an identical result since in the second one the values are specified
as a list object::

    Channel.from( 1, 3, 5, 7, 9 )
    Channel.from( [1, 3, 5, 7, 9] )


But when more than an argument is provided, they are always managed as `single` emission. Thus, the the following example
creates a channel emitting three, entries each of which is a list containing two elements::

    Channel.from( [1, 2], [5,6], [7,9] )



.. _channel-just:

just
-------

This method create a dataflow `variable`, that is a channel to which is possible to bind at most one entry. An optional,
not ``null`` value can be specified as parameters, which is bound to the newly created channel. For example::


    expl1 = Channel.just()
    expl2 = Channel.just( 'Hello there' )
    expl3 = Channel.just( [1,2,3,4,5] )


The first example creates am 'empty' variable. The second, creates a channel and bind a string to it. Finally, the last one
creates a channel and bind a list object to it.

.. _channel-path:

fromPath
--------

You can create a channel emitting one or more file paths by using the ``fromPath`` method and specifying a path string
as argument. For example::

    myFileChannel = Channel.fromPath( '/data/some/bigfile.txt' )

The above line creates a channel and bind to it an item of type
`Path <http://docs.oracle.com/javase/8/docs/api/java/nio/file/Path.html>`_ referring the specified file.

.. note:: It does not check the file existence.

Whenever the ``fromPath`` arguments contains a ``*`` or ``?`` wildcard characters it is interpreted as a `glob` path matcher.
For example::

    myFileChannel = Channel.fromPath( '/data/big/*.txt' )


Creates a channel and sends over it as many ``Path`` items as many are the files with ``txt`` extension in the ``/data/big`` folder.

.. tip:: Two asterisks, i.e. ``**``, works like ``*`` but crosses directory boundaries.
  This syntax is generally used for matching complete paths.

For example::

    files = Channel.fromPath( 'data/**.fa' )
    moreFiles = Channel.fromPath( 'data/**/.fa' )

The first returns a channel emitting the files ending with the suffix ``.fa`` in the ``data`` folder `and` recursively
in all the  its sub-folders. While the second only the files with the same suffix in `any` sub-folder in the ``data`` path.

.. note:: As in Linux BASH the ``*`` wildcard does not match against hidden files (i.e. files which name starts a ``.`` character).

In order the include hidden files you need to start your pattern with a period character. For example::

    expl1 = Channel.fromPath( '/path/.*' )
    expl2 = Channel.fromPath( '/path/.*.fa' )
    expl3 = Channel.fromPath( '/path/{.*,*}' )


The first example returns all hidden files in the path the specified path. The second one, returns all hidden files
ending with ``.fa``. Finally, the last example returns all files (hidden and non-hidden) in that path


Learn more about `glob` patterns at `this link <http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob>`_

.. _channel-watch:

watchPath
-----------

The ``watchPath`` factory method watches a folder for one or more files matching a specified pattern. As soon as a
there is a file that meets the specified condition, this file is emitted over the channel returned by the ``watchPath`` method.
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

You can specified one more of these events by using a comma separated string, as shown below::

     Channel
        .watchPath( '/path/*.fa', 'create,modify' )
        .subscribe { println "File created or modified: $it" }


.. warning:: The ``watchPath`` factory wait endlessly for files that matches the specified pattern and event(s).
  Thus, whenever you use it in your script, the resulting pipeline will never finish.

See also: `fromPath`_ factory method

Learn more about `glob` patterns at `this link <http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob>`_


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
=================


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


onNext, onCompleted, and onError
--------------------------------

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




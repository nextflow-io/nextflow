.. _channel-page:

**********
Channels
**********

Nextflow is based on the Dataflow programming model in which processes communicate through channels.

A `channel` is a non-blocking unidirectional FIFO queue which connects two processes.

The channel has the property that sending a message is an `asynchronous` operation which completes immediately,
without any interaction with or waiting for the receiving process. While receiving data is blocking
operation which stops the receiving process until the message has arrived.



Channel factory
===============

Channels may be created implicitly by the process output(s) declaration or explicitly using the following channel
factory methods.


create( )
---------

Creates a new `channel` by using the method ``create( )`` method, as showed below::

    myChannel = Channel.create()


from( )
-------

You can create transform any object that supports ``Iterable`` into a channel by using the ``from( )`` method.
Each iterable item is send over the newly created channel. For example::

    // creating a channel by a generic iterable object
    myChannel = Channel.from( myIterableObj )

    // creating a channel by a list
    myChannel = Channel.from( [1,2,3,4,5] )


Square brackets can be omitted when passing a list a parameter, thus the following example it is identical to the previous one::

     // square brackets can be omitted
     myChannel = Channel.from( 1,2,3,4,5 )



just( )
-------

This method create a dataflow `variable`, that is a channel to which is possible to bind at most entry. An optional,
not ``null`` value can be specified as parameters, which is bound to the newly created channel. For example::

    // creates am 'empty' variable
    myChannel = Channel.just()

    // creates a channel and bind a string to it
    myChannel = Channel.just( 'Hello there' )


    // creates a channel and bind a list object to it
    myChannel = Channel.just( [1,2,3,4,5] )



path( )
--------

You can create a channel emitting one or more file paths by using the ``path( )`` method and specifying a path matcher
as argument. For example::

    myFileChannel = Channel.path( '/data/some/bigfile.txt' )

The above line creates a channel and bind to it an item of type ``Path`` referring the specified file.

.. note:: It does not check the file existence.

Whenever the ``path`` arguments contains a ``*`` or ``?`` wildcard characters it is interpreted as a `glob` path matcher.
For example::

    myFileChannel = Channel.path( '/data/big/*.txt' )


Creates a channel and sends over it as many ``Path`` items as many are the files with ``txt`` extension in the ``/data/big`` folder.

.. note:: As in Linux BASH the ``*`` wildcard does not match against hidden files (i.e. files which name starts a ``.`` character).

In order the include hidden files you need to start your pattern with a period character. For example::

    // returns all hidden files in the path
    myFileChannel = Channel.path( '/path/.*' )

    // returns all hidden files ending with .fa
    myFileChannel = Channel.path( '/path/.*.fa' )

    // returns all files (hidden and non-hidden)
    myFileChannel = Channel.path( '/path/{.*,*}' )



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




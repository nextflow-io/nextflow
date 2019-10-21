.. _dsl2-page:

******
DSL 2
******

Nextflow implements an experimental syntax that implements new features and enhancements that
simplifies the implementation of data analysis applications.

To enable this feature you need to defined the following directive at the beginning of
your workflow script::

    nextflow.preview.dsl=2


.. warning:: THIS IS AN EXPERIMENT FEATURE UNDER DEVELOPMENT. SYNTAX MAY CHANGE IN FUTURE RELEASE.


Function
========

Nextflow allows the definition of custom function in the workflow script using the following syntax::

    def <function name> ( arg1, arg, .. ) {
        <function body>
    }

For example::

    def foo() {
        'Hello world'
    }

    def bar(alpha, omega) {
        alpha + omega
    }


The above snippet defines two simple functions, that can be invoked in the workflow script as ``foo()`` which
returns the ``Hello world`` string and ``bar(10,20)`` which return the sum of two parameters.

.. tip:: Functions implicitly return the result of the last evaluated statement.

The keyword ``return`` can be used to explicitly exit from a function returning the specified value.
for example::

    def fib( x ) {
        if( x <= 1 )
            return x
        else
            fib(x-1) + fib(x-2)
    }

Process
=======

Process definition
------------------

The new DSL separates the definition of a process by its invocation. The process definition follows the usual
for syntax as described in the :ref:`process documentation <process-page>`. The only difference is that the
``from`` and ``into`` channel declaration has to be omitted.

Then a process can be invoked as a function in the ``workflow`` scope, passing the expected
input channels as parameters as it if were a custom function. For example::

    nextflow.preview.dsl=2

    process foo {
        output:
          path 'foo.txt'
        script:
          """
          your_command > foo.txt
          """
    }

     process bar {
        input:
          path x
        output:
          file 'bar.txt'
        script:
          """
          another_command $x > bar.txt
          """
    }

    workflow {
        data = Channel.fromPath('/some/path/*.txt')
        foo()
        bar(data)
    }


.. warning:: A process component can be invoked only once in the same workflow context.


Process composition
-------------------

Processes having matching *input-output* declaration can be composed so that the output
of the first process is passed as input to the following process. Taking in consideration
the previous process definition, it's possible to write the following::

    workflow {
        bar(foo())
    }


Process outputs
---------------

A process output can also be accessed using the ``out`` attribute for the respective
process object. For example::

    workflow {
        foo()
        bar(foo.out)
        bar.out.view()
    }


When a process defines two or more output channels, each of them can be accessed
using the array element operator e.g. ``out[0]``, ``out[1]``, etc. or using
*named outputs* (see below).

Process named output
--------------------

The process output definition allow the use of the ``emit`` option to define a name identifier
that can be used to reference the channel in external scope. For example::

    process foo {
      output:
        path '*.bam' emit: samples_bam

      '''
      your_command --here
      '''
    }
    
    workflow {
        foo()
        foo.out.samples_bam.view()
    }


Workflow
========

Workflow definition
--------------------

The ``workflow`` keyword allows the definition of sub-workflow components that enclose the
invocation of one or more processes and operators::

    workflow my_pipeline {
        foo()
        bar( foo.out.collect() )
    }


For example, the above snippet defines a workflow component, named ``my_pipeline``, that can be invoked from
another workflow component definition as any other function or process i.e. ``my_pipeline()``.


Workflow parameters
---------------------

A workflow component can access any variable and parameter defined in the outer scope::

        params.data = '/some/data/file'

        workflow my_pipeline {
            if( params.data )
                bar(params.data)
            else
                bar(foo())
        }


Workflow inputs
---------------

A workflow component can declare one or more input channels using the ``get`` keyword. For example::

        workflow my_pipeline {
            get: data
            main:
            foo(data)
            bar(foo.out)
        }

.. warning:: When the ``get`` is used the beginning of the workflow body needs to be identified with the
  ``main`` keyword.

Then, the input can be specified a argument on the workflow invocation statement::

    workflow {
        my_pipeline( Channel.from('/some/data') )
    }

.. note:: Workflow inputs are by definition *channel* data structure. If a basic data type is provided
  instead, ie. number, string, list, etc. it's implicitly converted to a :ref:`channel value <channel-type-value>` (ie. non-consumable).


Workflow outputs
----------------

A workflow component can declare one or more out channels using the ``emit`` keyword. For example::

        workflow my_pipeline {
            main:
              foo(data)
              bar(foo.out)
            emit:
              bar.out
        }

Then, the result of the ``my_pipeline`` execution can be accessed using the ``out`` property ie.
``my_pipeline.out``. When is declared more than one output channels, use the array bracket notation
to access each output component as described for the `Process outputs`_ definition.

Alternatively, the output channel can be accessed using the identifier name to which it's assigned
in the ``emit`` declaration::

         workflow my_pipeline {
            main:
              foo(data)
              bar(foo.out)
            emit:
              my_data = bar.out
        }

Then, the result of the above snippet can accessed using the ``my_pipeline.out.my_data``.


Implicit workflow
-----------------

A workflow definition which does not declare any name is assumed to be the main workflow and it's
implicitly executed. Therefore it's the entry point of the workflow application.

.. note:: Implicit workflow definition is ignored when a script is included as module. This
  allow the writing a workflow script that can be used either as a library module and as
  application script. 

.. tip:: An alternative workflow entry can be specifying using the ``-entry`` command line option.


Workflow publish
----------------

The ``publish`` clause in the workflow declaration allow the definition of one or more output channels
whose content needs to be copied to a storage location of your choice. For example::

    workflow {
        main:
          foo()
          bar()
        publish:
          foo.out to: '/some/data/dir'
          bar.out to: '/some/other/dir'
    }

.. note:: The ``publish`` clause is only allowed in a *implicit workflow* definition.

The ``publish`` clause works in similar manner to the process :ref:`process-publishDir` directive,
however the former allow you to specify a different target directory for each published channel.

.. tip:: In the ``publish`` clause can be specified the same options as for the :ref:`process-publishDir`
  directive i.e. ``mode``, ``overwrite``, ``enabled``, etc.


Workflow composition
--------------------

Workflow defined in your script or import by a module inclusion can be invoked and composed
as any other process in your application.

::

    workflow flow1 {
        get: data
        main:
            foo(data)
            bar(foo.out)
        emit:
            bar.out
    }

    workflow flow2 {
        get: data
        main:
            foo(data)
            baz(foo.out)
        emit:
            baz.out
    }

    workflow {
        get: data
        main:
          flow1(data)
          flow2(flow1.out)
    }


.. note::
    Nested workflow execution determines an implicit scope. Therefore the same process can be
    invoked in two different workflow scopes, like for example ``foo`` in the above snippet that
    is used either in ``flow1`` and ``flow2``. The workflow execution path along with the
    process names defines the process *fully qualified name* that is used to distinguish the
    two different process invocation i.e. ``flow1:foo`` and ``flow2:foo`` in the above example.

.. tip::
    The process fully qualified name can be used as a valid :ref:`process selector <config-process-selectors>` in the
    ``nextflow.config`` file and it has priority over the process simple name.


Modules
=======

The new DSL allows the definition module scripts that
can be included and shared across workflow applications.

A module can contain the definition of function, process and workflow definitions
as described in the above sections.

Modules include
---------------

A module script can be included from another Nextflow script using the ``include`` keyword.
Then it's possible to reference the components (eg. functions, processes and workflow ) defined in the module
from the including script.

For example::

    nextflow.preview.dsl=2
    include './modules/my-module'

    workflow {
        data = Channel.fromPath('/some/data/*.txt')
        my_pipeline(data)
    }


Nextflow implicitly looks for the script file ``modules/my-module.nf`` resolving the path
against the *including* script location.

.. note:: Relative paths must begin with the ``./`` prefix.

Selective inclusion
-------------------

The module inclusion implicitly imports all the components defined in the module script.
It's possible to selective include only a specific component by its name using the
inclusion extended syntax as shown below::

    nextflow.preview.dsl=2
    include my_pipeline from './modules/my-module'

    workflow {
        data = Channel.fromPath('/some/data/*.txt')
        my_pipeline(data)
    }


Module aliases
--------------

When including a module component it's possible to specify a name *alias*.
This allows the inclusion and the invocation of the same component multiple times
in your script using different names. For example::

    nextflow.preview.dsl=2

    include foo from './modules/my-module'
    include foo as bar from './modules/my-module'

    workflow {
        foo(some_data)
        bar(other_data)
    }

.. warning:: The process configuration defined in the ``nextflow.config`` file will
  be resolved against the alias name. Following the above example the process ``bar``
  ignores any process specific settings eventually defined for the process with name ``foo``.


Module parameters
-----------------

A module script can define one or more parameters as any other Nextflow script::

    params.foo = 'hello'
    params.bar = 'world'

    def sayHello() {
        "$params.foo $params.bar"
    }


Then, parameters can be specified when the module is imported with the ``include`` statement::


    nextflow.preview.dsl=2

    include './modules/my-module' params(foo: 'Hola', bar: 'mundo')


.. tip::
    All parameters in the current context can be passed to a module using the ``params(params)`` idiom.


Channel forking
===============

Using the new DSL Nextflow channels are automatically forked when connecting two or more consumers.

For example::

    Channel
        .from('Hello','Hola','Ciao')
        .set{ cheers }

    cheers
        .map{ it.toUpperCase() }
        .view()

    cheers
        .map{ it.reverse() }
        .view()


The same is valid for the result (channel) of a process execution. Therefore a process output can be used by
two or more processes without the need to fork them using the :ref:`operator-into` operator, making the
writing of workflow script much fluent and readable.


Pipes
=====

The *pipe* operator
-------------------

Nextflow processes and operators can be composed using the ``|`` *pipe* operator. For example::

    process foo {
        input: val data
        output: val result
        exec:
        result = "$data world"
    }

    workflow {
       Channel.from('Hello','Hola','Ciao') | foo | map { it.toUpperCase() } | view
    }



The above snippet defines a process named ``foo`` then invoke it passing the content of the
``data`` channel. The result is piped to the :ref:`operator-map` operator which converts each string
to uppercase and finally, the last :ref:`operator-view` prints it.


The *and* operator
------------------

The ``&`` *and* operator allow the feed of two or more processes with the content of the same
channel(s) e.g.::

    process foo {
      input: val data
      output: val result
      exec:
        result = "$data world"
    }

    process bar {
        input: val data
        output: val result
        exec:
          result = data.toUpperCase()
    }

    workflow {
       Channel.from('Hello') | map { it.reverse() } | (foo & bar) | mix | view
    }


In the above snippet the channel emitting the ``Hello`` is piped with the :ref:`operator-map`
which reverse the string value. Then, the result is passed to either ``foo`` and ``bar``
processes which are executed in parallel. The result is pair of channels which content
is merged into a single channel using the :ref:`operator-mix`. Finally the result is printed
using the :ref:`operator-view`.

.. tip:: The break-line operator ``\`` can be use to split long pipes concatenation
  over multiple lines.


The above snippet can be written as shown below::

    workflow {
       Channel.from('Hello') \
         | map { it.reverse() } \
         | (foo & bar) \
         | mix \
         | view
    }


Deprecated methods and operators
================================

The following methods are not allowed any more when using Nextflow DSL 2:

* :ref:`channel-create`
* :ref:`channel-bind1`
* :ref:`channel-bind2`
* :ref:`operator-choice`
* :ref:`operator-close`
* :ref:`operator-countby`
* route
* :ref:`operator-separate`
* :ref:`operator-into`
* :ref:`operator-merge`

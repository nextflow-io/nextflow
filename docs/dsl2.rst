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


The above snippet defines two simple functions, that can be invoked in the workflow script as `foo()` which
returns the ``Hello world`` string and ``bar(10,20)`` which return the sum of two parameters.

.. tip:: Functions implicitly return the result of the function last evaluated statement.

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

Then processes can be invoked as a function, passing the expected input channels as parameters.

For example::

    nextflow.preview.dsl=2

    process foo {
        output:
          file 'foo.txt'
        script:
          """
          your_command > foo.txt
          """
    }

     process bar {
        input:
          file x
        output:
          file 'bar.txt'
        script:
          """
          another_command $x > bar.txt
          """
    }

    data = Channel.fromPath('/some/path/*.txt')
    foo()
    bar(data)


Process invocation
------------------


Process composition
-------------------

Processes having matching input-output declaration can be composed so that the output
of the first process is passed as input to the following process. Take in consideration
the previous process definition, it's possible to write the following::

    bar(foo())

Process outputs
---------------

A process output can also be accessed using the ``output`` attribute for the respective
process object. For example::

    foo()
    bar( foo.output )
    bar.output.println()


When a process defines two or more output channels, each of them can be accessed
using the array element operator e.g. ``output[0]``, etc or using the ``first``, ``second``, etc
sub-properties e.g. ``output.first``.

.. Workflow
.. ========
..
.. Workflow definition
.. --------------------
..
.. The ``workflow`` keyword allows the definition of sub-workflow components that enclose the
.. invocation of two or more processes or operators. For example::
..
..     workflow my_pipeline {
..         foo()
..         bar( foo.output.collect() )
..     }
..
..
.. Once defined it can be invoked from another (sub) workflow component definition.
..
.. Workflow parameters
.. -------------------
..
.. A workflow component can be define one or more parameter in a similar manner as for a function
.. definition. For example::
..
..         workflow my_pipeline( data )  {
..             foo()
..             bar( data.mix( foo.output ) )
..         }
..
.. The result channel of the last evaluated process is implicitly returned as the workflow output.
..
..
.. Main workflow
.. -------------
..
.. A workflow definition which does not define any name is assumed to be the main workflow and it's
.. implicitly executed. Therefore it's the entry point of the workflow application.

Modules
=======

The new DSL allows the definition module scripts that
can be included and shared across workflow applications.

A module can contain the definition of function, process and workflow definitions
as described above.

Modules include
---------------

A module script can be included from another Nextflow script using the ``include`` keyword.
Then it's possible to reference of components (eg. functions, processes and workflow ) defined in the module
from the importing script.

For example::

    nextflow.preview.dsl=2
    include 'modules/libx'

    data = Channel.fromPath('/some/data/*.txt')
    my_pipeline(data)

Nextflow implicitly looks for the module script ``modules/libx.nf`` resolving the path
against the main script location.

Selective inclusion
-------------------

The module inclusion implicitly imports all the components defined in the module script.
It's possible to selective include only a specific component by its name using the
inclusion extended syntax as shown below::

    nextflow.preview.dsl=2
    include my_pipeline from 'modules/libx'

    data = Channel.fromPath('/some/data/*.txt')
    my_pipeline(data)

The module component can be included using a name alias as shown below::


    nextflow.preview.dsl=2
    include my_pipeline as my_tool from 'modules/libx'

    data = Channel.fromPath('/some/data/*.txt')
    my_tool(data)

Module aliases
--------------

When including a module component it's possible to specify a name alias.
This allows the import and the invocation of the same component multiple times
in your script using different names. For example::

    nextflow.preview.dsl=2

    include foo from 'modules/my-library'
    include for as bar from 'modules/my-library'

    foo(some_data)
    bar(other_data)


Module parameters
-----------------

A module script can define one or more parameters as any other Nextflow script.::

    params.foo = 'hello'
    params.bar = 'world'

    def sayHello() {
        "$params.foo $params.bar"
    }


Then, parameters can be specified when the module is imported with the ``include`` statement::


    nextflow.preview.dsl=2

    include 'modules/library.nf' params(foo: 'Hola', bar: 'mundo')



Channel forking
===============

Using the new DSL Nextflow channels are automatically forked when connecting two or more consumers.
This means that, for example, a process output can be used by two or more processes without the
need to fork them using the :ref:`operator-into` operator, making the writing of workflow script
much fluent and readable.

Pipes
=====

Nextflow processes and operators can be composed using the ``|`` *pipe* operator. For example::

      process foo {
          input: val data
          output: val result
          exec:
            result = "$data mundo"
      }

      Channel.from('Hello','world') | foo


The above snippet defines a process named ``foo`` then invoke it passing the content of the
``data`` channel.

The ``&`` *and* operator allow the feed of two or more processes with the content of the same
channel e.g.::

    process foo {
      input: val data
      output: val result
      exec:
        result = "$data mundo"
    }

    process bar {
        input: val data
        output: val result
        exec:
          result = data.toUpperCase()
    }


    Channel.from('Hello') | map { it.reverse() } | (foo & bar)


Deprecated methods and operators
================================

The following methods are not allowed any more when using Nextflow DSL 2:

* :ref:`channel-create`
* :ref:`channel-bind1`
* :ref:`channel-bind2`
* :ref:`operator-close`
* :ref:`operator-countby`
* :ref:`operator-route`
* :ref:`operator-separate`
* :ref:`operator-into`
* :ref:`operator-merge`

.. _dsl2-page:

******
DSL 2
******

Nextflow provides a syntax extension that implements that allow the definition of module libraries and
simplifies the writing of complex data analysis pipelines.

To enable this feature you need to defined the following directive at the beginning of
your workflow script::

    nextflow.enable.dsl=2


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

The new DSL separates the definition of a process from its invocation. The process definition follows the usual 
syntax as described in the :ref:`process documentation <process-page>`. The only difference is that the
``from`` and ``into`` channel declaration has to be omitted.

Then a process can be invoked as a function in the ``workflow`` scope, passing the expected
input channels as parameters as it if were a custom function. For example::

    nextflow.enable.dsl=2

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
          path 'bar.txt'
        script:
          """
          another_command $x > bar.txt
          """
    }

    workflow {
        data = channel.fromPath('/some/path/*.txt')
        foo()
        bar(data)
    }


.. warning::
  A process component can be invoked only once in the same workflow context.


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

The process output definition allows the use of the ``emit`` option to define a name identifier
that can be used to reference the channel in the external scope. For example::

    process foo {
      output:
        path '*.bam', emit: samples_bam

      '''
      your_command --here
      '''
    }
    
    workflow {
        foo()
        foo.out.samples_bam.view()
    }
    
Process named stdout
--------------------

The process can name stdout using the ``emit`` option::

    process sayHello {
        input:
            val cheers
        output:
            stdout emit: verbiage
        script:
        """
        echo -n $cheers
        """
    }

    workflow {
        things = channel.of('Hello world!', 'Yo, dude!', 'Duck!')
        sayHello(things)
        sayHello.out.verbiage.view()
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

A workflow component can declare one or more input channels using the ``take`` keyword. For example::

        workflow my_pipeline {
            take: data
            main:
            foo(data)
            bar(foo.out)
        }

.. warning::
  When the ``take`` keyword is used, the beginning of the workflow body needs to be identified with the
  ``main`` keyword.

Then, the input can be specified as an argument in the workflow invocation statement::

    workflow {
        my_pipeline( channel.from('/some/data') )
    }

.. note::
  Workflow inputs are by definition *channel* data structures. If a basic data type is provided
  instead, ie. number, string, list, etc. it's implicitly converted to a :ref:`channel value <channel-type-value>`
  (ie. non-consumable).


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
``my_pipeline.out``. When there are multiple output channels declared, use the array bracket notation
to access each output component as described for the `Process outputs`_ definition.

Alternatively, the output channel can be accessed using the identifier name which it's assigned to
in the ``emit`` declaration::

     workflow my_pipeline {
        main:
          foo(data)
          bar(foo.out)
        emit:
          my_data = bar.out
     }

Then, the result of the above snippet can accessed using ``my_pipeline.out.my_data``.


Implicit workflow
-----------------

A workflow definition which does not declare any name is assumed to be the main workflow and it's
implicitly executed. Therefore it's the entry point of the workflow application.

.. note::
  Implicit workflow definition is ignored when a script is included as module. This
  allows the writing of a workflow script that can be used either as a library module and as
  application script. 

.. tip::
  An alternative workflow entry can be specified using the ``-entry`` command line option.


Workflow composition
--------------------

Workflows defined in your script or imported by a module inclusion can be invoked and composed
as any other process in your application.

::

    workflow flow1 {
        take: data
        main:
            foo(data)
            bar(foo.out)
        emit:
            bar.out
    }

    workflow flow2 {
        take: data
        main:
            foo(data)
            baz(foo.out)
        emit:
            baz.out
    }

    workflow {
        take: data
        main:
          flow1(data)
          flow2(flow1.out)
    }


.. note::
    Nested workflow execution determines an implicit scope. Therefore the same process can be
    invoked in two different workflow scopes, like for example ``foo`` in the above snippet that
    is used either in ``flow1`` and ``flow2``. The workflow execution path along with the
    process names defines the process *fully qualified name* that is used to distinguish the
    two different process invocations i.e. ``flow1:foo`` and ``flow2:foo`` in the above example.

.. tip::
    The process fully qualified name can be used as a valid :ref:`process selector <config-process-selectors>` in the
    ``nextflow.config`` file and it has priority over the process simple name.


Modules
=======

The new DSL allows the definition module scripts that
can be included and shared across workflow applications.

A module can contain the definition of a function, process and workflow definitions
as described in the above sections.

Modules include
---------------

A component defined in a module script can be imported into another Nextflow script using the ``include`` keyword.

For example::

    include { foo } from './some/module'

    workflow {
        data = channel.fromPath('/some/data/*.txt')
        foo(data)
    }

The above snippets includes a process with name ``foo`` defined in the module script in the main
execution context, as such it can be invoked in the ``workflow`` scope.

Nextflow implicitly looks for the script file ``./some/module.nf`` resolving the path
against the *including* script location.

.. note:: Relative paths must begin with the ``./`` prefix.

Multiple inclusions
-------------------

A Nextflow script allows the inclusion of any number of modules. When multiple
components need to be included from the same module script, the component names can be
specified in the same inclusion using the curly brackets notation as shown below::

    include { foo; bar } from './some/module'

    workflow {
        data = channel.fromPath('/some/data/*.txt')
        foo(data)
        bar(data)
    }


Module aliases
--------------

When including a module component it's possible to specify a name *alias*.
This allows the inclusion and the invocation of the same component multiple times
in your script using different names. For example::

    include { foo } from './some/module'
    include { foo as bar } from './other/module'

    workflow {
        foo(some_data)
        bar(other_data)
    }

The same is possible when including multiple components from the same module script as shown below::

    include { foo; foo as bar } from './some/module'

    workflow {
        foo(some_data)
        bar(other_data)
    }


Module parameters
-----------------

A module script can define one or more parameters using the same syntax of a Nextflow workflow script::

    params.foo = 'Hello'
    params.bar = 'world!'

    def sayHello() {
        println "$params.foo $params.bar"
    }


Then, parameters are inherited from the including context. For example::

    params.foo = 'Hola'
    params.bar = 'Mundo'

    include {sayHello} from './some/module'

    workflow {
        sayHello()
    }

The above snippet prints::

    Hola Mundo


.. note::
  The module inherits the parameters define *before* the ``include`` statement, therefore any further
  parameter set later is ignored.

.. tip::
  Define all pipeline parameters at the beginning of the script *before*
  any ``include`` declaration.

The option ``addParams`` can be used to extend the module parameters without affecting the external
scope. For example::


    include {sayHello} from './some/module' addParams(foo: 'Ciao')

    workflow {
        sayHello()
    }


The above snippet prints::

    Ciao world!


Finally the include option ``params`` allows the specification of one or more parameters without
inheriting any value from the external environment. 

.. _module-templates:

Module templates
-----------------
The module script can be defined in an external :ref:`template <process-template>` file. With DSL2 the template file
can be placed under the ``templates`` directory where the module script is located.

For example, let's suppose to have a project L with a module script defining 2 processes (P1 and P2) and both use templates.
The template files can be made available under the local ``templates`` directory::

	Project L
		|-myModules.nf
		|-templates
			|-P1-template.sh
			|-P2-template.sh

Then, we have a second project A with a workflow that includes P1 and P2::

	Pipeline A
		|-main.nf

Finally, we have a third project B with a workflow that includes again P1 and P2::

	Pipeline B
		|-main.nf

With the possibility to keep the template files inside the project L, A and B can use the modules defined in L without any changes.
A future prject C would do the same, just cloning L (if not available on the system) and including its module script.

Beside promoting sharing modules across pipelines, there are several advantages in keeping the module template under the script path::
1 - module components are *self-contained*
2 - module components can be tested independently from the pipeline(s) importing them
3 - it is possible to create libraries of module components

Ultimately, having multiple template locations allows a more structured organization within the same project. If a project
has several module components, and all them use templates, the project could group module scripts and their templates as needed. For example::

	baseDir
		|-main.nf
		|-Phase0-Modules
			|-mymodules1.nf
			|-mymodules2.nf
			|-templates
				|-P1-template.sh
				|-P2-template.sh
		|-Phase1-Modules
			|-mymodules3.nf
			|-mymodules4.nf
			|-templates
				|-P3-template.sh
				|-P4-template.sh
		|-Phase2-Modules
			|-mymodules5.nf
			|-mymodules6.nf
			|-templates
				|-P5-template.sh
				|-P6-template.sh
				|-P7-template.sh


Channel forking
===============

Using the new DSL, Nextflow channels are automatically forked when connecting two or more consumers.

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
writing of workflow scripts more fluent and readable.


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
       channel.from('Hello','Hola','Ciao') | foo | map { it.toUpperCase() } | view
    }



The above snippet defines a process named ``foo`` then invoke it passing the content of the
``data`` channel. The result is piped to the :ref:`operator-map` operator which converts each string
to uppercase and finally, the last :ref:`operator-view` operator prints it.


The *and* operator
------------------

The ``&`` *and* operator allows feeding of two or more processes with the content of the same
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
       channel.from('Hello') | map { it.reverse() } | (foo & bar) | mix | view
    }


In the above snippet the channel emitting the ``Hello`` is piped with the :ref:`operator-map`
which reverses the string value. Then, the result is passed to either ``foo`` and ``bar``
processes which are executed in parallel. The result is pair of channels whose content
is merged into a single channel using the :ref:`operator-mix` operator. Finally the result is printed
using the :ref:`operator-view` operator.

.. tip::
  The break-line operator ``\`` can be used to split long pipe concatenations
  over multiple lines.


The above snippet can be written as shown below::

    workflow {
       channel.from('Hello') \
         | map { it.reverse() } \
         | (foo & bar) \
         | mix \
         | view
    }



DSL2 migration notes
=====================

* DSL2 final version is activated using the declaration ``nextflow.enable.dsl=2`` in place of ``nextflow.preview.dsl=2``.
* Process inputs of type ``set`` have to be replaced with :ref:`tuple <process-input-tuple>`.
* Process outputs of type ``set`` have to be replaced with :ref:`tuple <process-out-tuple>`.
* Process output option ``mode flatten`` is not available anymore. Replace it using the :ref:`operator-flatten` operator on the corresponding output channel.
* Anonymous and unwrapped includes are not supported anymore. Replace it with a explicit module inclusion. For example::

        include './some/library'
        include bar from './other/library'

        workflow {
          foo()
          bar()
        }

  Should be replaced with::

        include { foo } from './some/library'
        include { bar } from './other/library'

        workflow {
          foo()
          bar()
        }
        
* The use of unqualified value and file elements into input tuples is not allowed anymore. Replace them with a corresponding
  ``val`` or ``path`` qualifier::

        process foo {
        input:
          tuple X, 'some-file.bam'
         script:
           '''
           your_command
           '''
        }

  Use::

        process foo {
        input:
          tuple val(X), path('some-file.bam')
         script:
           '''
           your_command --in $X some-file.bam
           '''
        }


* The use of unqualified value and file elements into output tuples is not allowed anymore. Replace them with a corresponding
  ``val`` or ``path`` qualifier::


        process foo {
        output:
          tuple X, 'some-file.bam'

        script:
           X = 'some value'
           '''
           your_command > some-file.bam
           '''
        }

  Use::

        process foo {
        output:
          tuple val(X), path('some-file.bam')

        script:
           X = 'some value'
           '''
           your_command > some-file.bam
           '''
        }


* Operator :ref:`channel-bind1` has been deprecated by DSL2 syntax
* Operator :ref:`channel-bind2` has been deprecated by DSL2 syntax.
* Operator :ref:`operator-choice` has been deprecated by DSL2 syntax. Use :ref:`operator-branch` instead.
* Operator :ref:`operator-close` has been deprecated by DSL2 syntax.
* Operator :ref:`channel-create` has been deprecated by DSL2 syntax.
* Operator ``countBy`` has been deprecated by DSL2 syntax.
* Operator :ref:`operator-into` has been deprecated by DSL2 syntax since it's not needed anymore.
* Operator ``fork`` has been renamed to :ref:`operator-multimap`.
* Operator ``groupBy`` has been deprecated by DSL2 syntax. Replace it with :ref:`operator-grouptuple`
* Operator ``print`` and ``println`` have been deprecated by DSL2 syntax. Use :ref:`operator-view` instead.
* Operator :ref:`operator-merge` has been deprecated by DSL2 syntax. Use :ref:`operator-join` instead.
* Operator :ref:`operator-separate` has been deprecated by DSL2 syntax.
* Operator :ref:`operator-spread` has been deprecated with DSL2 syntax. Replace it with :ref:`operator-combine`.
* Operator route has been deprecated by DSL2 syntax.


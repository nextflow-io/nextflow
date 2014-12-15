.. _config-page:

*******************
Configuration
*******************

Configuration file
===================

When a pipeline script is launched Nextflow looks for a file named ``nextflow.config`` in the current directory and
in the script base directory (if it is not the same as the current directory). Finally it checks for the file
``$HOME/.nextflow/config``.

When more than one on the above files exist they are merged, so that the settings in the first override the same ones
that may appear in the second one, and so on.

The default config file search mechanism can be extended proving an extra configuration file by using the command line
option ``-c <config file>``.

.. note:: It's worth noting that by doing this, the files ``nextflow.config`` and ``$HOME/.nextflow/config`` are not
  ignored and they are merged as explained above.

.. tip:: If you want to ignore any default configuration files and use only the custom one use the command line option
  ``-C <config file>``.

Config syntax
--------------

A Nextflow configuration file is a simple text file containing a set of properties defined using the syntax::

  name = value

Please note, string values need to be wrapped in quotation characters while numbers and boolean values (``true``, ``false``) do not.
Also note that values are typed, meaning for example that, ``1`` is different from ``'1'``, since the first is interpreted
as the number one, while the latter is interpreted as a string value.


Config Variables
--------------------

Configuration properties can be used as variables in the configuration file itself, by using the usual
``$propertyName`` or ``${expression}`` syntax.


For example::

     propertyOne = 'world'
     anotherProp = "Hello $propertyOne"
     customPath = "$PATH:/my/app/folder"

Please note, the usual rules for :ref:`string-interpolation` are applied, thus a string containing a variable
reference must be wrapped in double-quote chars instead of single-quote chars.

The same mechanism allows you to access environment variables defined in the hosting system. Any variable whose name is
not defined in the Nextflow configuration file(s) is supposed to be a reference to an environment variable with that name.
So, in the above example the property ``customPath`` is defined as the current system ``PATH`` to which
the string ``/my/app/folder`` is appended.

.. warning:: If you are accessing an environment variable that may not exist in the system, your property may contain
    an undefined value. You can avoid this by using a conditional expression in your property definition as shown below.

::

    mySafeProperty = "${MY_FANCY_VARIABLE?:''}"


Config comments
------------------

Configuration files use the same conventions for comments used by the Groovy or Java programming languages. Thus, use ``//`` to comment
a single line or ``/*`` .. ``*/`` to comment a block on multiple lines.


Config scopes
==============

Configuration settings can be organized in different scopes by dot prefixing the property names with a scope
identifier or grouping the properties in the same scope using the curly brackets notation. This is shown in the
following example::

   alpha.x  = 1
   alpha.y  = 'string value..'

   beta {
        p = 2
        q = 'another string ..'
   }



Scope 'env'
-----------

The ``env`` scope allows you to define one or more environment variables that will be exported to the system environment
where pipeline processes need to be executed.

Simply prefix your variable names with the ``env`` scope or surround them by curly brackets, as shown below::

   env.ALPHA = 'some value'
   env.BETA = "$HOME/some/path"

   env {
        DELTA = 'one more'
        GAMMA = "/my/path:$PATH"
   }




Scope 'params'
----------------

The ``params`` scope allows you to define parameters that will be accessible in the pipeline script. Simply prefix the
parameter names with the ``params`` scope or surround them by curly brackets, as shown below::

     params.custom_param = 123
     params.another_param = 'string value .. '

     params {

        alpha_1 = true
        beta_2 = 'another string ..'

     }



.. _config-process:

Scope 'process'
----------------

The ``process`` configuration scope allows you to provide the default configuration for the processes in your pipeline.

You can specify here any property described in the :ref:`process directive<process-directives>` and the executor sections.
For examples::

  process {
    executor='sge'
    queue='long'
    clusterOptions = '-pe smp 10 -l virtual_free=64G,h_rt=30:00:00'
  }


By using this configuration all processes in your pipeline will be executed through the SGE cluster, with the specified
settings.

It is possible to set the properties for a specific process in your pipeline by prefixing the process name with
the symbol ``$`` and using it as special scope identifier. For example::

  process.queue = 'short'
  process.$hello.queue = 'long'


The above configuration example sets the ``queue`` property to ``'short'`` as default value for all processes in your
pipeline, but the process ``hello`` for which the ``queue`` property is set to ``'long'``.

When using the curly brackets notation, the above can be written as shown below::

  process {
    queue = 'short'

    $hello {
        queue = 'long'
    }
  }



Scope 'executor'
------------------

The ``executor`` configuration scope allows you to set the optional executor settings, listed in the following table.

===================== =====================
Name                  Description
===================== =====================
name                  The name of the executor to be used e.g. ``local``, ``sge``, etc
queueSize             The number of tasks the executor will handle in a parallel manner.
pollInterval          Determines how often a poll occurs to check for a process termination.
dumpInterval          Determines how often the executor status is written in the application log file (default: ``5min``)
queueStatInterval     Determines how often the queue status is fetched from the cluster system. This setting is used only by grid executors (default: ``1min``)
exitReadTimeout       Determines how long the executor waits before return an error status when a process is terminated but the `exit` file does not exist or it is empty. This setting is used only by grid executors (default: ``90sec``)
===================== =====================


The executor settings can be defined as shown below::

    executor {
        name = 'sge'
        queueSize = 200
        pollInterval = '30 sec'
    }


When using two (or more) different executors in your pipeline, you can specify their settings separately by prefixing
the executor name with the symbol ``$`` and using it as special scope identifier. For example::

  executor {
    $sge {
        queueSize = 100
        pollInterval = '30sec'
    }

    $local {
        queueSize = 20
        pollInterval = '1sec'
    }
  }

The above configuration example can be rewritten using the dot notation as shown below::

  executor.$sge.queueSize = 100
  executor.$sge.pollInterval = '30sec'
  executor.$local.queueSize = 10
  executor.$local.pollInterval = '1sec'



.. _config-docker:

Scope 'docker'
---------------

The ``docker`` configuration scope controls how `Docker <http://www.docker.io>`_ containers are executed by Nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             Turn this flag to ``true`` to enable Docker execution (default: ``false``)
sudo                Executes Docker run command as ``sudo`` (default: ``false``)
tty                 Allocate a pseudo-tty (default: ``false``)
temp                Mounts a path of your choice as the ``/tmp`` directory in the container. Use the special value ``auto`` to create a temporary directory each time a container is created
remove              Clean-up the container after the execution (default: ``true``). For details see: http://docs.docker.com/reference/run/#clean-up-rm
runOptions          This attribute can be used to provide any extra command line options supported by the ``docker run`` command. For details see: http://docs.docker.com/reference/run
registry            The registry from where Docker images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. ``http://``
================== ================

The above options can be used by prefixing them with the ``docker`` scope or surrounding them by curly
brackets, as shown below::

    process.container = 'nextflow/examples'

    docker {
        enabled = true
        temp = 'auto'
    }



Read :ref:`docker-page` page to lean more how use Docker containers with Nextflow.


.. _config-manifest:

Scope `manifest`
-----------------

The ``manifest`` configuration scope allows you to define some meta-data information needed when publishing your
pipeline on GitHub or BitBucket.

The following settings are available:

================== ================
Name                Description
================== ================
homePage            Project home page URL
description         Free text describing the pipeline project
mainScript          Pipeline main script (default: ``main.nf``)
defaultBranch       Git repository default branch (default: ``master``)
================== ================

The above options can be used by prefixing them with the ``manifest`` scope or surrounding them by curly
brackets. For example::

    manifest {
        homePage = 'http://foo.com'
        description = 'Pipeline does this and that'
        mainScript = 'foo.nf'
    }


To learn how to publish your pipeline on GitHub or BitBucket code repositories read :ref:`sharing-page`
documentation page.

.. _config-trace:

Scope `trace`
--------------

The ``trace`` scope allows you to control the layout of the execution trace file generated by nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             When ``true`` turns on the generation of the execution trace report file (default: ``false``).
fields              Comma separated list of fields to be included in the report. The available fields are listed at :ref:`this page <trace-fields>`
file                Trace file name (default: ``trace.csv``).
sep                 Character used to separate values in each row (default: ``\t``).
raw                 When ``true`` turns on raw number report generation i.e. date and time are reported as milliseconds and memory as number of bytes
================== ================

The above options can be used by prefixing them with the ``trace`` scope or surrounding them by curly
brackets. For example::

    trace {
        enabled = true
        file = 'pipeline_trace.csv'
        fields = 'task_id,name,status,exit_code,run_time,%cpu,rss'
    }


To learn more about the execution report that can be generated by Nextflow read :ref:`trace-report` documentation page.


Environment variables
======================

The following environment variables control the configuration of the Nextflow runtime and
the Java virtual machine used by it.

================== ================
Name                Description
================== ================
NXF_HOME            Nextflow home directory (default: ``$HOME/.nextflow``).
NXF_VER             Defines what version of Nextflow to use.
NXF_ORG             Default `organization` prefix when looking for a hosted repository (default: ``nextflow-io``).
NXF_GRAB            Provides extra runtime dependencies downloaded from a Maven repository service.
NXF_OPTS            Provides extra options for the Java and Nextflow runtime. It must be a blank separated list of ``-Dkey[=value]`` properties.
NXF_CLASSPATH       Allows to extend the Java runtime classpath with extra jar files or class folders.
NXF_DRMAA           Defines the Java DRMAA binding library to be used. It can be specified as a jar file location or a Maven dependency.
NXF_ASSETS          Defined the directory where downloaded pipeline repositories are stored (default: `$NXF_HOME/assets`)
NXF_PID_FILE        Name of the file where the process PID is saved when Nextflow is launched in background.
NXF_WORK            Directory where working files are stored (usually your *scratch* directory)
NXF_TEMP            Directory where temporary files are stored
JAVA_HOME           Path location of the Java VM installation used to run Nextflow.
JAVA_CMD            Path location of the Java binary command used to launch  Nextflow.
HTTP_PROXY          Defines the HTTP proxy server
HTTPS_PROXY         Defines the HTTPS proxy server
================== ================

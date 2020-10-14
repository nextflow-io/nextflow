.. _cli-page:

*****************************
Command line interface (CLI)
*****************************

`Nextflow` provides a robust command line interface for the management and 
execution pipelines. The top-level interface consists of two aspects, 
*options** and **commands**.

Here's what you'll see at the top-level upon invoking the ``nextflow`` CLI. ::


    $ nextflow
    Usage: nextflow [options] COMMAND [arg...]



.. _cli-options:

Options
============

The top-level ``options`` are meant to be invoked in relation to the core 
**nextflow runtime** and are applied to all commands. For options 
specific to any ``command``, please refer the `CLI commands` section.

An overview of the `Nextflow` top-level options. ::


    $ nextflow
    Usage: nextflow [options] COMMAND [arg...]

    Options:
    -C
        Use the specified configuration file(s) overriding any defaults
    -D
        Set JVM properties
    -bg
        Execute nextflow in background
    -c, -config
        Add the specified file to configuration set
    -d, -dockerize
        Launch nextflow via Docker (experimental)
    -h
        Print this help
    -log
        Set nextflow log file path
    -q, -quiet
        Do not print information messages
    -syslog
        Send logs to syslog server (eg. localhost:514)
    -v, -version
        Print the program version

    Commands...

---------------------------
Hard configuration override
---------------------------


**Description**

Use the specified configuration file(s) overriding any defaults.


**Usage** ::

   $ nextflow -C nxf.config COMMAND [arg...]


**Extended description**

The ``-C`` option is used to override *all* settings specified in the default config file. 
For soft override, please refer the ``c`` option.


**Examples**


- Override **any** previous config with a custom ``config`` file. ::
    
  $ nextflow -C nxf.config run nextflow-io/hello

--------------------
JVM properties
--------------------


**Description**

Set JVM properties.

**Usage**
::

   $ nextflow -Dkey=value COMMAND [arg...]

**Extended description**

The Java Virtual Machine (JVM) platform used by ``nextflow`` allows for the addition of 
system level properties and values during the startup that can then be accessed during 
runtime. For specifying other JVM level options, please refer to the `environment variables section <https://www.nextflow.io/docs/latest/config.html#environment-variables>`__.


**Examples**

Add ``JVM properties`` to the invoked pipeline. ::
    
    $ nextflow -Dkey=value run nextflow-io/hello


-----------------------------
Execution as a background job
-----------------------------


**Description**

Execute ``nextflow`` in the background.


**Usage**
::

   $ nextflow -bg COMMAND [arg...]

**Extended description**

The ``-bg`` option is used to invoke the nextflow execution in the background and allows 
the user to continue interacting with the terminal. This option is similar to ``nohup`` in 
behavior.


**Examples**

Invoke any execution as a background job. ::
    
    $ nextflow -bg run nextflow-io/hello 



---------------------------
Soft configuration override
---------------------------


**Description**

Add the specified file to configuration set.

**Usage**

::

   $ nextflow -c nxf.config COMMAND [arg...]


**Extended description**

The ``-c`` option is used to append a new configuration to the default configuration. 
The ``-c`` option allows us to update the config in an additive manner. For **hard override** 
please refer the ``-C`` option.


**Examples**

Update **some** fields of the default config for any pipeline. ::

  $ nextflow -c nxf.config run nextflow-io/hello



-----------------------
Docker driven execution
-----------------------


**Description**

Launch ``nextflow`` via Docker (experimental).


**Usage**
::

   $ nextflow -dockerize COMMAND [arg...]


**Extended description**

The ``-dockerize`` option is used to invoke **nextflow runtime** as a docker container 
itself. For invoking a pipeline with the ``docker`` profile or executor, please 
to refer the ``-with-docker`` options the ``run`` and ``kuberun`` commands.



**Examples**

Invoke ``nextflow`` as a docker container to execute a pipeline. ::

   $ nextflow -dockerize run nextflow-io/hello




--------------------
Help
--------------------


**Description**

Print the help message.


**Usage**
::

   $ nextflow -h

**Extended description**

The ``-h`` option prints out the overview of the CLI interface and enumerates the top-level *options* 
and *commands*.


--------------------
Execution logs
--------------------


**Description**

Sets the path of the nextflow log file.


**Usage**
::

   $ nextflow -log custom.log COMMAND [arg...]


**Extended description**

The ``log`` option takes a path of the new log file which to be used instead of the 
default ``.nextflow.log`` for storing execution logs.


**Examples**

Save all execution logs to the custom ``nxf.log`` file. ::

   $ nextflow -log nxf.log run nextflow-io/hello



--------------------
Quiet execution
--------------------


**Description**

Disable the printing of information to the terminal.

**Usage**
::

    $ nextflow -q COMMAND [arg...]

**Extended description**

The ``-q`` option suppresses the banner, process related info and exits once the 
execution is completed. Please note that it does not affect any explicit print 
statement within a pipeline.


**Examples**

Invoke the pipeline execution without the banner and pipeline information. ::

   $ nextflow -q run nextflow-io/hello



---------------------------
Logging to a syslog server
---------------------------


**Description**

Send logs to syslog server.

**Usage**
::

    $ nextflow -syslog localhost:1234 COMMAND [arg...]


**Extended description**

The ``-syslog`` option is used to send logs to a ``syslog`` logging server at the specified endpoint.


**Examples**

Send the logs to a ``syslog`` server at specific endpoint. ::

    $ nextflow -syslog localhost:1234 run nextflow-io/hello





--------------------
Version
--------------------


**Description**

Print the ``nextflow`` version information.

**Usage**

::

    $ nextflow -v


**Extended description**

The ``-v`` option prints out information about *Nextflow* such as the version and build. 
The ``-version`` option in addition prints out the citation reference and official website.

**Examples**

- The short version. ::

      $ nextflow -v
      nextflow version 20.07.1.5412


- The full version info with citation and website link. ::

      $ nextflow -version
      N E X T F L O W
      version 20.07.1 build 5412
      created 24-07-2020 15:18 UTC (20:48 IDT)
      cite doi:10.1038/nbt.3820
      http://nextflow.io


.. _cli-commands:

Commands
============

An overview of the Nextflow top-level commands. ::


    $ nextflow

    Usage: nextflow [options] COMMAND [arg...]
    
    Options...

    Commands:
    clean         Clean up project cache and work directories
    clone         Clone a project into a folder
    config        Print a project configuration
    console       Launch Nextflow interactive console
    drop          Delete the local copy of a project
    help          Print the usage help for a command
    info          Print project and system runtime information
    kuberun       Execute a workflow in a Kubernetes cluster (experimental)
    list          List all downloaded projects
    log           Print executions log and runtime info
    pull          Download or update a project
    run           Execute a pipeline project
    self-update   Update nextflow runtime to the latest available version
    view          View project script file(s)

--------------------
clean
--------------------


**Description**

Clean up *cache* and *work* directories.

**Usage**


::

    $ nextflow clean [run_name|session_id] [options]


**Extended description**

Upon invocation within a directory, ``nextflow`` creates a project specific ``.nextflow.log`` 
file, ``.nextflow`` cache directory as well as a ``work`` directory. The ``clean`` command is 
designed to facilitate removal of these files from previous executions. 
A list of of run names and session ids can be generated by invoking ``nextflow log -q``.


**Options**


+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -after                    |            | Clean up runs executed *after* the specified one.                              |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -before                   |            | Clean up runs executed *before* the specified one.                             |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -but                      |            | Clean up all runs *except* the specified one.                                  |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -dry-run, -n              |   false    | Print names of files to be removed without deleting them.                      | 
+---------------------------+------------+--------------------------------------------------------------------------------+
| -force, -f                |   false    | Force clean command.                                                           |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -help, -h                 |   false    | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -keep-logs, -k            |   false    | Removes only temporary files but retains execution log entries and metadata.   |                                           
+---------------------------+------------+--------------------------------------------------------------------------------+
| -quiet, -q                |   false    | Do not print names of files removed.                                           |
+---------------------------+------------+--------------------------------------------------------------------------------+



**Examples**

Dry run to remove work directories for the run name `boring_euler`.::

   $ nextflow clean boring_euler -n

   Would remove work/92/c1a9cd9a96e0531d81ca69f5dc3bb7
   Would remove work/3f/70944c7a549b6221e1ccc7b4b21b62
   Would remove work/0e/2ebdba85f76f6068b21a1bcbf10cab

Remove work directories for the run name `boring_euler`. ::

   $ nextflow clean boring_euler -f

   Removed work/92/c1a9cd9a96e0531d81ca69f5dc3bb7
   Removed work/3f/70944c7a549b6221e1ccc7b4b21b62
   Removed work/0e/2ebdba85f76f6068b21a1bcbf10cab


Remove the execution entries *except* for a specific execution. ::

    $ nextflow clean -but tiny_leavitt -f

    Removed work/1f/f1ea9158fb23b53d5083953121d6b6
    Removed work/bf/334115deec60929dc18edf0010032a
    Removed work/a3/06521d75da296d4dd7f4f8caaddad8

Dry run to remove the execution data *before* a specific execution. ::

   $ nextflow clean -before tiny_leavitt -n

   Would remove work/5d/ad76f7b7ab3500cf616814ef644b61
   Would remove work/c4/69a82b080a477612ba8d8e4c27b579
   Would remove work/be/a4fa2aa38f76fd324958c81c2e4603
   Would remove work/54/39116773891c47a91e3c1733aad4de


Dry run to remove the execution data *after* a specific execution. ::

   $ nextflow clean -after focused_payne -n

   Would remove work/1f/f1ea9158fb23b53d5083953121d6b6
   Would remove work/bf/334115deec60929dc18edf0010032a
   Would remove work/a3/06521d75da296d4dd7f4f8caaddad8


- Dry run to remove the temporary execution data for a specific execution, while keeping the log files. ::

   $ nextflow clean -keep-logs tiny_leavitt -n

   Would remove temp files from work/1f/f1ea9158fb23b53d5083953121d6b6
   Would remove temp files from work/bf/334115deec60929dc18edf0010032a
   Would remove temp files from work/a3/06521d75da296d4dd7f4f8caaddad8


--------------------
clone         
--------------------


**Description**

Clone a remote project into a folder.


**Usage**


::

    $ nextflow clone [options] [project]


**Extended description**


The ``clone`` command downloads a pipeline from a Git-hosting platform into the *current directory* 
and modifies it accordingly. For downloading a pipeline into the global cache ``~/.nextflow/assets`` , 
please refer to the ``nextflow pull`` command.

**Options**


+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -hub                      |  github    | Service hub where the project is hosted. Options: ``gitlab`` or ``bitbucket``  |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -r                        |  master    | Revision to clone - It can be a git ``branch``, ``tag`` or ``revision number`` |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -user                     |            | Private repository user name                                                   |
+---------------------------+------------+--------------------------------------------------------------------------------+




**Examples**


Clone the latest revision of a pipeline. ::

    $ nextflow clone nextflow-io/hello
    nextflow-io/hello cloned to: hello


Clone a specific revision of a pipeline. ::

    $ nextflow clone nextflow-io/hello -r v1.1
    nextflow-io/hello cloned to: hello



--------------------
config        
--------------------


**Description**

Print the resolved pipeline configuration.

**Usage**


::

    $ nextflow config [options]


**Extended description**


The ``config`` command is used for printing the project's configuration i.e. the ``nextflow.config`` 
and is especially useful for understanding the resolved profiles and parameters that Nextflow will use 
run a pipeline. For in-depth information, please refer `config-profiles section <https://www.nextflow.io/docs/latest/config.html#config-profiles>`_.

**Options**


+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -flat                     |  false     | Print config using flat notation.                                              |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -profile                  |            | Choose a configuration profile.                                                |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -properties               |  false     | Print config using Java properties notation.                                   |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -show-profiles, -a        |  false     | Show all configuration profiles.                                               |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -sort                     |  false     | Sort config attributes.                                                        |
+---------------------------+------------+--------------------------------------------------------------------------------+




**Examples**


Print out the inferred config using a the default group key-value notation. ::

   $ nextflow config

   docker {
      enabled = true
   }

   process {
      executor = 'local'
   }

Print out the config using a flat notation. ::

   $ nextflow config -flat

   docker.enabled = true
   process.executor = 'local'


Print out the config using the Java properties notation. ::

   $ nextflow config -properties

   docker.enabled = true
   process.executor = local


Print out all profiles from the project's configuration. ::

   $ nextflow config -show-profiles

   docker {
      enabled = true
   }

   profiles {
      standard {
         process {
            executor = 'local'
         }
      }
      cloud {
         process {
            executor = 'cirrus'
            container = 'cbcrg/imagex'
         }
      }
   }

--------------------
console       
--------------------


**Description**

Launch the *Nextflow* interactive console.


**Usage**


::

    $ nextflow console



**Extended description**

The ``console`` command is a wrapper over the Groovy *console* and provides a Graphic User 
Interface (GUI) and an interactive REPL (Read-Eval-Print-Loop) for quick experimentation.


**Options**

None available


**Examples**


Launch the ``console`` GUI. ::

  $ nextflow console


--------------------
drop          
--------------------


**Description**

Delete the local copy of a project.


**Usage**


::

    $ nextflow drop [options] [project]


**Extended description**


The ``drop`` command is used to remove the projects which have been downloaded into the 
global cache. Please refer the ``list`` command for generating a list of downloaded pipelines.

**Options**


+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -f                        |  false     | Delete the repository without taking care of local changes.                    |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+



**Examples**


Drop the ``nextflow-io/hello`` project. ::

  $ nextflow drop nextflow-io/hello


Forcefully drop the ``nextflow-io/hello`` pipeline, ignoring any local changes. ::

  $ nextflow drop nextflow-io/hello -f


--------------------
help          
--------------------


**Description**

Print the top-level help or specific help for a command.


**Usage**


::

    $ nextflow help [options] [command]


**Extended description**

The ``help`` command prints out the overview of the CLI interface and enumerates the top-level 
*options* and *commands*. Note that this command is equivalent to simply invoking ``nextflow`` 
at the command line.

**Options**


+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+



**Examples**

Invoke the ``help`` option for the ``drop`` command. ::

     $ nextflow help drop
 
     Delete the local copy of a project
     Usage: drop [options] name of the project to drop
        Options:
          -f
               Delete the repository without taking care of local changes
               Default: false
          -h, -help
               Print the command usage
               Default: false


--------------------
info          
--------------------


**Description**

Print project or system runtime information.


**Usage**


::

    $ nextflow info [options] [project]



**Extended description**


The ``info`` command prints out the nextflow runtime information about the hardware as 
well as the software versions of the ``Nextflow version and build``, ``Operating System`` 
and ``Groovy and Java runtime``. It can also be used to display information about a 
specific project.

**Options**


+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -check-updates, -u        |  false     | Check for remote updates.                                                      |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -d                        |  false     | Show detailed information.                                                     |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -o                        |  text      | Output format, either ``text``, ``json`` or ``yaml``.                          |
+---------------------------+------------+--------------------------------------------------------------------------------+



**Examples**

Display `nextflow` run-time and system info::

    $ nextflow info

      Version: 20.07.1 build 5412
      Created: 24-07-2020 15:18 UTC (20:48 IDT)
      System: Mac OS X 10.15.6
      Runtime: Groovy 2.5.11 on OpenJDK 64-Bit Server VM 1.8.0_192-b01
      Encoding: UTF-8 (UTF-8)

Display information about a specific project::

    $ nextflow info nextflow-io/hello

      project name: nextflow-io/hello
      repository  : https://github.com/nextflow-io/hello
      local path  : /Users/evanfloden/.nextflow/assets/nextflow-io/hello
      main script : main.nf
      revisions   : 
      * master (default)
        mybranch
        testing
        v1.1 [t]
        v1.2 [t]


--------------------
kuberun       
--------------------


**Description**

Deploy Nextflow into a Kubernetes cluster (experimental)


**Usage**

::

    $ nextflow kuberun [options] [project]


**Extended description**

The ``kuberun`` command builds upon the ``run`` command and offers a deep integration with 
the Kubernetes execution environment. This command deploys the Nextflow runtime as a Kubernetes 
pod and assumes that you've already installed the ``kubectl`` CLI. The ``kuberun`` command 
does not allow the execution of **local** Nextflow scripts. For more information please refer 
the `Kubernetes executor section <https://www.nextflow.io/docs/latest/config/kubernetes.html>`__.


**Options**


+---------------------------+-------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default     | Description                                                                    |
+===========================+=============+================================================================================+
| -E                        | false       | Exports all current system environment.                                        |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -ansi-log                 |             | Enable/disable ANSI console logging.                                           |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -bucket-dir               |             | Remote bucket where intermediate result files are stored.                      |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -cache                    |             | Enable/disable processes caching.                                              |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -dsl2                     | false       | Execute the workflow using DSL2 syntax.                                        |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -dump-channels            |             | Dump channels for debugging purpose.                                           |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -dump-hashes              | false       | Dump task hash keys for debugging purpose.                                     |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -e.                       | {}          | Add the specified variable to execution environment. Syntax: ``-e.key=value``  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -entry                    |             | Entry workflow name to be executed.                                            |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -h, -help                 | false       | Print the command usage.                                                       |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -hub                      | github      | Service hub where the project is hosted. Options: ``gitlab`` or ``bitbucket``  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -latest                   | false       | Pull latest changes before run.                                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -lib                      |             | Library extension path.                                                        |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -name                     |             | Assign a mnemonic name to the a pipeline run.                                  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -n, -namespace            |             | Specify the K8s namespace to use.                                              |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -offline                  | false       | Do not check for remote project updates.                                       |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -params-file              |             | Load script parameters from a JSON/YAML file.                                  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -pod-image                |             | Specify the container image for the Nextflow pod.                              |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -process.                 | {}          | Set process options. Syntax ``-process.key=value``                             |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -profile                  |             | Choose a configuration profile.                                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -qs, -queue-size          |             | Max number of processes that can be executed in parallel by each executor.     |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -resume                   |             | Execute the script using the cached results, useful to continue executions that|
|                           |             | was stopped by an error.                                                       |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -r, -revision             |             | Revision of the project to run (either a git branch, tag or commit SHA number) |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -test                     |             | Test a script function with the name specified.                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -user                     |             | Private repository user name.                                                  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -v, -volume-mount         |             | Volume claim mounts eg. ``my-pvc:/mnt/path``                                   |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-conda               |             | Use the specified Conda environment package or                                 |
|                           |             | file (must end with ``.yml|.yaml``)                                            |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-dag                 | dag.dot     | Create pipeline DAG file.                                                      |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-docker              |             | Enable process execution in a Docker container.                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -N, -with-notification    |             | Send a notification email on workflow completion to the specified recipients.  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-podman              |             | Enable process execution in a Podman container.                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-report              | report.html | Create processes execution html report.                                        |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-singularity         |             | Enable process execution in a Singularity container.                           |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-timeline            |timeline.html| Create processes execution timeline file.                                      |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-tower               |             | Monitor workflow execution with Seqera Tower service.                          |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-trace               | trace.txt   | Create processes execution tracing file.                                       |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-weblog              |             | Send workflow status messages via HTTP to target URL.                          |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -without-docker           | false       | Disable process execution with Docker.                                         |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -without-podman           |             | Disable process execution in a Podman container.                               |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -w, -work-dir             | work        | Directory where intermediate result files are stored.                          |
+---------------------------+-------------+--------------------------------------------------------------------------------+



**Examples**

Execute a pipeline into a Kubernetes cluster. ::

     $ nextflow kuberun nextflow-io/hello 



--------------------
list          
--------------------


**Description**

List all downloaded projects.


**Usage**

::

    $ nextflow list [options]



**Extended description**


The ``list`` commands prints a list of the projects which are already downloaded into the global cache ``~/.nextflow/assets``.


**Options**


+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+



**Examples**

List the downloaded pipelines. ::

    $ nextflow list

    nextflow-io/hello
    nextflow-hub/fastqc


--------------------
log           
--------------------


**Description**

Print the execution history and log information.


**Usage**


::

    $ nextflow log [options] [run_name | session_id]




**Extended description**

The ``log`` command is used to query the execution metadata associated with pipelines executed 
by Nextflow. The list of executed pipelines can be generated by issuing ``nextflow log`` at the terminal. 
Instead of run name, it's also possible to use a session id. Moreover, this command contains multiple options 
to facilitate the queries and is especially useful while debugging a pipeline and while inspecting pipeline 
execution metadata.


**Options**



+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -after                    |            | Show log entries for runs executed *after* the specified one.                  |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -before                   |            | Show log entries for runs executed *before* the specified one.                 |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -but                      |            | Show log entries for runs executed *but* the specified one.                    |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -filter, -F               |            | Filter log entires by a custom expression                                      |
|                           |            | e.g. ``process =~ /foo.*/ && status == 'COMPLETED'``                           |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -list-fields, -l          |  false     | Show all available fields.                                                     |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -quiet                    |  false     | Show only run names.                                                           |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -s                        |            | Character used to separate column values                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -template, -t             |            | Text template used to each record in the log.                                  |
+---------------------------+------------+--------------------------------------------------------------------------------+






**Examples**


Listing the execution logs of previous invocations of all pipelines in a project. ::

    $ nextflow log

    TIMESTAMP          	DURATION	RUN NAME     	STATUS	REVISION ID	SESSION ID                          	COMMAND
    2020-10-07 11:52:24	2.1s    	focused_payne	OK    	96eb04d6a4 	af6adaaa-ad4f-48a2-9f6a-b121e789adf5	nextflow run nextflow-io/hello -r master
    2020-10-07 11:53:00	3.1s    	tiny_leavitt 	OK    	e3b475a61b 	4d3b95c5-4385-42b6-b430-c865a70d56a4	nextflow run ./tutorial.nf
    2020-10-07 11:53:29	2.5s    	boring_euler 	OK    	e3b475a61b 	a6276975-7173-4208-ae09-ab9d6dce8737	nextflow run tutorial.nf


Listing only the *run names* of the execution logs of all pipelines invocations in a project. ::

    $ nextflow log -quiet

    focused_payne
    tiny_leavitt
    boring_euler

List the execution entries *only* a specific execution. ::

   $ nextflow log tiny_leavitt

   work/1f/f1ea9158fb23b53d5083953121d6b6
   work/bf/334115deec60929dc18edf0010032a
   work/a3/06521d75da296d4dd7f4f8caaddad8


List the execution entries *after* a specific execution. ::

    $ nextflow log -after tiny_leavitt

    work/92/c1a9cd9a96e0531d81ca69f5dc3bb7
    work/3f/70944c7a549b6221e1ccc7b4b21b62
    work/0e/2ebdba85f76f6068b21a1bcbf10cab

List the execution entries *before* a specific execution. ::

    $ nextflow log -before tiny_leavitt

    work/5d/ad76f7b7ab3500cf616814ef644b61
    work/c4/69a82b080a477612ba8d8e4c27b579
    work/be/a4fa2aa38f76fd324958c81c2e4603
    work/54/39116773891c47a91e3c1733aad4de

List the execution entries *except* for a specific execution. ::

   $ nextflow log -but tiny_leavitt

    work/5d/ad76f7b7ab3500cf616814ef644b61
    work/c4/69a82b080a477612ba8d8e4c27b579
    work/be/a4fa2aa38f76fd324958c81c2e4603
    work/54/39116773891c47a91e3c1733aad4de

Filter specific fields from the execution log of a process. ::

    $ nextflow log tiny_leavitt -f 'process,exit,hash,duration'

    splitLetters	0	1f/f1ea91	112ms
    convertToUpper	0	bf/334115	144ms
    convertToUpper	0	a3/06521d	139ms

Filter fields from the execution log of a process based on a criteria. ::

    $ nextflow log tiny_leavitt -F 'process =~ /splitLetters/'

    work/1f/f1ea9158fb23b53d5083953121d6b6


--------------------
pull          
--------------------


**Description**

Download or update a project.


**Usage**


::

    $ nextflow pull [options] [project]


**Extended description**


The ``pull`` command downloads a pipeline from a Git-hosting platform into the global cache ``~/.nextflow/assets`` 
and modifies it accordingly. For downloading a pipeline into a local directory, please refer to the ``nextflow clone`` command.


**Options**


+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -all                      |  false     | Update all downloaded projects.                                                |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -hub                      |  github    | Service hub where the project is hosted. Options: ``gitlab`` or ``bitbucket``  |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -r                        |            | Revision to run (either a git ``branch``, ``tag`` or commit ``SHA`` number).   |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -user                     |            | Private repository user name                                                   |
+---------------------------+------------+--------------------------------------------------------------------------------+



**Examples**


Download a new pipeline or pull the latest revision for a specific project. ::

    $ nextflow pull nextflow-io/hello

    Checking nextflow-io/hello ...
    done - revision: 96eb04d6a4 [master]

Pull the latest revision for all downloaded projects. ::

    $ nextflow pull -all

    Checking nextflow-io/hello ...
    done - revision: 96eb04d6a4 [master]
    Checking nextflow-hub/fastqc ...
    done - revision: 087659b18e [master]

Download a specific revision of a new project or pull the latest revision for a specific project. ::

    $ nextflow pull nextflow-io/hello -r v1.1

    Checking nextflow-io/hello ...
    checkout-out at AnyObjectId[1c3e9e7404127514d69369cd87f8036830f5cf64] - revision: 1c3e9e7404 [v1.1]


--------------------
run           
--------------------


**Description**

Execute a pipeline.


**Usage**

::

    $ nextflow run [options] [project]


**Extended description**


The ``run`` command is used to initiate the execution of ``nextflow`` script or 
downloaded pipeline. Along with serving the purpose of script execution, this command 
facilitates rapid iterations, inspections of any pipeline as well as debugging.


**Options**


+---------------------------+-------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default     | Description                                                                    |
+===========================+=============+================================================================================+
| -E                        |  false      | Exports all current system environment.                                        |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -ansi-log                 |             | Enable/disable ANSI console logging.                                           |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -bucket-dir               |             | Remote bucket where intermediate result files are stored.                      |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -cache                    |             | Enable/disable processes caching.                                              |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -dsl2                     | false       | Execute the workflow using DSL2 syntax.                                        |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -dump-channels            |             | Dump channels for debugging purpose.                                           |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -dump-hashes              | false       | Dump task hash keys for debugging purpose.                                     |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -e.                       | {}          | Add the specified variable to execution environment. Syntax: ``-e.key=value``  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -entry                    |             | Entry workflow name to be executed.                                            |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -h, -help                 | false       | Print the command usage.                                                       |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -hub                      | github      | Service hub where the project is hosted. Options: ``gitlab`` or ``bitbucket``  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -latest                   | false       | Pull latest changes before run.                                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -lib                      |             | Library extension path.                                                        |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -name                     |             | Assign a mnemonic name to the a pipeline run.                                  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -offline                  | false       | Do not check for remote project updates.                                       |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -params-file              |             | Load script parameters from a JSON/YAML file.                                  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -process.                 | {}          | Set process options. Syntax ``-process.key=value``                             |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -profile                  |             | Choose a configuration profile.                                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -qs, -queue-size          |             | Max number of processes that can be executed in parallel by each executor.     |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -resume                   |             | Execute the script using the cached results, useful to continue executions that|
|                           |             | was stopped by an error.                                                       |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -r, -revision             |             | Revision of the project to run                                                 |
|                           |             | (either a git ``branch``, ``tag`` or commit ``SHA`` number).                   |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -test                     |             | Test a script function with the name specified.                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -user                     |             | Private repository user name.                                                  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-conda               |             | Use the specified Conda environment package or                                 |
|                           |             | file (must end with ``.yml|.yaml``)                                            |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-dag                 | dag.dot     | Create pipeline DAG file.                                                      |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-docker              |             | Enable process execution in a Docker container.                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -N, -with-notification    |             | Send a notification email on workflow completion to the specified recipients.  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-podman              |             | Enable process execution in a Podman container.                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-report              | report.html | Create processes execution html report.                                        |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-singularity         |             | Enable process execution in a Singularity container.                           |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-timeline            |timeline.html| Create processes execution timeline file.                                      |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-tower               |             | Monitor workflow execution with Seqera Tower service.                          |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-trace               | trace.txt   | Create processes execution tracing file.                                       |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-weblog              |             | Send workflow status messages via HTTP to target URL.                          |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -without-docker           | false       | Disable process execution with Docker.                                         |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -without-podman           |             | Disable process execution in a Podman container.                               |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -w, -work-dir             | work        | Directory where intermediate result files are stored.                          |
+---------------------------+-------------+--------------------------------------------------------------------------------+


**Examples**

- Run a specific revision of a downloaded pipeline. ::

    $ nextflow run nextflow-io/hello -r v1.1

    N E X T F L O W  ~  version 20.07.1
    Launching `nextflow-io/hello` [grave_cajal] - revision: 1c3e9e7404 [v1.1]


- Choose a ``profile`` for running the project. Assumes that a profile named `docker` has already been defined in the config file. ::

    $ nextflow run main.nf -profile docker


- Invoke the pipeline execution and generate the summary HTML report. For more information on the metrics, please refer the `Tracing & visualization section <https://www.nextflow.io/docs/latest/tracing.html>`__. ::

    $ nextflow run main.nf -with-report


- Invoke the nextflow pipeline execution with a custom queue size. By default, the value of **queue-size** is the same as the number of available CPUs. ::

    $ nextflow run nextflow-io/hello -qs 4


- Execute the pipeline with ``dsl2`` .::

    $ nextflow run nextflow-io/hello -dsl2


- Invoke the pipeline with a specific workflow as the entry-point, this option is meant to be used with ``dsl2``. For more information on ``dsl2``, please refer the `DSL2 section <https://www.nextflow.io/docs/latest/dsl2.html>`__.. ::

   $ nextflow run main.nf -entry workflow_A


- Invoke the nextflow pipeline execution with the integrated monitoring dashboard Tower. For more information, please refer to the tower.nf website <https://www.tower.nf>`__. ::

    $ nextflow run nextflow-io/hello -with-tower
 


--------------------
self-update   
--------------------


**Description**

Update the nextflow runtime to the latest available version.


**Usage**

::

    $ nextflow self-update


**Extended description**

The ``self-update`` command directs the ``nextflow`` cli to update itself to the latest stable release.


**Examples**

Update Nextflow. ::

    $ nextflow self-update

          N E X T F L O W
          version 20.07.1 build 5412
          created 24-07-2020 15:18 UTC (20:48 IDT)
          cite doi:10.1038/nbt.3820
          http://nextflow.io


    Nextflow installation completed. Please note:
    - the executable file `nextflow` has been created in the folder: /usr/local/bin



--------------------
view          
--------------------


**Description**

View a projects script file(s).


**Usage**

::

    $ nextflow view [options] [project]



**Extended description**


The ``view`` command is used to inspect the pipelines which are already stored in the global nextflow cache. 
For downloading a pipeline into the global cache ``~/.nextflow/assets``, please refer to the ``pull`` command.

**Options**


+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -l                        |  false     | List repository content.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -q                        |  false     | Hide header line.                                                              |
+---------------------------+------------+--------------------------------------------------------------------------------+


**Examples**


Viewing the contents of a downloaded pipeline. ::

   $ nextflow view nextflow-io/hello

   == content of file: .nextflow/assets/nextflow-io/hello/main.nf
   #!/usr/bin/env nextflow

   cheers = Channel.from 'Bonjour', 'Ciao', 'Hello', 'Hola'

    process sayHello {
      echo true
      input:
        val x from cheers
      script:
        """
        echo '$x world!'
        """
    }


Listing the folder structure of the downloaded pipeline. ::

   $ nextflow view -l nextflow-io/hello

   == content of path: .nextflow/assets/nextflow-io/hello
   LICENSE
   README.md
   nextflow.config
   .gitignore
   circle.yml
   foo.nf
   .git
   .travis.yml
   main.nf


Viewing the contents of a downloaded pipeline without the header ``== contents of file ...``. ::

   $ nextflow view -q nextflow-io/hello

   #!/usr/bin/env nextflow

   cheers = Channel.from 'Bonjour', 'Ciao', 'Hello', 'Hola'

    process sayHello {
      echo true
      input:
        val x from cheers
      script:
        """
        echo '$x world!'
        """
    }

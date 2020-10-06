.. _cli-page:

******************
Command line interface
******************

`Nextflow` provides a robust command line interface allowing you to manage the execution of a pipeline. The top-level interface consists of two aspects, ``options`` and ``commands``.

Here's what you'll see at the top-level upon invoking ``nextflow`` CLI. ::


    $ nextflow
    Usage: nextflow [options] COMMAND [arg...]



.. _cli-options:
CLI options
============

The top-level ``options`` are meant to be invoked in relation to the core **nextflow runtime**. For child options specific to any ``command``, please refer the `CLI commands` section.

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

--------------------
Hard configuration override
--------------------


**Description**

Use the specified configuration file(s) overriding any defaults


**Usage**
::
   $ nextflow -C nxf.config COMMAND [arg...]

**Extended description**

The ``C`` option is used to override *all* settings specified in the default config file. For soft override, please refer the ``c`` option.


**Examples**


- Override **any** previous config with a custom ``config`` file. ::
    
  $ nextflow -C nxf.config run nextflow-io/hello

--------------------
JVM properties
--------------------


**Description**

Set JVM properties

**Usage**
::

   $ nextflow -Dkey=value COMMAND [arg...]

**Extended description**

The JVM platform, which is utilized by ``nextflow`` as its runtime, allows adding system level properties/values during the startup, which could be accessed during runtime. For specifying other JVM level options, please refer the
`environment variables section <https://www.nextflow.io/docs/latest/config.html#environment-variables>`__.


**Examples**

- Add ``JVM properties`` to the invoked pipeline. For more information on the ``-r`` child-option, please refer its parent, ``run`` command. ::
    
    $ nextflow -Dkey=value run nextflow-io/hello -r jvm-properties


--------------------
Execution as a background job
--------------------


**Description**

Execute nextflow in background


**Usage**
::

   $ nextflow -bg COMMAND [arg...]

**Extended description**

The ``bg`` option is used to invoke nextflow execution in the background and allows the user to continue interacting with the terminal. This option is similar to ``nohup`` in behavior.


**Examples**

- Invoke any execution as a background job. ::
    
    $ nextflow -bg run nextflow-io/hello 



--------------------
Soft configuration override
--------------------


**Description**

Add the specified file to configuration set

**Usage**

::

   $ nextflow -c nxf.config COMMAND [arg...]


**Extended description**

The ``c`` option is used to append or derive new configuration from the default configuration. The ``c`` option allows us to update the config in an additive manner, for **hard override** please refer the ``C`` option.



**Examples**

- Update **some** fields of the default config for any pipeline. ::

  $ nextflow -c nxf.config run nextflow-io/hello



--------------------
Docker driven execution
--------------------


**Description**

Launch nextflow via Docker (experimental)


**Usage**
::

   $ nextflow -dockerize COMMAND [arg...]


**Extended description**

The ``dockerize`` option is used to invoke **nextflow runtime** as a docker container itself. For invoking a pipeline with the ``docker`` profile or executor, please refer ``with-docker`` option of the ``run`` and ``kuberun`` commands.



**Examples**

- Invoke ``nextflow`` as a docker container to execute a pipeline. ::

   $ nextflow -dockerize run nextflow-io/hello




--------------------
Help
--------------------


**Description**

Print this help


**Usage**
::

   $ nextflow -h

**Extended description**

The ``h`` option prints out the overview of the CLI interface and enumerates the top-level *options* and *commands*.


--------------------
Execution logs
--------------------


**Description**

Set nextflow log file path


**Usage**
::

   $ nextflow -log custom.log COMMAND [arg...]


**Extended description**

The ``log`` option takes path of the new log file which would be used instead of the default ``.nextflow.log`` for storing execution logs for the pipeline.


**Examples**

- Save all execution logs to the custom ``nxf.log`` file. ::

   $ nextflow -log nxf.log run nextflow-io/hello



--------------------
Quiet execution
--------------------


**Description**

Do not print information messages

**Usage**
::

    $ nextflow -q COMMAND [arg...]

**Extended description**

The ``q`` option suppresses the banner, process related info and exits once the execution is completed. Please note that it does not affect any explicit print statement which is in the pipeline.


**Examples**

- Invoke the pipeline execution without the banner and pipeline informaton. ::

   $ nextflow -q run nextflow-io/hello





--------------------
Logging to a syslog server
--------------------


**Description**

Send logs to syslog server (eg. localhost:514)

**Usage**
::

    $ nextflow -syslog localhost:1234 COMMAND [arg...]


**Extended description**

The ``syslog`` option is used to send logs to a ``syslog`` logging server at the specified endpoint.


**Examples**

- Send the logs to a ``syslog`` server at a specific endpointt. ::

    $ nextflow -syslog localhost:1234 run nextflow-io/hello





--------------------
Version
--------------------


**Description**

Print the program version

**Usage**

::

    $ nextflow -v


**Extended description**

The ``v`` option prints out information about *Nextflow* such as ``Nextflow version and build``. The ``version`` variant in addition prints out the ``citation doi`` and ``official website`` as well.



- The short version info with citation and website link. ::

      nextflow version 20.07.1.5412


- The full version info with citation and website link. ::

      N E X T F L O W
      version 20.07.1 build 5412
      created 24-07-2020 15:18 UTC (20:48 IDT)
      cite doi:10.1038/nbt.3820
      http://nextflow.io


.. _cli-commands:
CLI commands
============

The top-level ``commands`` are meant to be invoked to initiate and inspect a **nextflow process**. For top-level options specific to **nextflow runtime**, please refer the `CLI options` section.

An overview of the `Nextflow` top-level command. ::


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

Clean up project specific *cache* and *work* directories

**Usage**


::

    $ nextflow clean RUN_NAME [options]


**Extended description**

Upon invocation within a directory, ``nextflow`` creates a project specific ``.nextflow.log`` file, ``.nextflow`` cache directory as well as a ``work`` directory. The ``clean`` option is designed to facilitate rapid iteration without the clutter introduced by previous executions. A list of ``RUN_NAME`` can be generated by invoking ``nextflow log -q``.



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

- Check what would be removed upon invocation of ``clean`` command using the ``dry-run`` option. ::

   $ nextflow clean boring_euler -n

   Would remove /Users/eklavya/projects/code/nextflow/_resources/work/92/c1a9cd9a96e0531d81ca69f5dc3bb7
   Would remove /Users/eklavya/projects/code/nextflow/_resources/work/3f/70944c7a549b6221e1ccc7b4b21b62
   Would remove /Users/eklavya/projects/code/nextflow/_resources/work/0e/2ebdba85f76f6068b21a1bcbf10cab

- Remove the execution data for a specific execution. ::

   $ nextflow clean boring_euler -f

   Removed /Users/eklavya/projects/code/nextflow/_resources/work/92/c1a9cd9a96e0531d81ca69f5dc3bb7
   Removed /Users/eklavya/projects/code/nextflow/_resources/work/3f/70944c7a549b6221e1ccc7b4b21b62
   Removed /Users/eklavya/projects/code/nextflow/_resources/work/0e/2ebdba85f76f6068b21a1bcbf10cab


- Remove the execution entries *except* for a specific execution. ::

    $ nextflow clean -but tiny_leavitt -f

    Removed /Users/eklavya/projects/code/nextflow/_resources/work/1f/f1ea9158fb23b53d5083953121d6b6
    Removed /Users/eklavya/projects/code/nextflow/_resources/work/bf/334115deec60929dc18edf0010032a
    Removed /Users/eklavya/projects/code/nextflow/_resources/work/a3/06521d75da296d4dd7f4f8caaddad8

- Dry run to remove the execution data *before* a specific execution. ::

   $ nextflow clean -before tiny_leavitt -n

   Would remove /Users/eklavya/projects/code/nextflow/_resources/work/5d/ad76f7b7ab3500cf616814ef644b61
   Would remove /Users/eklavya/projects/code/nextflow/_resources/work/c4/69a82b080a477612ba8d8e4c27b579
   Would remove /Users/eklavya/projects/code/nextflow/_resources/work/be/a4fa2aa38f76fd324958c81c2e4603
   Would remove /Users/eklavya/projects/code/nextflow/_resources/work/54/39116773891c47a91e3c1733aad4de


- Dry run to remove the execution data *after* a specific execution. ::

   $ nextflow clean -after focused_payne -n

   Would remove /Users/eklavya/projects/code/nextflow/_resources/work/1f/f1ea9158fb23b53d5083953121d6b6
   Would remove /Users/eklavya/projects/code/nextflow/_resources/work/bf/334115deec60929dc18edf0010032a
   Would remove /Users/eklavya/projects/code/nextflow/_resources/work/a3/06521d75da296d4dd7f4f8caaddad8


- Dry run to remove the temporary execution data for a specific execution, while saving the log files. ::

   $ nextflow clean -keep-logs tiny_leavitt -n

   Would remove temp files from /Users/eklavya/projects/code/nextflow/_resources/work/1f/f1ea9158fb23b53d5083953121d6b6
   Would remove temp files from /Users/eklavya/projects/code/nextflow/_resources/work/bf/334115deec60929dc18edf0010032a
   Would remove temp files from /Users/eklavya/projects/code/nextflow/_resources/work/a3/06521d75da296d4dd7f4f8caaddad8

--------------------
clone         
--------------------


**Description**

Clone a project into a folder



**Usage**


::

    $ nextflow clone [options]



**Extended description**


The ``clone`` command faciliatates collaboration by allowing the users to download any existing pipeline from the specified ``-hub`` into the *current directory* and modify it accordingly. For downloading a pipeline into the global cache ``~/.nextflow/assets`` , please refer ``pull`` command.

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


- Clone the latest revision of a pipeline. ::

    $ nextflow clone nextflow-io/hello
    nextflow-io/hello cloned to: hello


- Clone a specific revision of a pipeline. ::

    $ nextflow clone nextflow-io/hello -r mybranch #OR v1.1
    nextflow-io/hello cloned to: hello




--------------------
config        
--------------------


**Description**

Print a project configuration

**Usage**


::

    $ nextflow config [options]


**Extended description**


The ``config`` command is used for printing the project's configuration i.e. the ``nextflow.config`` and is especially useful for accomodating alternative executors, profiles, tools and parameters. For in-depth information, please refer `config-profiles section <https://www.nextflow.io/docs/latest/config.html#config-profiles>`_.

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


- Print out the inferred config using a the default group key-value notation. ::

   $ nextflow config

   docker {
      enabled = true
   }

   process {
      executor = 'local'
   }

- Print out the config using a flat notation. ::

   $ nextflow config -flat

   docker.enabled = true
   process.executor = 'local'


- Print out the config using the Java properties notation. ::

   $ nextflow config -properties

   docker.enabled = true
   process.executor = local


- Print out all profiles from the project's configuration. ::

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

Launch *Nextflow* interactive console


**Usage**


::

    $ nextflow console



**Extended description**

The ``console`` command is a wrapper over the Groovy *console* and provides a Graphic User Interface (GUI) and an interactive REPL (Read-Eval-Print-Loop) for quick experimentation.


**Options**

None available


**Examples**


- Launch the ``console`` GUI. ::

  $ nextflow console


--------------------
drop          
--------------------


**Description**

Delete the local copy of a project


**Usage**


::

    $ nextflow drop [options]




**Extended description**


The ``drop`` command is used to remove the piplines which have already been downloaded into the global cache. Please refer the ``list`` command for generating a list of downloaded pipelines.

**Options**


+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -f                        |  false     | Delete the repository without taking care of local changes.                    |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+



**Examples**


- Drop the downloaded ``nextflow-io/hello`` pipeline. ::

  $ nextflow drop nextflow-io/hello


- Forcefully drop the ``nextflow-io/hello`` pipeline, ignoring any local changes. ::

  $ nextflow drop nextflow-io/hello -f


--------------------
help          
--------------------


**Description**

Print the usage help for a command


**Usage**


::

    $ nextflow help [options]


**Extended description**

The ``help`` command prints out the overview of the CLI interface and enumerates the top-level *options* and *commands*. Note that, this command is equivalent to simply issuing ``nextflow`` at the command line.

**Options**


+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+



**Examples**

- Invoke the ``help`` option from the command line to see an overview of top-level commands and options. ::

    $ nextflow help

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
info          
--------------------


**Description**

Print project and system runtime information



**Usage**


::

    $ nextflow info [options]



**Extended description**


The ``info`` command prints out the nextflow runtime information about the hardware as well as the software versions of the ``Nextflow version and build``, ``Operating System`` and ``Groovy and Java runtime``.

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

- ::

    $ nextflow info

      Version: 20.07.1 build 5412
      Created: 24-07-2020 15:18 UTC (20:48 IDT)
      System: Mac OS X 10.15.6
      Runtime: Groovy 2.5.11 on OpenJDK 64-Bit Server VM 1.8.0_192-b01
      Encoding: UTF-8 (UTF-8)

- ::

    $ nextflow info -check-updates

      Version: 20.07.1 build 5412
      Created: 24-07-2020 15:18 UTC (20:48 IDT)
      System: Mac OS X 10.15.6
      Runtime: Groovy 2.5.11 on OpenJDK 64-Bit Server VM 1.8.0_192-b01
      Encoding: UTF-8 (UTF-8)

- Print out the detaild information from the ``nextflow`` runtime and system environment. ::

    $ nextflow info -d

      Version: 20.07.1 build 5412
      Created: 24-07-2020 15:18 UTC (20:48 IDT)
      System: Mac OS X 10.15.6
      Runtime: Groovy 2.5.11 on OpenJDK 64-Bit Server VM 1.8.0_192-b01
      Encoding: UTF-8 (UTF-8)
      File systems: file, jar, ftp, http, https, s3
      JVM opts:
        -XX:+TieredCompilation
        -XX:TieredStopAtLevel=1
        -Djava.io.tmpdir=/var/folders/zp/63677vtx23d_b2_nd7mm92040000gn/T/
        -Djava.library.path=
          /Users/eklavya/Library/Java/Extensions
          /Library/Java/Extensions
          /Network/Library/Java/Extensions
          /System/Library/Java/Extensions
          /usr/lib/java
        -Dfile.encoding=UTF-8
        -Djava.awt.headless=true
      Capsule:
        capsule.app=nextflow_20.07.1
        capsule.jar=/Users/eklavya/.nextflow/framework/20.07.1/nextflow-20.07.1-one.jar
      Environment:
        NXF_CLI=/usr/local/bin/nextflow info -d
        NXF_HOME=/Users/eklavya/.nextflow
        NXF_ORG=nextflow-io
      Class-path:
        /Users/eklavya/.nextflow/framework/20.07.1/nextflow-20.07.1-one.jar
        /Users/eklavya/.nextflow/capsule/deps/io/nextflow/nf-tower/20.07.1/nf-tower-20.07.1.jar
        /Users/eklavya/.nextflow/capsule/deps/io/nextflow/nf-amazon/20.07.1/nf-amazon-20.07.1.jar
        /Users/eklavya/.nextflow/capsule/deps/io/nextflow/nextflow/20.07.1/nextflow-20.07.1.jar
        /Users/eklavya/.nextflow/capsule/deps/org/apache/ivy/ivy/2.3.0/ivy-2.3.0.jar
        /Users/eklavya/.nextflow/capsule/deps/io/nextflow/nf-httpfs/20.07.1/nf-httpfs-20.07.1.jar
        /Users/eklavya/.nextflow/capsule/deps/io/nextflow/nf-commons/20.07.1/nf-commons-20.07.1.jar
        /Users/eklavya/.nextflow/capsule/deps/org/codehaus/groovy/groovy-nio/2.5.11/groovy-nio-2.5.11.jar
        /Users/eklavya/.nextflow/capsule/deps/org/codehaus/groovy/groovy-templates/2.5.11/groovy-templates-2.5.11.jar
        /Users/eklavya/.nextflow/capsule/deps/org/codehaus/groovy/groovy-xml/2.5.11/groovy-xml-2.5.11.jar
        /Users/eklavya/.nextflow/capsule/deps/org/codehaus/groovy/groovy-json/2.5.11/groovy-json-2.5.11.jar
        /Users/eklavya/.nextflow/capsule/deps/org/codehaus/groovy/groovy/2.5.11/groovy-2.5.11.jar
        /Users/eklavya/.nextflow/capsule/deps/com/amazonaws/aws-java-sdk-s3/1.11.542/aws-java-sdk-s3-1.11.542.jar
        /Users/eklavya/.nextflow/capsule/deps/com/amazonaws/aws-java-sdk-ec2/1.11.542/aws-java-sdk-ec2-1.11.542.jar
        /Users/eklavya/.nextflow/capsule/deps/com/amazonaws/aws-java-sdk-batch/1.11.542/aws-java-sdk-batch-1.11.542.jar
        /Users/eklavya/.nextflow/capsule/deps/com/amazonaws/aws-java-sdk-iam/1.11.542/aws-java-sdk-iam-1.11.542.jar
        /Users/eklavya/.nextflow/capsule/deps/com/amazonaws/aws-java-sdk-ecs/1.11.542/aws-java-sdk-ecs-1.11.542.jar
        /Users/eklavya/.nextflow/capsule/deps/com/amazonaws/aws-java-sdk-kms/1.11.542/aws-java-sdk-kms-1.11.542.jar
        /Users/eklavya/.nextflow/capsule/deps/com/amazonaws/aws-java-sdk-core/1.11.542/aws-java-sdk-core-1.11.542.jar
        /Users/eklavya/.nextflow/capsule/deps/org/apache/httpcomponents/httpclient/4.5.5/httpclient-4.5.5.jar
        /Users/eklavya/.nextflow/capsule/deps/org/slf4j/jcl-over-slf4j/1.7.25/jcl-over-slf4j-1.7.25.jar
        /Users/eklavya/.nextflow/capsule/deps/org/slf4j/jul-to-slf4j/1.7.25/jul-to-slf4j-1.7.25.jar
        /Users/eklavya/.nextflow/capsule/deps/org/slf4j/log4j-over-slf4j/1.7.25/log4j-over-slf4j-1.7.25.jar
        /Users/eklavya/.nextflow/capsule/deps/ch/qos/logback/logback-classic/1.1.11/logback-classic-1.1.11.jar
        /Users/eklavya/.nextflow/capsule/deps/ch/qos/logback/logback-core/1.1.11/logback-core-1.1.11.jar
        /Users/eklavya/.nextflow/capsule/deps/org/codehaus/gpars/gpars/1.2.1/gpars-1.2.1.jar
        /Users/eklavya/.nextflow/capsule/deps/ch/grengine/grengine/1.3.0/grengine-1.3.0.jar
        /Users/eklavya/.nextflow/capsule/deps/commons-lang/commons-lang/2.6/commons-lang-2.6.jar
        /Users/eklavya/.nextflow/capsule/deps/commons-codec/commons-codec/1.10/commons-codec-1.10.jar
        /Users/eklavya/.nextflow/capsule/deps/com/beust/jcommander/1.35/jcommander-1.35.jar
        /Users/eklavya/.nextflow/capsule/deps/com/esotericsoftware/kryo/kryo/2.24.0/kryo-2.24.0.jar
        /Users/eklavya/.nextflow/capsule/deps/org/iq80/leveldb/leveldb/0.12/leveldb-0.12.jar
        /Users/eklavya/.nextflow/capsule/deps/org/eclipse/jgit/org.eclipse.jgit/5.2.1.201812262042-r/org.eclipse.jgit-5.2.1.201812262042-r.jar
        /Users/eklavya/.nextflow/capsule/deps/javax/mail/mail/1.4.7/mail-1.4.7.jar
        /Users/eklavya/.nextflow/capsule/deps/javax/activation/activation/1.1.1/activation-1.1.1.jar
        /Users/eklavya/.nextflow/capsule/deps/org/yaml/snakeyaml/1.18/snakeyaml-1.18.jar
        /Users/eklavya/.nextflow/capsule/deps/org/jsoup/jsoup/1.11.2/jsoup-1.11.2.jar
        /Users/eklavya/.nextflow/capsule/deps/jline/jline/2.9/jline-2.9.jar
        /Users/eklavya/.nextflow/capsule/deps/io/nextflow/nxf-s3fs/1.0.8/nxf-s3fs-1.0.8.jar
        /Users/eklavya/.nextflow/capsule/deps/com/google/guava/guava/21.0/guava-21.0.jar
        /Users/eklavya/.nextflow/capsule/deps/org/slf4j/slf4j-api/1.7.25/slf4j-api-1.7.25.jar
        /Users/eklavya/.nextflow/capsule/deps/org/multiverse/multiverse-core/0.7.0/multiverse-core-0.7.0.jar
        /Users/eklavya/.nextflow/capsule/deps/org/codehaus/jsr166-mirror/jsr166y/1.7.0/jsr166y-1.7.0.jar
        /Users/eklavya/.nextflow/capsule/deps/org/objenesis/objenesis/2.1/objenesis-2.1.jar
        /Users/eklavya/.nextflow/capsule/deps/org/iq80/leveldb/leveldb-api/0.12/leveldb-api-0.12.jar
        /Users/eklavya/.nextflow/capsule/deps/com/jcraft/jsch/0.1.54/jsch-0.1.54.jar
        /Users/eklavya/.nextflow/capsule/deps/com/jcraft/jzlib/1.1.1/jzlib-1.1.1.jar
        /Users/eklavya/.nextflow/capsule/deps/com/googlecode/javaewah/JavaEWAH/1.1.6/JavaEWAH-1.1.6.jar
        /Users/eklavya/.nextflow/capsule/deps/com/amazonaws/jmespath-java/1.11.542/jmespath-java-1.11.542.jar
        /Users/eklavya/.nextflow/capsule/deps/software/amazon/ion/ion-java/1.0.2/ion-java-1.0.2.jar
        /Users/eklavya/.nextflow/capsule/deps/com/fasterxml/jackson/core/jackson-databind/2.6.7.2/jackson-databind-2.6.7.2.jar
        /Users/eklavya/.nextflow/capsule/deps/com/fasterxml/jackson/dataformat/jackson-dataformat-cbor/2.6.7/jackson-dataformat-cbor-2.6.7.jar
        /Users/eklavya/.nextflow/capsule/deps/joda-time/joda-time/2.8.1/joda-time-2.8.1.jar
        /Users/eklavya/.nextflow/capsule/deps/org/apache/httpcomponents/httpcore/4.4.9/httpcore-4.4.9.jar
        /Users/eklavya/.nextflow/capsule/deps/com/fasterxml/jackson/core/jackson-annotations/2.6.0/jackson-annotations-2.6.0.jar
        /Users/eklavya/.nextflow/capsule/deps/com/fasterxml/jackson/core/jackson-core/2.6.7/jackson-core-2.6.7.jar




--------------------
kuberun       
--------------------


**Description**

Execute a workflow in a Kubernetes cluster (experimental)


**Usage**

::

    $ nextflow kuberun [options]


**Extended description**

The ``kuberun`` command builds upon the ``run`` command and offers a deep integration with the ``Kubernetes`` execution environment. This ``command`` assumes that you've already installed the ``kubectl`` CLI. Also, please note that currently, the ``kuberun`` command does not allow the execution of **local** Nextflow scripts. For more information please refer the `Kubernetes executor section <https://www.nextflow.io/docs/latest/config/kubernetes.html>`__.

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
Please note that many of the options for ``kuberun`` are similar to ``run`` command, please also refer the examples from that section.



--------------------
list          
--------------------


**Description**

List all downloaded projects


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

- List the downloaded pipelines. ::

    $ nextflow list

    nextflow-io/hello
    nextflow-hub/fastqc


--------------------
log           
--------------------


**Description**

Print executions log and runtime info.


**Usage**


::

    $ nextflow log RUN_NAME [options]




**Extended description**

The ``log`` command is used to query the execution metadata associated with every pipeline executed by *Nextflow*, the list of executed pipelines can be generated by issuing ``nextflow log`` at the terminal. Instead of ``RUN_NAME``, it's also possible to use ``SESSION_ID`` Moreover, this command contains multiple options to facilitate the queries and is especially useful while debugging a pipeline and while inspecting the pipelines' execution metadata.


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


- Listing the execution logs of previous invocations of all pipelines in a project. ::

    $ nextflow log

    TIMESTAMP          	DURATION	RUN NAME     	STATUS	REVISION ID	SESSION ID                          	COMMAND
    2020-10-07 11:52:24	2.1s    	focused_payne	OK    	96eb04d6a4 	af6adaaa-ad4f-48a2-9f6a-b121e789adf5	nextflow run nextflow-io/hello -r master
    2020-10-07 11:53:00	3.1s    	tiny_leavitt 	OK    	e3b475a61b 	4d3b95c5-4385-42b6-b430-c865a70d56a4	nextflow run ./tutorial.nf
    2020-10-07 11:53:29	2.5s    	boring_euler 	OK    	e3b475a61b 	a6276975-7173-4208-ae09-ab9d6dce8737	nextflow run tutorial.nf


- Listing only the *run names* of the execution logs of all pipelines invocations in a project. ::

    $ nextflow log -quiet

    focused_payne
    tiny_leavitt
    boring_euler

- List the execution entries *only* a specific execution. ::

   $ nextflow log tiny_leavitt

   /Users/eklavya/projects/code/nextflow/_resources/work/1f/f1ea9158fb23b53d5083953121d6b6
   /Users/eklavya/projects/code/nextflow/_resources/work/bf/334115deec60929dc18edf0010032a
   /Users/eklavya/projects/code/nextflow/_resources/work/a3/06521d75da296d4dd7f4f8caaddad8


- List the execution entries *after* a specific execution. ::

    $ nextflow log -after tiny_leavitt

    /Users/eklavya/projects/code/nextflow/_resources/work/92/c1a9cd9a96e0531d81ca69f5dc3bb7
    /Users/eklavya/projects/code/nextflow/_resources/work/3f/70944c7a549b6221e1ccc7b4b21b62
    /Users/eklavya/projects/code/nextflow/_resources/work/0e/2ebdba85f76f6068b21a1bcbf10cab

- List the execution entries *before* a specific execution. ::

    $ nextflow log -before tiny_leavitt

    /Users/eklavya/projects/code/nextflow/_resources/work/5d/ad76f7b7ab3500cf616814ef644b61
    /Users/eklavya/projects/code/nextflow/_resources/work/c4/69a82b080a477612ba8d8e4c27b579
    /Users/eklavya/projects/code/nextflow/_resources/work/be/a4fa2aa38f76fd324958c81c2e4603
    /Users/eklavya/projects/code/nextflow/_resources/work/54/39116773891c47a91e3c1733aad4de

- List the execution entries *except* for a specific execution. ::

   $ nextflow log -but tiny_leavitt

    /Users/eklavya/projects/code/nextflow/_resources/work/5d/ad76f7b7ab3500cf616814ef644b61
    /Users/eklavya/projects/code/nextflow/_resources/work/c4/69a82b080a477612ba8d8e4c27b579
    /Users/eklavya/projects/code/nextflow/_resources/work/be/a4fa2aa38f76fd324958c81c2e4603
    /Users/eklavya/projects/code/nextflow/_resources/work/54/39116773891c47a91e3c1733aad4de

- Filter specific fields from the execution log of a process. ::

    $ nextflow log tiny_leavitt -f 'process,exit,hash,duration'

    splitLetters	0	1f/f1ea91	112ms
    convertToUpper	0	bf/334115	144ms
    convertToUpper	0	a3/06521d	139ms

- Filter fields from the execution log of a process based on a criteria. ::

    $ nextflow log tiny_leavitt -F 'process =~ /splitLetters/'

    /Users/eklavya/projects/code/nextflow/_resources/work/1f/f1ea9158fb23b53d5083953121d6b6

--------------------
pull          
--------------------


**Description**

Download or update a project.


**Usage**


::

    $ nextflow pull REPO_NAME [options]




**Extended description**


The ``pull`` command faciliatates collaboration by allowing the users to download any existing pipeline from the specified ``-hub`` and execute it using the ``run`` command. For downloading a pipeline into the project directory, please refer the ``clone`` command.


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


- Download a new pipeline or pull the latest revision for a specific downloaded pipelines. ::

    $ nextflow pull nextflow-io/hello

    Checking nextflow-io/hello ...
    done - revision: 96eb04d6a4 [master]

- Pull the latest revision for all downloaded pipelines. ::

    $ nextflow pull -all

    Checking nextflow-io/hello ...
    done - revision: 96eb04d6a4 [master]
    Checking nextflow-hub/fastqc ...
    done - revision: 087659b18e [master]

- Download a specific revision of a new pipeline or pull the latest revision for a specific downloaded pipelines. ::

    $ nextflow pull nextflow-io/hello -r mybranch #OR v1.1

    Checking nextflow-io/hello ...
    checkout-out at AnyObjectId[1c3e9e7404127514d69369cd87f8036830f5cf64] - revision: 1c3e9e7404 [mybranch]


--------------------
run           
--------------------


**Description**

Execute a pipeline project


**Usage**

::

    $ nextflow run [options]



**Extended description**


The ``run`` command is used to initiate the execution of ``nextflow`` script or downloaded pipeline. Along with serving the purpose of script execution, this command facilitates rapid iterations, inspections of any pipeline as well as debugging them via numerous options.

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

- Execute the pipeline with ``dsl2`` .::

    $ nextflow run nextflow-io/hello -dsl2

- Invoke the pipeline with a specific workflow as the entrypoint, this option is meant to be used with ``dsl2``. For more information on ``dsl2``, please refer the `DSL2 section <https://www.nextflow.io/docs/latest/dsl2.html>`__.. ::

   $ nextflow run main.nf -entry workflow_A

- Run a specific revision of a downloaded pipeline. ::

    $ nextflow run nextflow-io/hello -r mybranch #OR v1.1

    N E X T F L O W  ~  version 20.07.1
    Launching `nextflow-io/hello` [grave_cajal] - revision: 1c3e9e7404 [mybranch]
    executor ...

- Choose a ``profile`` for running the project. This option assumes that the ``profiles`` have already been defined in the ``config`` file. ::

    $ nextflow run main.nf -profile docker

- Invoke the pipeline execution and generate the summary ``HTML report``. For more information on the metrics, please refer the `Tracing & visualization section <https://www.nextflow.io/docs/latest/tracing.html>`__. ::

    $ nextflow run main.nf -with-report


- Invoke the nextflow pipeline execution with a custom queue size. By default, the value of **queue-size** is the same as the number of available CPUs. ::

    $ nextflow run nextflow-io/hello -qs 4


- Invoke the nextflow pipeline execution with integrated monitoring dashboard with ``tower.nf``. For more information, please refer `the tower.nf website <https://www.tower.nf>`__. ::

    $ nextflow run nextflow-io/hello -with-tower
 
--------------------
self-update   
--------------------


**Description**

Update nextflow runtime to the latest available version.


**Usage**

::

    $ nextflow self-update

**Extended description**


The ``self-update`` command directs the ``nextflow`` cli to update itself to the latest stable release.

**Examples**

::

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

View project script file(s).


**Usage**

::

    $ nextflow view [options]



**Extended description**


The ``view`` command is used to inspect the pipelines which are already stored in the global nextflow cache. For downloading a pipeline into the global cache ``~/.nextflow/assets`` , please refer ``pull`` command.

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


- Viewing the contents of a downloaded pipeline. ::

   $ nextflow view nextflow-io/hello

   == content of file: /Users/eklavya/.nextflow/assets/nextflow-io/hello/main.nf
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

- Listing the folder structure of the downloaded pipeline. ::

   $ nextflow view -l nextflow-io/hello

   == content of path: /Users/eklavya/.nextflow/assets/nextflow-io/hello
   LICENSE
   README.md
   nextflow.config
   .gitignore
   circle.yml
   foo.nf
   .git
   .travis.yml
   main.nf

- Viewing the contents of a downloaded pipeline without the header ``== contents of file ...``. ::


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

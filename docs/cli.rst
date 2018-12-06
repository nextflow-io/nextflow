.. _cli-page:

**********************
Command Line Interface
**********************

Nextflow has different command line arguments that can be used to modify the execution of workflows or interact with existing runs.  

Top options
===========

Nextflow Docker
---------------

Nextflow can be launched using a provided Docker container. This container has been modified to allow additional sibling containers to be launched from within the Nextflow container. Running Nextflow in this manner::

    nextflow -d run main.nf

clean
=====

Clean up project cache and work directories.

**Usage** 
::

    nextflow clean [options] [run name or session id]

**Options**

=============== =============
Name            Description
=============== =============
-after          Clean up runs executed after the specified one
-before         Clean up runs executed before the specified one
-but            Clean up all runs except the specified one
-n, -dry-run    Print names of file to be removed without deleting them
-f, -force      Force clean command
-k, -keep-logs  Removes only temporary files but retains execution log entries and metadata
-q, -quiet      Do not print names of files removed
=============== =============

**Examples**
::

    $ nextflow clean dreamy_swartz

Delete temporary data created by run with name ``dreamy_swartz``.

::

    $ nextflow -after dreamy_swartz

Delete temporary data created by runs executed *after* the one with name ``dreamy_swartz``. 

::

    $ nextflow -before dreamy_swartz

Delete temporary data created by runs executed *before* the one with name ``dreamy_swartz``. 

::

    nextflow -but dreamy_swartz

Delete temporary produced by any run except the one with name ``dreamy_swartz``. 


clone
=====

Clone a project into a folder.

**Usage** 
::

    nextflow clone [options] <project name> [target dir]

**Options**

=============== =============
Name            Description
=============== =============
-hub            Service hub where the project is hosted: ``github``, ``gitlab`` or ``bitbucket`` (default: ``github``).
-r              Revision to clone - It can be a git branch, tag or revision number.
-user           Private repository user name.
=============== =============

**Examples**
::

    nextflow clone nextflow-io/rnaseq-nf foo

Clone the workflow project ``nextflow-io/rnaseq-nf`` hosted on GitHub into the 
local directory ``foo``.

cloud
=====

Manage Nextflow clusters in the cloud


config
======

Show a project configuration

**Usage**
::

    nextflow config [options] [script or project]


======================= =============
Name                    Description
======================= =============
-flat                   Print config using flat notation
-profile                Choose the configuration profile
-properties             Prints config using Java properties notation
-a, -show-profiles      Show all configuration profiles
-sort                   Sort config attributes
======================= =============

**Examples**
::

    nextflow config nextflow-io/rnaseq-nf

drop
====

Delete the local copy of a project

**Usage**
::

    nextflow drop [options] <project name>

**Options**

=============== =============
Name            Description
=============== =============
-force          Delete the repository without taking care of local changes
=============== =============


help
====

Print the usage help for a command



info
====

Print a project and system runtime information

**Usage**
::

    nextflow info [options] [project name]

When no project name is specified the ``info`` command print
the Nextflow runtime information::

    $ nextflow info
      Version: 18.11.0-edge build 5016
      Modified: 12-11-2018 17:04 UTC (18:04 CEST)
      System: Mac OS X 10.14
      Runtime: Groovy 2.5.4 on Java HotSpot(TM) 64-Bit Server VM 1.8.0_161-b12
      Encoding: UTF-8 (UTF-8)


When it's specified the name of a project previously downloaded either with the command ``run`` or ``pull``,
the command ``info`` prints the project information::
    
    $ nextflow info nextflow-io/rnaseq-nf
     project name: nextflow-io/rnaseq-nf
     repository  : https://github.com/nextflow-io/rnaseq-nf
     local path  : /Users/pditommaso/.nextflow/assets/nextflow-io/rnaseq-nf
     main script : main.nf
     description : Proof of concept of a RNA-seq pipeline implemented with Nextflow
     author      : Paolo Di Tommaso
     revisions   :
     * master (default)
       dev
       hybrid
       k8s-demo
       v1.0 [t]
       v1.1 [t]

**Options**

=============== =============
Name            Description
=============== =============
-d              Show detailed information
=============== =============


list
====

List all projects downloaded::

    $ nextflow list
    nextflow-io/rnaseq-nf
    nextflow-io/hello
    CRG-CNAG/CalliNGS-NF



log
===

Print executions log and runtime info


The logs from previous Nextflow runs can be viewed using::

    nextflow log

You can specify a run name or use the special name ``last``::

    nextflow log last

Listing the available fields that can be accessed from the log::
    
    nextflow log <run_name> -l

Listing the values for specified fields::

    nextflow log <run_name> -f process,container,status

You can specify a seperator for the outputs, such as using a comma here::

    nextflow log <run_name> -f process,container -s ,

You can filter the results based on any of the fields with a expression::

    nextflow log <run_name> -f 'container,process' -F 'container == "ubuntu"'

The filter can also include regular expressions::

    nextflow log last -f "process,container" -F 'process =~ /bar.*/ && container =~ /biocontainers.*/'

You can specify either a string or file template to be used with the log this allows more clear formatting of your log outputs::

    nextflow log <run_name> -t 'container: $container\nprocess: $process'

or saved in a file ``template.md``::

    process: $process
    container: $container

the template file can then by specified::
    
    nextflow log <run_name> -t template.md


pull
====

Download a workflow from the source code management platform given the project name::

    nextflow pull <project name>

The project name should is composed by the `owner` name (or organisation name) and the `repository` name
separated by a ``/`` character e.g. ``nextflow-io/rnaseq-nf``. When the owner names is omitted the
value defined by environment variable ``NXF_ORG`` is used (default: ``nextflow-io``)::

    nextflow pull nextflow-io/rnaseq-nf

By default Nextflow checks for the specific project on `GitHub <https://github.com>`_ hosting platform.
Use the ``-hub`` option to use a different platform. Alternatively the full URL can be specified::

    nextflow pull https://github.com/nextflow-io/rnaseq-nf


**Options**

=============== =============
Name            Description
=============== =============
-all            Update all downloaded projects
-hub            Service hub where the project is hosted: ``github``, ``gitlab`` or ``bitbucket`` (default: ``github``).
-r              Revision to clone - It can be a git branch, tag or revision number.
-user           Private repository user name.
=============== =============


run
===

Launch the execution of a workflow script or project. If the workflow has not been downloaded before, Nextflow will automatically pull the workflow from the respective hub platform and store it in `~/.nextflow/assets/URI`. 

**Options**

=============== =============
Name            Description
=============== =============
-E              Exports all current system environment
-bucket_dir     Remote bucket where intermediate result files are stored
-cache          Enable/disable processes caching
-dump-channels  Dump channels for debugging purposes
-dump-hashes    Dump task hash keys for debugging purposes
-e.key          Add the specified variable <key> to execution environment
-h              Print the command usage
-hub            Service hub where the project is hosted
-latest         Pull latest changes before run
-lib            Library extension path
-name           Assign a mnemonic name to the a pipeline run
-offline        Do not check for remote project updates
-params-file    Load script parameters from a JSON/YAML file
-process.       Set process options, Syntax: -process.key=value        
-profile        Choose a configuration profile
-qs/-queue-size Max number of processes that can be executed in parallel by each executor
-resume         Execute the script using the cached results, useful to continue executions that was stopped by an error
-r              Revision of the project to run (either a git branch, tag or commit SHA number)
-test           Test a script function with the name specified
-user           Private repository user name
-with-conda     Use the specified Conda environment package or file (must end with .yml|.yaml suffix)
-with-dag       Create pipeline DAG file
-with-docker    Enable process execution in a Docker container
-N              Send a notification email on workflow completion to the specified recipients
-with-report    Create processes execution html report
-with-singularity   Enable process execution in a Singularity container
-with-timeline  Create processes execution timeline file
-with-trace     Create processes execution tracing file
-with-weblog    Send workflow status messages via HTTP to target URL
-without-docker Disable process execution with Docker
-w              Directory where intermediate result files are stored

**** params-file ****

Users can specify either a JSON or YAML based file containing parameters for a Nextflow workflow using the option ``-params-file``. The YAML format should follow a key - value pattern, e.g. 
    foo:
    - bar: "2.0"
    baz:
    - type: "test"
Similar to this, you could specify a JSON file with parameters:
    {
     "foo": "bar",
    "baz": "1.0"
    }
The above feature requires version 0.24.0 or higher.


self-update
===========

Update Nextflow runtime to the latest available version.

Example:: 

    nextflow self-update



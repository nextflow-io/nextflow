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

**Example**
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

drop
====

Delete the local copy of a project


help
====

Print the usage help for a command



info
====

Print a project and system runtime information



list
====

List all downloaded projects




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

Download or update a project



run
===

Launch the execution of workflow script or project

self-update
===========

Update Nextflow runtime to the latest available version.

Example:: 

    nextflow self-update



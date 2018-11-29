.. _cli-page:

**********************
Command Line Interface
**********************

Command line
============

Nextflow has different command line arguments that can be used to modify the execution of workflows or interact with existing runs.  


Nextflow Docker
---------------

Nextflow can be launched using a provided Docker container. This container has been modified to allow additional sibling containers to be launched from within the Nextflow container. Running Nextflow in this manner::

    nextflow -d run main.nf

clean
-----

Clean up project cache and work directories



clone
-----

Clone a project into a folder


cloud
-----

Manage Nextflow clusters in the cloud


config
------

Show a project configuration

drop
----

Delete the local copy of a project


help
----

Print the usage help for a command



info
----

Print a project and system runtime information



list
----

List all downloaded projects




log
---

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
----

Download or update a project



run
---

Launch the execution of workflow script or project

self-update
-----------

Update nextflow runtime to the latest available version


.. _cli-page:

*************
Command Line Interface
*************

Command line
==================

Nextflow has different command line arguments that can be used to modify the execution of workflows or interact with existing runs.  


Nextflow Docker
--------------

Nextflow can be launched using a provided Docker container. This container has been modified to allow additional sibling containers to be launched from within the Nextflow container. Running Nextflow in this manner::

    nextflow -d run main.nf

Nextflow log
--------------

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



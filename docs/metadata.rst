.. _metadata-page:

***********************
Workflow introspection
***********************

.. _metadata-workflow:

Runtime metadata
----------------

The implicit ``workflow`` object allows you to access some workflow and runtime metadata in your Nextflow scripts.
For example::

    println "Project : $workflow.projectDir"
    println "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"
    println "Cmd line: $workflow.commandLine"
    println "Manifest's pipeline version: $workflow.manifest.version"


.. tip:: To shortcut the access to multiple ``workflow`` properties you can use the Groovy
    `with <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/lang/Object.html#with(groovy.lang.Closure)>`_ method.


The following table lists the properties that can be accessed on the ``workflow`` object:

=========================== ===========================
Name                        Description
=========================== ===========================
scriptId                    Project main script unique hash ID.
scriptName                  Project main script file name.
scriptFile                  Project main script file path.
repository                  Project repository Git remote URL.
commitId                    Git commit ID of the executed workflow repository.
revision                    Git branch/tag of the executed workflow repository.
projectDir                  Directory where the workflow project is stored in the computer.
launchDir                   Directory where the workflow execution has been launched.
workDir                     Workflow working directory.
homeDir                     User system home directory.
userName                    User system account name.
configFiles                 Configuration files used for the workflow execution.
container                   Docker image used to run workflow tasks. When more than one image is used
                            it returns a map object containing `[process name, image name]` pair entries.
containerEngine             Returns the name of the container engine (e.g. docker or singularity) or null
                            if no container engine is enabled. 
commandLine                 Command line as entered by the user to launch the workflow execution.
profile                     Used configuration profile.
runName                     Mnemonic name assigned to this execution instance.
sessionId                   Unique identifier (UUID) associated to current execution.
resume                      Returns ``true`` whenever the current instance is resumed from a previous execution.
start                       Timestamp of workflow at execution start.
manifest                    Entries of the workflow manifest.
:sup:`✝` complete           Timestamp of workflow when execution is completed.
:sup:`✝` duration           Time elapsed to complete workflow execution.
:sup:`*` success            Reports if the execution completed successfully.
:sup:`*` exitStatus         Exit status of the task that caused the workflow execution to fail.
:sup:`*` errorMessage       Error message of the task that caused the workflow execution to fail.
:sup:`*` errorReport        Detailed error of the task that caused the workflow execution to fail.
=========================== ===========================

| Properties marked with a `✝` are accessible only in the workflow completion handler.
| Properties marked with a `*` are accessible only in the workflow completion and error handlers. See the `Completion handler`_ section for details.
|

.. _metadata-nextflow:

Nextflow metadata
-----------------

The implicit ``nextflow`` object allows you to access the metadata information of the Nextflow runtime.

=========================== ===========================
Name                        Description
=========================== ===========================
nextflow.version            Nextflow runtime version number.
nextflow.build              Nextflow runtime build number.
nextflow.timestamp          Nextflow runtime compile timestamp.
=========================== ===========================

The method ``nextflow.version.matches`` allows you to check if the Nextflow runtime satisfies the version
requirement eventually needed by your workflow script. The required version string can be prefixed with the usual
comparison operators eg ``>``, ``>=``, ``=``, etc. or postfixed with the ``+`` operator to specify a minimal version
requirement. For example::

    if( !nextflow.version.matches('0.22+') ) {
        println "This workflow requires Nextflow version 0.22 or greater -- You are running version $nextflow.version"
        exit 1
    }


.. _metadata-completion-handler:

Completion handler
------------------

Due to the asynchronous nature of Nextflow the termination of a script does not correspond to the termination
of the running workflow. Thus some information, only available on execution completion, needs to be accessed by
using an asynchronous handler.

The ``onComplete`` event handler is invoked by the framework when the workflow execution is completed. It allows one
to access the workflow termination status and other useful information. For example::

    workflow.onComplete {
        println "Pipeline completed at: $workflow.complete"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    }

If you want an e-mail notification on completion, check :ref:`mail-page`.

.. _metadata-error-handler:

Error handler
-------------

The ``onError`` event handler is invoked by Nextflow when a runtime or process error caused the pipeline execution to stop.
For example::

    workflow.onError {
        println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
    }

.. note:: Both the ``onError`` and ``onComplete`` handlers are invoked when an error condition is encountered.
    However the first is called as soon as the error is raised, while the second just before the pipeline execution
    is going terminate. When using the ``finish`` :ref:`process-page-error-strategy`, between the two there could be
    a significant time gap depending by the time required to complete any pending job.


Decoupling metadata
-----------------------

The workflow event handlers can be defined also in the ``nextflow.config`` file. This is useful to
decouple the handling of pipeline events from the main script logic.

When the event handlers are included in a configuration file the only difference is that the ``onComplete`` and
the ``onError`` closures have to be defined by using the assignment operator as shown below::

    workflow.onComplete = {
        // any workflow property can be used here
        println "Pipeline complete"
        println "Command line: $workflow.commandLine"
    }


    workflow.onError = {
        println "Oops .. something when wrong"
    }


.. note:: It is possible to define a workflow event handlers both in the pipeline script **and** in the
  configuration file.


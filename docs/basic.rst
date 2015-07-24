***************
Basic concepts
***************


`Nextflow` is a reactive workflow framework and a programming `DSL <http://en.wikipedia.org/wiki/Domain-specific_language>`_
that eases the writing of computational pipelines with complex data.

It is designed around the idea that the Linux platform is the lingua franca of data science. Linux provides many
simple command line and scripting tools, which by themselves are powerful, but when chained together facilitate complex
data manipulations.

`Nextflow` extends this approach adding to it the ability to define complex program interactions and an high-level
parallel computational environment based on the `dataflow` programming model.


Processes and channels
----------------------

In practice a Nextflow pipeline script is made by joining together many different processes.
Each process can be written in any scripting language that can be executed by the Linux platform (BASH, Perl, Ruby, Python, etc).

Processes are executed independently and are isolated from each other i.e. they do not share a common (writable) state.
The only way they can communicate is by using asynchronous FIFO queues, called `channels` in the Nextflow lingo.

Any process can define one or more channels as `input` and `output`. The interaction between these processes,
and ultimately the pipeline execution flow itself, is implicitly defined by these input and output declarations.

A Nextflow script looks like the following example::

    params.query = "$HOME/sample.fa"
    params.db = "$HOME/tools/blast-db/pdb/pdb"

    db = file(params.db)
    query = file(params.query)

    process blastSearch {
        input:
        file query

        output:
        file top_hits

        """
        blastp -db $db -query $query -outfmt 6 > blast_result
        cat blast_result | head -n 10 | cut -f 2 > top_hits
        """
    }


    process extractTopHits {
        input:
        file top_hits

        output:
        file sequences

        """
        blastdbcmd -db ${db} -entry_batch $top_hits > sequences
        """
    }



In the above example two processes are defined. The execution order is not set by the fact that ``blastSearch`` process comes
before the ``extractTopHits`` in the script (it could also be written the other way around) but because the first defines
the channel ``top_hits`` in its output declarations while the ``extractTopHits`` process defines the same channel in its
input declaration, and thus establishing a communication link from the `blastSearch` task towards the `extractTopHits` task.

.. TODO describe that both processes are launched at the same time

Read :ref:`Channel <channel-page>` and :ref:`Process <process-page>` sections to learn more about these features.


Execution abstraction
---------------------

While a process defines `what` command or script has to be executed, the `executor` has the role to determine `how`
that script is actually run in the target system.

If not specified otherwise processes are executed in the local computer. The local executor is very useful for pipeline
development and test purposes, but for real world computational pipelines, an HPC or cloud platform is required.

In other words, `Nextflow` provides an abstraction between the pipeline's functional logic and the underlying execution system.
Thus it is possible to write a pipeline once and have it running on your computer, a grid platform or the cloud seamlessly,
without modifying it, by simply defining the target execution platform in the configuration file.

The following HPC and cloud platforms are supported:

* `Open grid engine <http://gridscheduler.sourceforge.net/>`_
* `Univa grid engine <http://www.univa.com/>`_
* `Platform LSF <http://www.ibm.com/systems/technicalcomputing/platformcomputing/products/lsf/>`_
* `Linux SLURM <https://computing.llnl.gov/linux/slurm/>`_
* `PBS Works <http://www.pbsworks.com/gridengine/>`_
* `Torque <http://www.adaptivecomputing.com/products/open-source/torque/>`_
* `ClusterK <http://clusterk.com/>`_
* `DNAnexus <http://www.dnanexus.com>`_


Read :ref:`executor-page` section to learn more about Nextflow executors.


Scripting language
------------------

Although `Nextflow` is designed to be used with a minimal learning curve, without having to study
a new programming language and using your current skills, it also provides a powerful DSL scripting language.

Nextflow scripting is an extension of the `Groovy programming language <http://en.wikipedia.org/wiki/Groovy_(programming_language)>`_
which in turn is a super-set of the Java programming language, thus if you have some knowledge of these languages,
or even only some confidence with the `C/C++` syntax, you will be comfortable using it.

Read :ref:`pipeline-page` section to learn about Nextflow scripting language.



.. TODO Running pipeline


.. TODO Pipeline parameters


Configuration options
---------------------

Pipeline configuration properties are defined in a file named ``nextflow.config`` in the pipeline execution directory. 

It allows to set what executor to use, the processes environment variables, pipeline parameters etc. 

A basic configuration file looks like the following example::

	process { 
	  executor='sge'
	  queue = 'cn-el6' 
	}

	env {
	  PATH="$PWD/bowtie2:$PWD/tophat2:$PATH"
	}

Read :ref:`config-page` section to learn more about the Nextflow configuration file and settings.




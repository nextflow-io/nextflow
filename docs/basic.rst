***************
Basic concepts
***************


`Nextflow` is a reactive workflow framework and a programming `DSL <http://en.wikipedia.org/wiki/Domain-specific_language>`_
that eases the writing of data-intensive computational pipelines.

It is designed around the idea that the Linux platform is the lingua franca of data science. Linux provides many
simple but powerful command-line and scripting tools that, when chained together, facilitate complex
data manipulations.

`Nextflow` extends this approach, adding the ability to define complex program interactions and a high-level
parallel computational environment based on the `dataflow` programming model.


Processes and channels
----------------------

In practice a Nextflow pipeline script is made by joining together different processes.
Each process can be written in any scripting language that can be executed by the Linux platform (Bash, Perl, Ruby, Python, etc.).

Processes are executed independently and are isolated from each other, i.e. they do not share a common (writable) state.
The only way they can communicate is via asynchronous FIFO queues, called `channels` in Nextflow.

Any process can define one or more channels as `input` and `output`. The interaction between these processes,
and ultimately the pipeline execution flow itself, is implicitly defined by these input and output declarations.

A Nextflow script looks like this::

    // Script parameters
    params.query = "/some/data/sample.fa"
    params.db = "/some/path/pdb"

    db = file(params.db)
    query_ch = Channel.fromPath(params.query)

    process blastSearch {
        input:
        file query from query_ch

        output:
        file "top_hits.txt" into top_hits_ch

        """
        blastp -db $db -query $query -outfmt 6 > blast_result
        cat blast_result | head -n 10 | cut -f 2 > top_hits.txt
        """
    }

    process extractTopHits {
        input:
        file top_hits from top_hits_ch

        output:
        file "sequences.txt" into sequences_ch

        """
        blastdbcmd -db $db -entry_batch $top_hits > sequences.txt
        """
    }



The above example defines two processes. Their execution order is not determined by the fact that the ``blastSearch``
process comes before ``extractTopHits`` in the script (it could also be written the other way around).

Instead, because the first process defines the channel ``top_hits_ch`` in its output declarations, and the
process ``extractTopHits`` defines the channel in its input declaration, a communication link is established.

This linking via the channels means that `extractTopHits` is waiting for the output of `blastSearch`, and then
runs `reactively` when the channel has contents.

.. TODO describe that both processes are launched at the same time

Read the :ref:`Channel <channel-page>` and :ref:`Process <process-page>` sections to learn more about these features.


Execution abstraction
---------------------

While a process defines `what` command or script has to be executed, the `executor` determines `how`
that script is actually run on the target system.

If not otherwise specified, processes are executed on the local computer. The local executor is very useful for pipeline
development and testing purposes, but for real world computational pipelines an HPC or cloud platform is often required.

In other words, `Nextflow` provides an abstraction between the pipeline's functional logic and the underlying execution system.
Thus it is possible to write a pipeline once and to seamlessly run it on your computer, a grid platform, or the cloud,
without modifying it, by simply defining the target execution platform in the configuration file.

The following batch schedulers are supported:

* `Open grid engine <http://gridscheduler.sourceforge.net/>`_
* `Univa grid engine <http://www.univa.com/>`_
* `Platform LSF <http://www.ibm.com/systems/technicalcomputing/platformcomputing/products/lsf/>`_
* `Linux SLURM <https://computing.llnl.gov/linux/slurm/>`_
* `PBS Works <http://www.pbsworks.com/gridengine/>`_
* `Torque <http://www.adaptivecomputing.com/products/open-source/torque/>`_
* `HTCondor <https://research.cs.wisc.edu/htcondor/>`_


The following cloud platforms are supported:

* `Amazon Web Services (AWS) <https://aws.amazon.com/>`_
* `Google Cloud Platform (GCP) <https://cloud.google.com/>`_
* `Kubernetes <https://kubernetes.io/>`_

Read the :ref:`executor-page` to learn more about the Nextflow executors.


Scripting language
------------------

`Nextflow` is designed to have a minimal learning curve, without having to pick up
a new programming language. In most cases, users can utilise their current skills to develop
Nextflow workflows. However, it also provides a powerful scripting DSL.

Nextflow scripting is an extension of the `Groovy programming language <http://en.wikipedia.org/wiki/Groovy_(programming_language)>`_,
which in turn is a super-set of the Java programming language. Groovy can be considered as Python for Java
in that is simplifies the writing of code and is more approachable.

Read the :ref:`script-page` section to learn about the Nextflow scripting language.


.. TODO Running pipeline


.. TODO Pipeline parameters


Configuration options
---------------------

Pipeline configuration properties are defined in a file named ``nextflow.config`` in the pipeline execution directory. 

This file can be used to define which executor to use, the process's environment variables, pipeline parameters etc. 

A basic configuration file might look like this::

	process { 
	  executor='sge'
	  queue = 'cn-el6' 
	}


Read the :ref:`config-page` section to learn more about the Nextflow configuration file and settings.




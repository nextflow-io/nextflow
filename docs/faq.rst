.. _faq-page:

***
FAQ
***

How do I process multiple input files in parallel?
--------------------------------------------------

Q: *I have a collection of input files (e.g. carrots.fa, onions.fa, broccoli.fa). How can I specify that a process is performed on each input file in a parallel manner?*

A: The idea here is to create a ``channel`` that will trigger a process
execution for each of your files. First define a parameter that specifies where
the input files are:

::

    params.input = "data/*.fa"

Each of the files in the data directory can be made into a channel with:

::

    vegetable_datasets = Channel.fromPath(params.input)

From here, each time the variable ``vegetable_datasets`` is called as an
input to a process, the process will be performed on each of the files
in the vegetable datasets. For example, each input file may contain a
collection of unaligned sequences. We can specify a process to align
them as follows:

::

    process clustalw2_align {
        input:
        file vegetable_fasta from vegetable_datasets

        output:
        file "${vegetable_fasta.baseName}.aln" into vegetable_alns

        script:
        """
        clustalw2 -INFILE=${vegetable_fasta}
        """
    }

This would result in the alignment of the three vegetable fasta files
into ``carrots.aln``, ``onions.aln`` and ``broccoli.aln``.

These aligned files are now in the channel ``vegetable_alns`` and can be
used as input for a further process.

How do I get a unique ID based on the file name?
------------------------------------------------

*Q: How do I get a unique identifier based on a dataset file names (e.g. broccoli from broccoli.fa) and have the results going to a specific folder (e.g. results/broccoli/)?*

A: First we can specify a results directory as shown below:

::

    results_path = $PWD/results

The best way to manage this is to have the channel emit a tuple
containing both the file base name (``broccoli``) and the full file path
(``data/broccoli.fa``):

::

    datasets = Channel
                    .fromPath(params.input)
                    .map { file -> tuple(file.baseName, file) }

And in the process we can then reference these variables (``datasetID``
and ``datasetFile``):

::

    process clustalw2_align {
        publishDir "$results_path/$datasetID"

        input:
        set datasetID, file(datasetFile) from datasets

        output:
        set datasetID, file("${datasetID}.aln") into aligned_files

        script:
        """
        clustalw2 -INFILE=${datasetFile} -OUTFILE=${datasetID}.aln
        """
    }

In our example above would now have the folder ``broccoli`` in the results directory which would
contain the file ``broccoli.aln``.

If the input file has multiple extensions (e.g. ``brocolli.tar.gz``), you will want to use
``file.simpleName`` instead, to strip all of them (available since Nextflow 0.25+).


How do I use the same channel multiple times?
---------------------------------------------

*Q: Can a channel be used in two input statements? For example, I want carrots.fa to be aligned by both ClustalW and T-Coffee.*


A: A channel can be consumed only by one process or operator (except if channel only ever contains one item). You must
duplicate a channel before calling it as an input in different processes.
First we create the channel emitting the input files:

::

    vegetable_datasets = Channel.fromPath(params.input)

Next we can split it into two channels by using the :ref:`operator-into` operator:

::

    vegetable_datasets.into { datasets_clustalw; datasets_tcoffee }

Then we can define a process for aligning the datasets with *ClustalW*:

::

    process clustalw2_align {
        input:
        file vegetable_fasta from datasets_clustalw
        
        output:
        file "${vegetable_fasta.baseName}.aln" into clustalw_alns

        script:
        """
        clustalw2 -INFILE=${vegetable_fasta}
        """
    }

And a process for aligning the datasets with *T-Coffee*:

::

    process tcoffee_align {
        input:
        file vegetable_fasta from datasets_tcoffee
        
        output:
        file "${vegetable_fasta.baseName}.aln" into tcoffee_alns

        script:
        """
        t_coffee ${vegetable_fasta}
        """
    }

The upside of splitting the channels is that given our three unaligned
fasta files (``broccoli.fa``, ``onion.fa`` and ``carrots.fa``) six
alignment processes (three x ClustalW) + (three x T-Coffee) will be
executed as parallel processes.


How do I invoke custom scripts and tools?
-----------------------------------------

*Q: I have executables in my code, how should I call them in Nextflow?*

A: Nextflow will automatically add the directory ``bin`` into the ``PATH``
environmental variable. So therefore any executable in the ``bin``
folder of a Nextflow pipeline can be called without the need to
reference the full path.

For example, we may wish to reformat our *ClustalW* alignments from
Question 3 into *PHYLIP* format. We will use the handy tool
``esl-reformat`` for this task.

First we place copy (or create a symlink to) the ``esl-reformat``
executable to the project's bin folder. From above we see the *ClustalW*
alignments are in the channel ``clustalw_alns``:

::

    process phylip_reformat {
        input:
        file clustalw_alignment from clustalw_alns
        
        output:
        file "${clustalw_alignment.baseName}.phy" to clustalw_phylips

        script:
        """
        esl-reformat phylip ${clustalw_alignment} ${clustalw_alignment.baseName}.phy
        """
    }


    process generate_bootstrap_replicates {
        input:
        file clustalw_phylip from clustalw_phylips
        output:
        file "${clustalw_alignment.baseName}.phy" to clustalw_phylips

        script:
        """
        esl-reformat phylip ${clustalw_alignment} ${clustalw_alignment.baseName}.phy
        """
    }

How do I iterate over a process n times?
-----------------------------------------

To perform a process *n* times, we can specify the input to be
``each x from y..z``. For example:

::

    bootstrapReplicates=100

    process bootstrapReplicateTrees {
        publishDir "$results_path/$datasetID/bootstrapsReplicateTrees"

        input:
        each x from 1..bootstrapReplicates
        set val(datasetID), file(ClustalwPhylips)

        output:
        file "bootstrapTree_${x}.nwk" into bootstrapReplicateTrees

        script:
        // Generate Bootstrap Trees

        """
        raxmlHPC -m PROTGAMMAJTT -n tmpPhylip${x} -s tmpPhylip${x}
        mv "RAxML_bestTree.tmpPhylip${x}" bootstrapTree_${x}.nwk
        """
    }


How do I iterate over nth files from within a process?
------------------------------------------------------

*Q: For example, I have 100 files emitted by a channel. I wish to perform one process where I iterate over each file inside the process.*

A: The idea here to transform a channel emitting multiple items into a channel
that will collect all files into a list object and produce that list as a single emission. We do this using the ``collect()`` operator. The process script would then be able to iterate over
the files by using a simple for-loop.

This is also useful if all the items of a channel are required to be in the work directory.

::

    process concatenateBootstrapReplicates {
        publishDir "$results_path/$datasetID/concatenate"

        input:
        file bootstrapTreeList from bootstrapReplicateTrees.collect()

        output:
        file "concatenatedBootstrapTrees.nwk"

        // Concatenate Bootstrap Trees
        script:
        """
        for treeFile in ${bootstrapTreeList}
        do
            cat \$treeFile >> concatenatedBootstrapTrees.nwk
        done

        """
    }

How do I use a specific version of Nextflow?
------------------------------------------------------

*Q: I need to specify a version of Nextflow to use, or I need to pull a snapshot release.*

A: Sometimes it is necessary to use a different version of Nextflow for a specific feature or testing purposes. Nextflow is able to automatically pull versions when the ``NXF_VER`` environment variable is defined on the commandline. 

::

    NXF_VER=0.28.0 nextflow run main.nf


.. _faq-no-resource-enforcement:

How does nextflow allocate resources to processes?
------------------------------------------------------

*Q: How does nextflow deal with resource allocation? Can I configure nextflow to make sure that my workflow doesn't overload the resources available to it?*

A: Nextflow sets processes to run and collects the results, it does **not** control what they request. Nextflow does **not** enforce resource allocation or provide resource isolation (nor do containers). For example, if your process declares, for example, 2 cpus and the task uses 20, on a compute node with 8 cpus - then, you will overload the running node. 

Nextflow configuration files allow you to set expectations for resource uses, even though it does not enforce these. For example, to use a local executor expecting 4 cpus and 16GB memory, add to your configuration file: 

```
executor {
    name   = 'local'
    cpus   = 4
    memory = '16GB'
}
```

This tells nextflow to run up to 4 *processes* at once, however each of the processes that nextflow starts may (without telling nextflow) itself request multiple threads.

Multi-thread applications usually have a command line option to specify the maximum number of cpus that may be used. For example, [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) has the `-p/--threads` option, [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) has the `-j/--cores` option, and so on. For a specific workflow, you can specify the number of threads for your process to use, and match that to the process-specific entry in your nextflow config file. Unless you explicitly specify the task cpus for a multi-thread application, it is possible that each of them will try to use all the available cpus. In this case, the total number of cpus could exceed the number specified by the executor. Essentially, nextflow only "knows" how many tasks to run, and some of your tasks, if multi-threaded, may request more resources than nextflow expects. 

Also, it is not common for (bioinformatics) applications to allow you to specify maximum memory usage.

In order to run a workflow on a large dataset without overloading requested resources, it can be helpful to run a workflow on a small "test" dataset and track cpu/memory usage from there. Then, afterwards, use that information to adjust which resources are requested for the full-size dataset.

This is likely to be a more severe problem using a [local executor](/docs/latest/executor.html#local). Using a nextflow executor designed to run with a job submission system on a cluster (SGE/LSF/etc) allows each process to be despatched, and resources allocated, using the job submission system.

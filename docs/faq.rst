.. _faq-page:

***
FAQ
***

How do I process multiple input files in parallel?
--------------------------------------------------

Q: *I have a collection of input files (e.g. carrots.fa, onions.fa, broccoli.fa). How can I specify that a process is performed on each input file in a parallel manner?*

A: The idea here is to create a ``channel`` that will trigger a process
execution for each of your files. First define at the top of your script file
a parameter that specifies where the input files are:

::

    params.input = "data/*.fa"

Each of the files in the data directory can be made into a channel with:

::

    vegetable_datasets = Channel.fromPath(params.input)

This should be placed within the workflow block. From here, each time
the variable ``vegetable_datasets`` is called as an input to a process,
the process will be performed on each of the files in the vegetable datasets.
For example, each input file may contain a collection of unaligned sequences.
We can specify a process to align them as follows:

::

    process clustalw2_align {
        input:
        path vegetable_fasta

        output:
        path "${vegetable_fasta.baseName}.aln"

        script:
        """
        clustalw2 -INFILE=${vegetable_fasta}
        """
    }

And the workflow block with the logic as follows:

::

    workflow {
        vegetable_datasets = Channel.fromPath(params.input)
        vegetable_alns = clustalw2_align(vegetable_datasets)
    }

This would result in the alignment of the three vegetable fasta files
into ``carrots.aln``, ``onions.aln`` and ``broccoli.aln``.

These aligned files are now in the channel ``vegetable_alns`` and can be
used as input for a further process.

How do I get a unique ID based on the file name?
------------------------------------------------

*Q: How do I get a unique identifier based on a dataset file names (e.g. broccoli from broccoli.fa) and have the results going to a specific folder (e.g. results/broccoli/)?*

A: First we can specify a results directory as shown below, at the top of your
script file (along with the parameter definition for input files):

::

    results_path = $PWD/results
    params.input = "data/*.fa"

The best way to manage this is to have the channel emit a tuple
containing both the file base name (``broccoli``) and the full file path
(``data/broccoli.fa``):

::

    datasets = Channel
        .fromPath(params.input)
        .map { file -> tuple(file.baseName, file) }

In the process we can then reference these variables (``datasetID``
and ``datasetFile``):

::

    process clustalw2_align {
        publishDir "$results_path/$datasetID"

        input:
        tuple val(datasetID), path(datasetFile)

        output:
        tuple val(datasetID), path("${datasetID}.aln")

        script:
        """
        clustalw2 -INFILE=${datasetFile} -OUTFILE=${datasetID}.aln
        """
    }

With the logic written in the workflow block, as follows:

::

    workflow {
        datasets = Channel
            .fromPath(params.input)
            .map { file -> tuple(file.baseName, file) }
        aligned_files = clustalw2_align(datasets)
    }

In our example above, there would now be the folder ``broccoli`` in the results
directory which would contain the file ``broccoli.aln``.

If the input file has multiple extensions (e.g. ``broccoli.tar.gz``), you will
want to use ``file.simpleName`` instead, to strip all of them.


How do I invoke custom scripts and tools?
-----------------------------------------

*Q: I have executables in my code, how should I call them in Nextflow?*

A: Nextflow will automatically add the directory ``bin`` into the ``PATH``
environmental variable. So therefore any executable in the ``bin``
folder of a Nextflow pipeline can be called without the need to
reference the full path, as long as it has execution permissions (`chmod +x`).

For example, we may wish to reformat our *ClustalW* alignments into
*PHYLIP* format. We will use the handy tool ``esl-reformat`` for this task.

First we place a copy (or create a symlink) of the ``esl-reformat``
executable to the project's bin folder. From above we see the *ClustalW*
alignments are in the channel ``aligned_files``:

::

    process phylip_reformat {
        input:
        path clustalw_alignment

        output:
        path "${clustalw_alignment.baseName}.phy"

        script:
        """
        esl-reformat phylip ${clustalw_alignment} ${clustalw_alignment.baseName}.phy
        """
    }

    workflow {
        clustalw_phylips = phylip_reformat(aligned_files)
    }

How do I iterate over a process n times?
-----------------------------------------

To perform a process *n* times, we can use the ``each`` process directive.
For example:

::

    bootstrapReplicates=100

    process bootstrapReplicateTrees {
        publishDir "$results_path/$datasetID/bootstrapsReplicateTrees"

        input:
        each x
        tuple val(datasetID), path(ClustalwPhylips)

        output:
        path "bootstrapTree_${x}.nwk"

        script:
        // Generate Bootstrap Trees
        """
        raxmlHPC -m PROTGAMMAJTT -n tmpPhylip${x} -s tmpPhylip${x}
        mv "RAxML_bestTree.tmpPhylip${x}" bootstrapTree_${x}.nwk
        """
    }

    workflow {
        Channel
            .of(1..bootstrapReplicates)
            .set { x }
        datasets = Channel
            .fromPath(params.input)
            .map { file -> tuple(file.baseName, file) }
        clustalw2_align(x, datasets)

    }


How do I iterate over nth files from within a process?
------------------------------------------------------

*Q: For example, I have 100 files emitted by a channel. I wish to perform one process where I iterate over each file inside the process.*

A: The idea here is to transform a channel emitting multiple items into a
channel that will collect all files into a list object and produce that list
as a single emission. We do this using the ``collect()`` operator. The process
script would then be able to iterate over the files by using a simple for-loop.

This is also useful if all the items of a channel are required to be in the
work directory.

::

    process concatenateBootstrapReplicates {
        publishDir "$results_path/$datasetID/concatenate"

        input:
        path bootstrapTreeList

        output:
        path "concatenatedBootstrapTrees.nwk"

        // Concatenate Bootstrap Trees
        script:
        """
        for treeFile in ${bootstrapTreeList}
        do
            cat \$treeFile >> concatenatedBootstrapTrees.nwk
        done
        """
    }

    workflow {
        concatenateBootstrapReplicates(bootstrapReplicateTrees.collect())
    }

How do I use a specific version of Nextflow?
------------------------------------------------------

*Q: I need to specify a version of Nextflow to use, or I need to pull a snapshot release.*

A: Sometimes it is necessary to use a different version of Nextflow for a
specific feature or testing purposes. Nextflow is able to automatically pull
versions when the ``NXF_VER`` environment variable is defined on the
command line.

::

    NXF_VER=21.04 nextflow run main.nf

(overview-page)=

# Overview

## Why Nextflow?

The rise of big data has made it increasingly necessary to be able to analyze and perform experiments on large datasets in a portable and reproducible manner. Parallelization and distributed computing are the best ways to tackle this challenge, but the tools commonly available to computational scientists often lack good support for these techniques, or they provide a model that fits poorly with the needs of computational scientists and often require knowledge of complex tools and APIs. Nextflow was created to address these challenges.

The Nextflow language is inspired by [the Unix philosophy](https://en.wikipedia.org/wiki/Unix_philosophy), in which many simple command line tools can be chained together into increasingly complex tasks. Similarly, a Nextflow script consists of composing many simple processes into increasingly complex pipelines. Each process executes a given tool or scripting language, and by specifying the process inputs and outputs, Nextflow coordinates the execution of tasks for you.

The Nextflow runtime integrates with many popular execution platforms (HPC schedulers, cloud providers) and software tools (Git, Docker, Conda), allowing you to fully describe a computational pipeline with all of its dependencies and run it in nearly any environment -- write once, run anywhere.

## Processes and dataflow

In practice a Nextflow pipeline script is made by joining together different processes. Each process can be written in any scripting language that can be executed by the Linux platform (Bash, Perl, Ruby, Python, etc.).

A process can define one or more *inputs* and *outputs*. Data flows from process to process through asynchronous dataflow structures, known as *channels* and *values* in Nextflow. The data dependencies between these processes implicitly determines the flow of execution.

A Nextflow script looks like this:

```nextflow
// Script parameters
params.query = "/some/data/sample.fa"
params.db = "/some/path/pdb"

process blast_search {
  input:
  path query
  path db

  output:
  path "top_hits.txt"

  script:
  """
  blastp -db $db -query $query -outfmt 6 > blast_result
  cat blast_result | head -n 10 | cut -f 2 > top_hits.txt
  """
}

process extract_top_hits {
  input:
  path top_hits
  path db

  output:
  path "sequences.txt"

  script:
  """
  blastdbcmd -db $db -entry_batch $top_hits > sequences.txt
  """
}

workflow {
  def query_ch = channel.fromPath(params.query)
  blast_search(query_ch, params.db)
  extract_top_hits(blast_search.out, params.db).view()
}
```

The above example defines two processes. Their execution order is not determined by the fact that the `blast_search` process comes before `extract_top_hits` in the script (it could also be written the other way around). Instead, execution order is determined by their _dependencies_ -- `extract_top_hits` depends on the output of `blast_search`, so `blast_search` will be executed first, and then `extract_top_hits`.

When the workflow is executed, it creates two processes (`blast_search` and `extract_top_hits`) connected by the channel `query_ch`. Each process executes a task and emits a value for each input that it receives. Whenever `blast_search` emits a value, `extract_top_hits` receives it through the `query_ch` channel.

See {ref}`process-page`, {ref}`dataflow-page`, and {ref}`workflow-page` to learn more about these features.

## Execution abstraction

While a process defines *what* command or script has to be executed, the *executor* determines *how* that script is actually run on the target system.

If not otherwise specified, processes are executed on the local computer. The local executor is very useful for pipeline development and testing purposes, but for real world computational pipelines an HPC or cloud platform is often required.

In other words, Nextflow provides an abstraction between the pipeline's functional logic and the underlying execution system. Thus it is possible to write a pipeline once and to seamlessly run it on your computer, a grid platform, or the cloud, without modifying it, by simply defining the target execution platform in the configuration file.

The following batch schedulers are supported:

- [Open Grid Engine](http://gridscheduler.sourceforge.net/)
- [Univa Grid Engine](http://www.univa.com/)
- [Platform LSF](http://www.ibm.com/systems/technicalcomputing/platformcomputing/products/lsf/)
- [SLURM](https://computing.llnl.gov/linux/slurm/)
- [Flux Framework](https://flux-framework.org/)
- [PBS](http://www.pbsworks.com/gridengine/)
- [Torque](http://www.adaptivecomputing.com/products/open-source/torque/)
- [HTCondor](https://research.cs.wisc.edu/htcondor/)

The following cloud platforms are supported:

- [Amazon Web Services (AWS)](https://aws.amazon.com/)
- [Microsoft Azure](https://azure.microsoft.com/)
- [Google Cloud Platform (GCP)](https://cloud.google.com/)
- [Kubernetes](https://kubernetes.io/)

Read the {ref}`executor-page` to learn more about the Nextflow executors.

## Scripting language

Nextflow is a workflow language based on [Java](https://en.wikipedia.org/wiki/Java_(programming_language)) and [Groovy](https://groovy-lang.org/). It is designed to simplify writing scalable and reproducible pipelines. In most cases, users can leverage their existing programming skills to develop Nextflow pipelines without the steep learning curve that usually comes with a new programming language.

See {ref}`script-page` for more information about the Nextflow scripting language.

## Configuration options

Pipeline configuration properties are defined in a file named `nextflow.config` in the pipeline execution directory.

This file can be used to define which executor to use, the process's environment variables, pipeline parameters etc.

A basic configuration file might look like this:

```groovy
process {
  executor = 'sge'
  queue = 'cn-el6'
}
```

Read the {ref}`config-page` section to learn more about the Nextflow configuration file and settings.

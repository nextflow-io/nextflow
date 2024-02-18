<!-- FROM basic.md -->

# Concepts

Nextflow is a reactive workflow framework and a programming [DSL](http://en.wikipedia.org/wiki/Domain-specific_language) that eases the writing of data-intensive computational pipelines.

It is designed around the idea that the Linux platform is the lingua franca of data science. Linux provides many simple but powerful command-line and scripting tools that, when chained together, facilitate complex data manipulations.

Nextflow extends this approach, adding the ability to define complex program interactions and a high-level parallel computational environment based on the *dataflow* programming model.

## Processes and channels

In practice a Nextflow pipeline script is made by joining together different processes. Each process can be written in any scripting language that can be executed by the Linux platform (Bash, Perl, Ruby, Python, etc.).

Processes are executed independently and are isolated from each other, i.e. they do not share a common (writable) state. The only way they can communicate is via asynchronous FIFO queues, called *channels* in Nextflow.

Any process can define one or more channels as *input* and *output*. The interaction between these processes, and ultimately the pipeline execution flow itself, is implicitly defined by these input and output declarations.

A Nextflow script looks like this:

```groovy
// Declare syntax version
nextflow.enable.dsl=2

// Script parameters
params.query = "/some/data/sample.fa"
params.db = "/some/path/pdb"

process blastSearch {
  input:
    path query
    path db
  output:
    path "top_hits.txt"

    """
    blastp -db $db -query $query -outfmt 6 > blast_result
    cat blast_result | head -n 10 | cut -f 2 > top_hits.txt
    """
}

process extractTopHits {
  input:
    path top_hits

  output:
    path "sequences.txt"

    """
    blastdbcmd -db $db -entry_batch $top_hits > sequences.txt
    """
}

workflow {
  def query_ch = Channel.fromPath(params.query)
  blastSearch(query_ch, params.db) | extractTopHits | view
}
```

The above example defines two processes. Their execution order is not determined by the fact that the `blastSearch` process comes before `extractTopHits` in the script (it could also be written the other way around). Instead, the pipe operator (`|`) in the workflow between `blastSearch` and `extractTopHits` forwards the outputs from one process to the inputs of the following one.

When the workflow is started, it will create two processes and one channel (`query_ch`) and it will link all of them. Both processes will be started at the same time and they will listen to their respective input channels. Whenever `blastSearch` emits a value, `extractTopHits` will receive it (i.e. `extractTopHits` consumes the channel in a *reactive* way).

Read the {ref}`Channel <channel-page>` and {ref}`Process <process-page>` sections to learn more about these features.

## Execution abstraction

While a process defines *what* command or script has to be executed, the *executor* determines *how* that script is actually run on the target system.

If not otherwise specified, processes are executed on the local computer. The local executor is very useful for pipeline development and testing purposes, but for real world computational pipelines an HPC or cloud platform is often required.

In other words, Nextflow provides an abstraction between the pipeline's functional logic and the underlying execution system. Thus it is possible to write a pipeline once and to seamlessly run it on your computer, a grid platform, or the cloud, without modifying it, by simply defining the target execution platform in the configuration file.

The following batch schedulers are supported:

- [Open grid engine](http://gridscheduler.sourceforge.net/)
- [Univa grid engine](http://www.univa.com/)
- [Platform LSF](http://www.ibm.com/systems/technicalcomputing/platformcomputing/products/lsf/)
- [Linux SLURM](https://computing.llnl.gov/linux/slurm/)
- [Flux Framework](https://flux-framework.org/)
- [PBS Works](http://www.pbsworks.com/gridengine/)
- [Torque](http://www.adaptivecomputing.com/products/open-source/torque/)
- [HTCondor](https://research.cs.wisc.edu/htcondor/)

The following cloud platforms are supported:

- [Amazon Web Services (AWS)](https://aws.amazon.com/)
- [Google Cloud Platform (GCP)](https://cloud.google.com/)
- [Kubernetes](https://kubernetes.io/)

Read the {ref}`executor-page` to learn more about the Nextflow executors.

## Scripting language

Nextflow is designed to have a minimal learning curve, without having to pick up a new programming language. In most cases, users can utilise their current skills to develop Nextflow workflows. However, it also provides a powerful scripting DSL.

Nextflow scripting is an extension of the [Groovy programming language](<http://en.wikipedia.org/wiki/Groovy_(programming_language)>), which in turn is a super-set of the Java programming language. Groovy can be considered as Python for Java in that it simplifies the writing of code and is more approachable.

Read the {ref}`script-page` section to learn about the Nextflow scripting language.

<!-- FROM getstarted.md -->

(getstarted-first)=

## Your first script

Copy the following example into your favorite text editor and save it to a file named `tutorial.nf`:

```{literalinclude} ../../snippets/your-first-script.nf
:language: groovy
```

:::{note}
For versions of Nextflow prior to `22.10.0`, you must explicitly enable DSL2 by adding `nextflow.enable.dsl=2` to the top of the script or by using the `-dsl2` command-line option.
:::

This script defines two processes. The first splits a string into 6-character chunks, writing each one to a file with the prefix `chunk_`, and the second receives these files and transforms their contents to uppercase letters. The resulting strings are emitted on the `result` channel and the final output is printed by the `view` operator.

Execute the script by entering the following command in your terminal:

```console
$ nextflow run tutorial.nf

N E X T F L O W  ~  version 22.10.0
executor >  local (3)
[69/c8ea4a] process > splitLetters   [100%] 1 of 1 ✔
[84/c8b7f1] process > convertToUpper [100%] 2 of 2 ✔
HELLO
WORLD!
```

You can see that the first process is executed once, and the second twice. Finally the result string is printed.

It's worth noting that the process `convertToUpper` is executed in parallel, so there's no guarantee that the instance processing the first split (the chunk `Hello`) will be executed before the one processing the second split (the chunk `world!`).

Thus, it is perfectly possible that you will get the final result printed out in a different order:

```
WORLD!
HELLO
```

:::{tip}
The hexadecimal string, e.g. `22/7548fa`, is the unique hash of a task, and the prefix of the directory where the task is executed. You can inspect a task's files by changing to the directory `$PWD/work` and using this string to find the specific task directory.
:::

(getstarted-resume)=

### Modify and resume

Nextflow keeps track of all the processes executed in your pipeline. If you modify some parts of your script, only the processes that are actually changed will be re-executed. The execution of the processes that are not changed will be skipped and the cached result used instead.

This helps a lot when testing or modifying part of your pipeline without having to re-execute it from scratch.

For the sake of this tutorial, modify the `convertToUpper` process in the previous example, replacing the process script with the string `rev $x`, so that the process looks like this:

```groovy
process convertToUpper {
  input:
    path x
  output:
    stdout

  """
  rev $x
  """
}
```

Then save the file with the same name, and execute it by adding the `-resume` option to the command line:

```bash
nextflow run tutorial.nf -resume
```

It will print output similar to this:

```
N E X T F L O W  ~  version 22.10.0
executor >  local (2)
[69/c8ea4a] process > splitLetters   [100%] 1 of 1, cached: 1 ✔
[d0/e94f07] process > convertToUpper [100%] 2 of 2 ✔
olleH
!dlrow
```

You will see that the execution of the process `splitLetters` is actually skipped (the process ID is the same), and its results are retrieved from the cache. The second process is executed as expected, printing the reversed strings.

:::{tip}
The pipeline results are cached by default in the directory `$PWD/work`. Depending on your script, this folder can take up a lot of disk space. It's a good idea to clean this folder periodically, as long as you know you won't need to resume any pipeline runs.
:::

For more information, see the {ref}`cache-resume-page` page.

### Pipeline parameters

Pipeline parameters are simply declared by prepending to a variable name the prefix `params`, separated by dot character. Their value can be specified on the command line by prefixing the parameter name with a double dash character, i.e. `--paramName`

For the sake of this tutorial, you can try to execute the previous example specifying a different input string parameter, as shown below:

```bash
nextflow run tutorial.nf --str 'Bonjour le monde'
```

The string specified on the command line will override the default value of the parameter. The output will look like this:

```
N E X T F L O W  ~  version 22.10.0
executor >  local (4)
[8b/16e7d7] process > splitLetters   [100%] 1 of 1 ✔
[eb/729772] process > convertToUpper [100%] 3 of 3 ✔
m el r
edno
uojnoB
```

:::{versionchanged} 20.11.0-edge
Any `.` (dot) character in a parameter name is interpreted as the delimiter of a nested scope. For example, `--foo.bar Hello` will be interpreted as `params.foo.bar`. If you want to have a parameter name that contains a `.` (dot) character, escape it using the back-slash character, e.g. `--foo\.bar Hello`.
:::

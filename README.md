Nextflow
========

*"Dataflow variables are spectacularly expressive in concurrent programming"*
<br>[Henri E. Bal , Jennifer G. Steiner , Andrew S. Tanenbaum](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.145.7873)

Rationale
---------

With the rise of big data, techniques to analyse and run experiments on large datasets are increasingly necessary.

Parallelization and distributed computing are the best ways to tackle this kind of problem, but the tools commonly available
to the bioinformatics community traditionally lack good support for these techniques, or provide a model that fits
badly with the specific requirements in the bioinformatics domain and, most of the time, require the knowledge
of complex tools or low-level APIs.

Nextflow framework is based on the *dataflow* programming model, which greatly simplifies writing parallel and distributed pipelines
without adding unnecessary complexity and letting you concentrate on the flow of data, i.e. the functional logic of the application/algorithm.

It doesn't aim to be another pipeline scripting language yet, but it is built around the idea that the Linux platform 
is the *lingua franca* of data science, since it provides many simple command line and scripting tools, which by themselves 
are powerful, but when chained together facilitate complex data manipulations. 

In practice, this means that a Nextflow script is defined by composing  many different processes. 
Each process can be written in any scripting language that can be executed by the Linux platform (BASH, Perl, Ruby, Python, etc), 
to which is added the ability to coordinate and synchronize the processes execution by simply specifying their inputs and outputs.   

Quick start
-----------

Nextflow does not require any installation procedure, just download the distribution package by copying and pasting
this command in your terminal:

```
curl -fsSL get.nextflow.io | bash
```

It creates the ``nextflow`` executable file in the current directory. You may want to move it to a folder accessible from your ``$PATH``.

Create a file named `hello.nf` with the following content and copy it to the path where you downloaded the Nextflow package.

```nextflow
process sayHello {

    """
    printf 'Hello world! \n'
    """
}
```

Launch the above example by typing the following command on your terminal console:

```bash
./nextflow run -process.echo true hello.nf
```

Congratulations! You have just run your first program with Nextflow.


Something more useful
---------------------

Let's see a more real example: execute a BLAST search, get the top 10 hits, extract the found protein sequences and align them.

Copy the following example into a file named `pipeline.nf` .

```nextflow
params.query = "$HOME/sample.fa"
params.db = "$HOME/tools/blast-db/pdb/pdb"

process blast {
    output:
     file top_hits

    """
    blastp -query ${params.query} -db ${params.db} -outfmt 6 \
    | head -n 10 \
    | cut -f 2 > top_hits
    """
}

process extract {
    input:
     file top_hits
    output:
     file sequences

    "blastdbcmd -db ${params.db} -entry_batch $top_hits > sequences"
}

process align {
    input:
     file sequences
    echo true

    "t_coffee $sequences 2>&- | tee align_result"
}
```

The `input` and `output` declarations in each process, define what it is expecting to receive as input and what file(s)
are going to be produced as output.

Since the two variables `query` and `db` are prefixed by the `params` qualifier, their values can be overridden quickly
when the script is launched, by simply adding them on the Nextflow command line and prefixing them with the `--` characters.
For example:

```bash
./nextflow run pipeline.nf --db=/path/to/blast/db --query=/path/to/query.fasta
```


Mixing scripting languages
--------------------------

Processes in your pipeline can be written in any scripting language supported by the underlying Linux platform. To use a scripting
other than Linux BASH (e.g. Perl, Python, Ruby, R, etc), simply start your process script with the corresponding
<a href='http://en.wikipedia.org/wiki/Shebang_(Unix)' target='_bank'>shebang</a> declaration. For example:

```nextflow
process perlStuff {

    """
    #!/usr/bin/env perl

    print 'Hi there!' . '\n';
    """
}

process pyStuff {
    """
    #!/usr/bin/env python

    x = 'Hello'
    y = 'world!'
    print "%s - %s" % (x,y)
    """
}
```

Cluster Resource Managers support
---------------------------------

*Nextflow* provides an abstraction between the pipeline functional logic and the underlying processing system. 
Thus it is possible to write your pipeline once and have it running on your computer or a cluster resource
manager without modifying it. 

Currently the following clusters are supported: 
  
  + Open Grid Engine (SGE)
  + Univa Grid Engine
  + IBM Platform LSF
  + Linux SLURM
  + PBS/Torque
  + HTCondor (experimental)


By default processes are parallelized by spanning multiple threads in the machine where the pipeline is launched.

To submit the execution to a SGE cluster create a file named `nextflow.config`, in the directory
where the pipeline is going to be launched, with the following content: 

```nextflow
process {
  executor='sge'
  queue='<your execution queue>'
}
```

In doing that, processes will be executed as SGE jobs by using the `qsub` command, and so your pipeline will behave like any 
other SGE job script, with the benefit that *Nextflow* will automatically and transparently manage the processes 
synchronisation, file(s) staging/un-staging, etc.  

Alternatively the same declaration can be defined in the file `$HOME/.nextflow/config`, which is supposed to hold 
the global *Nextflow* configuration.

Cloud support
-------------

Nextflow provides out of the box support for the Amazon AWS cloud allowing you to setup a computing cluster,
deploy it and run your pipeline in the AWS infrastructure in a few commands.

The cloud configuration settings need to be specified in the `nextflow.config` file as shown below: 

    cloud {
          imageId = 'ami-43f49030'
          instanceType = 't2.micro'
          subnetId = 'subnet-05222a43'
          sharedStorageId = 'fs-1803efd1'
          spotPrice = 0.04 
    }

    aws {
        accessKey = 'xxx'
        secretKey = 'yyy'
        region = 'eu-west-1'
    }
    
    
Replace the settings in the above example with values of your choice. The attribute `sharedStorageId` is optional, 
when provided the [Amazon EFS](https://aws.amazon.com/efs/) file system is automatically mounted in the configured 
cloud environment. The `spotPrice` attribute allows you to use [EC2 Spot instances](https://aws.amazon.com/ec2/spot/) 
in place of regular on-request instances, bidding for the specified price.
 
The settings in the `aws` block can be omitted, in that case Nextflow will use the AWS credentials defined in 
your environment, using the standard AWS variables and configuration files. 

Once defined the configuration of your cloud environment, run the following command in the folder where the file 
`nextflow.config` was created: 

    nextflow cloud create my-cluster -c <num-of-nodes>
    
The string `my-cluster` identifies the cluster instance. Replace it with a name of your choice. 
Finally replace `num-of-nodes` with the actual number of instances that will made-up the cluster.
WARNING: you will be charged accordingly the type and the number of instances chosen.

Once the cluster deployment completes, SSH in the master node following the instruction that will be printed. Then 
you will be able to run your Nextflow pipeline as usual. 


Required dependencies
---------------------

* Compiler Java 8
* Runtime Java 8 or later

Build from source
-----------------

*Nextflow* is written in [Groovy](http://groovy-lang.org) (a scripting language for the JVM). A pre-compiled,
ready-to-run, package is available at the [Github releases page](https://github.com/nextflow-io/nextflow/releases),
thus it is not necessary to compile it in order to use it.

If you are interested in modifying the source code, or contributing to the project, it worth knowing that 
the build process is based on the [Gradle](http://www.gradle.org/) build automation system. 

You can compile *Nextflow* by typing the following command in the project home directory on your computer:

```bash
make compile
```

The very first time you run it, it will automatically download all the libraries required by the build process. 
It may take some minutes to complete.

When complete, execute the program by using the `launch.sh` script in the project directory.

The self-contained runnable Nextflow packages can be created by using the following command:

```bash
make pack
```

In order to install the compiled packages use the following command: 

```bash
make install
```

Then you will be able to run nextflow using the `nextflow` launcher script in the project root folder. 

Known compilation problems 
---------------------------

Nextflow required JDK 8 to be compiled. The Java compiler used by the build process can be choose by setting the
`JAVA_HOME` environment variable accordingly. 
 

If the compilation stops reporting the error: `java.lang.VerifyError: Bad <init> method call from inside of a branch`,
this is due to a bug affecting the following Java JDK:

- 1.8.0 update 11
- 1.8.0 update 20

Upgrade to a newer JDK to avoid to this issue. Alternatively a possible workaround is to define the following variable
in your environment:

```bash
_JAVA_OPTIONS='-Xverify:none'
```

Read more at these links:

- https://bugs.openjdk.java.net/browse/JDK-8051012
- https://jira.codehaus.org/browse/GROOVY-6951


IntelliJ IDEA
---------------

Nextflow development with [IntelliJ IDEA](https://www.jetbrains.com/idea/) requires the latest version of the IDE (2017.3 or higher).

If you have it installed in your computer, follow the steps below in order to use it with Nextflow:

1. Clone the Nextflow repository to a directory in your computer.
2. Open IntelliJ IDEA and choose "Import project" in the "File" menu bar.
3. Select the Nextflow project root directory in your computer and click "OK".
4. Then, choose the "Gradle" item in the "external module" list and click on "Next" button.
5. Confirm the default import options and click on "Finish" to finalize the project configuration.
6. When the import process complete, select the "Project structure" command in the "File" menu bar.
7. In the showed dialog click on the "Project" item in the list of the left, and make sure that
   the "Project SDK" choice on the right contains Java 8.


Documentation
--------------

Nextflow documentation is available at this link http://docs.nextflow.io

Contributing 
------------

Project contribution are more than welcome. See the [CONTRIBUTING](CONTRIBUTING.md) file for details.

Community
----------

You can post questions, or report problems by using the Nextflow [discussion forum](https://groups.google.com/forum/#!forum/nextflow)
or the [Nextflow channel on Gitter](https://gitter.im/nextflow-io/nextflow).


Build servers 
--------------

  * [Travis-CI](https://travis-ci.org/nextflow-io/nextflow)
  * [Groovy Joint build](http://ci.groovy-lang.org/project.html?projectId=JointBuilds_Nextflow&guest=1)

License
-------

The *Nextflow* framework is released under the Apache 2.0 license.

Citations
----------

If you use Nextflow for research purpose, please cite: 

P. Di Tommaso, et al. Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319 (2017) doi:[10.1038/nbt.3820](http://www.nature.com/nbt/journal/v35/n4/full/nbt.3820.html) 

Credits
-------

Nextflow is built on two great pieces of open source software, namely <a href='http://groovy-lang.org' target='_blank'>Groovy</a>
and <a href='http://www.gpars.org/' target='_blank'>Gpars</a>.

YourKit is kindly supporting this open source project with its full-featured Java Profiler.
Read more http://www.yourkit.com


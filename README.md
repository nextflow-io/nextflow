Nextflow
========

A *dataflow* oriented workflow framwork for bioinformatics pipeline

Rationale
---------

With the arise of big data, techniques to process and run experiments on large dataset are increasingly necessary.

Parallelization and distributed computing are the best way to tackle this problem, but the tools commonly available
to the bioinformaticians community traditionally lack of a good support of these techniques, or provide a model that fits
badly with the specific requirements in the bioinformatics domain, and most of the time, require the knowledge
of complex tools or low-level APIs.

Nextflow framework is based on the *dataflow* programming model, which greatly simplify writing parallel and distributed pipelines
without adding unnecessary complexity and letting to concentrate on the flow of data, i.e. the functional logic of the application/algorithm.



Quick start
-----------

Nextflow does not require any installation procedure, just download the <a href="http://dl.dropbox.com/u/376524/nextflow/nextflow">executable package here</a> and
save it somewhere on your computer.

Grants the execute permission to the donwload packages using the following command `chmod +x nextflow`, after that you are ready to use it.
You may try to execute the command `./nextflow -h` to show the program help.

Create a file named `hello.nf` with the following content and copy it
to the path where you downloaded the Nextflow executable.

    echo true

    task {
        '''
        #!/bin/sh
        printf 'Hello world! \n'
        '''
    }



Launch the above trivial example by typing the following command on your terminal console:

    ./nextflow -q hello.nf


Congratulations! You have just run your first task with Nextflow.


Something more usefull
----------------------

Let's see a more real example, execute BLAST search, get the top 10 hits, extract the found proteins sequences an align.
Copy the following example into a file named `pipeline.nf` .


    params.query = "$HOME/sample.fa"
    params.db = "$HOME/tools/blast-db/pdb/pdb"

    task ('blast') {
        output top_hits

        """
        blastp -db ${params.db} -query ${params.query} -outfmt 6 > blast_result
        cat blast_result | head -n 10 | cut -f 2 > top_hits
        """
    }

    task ('extract') {
        input top_hits
        output sequences

        "blastdbcmd -db ${params.db} -entry_batch $top_hits > sequences"
    }

    task ('align') {
        input sequences
        echo true

        "t_coffee $sequences 2>&- | tee align_result"
    }


The `input` and `output` declaration in each tasks, define what the task is expeting to receive as input and what file(s)
it is going to produce as output.

Since the two variables `query` and `db` are prefixed by the `params` qualifier, their values can be overriden quickly
when the script is launched, simply adding them on the Nextflow command line and prefixing them with the `--` characters.
For example:

    $ ./nextflow pipeline.nf --db=/path/to/blast/db --query=/path/to/query.fasta




Compile from sources
--------------------

The Nextflow build process is based on the Gradle build automation system. You can compile Nextflow by typing the
following command in the project home directory on your computer:

    $ ./gradlew compile

The very first time you run it, it will download automatically all the libraries required by the build process. It can requires
some minutes to complete.

When it completes, execute the program by using the `nextflow.sh` launch script in the project directory.

In order to create the self-contained executable package launch gradle specifying the *pack* task, using the follwing command:

    $ ./gradlew pack


Required dependencies
---------------------

Java 6 or higher


License
-------

The *Nextflow* framework source code is released under the GNU GPL3 License.


Credits
-------

Nextflow is built on two great pieces of open source software, namely <a href='groovy.codehaus.org' target='_blank'>Groovy</a>
and <a href='http://www.gpars.org/' target='_blank'>Gpars</a>

YourKit is kindly supporting this open source project with its full-featured Java Profiler.
Read more http://www.yourkit.com



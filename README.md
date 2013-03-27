NEXTFLOW
========

A workflow framework based on Dataflow programming model


Quick start
-----------

Download the <a href="http://dl.dropbox.com/u/376524/nextflow/nextflow.run">Nextflow executable package here</a> and
save it somewhere on your computer.

Create a file named `hello.nf` with the following content and copy it
to the path where you downloaded the Nextflow executable.

    echo true
    task {
        """
        echo 'Hello world'
        """
    }



If you are using a *nix system (Linux/MacOSX) grants the executable permission to the downloaded file with the command
`chmod +x ./nextflow.run`, after that you will able to run it using the below command:

    ./nextflow.run [program arguments]

If you are running a Windows OS, you will have to use the following syntax to run it:

    java -jar nextflow.run [program arguments]


Launch the above trivial example by typing the following command on your terminal console:

    ./nextflow.run hello.nf


Congratulations! You have just run your first task with Nextflow.


Build and run
-------------

The Nextflow build process is based on the Gradle build automation system. It can be download at the following link
http://gradle.org/


By having Gradle installed, you can compile Nextflow by typing the following command in the project home directory on your computer:

    $ gradle compile


Execute the program by using the `nextflow.sh` launch script in the project directory.


Required dependencies
---------------------

Java 6 or higher


License
-------

The *Nextflow* framework source code is released under the GNU GPL3 License.


Credits
-------

YourKit is kindly supporting this open source project with its full-featured Java Profiler.
Read more http://www.yourkit.com


Contact
-------
Paolo Di Tommaso - paolo (dot) ditommaso (at) gmail (dot) com
*******************
Get started
*******************

Requirements
============

`Nextflow` can be used on any POSIX compatible system (Linux, Solaris, OS X, etc.) and
requires `Java 7 <http://www.oracle.com/technetwork/java/javase/downloads/index.html>`_ to be installed.

Windows systems may be supported using a POSIX compatibility layer like `Cygwin <http://www.cygwin.com>`_ (unverified) or
in alternative installing it into a Linux VM using a virtualization software like `VirtualBox <http://www.virtualbox.org>`_
or `VMware <http://www.vmware.com/>`_.


Installation
============

Download the latest version by clicking the "`Nextflow`" green button on the `Github releases <https://github.com/paoloditommaso/nextflow/releases>`_ page.

`Nextflow` is distributed as a self-contained executable package, this means that it does not require any special installation procedure.

It only needs two easy steps:

#.  grant the execute permission to the downloaded file by using the command: ``chmod +x nextflow``
#.  save the ``nextflow`` file in a directory accessible by your ``$PATH`` variable
    (this is only required to avoid to remember and type the Nextflow full path each time you need to launch).



Verify it
==========

In order to verify that you have installed it successfully, just type ``nextflow`` and press enter on your
terminal command prompt. It will print the program help text, like showed below::


          N E X T F L O W
          Version 0.5.2 build 1109
          last modified 28-11-2013 18:33 UTC (19:33 CEST)
          http://nextflow.io

    Usage: nextflow [options] <script file> [script arguments]
      Options:
            --
           Set a parameter used by the workflow
           Syntax: --key=value
           Default: {}
        -E
           Exports all the current system environment
           Default: false
        -cache
           Enable/disable task(s) caching
           Default: true
        -c, -config
           Use the specified configuration file(s)
        -e
           Add the specified variable to execution environment
           Syntax: -ekey=value
           Default: {}
        -h, -help
           Show this help
           Default: false
        -history
           Show history of executed commands
           Default: false
        -lib
           Library extension path
        -log
           Define the application log file
           Default: .nextflow.log
        -process.
           Set default process options
           Syntax: -process.key=value
           Default: {}
        -qs, -queue-size
           The max number of task in execution queue
        -q, -quiet
           Do not print information messages
           Default: false
        -resume
           Execute the script using the cached results, useful to continue
           executions that stopped by an error
        -test
           Test the function with the name specified
        -v, -version
           Show the program version
           Default: false
        -w, -work-dir
           Directory where tasks results are stored
           Default: ./work



.. hint:: If you are a Linux hacker and you are wondering what contains the ``nextflow`` binary file,
    you may be interested to know that it is a Java executable JAR file that contains the application
    classes and libraries, moreover it is pre-prepended by a small BASH script fragment in order to make it self-runnable.

    You can also launch it as any other Java executable JAR, by using the command: ``java -jar <nextflow full path> [program options]``



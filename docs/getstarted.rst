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

`Nextflow` is distributed as a self-contained executable package, this means that it does not require any special installation procedure.

It only needs two easy steps:

#.  Download the executable package by copying and pasting the following command in your terminal
    window: ``wget -qO- get.nextflow.io | bash``.
    It will create the ``nextflow`` main executable in the current directory.

#.  Optionally, move the ``nextflow`` file in a directory accessible by your ``$PATH`` variable
    (this is only required to avoid to remember and type the Nextflow full path each time you need to run it).

.. tip:: In the case you don't have ``wget`` installed you can use the ``curl`` utility instead by entering
   the following command: ``curl -fsSL get.nextflow.io | bash``


Your first script
==================

By using your favourite text editor copy the following example, and save it to a file named ``tutorial.nf``. ::

    #!/usr/bin/env nextflow

    params.str = 'Hello world!'

    process splitLetters {

        output:
        file 'chunk_*' into letters mode flatten

        """
        printf '${params.str}' | split -b 6 - chunk_
        """
    }


    process massage {

        input:
        file x from letters

        output:
        stdout result

        """
        cat $x | tr '[a-z]' '[A-Z]'
        """
    }

    result.subscribe {
        println it.trim()
    }


Execute the script by entering the following command in your terminal::

   nextflow tutorial.nf

It will output something similar to the text shown below::

    N E X T F L O W  ~  version 0.7.0
    [warm up] executor > local
    [59/89680b] Submitted process > splitLetters (1)
    [42/78eb18] Submitted process > massage (2)
    [04/86d8e1] Submitted process > massage (1)
    WORLD!
    HELLO


You can see that the first process is execute once, it splits the string ``Hello world!`` in two files that are
received by the following process which is executed two times, transforming the content of the files to upper-case.

The resulting strings are emitted on the ``result`` channel and finally the output is printed out.


.. note:: The hexadecimal numbers like ``cb/52c058``


Pipeline parameters
--------------------


Modify and resume
------------------

For the sake of example, modify of the second process by replacing with script with the string ``rev $x``
so that the process will like the one shown below::

    process massage {

        input:
        file x from letters

        output:
        stdout result

        """
        cat $x | tr '[a-z]' '[A-Z]'
        """
    }

Then save the file with save name as before, and execute it by adding the ``-resume`` command line option::

    nextflow tutorial.nf -resume

It will print a similar output::

    N E X T F L O W  ~  version 0.7.0
    [warm up] executor > local
    [59/89680b] Cached process > splitLetters (1)
    [d0/7b79a3] Submitted process > massage (1)
    [b0/c99ef9] Submitted process > massage (2)
    olleH
    !dlrow






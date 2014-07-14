.. _perfanalysis-page:

*********************
Performance Analysis
*********************

Nextflow integrates the support for `Extrae`_ a performance tracing library. Trace files
created by Extrae can be analysed later with `Paraver`_ a visual performance analysis tool.

*Extrae* together with *Paraver* will allows you to analyse the execution performance of your
pipeline.

.. note:: Both *Extrae* and *Paraver* are tools developed by the `Barcelona Supercomputing Center`_.

How to use it
================

This feature currently depends on a custom version of Extrae 2.5.0 that needs to be installed in
the computer where the pipeline is going be executed.

Pre-requisite
---------------

* Install ``libxml2``
* Install ``binutils``

If you are using an Ubuntu Linux distribution these packages can be installed using the following
commands::

    sudo apt-get install libxml2-dev binutils-dev


It may change depending your linux distribution and the available package installer tool.

Installation
--------------

Download the Extrae 2.5.0 at this link http://www.nextflow.io/misc/extrae-2.5.0.tar.gz

Compile and install it by using the following command::

  ./configure \
   --without-mpi \
   --without-unwind \
   --without-dyninst \
   --without-papi \
   --without-java \
   --prefix=<extrae_install_dir>

  make
  make install


When the compilation process is completed define the following variables in your
environment profile file::

  export EXTRAE_HOME=<extrae_install_dir>
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${EXTRAE_HOME}/lib


Profile your pipeline
-----------------------

In order to use the Extrae tracing simply add the option `-with-extrae` to your Nextflow
launch command line, for example::

  nextflow pipeline.nf -with-extra


When the pipeline execution completes, other then the expected result files, it will produce the Extrae
trace file having the following name: ``TRACE.mpits``

Converts this file to the Paraver format by using the following command::

  ${EXTRAE_HOME}/bin/mpi2prv -task-view -f TRACE.mpits -o <your file name>.prv


Analysis with Paraver
-----------------------

If do not have Paraver installed, you need to download and install it in your computer.
You can download it from this page: http://www.bsc.es/performance_tools/downloads

Use the ``File > Load Trace`` command in the Paraver menu to load the trace file
(the file with ``.prv`` suffix).

To perform a basic analysis download the `configuration file available
at this link <http://www.nextflow.io/misc/nextflow_runtime_analysis.cfg>`_ and open it
by using the command ``File -> Load Configuration`` in the Paraver menu.

It defines some views useful for the observe the different tasks duration.

Please refers the `Paraver`_ documentation for information about it.


Advanced configuration
-------------------------

In order to fine control the Extra tracing you can provide a custom Extrae
configuration file by specifying the its location by using the environment
variable ``EXTRAE_CONFIG_FILE``.

Please read the `Extrae`_ documentation for more information about it.



.. _Barcelona Supercomputing Center: http://www.bsc.es
.. _Paraver: http://www.bsc.es/computer-sciences/performance-tools/paraver
.. _Extrae: http://www.bsc.es/computer-sciences/extrae
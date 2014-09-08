.. _perfanalysis-page:

*********************
Performance Analysis
*********************

Nextflow can produce an execution tracing report that provides some information for a basic performance analysis
of a pipeline execution. A more advanced analysis is possible by using the `Extrae` and `Paraver` tools integrated with Nextflow.


Execution report
===================

Nextflow creates a execution tracing file that contains some useful information about each process executed in your pipeline
script, including: submission time, start time, completion time, cpu and memory used.

.. warning:: This is an incubating feature. It may change in future Nextflow releases.


In order to create the execution trace file add the ``-with-trace`` command line option when launching the pipeline execution.
For example::

  nextflow run <pipeline name> -with-trace

It will create a file named ``trace.csv`` in the current directory. The content looks like the above example::

    task-id	hash	native-id	name	status	exit-status	submit	start	complete	wall-time	run-time	cpu	mem
    14	15/93fcd5	952	blast (2)	COMPLETED	0	8/26/2014 18:11:45	8/26/2014 18:11:45	8/26/2014 18:11:53	8s	7s	2.74	422.1 MB
    15	1d/5c2e2b	953	blast (3)	COMPLETED	0	8/26/2014 18:11:45	8/26/2014 18:11:45	8/26/2014 18:11:49	3s	3s	3.15	585.1 MB
    13	b4/e029e3	954	blast (1)	COMPLETED	0	8/26/2014 18:11:45	8/26/2014 18:11:45	8/26/2014 18:12:00	14s	14s	4.93	778.7 MB
    22	3b/c3932b	956	blast (4)	COMPLETED	0	8/26/2014 18:11:45	8/26/2014 18:11:45	8/26/2014 18:11:52	7s	7s	5.72	405.7 MB
    24	9f/290cc2	957	blast (6)	COMPLETED	0	8/26/2014 18:11:45	8/26/2014 18:11:45	8/26/2014 18:11:54	8s	8s	2.44	423 MB
    31	9e/d06eb2	964	exonerate (4)	COMPLETED	0	8/26/2014 18:12:16	8/26/2014 18:12:17	8/26/2014 18:12:22	6s	5s	5.07	126.4 MB
    34	94/b3d5f5	967	exonerate (7)	COMPLETED	0	8/26/2014 18:12:16	8/26/2014 18:12:17	8/26/2014 18:12:26	9s	8s	8.72	566.4 MB
    29	57/a9e2fd	962	exonerate (2)	COMPLETED	0	8/26/2014 18:12:16	8/26/2014 18:12:17	8/26/2014 18:12:32	15s	14s	13.88	92.5 MB
    30	18/3355ae	963	exonerate (3)	COMPLETED	0	8/26/2014 18:12:16	8/26/2014 18:12:17	8/26/2014 18:12:27	11s	9s	9.78	772.6 MB
    80	1d/9655d1	986	similarity (3)	COMPLETED	0	8/26/2014 18:14:31	8/26/2014 18:14:31	8/26/2014 18:14:31	541ms	415ms	0.18	0
    78	ac/3067fe	987	similarity (1)	COMPLETED	0	8/26/2014 18:14:31	8/26/2014 18:14:33	8/26/2014 18:14:34	2s	671ms	0.53	0
    79	c9/b2fc2a	988	similarity (2)	COMPLETED	0	8/26/2014 18:14:31	8/26/2014 18:14:33	8/26/2014 18:14:34	3s	1s	1.08	0
    81	b6/d5cc53	989	similarity (4)	COMPLETED	0	8/26/2014 18:14:31	8/26/2014 18:14:33	8/26/2014 18:14:34	2s	479ms	0.31	0
    124	02/18b20d	1032	matrix (1)	COMPLETED	0	8/26/2014 18:15:04	8/26/2014 18:15:05	8/26/2014 18:15:05	988ms	326ms	0.06	0


.. note:: Currently memory and cpu utilization are included only by using the :ref:`DRMAA executor <drmaa-executor>`.



Paraver integration
=====================


Nextflow integrates the support for `Extrae`_ a library for performance tracing. Trace files
created by Extrae can be analysed with `Paraver`_, a visual performance analysis tool.

*Extrae* together with *Paraver* will allows you to analyse the execution performance of your
pipeline.

.. note:: Both *Extrae* and *Paraver* are tools developed by the `Barcelona Supercomputing Center`_.


How to use it
---------------

This feature currently depends on a custom version of Extrae 2.5.0 that needs to be installed in
the computer where the pipeline is executed.

.. warning:: This is an incubating feature. It may change in future Nextflow releases.

Pre-requisite
---------------

* Install ``libxml2``
* Install ``binutils``

If you are using an Ubuntu Linux distribution these packages can be installed using the following
commands::

    sudo apt-get install libxml2-dev binutils-dev


It may change depending your Linux distribution and the available package installer tool.

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

In order to use the Extrae tracing simply add the option ``-with-extrae`` to your Nextflow
launch command line, for example::

  nextflow run <your pipeline> -with-extrae


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

It defines some views useful to observe the different tasks duration.

Please refers the `Paraver`_ documentation for information about it.


Advanced configuration
-------------------------

In order to fine control the Extra tracing you can provide a custom Extrae
configuration file by specifying its location by using the environment
variable ``EXTRAE_CONFIG_FILE``.

Read the `Extrae`_ documentation for more information about it.



.. _Barcelona Supercomputing Center: http://www.bsc.es
.. _Paraver: http://www.bsc.es/computer-sciences/performance-tools/paraver
.. _Extrae: http://www.bsc.es/computer-sciences/extrae
.. _ignite-page:

*************
Apache Ignite
*************


Nextflow can be be deployed in a *cluster* mode by using `Apache Ignite <https://ignite.apache.org/>`_, an in-memory data-grid
and clustering platform.

Apache Ignite is packaged with Nextflow itself, so you won't need to install it separately or configure other third party
software.

.. warning:: This is an incubating feature. It may change in future Nextflow releases.

.. _ignite-daemon:

Cluster daemon
---------------------

In order to setup a cluster you will need to run a cluster daemon on each node that made up your cluster.
If you want to support the :ref:`Docker integration <docker-page>` feature provided by Nextflow, the Docker engine has
to be installed and must run in each node.

In its simplest form just launch the Nextflow daemon in each cluster node as shown below::

    nextflow node -bg

The command line options ``-bg`` launches the node daemon in the background. The output is stored in the log file ``.node-nextflow.log``. The daemon
process ``PID`` is saved in the file ``.nextflow.pid`` in the same folder.


Multicast discovery
====================

By default, the Ignite daemon tries to discover the members in the cluster by using the network *multicast* discovery.
Note that NOT all networks support this feature (Amazon EC2 does not).

.. tip::  To check if multicast is available, `iperf <http://sourceforge.net/projects/iperf/>`_ is a useful tool which is available for Windows/\*NIX/OSX.
  To test multicast open a terminal one 2 machines within the network and run the following in the first terminal
  ``iperf -s -u -B 228.1.2.4 -i 1`` and ``iperf -c 228.1.2.4 -u -T 32 -t 3 -i 1`` in the other terminal.
  If data is being transferred then multicast is working.


Ignite uses the multicast group ``228.1.2.4`` and port ``47400`` by default. You can change these values, by using the
``cluster.join`` command line option, as shown below::

    nextflow node -bg -cluster.join multicast:224.2.2.3:55701



In the case that multicast discovery is not available in your network, you can try one of the following alternative methods:

Shared file system
====================

Simply provide a path shared across the cluster by a network file system, as shown below::

    nextflow node -bg -cluster.join path:/net/shared/cluster


The cluster members will use that path to discover each other.


IP address
============

Provide the list of pre-configured IP addresses on the daemon launch command line, for example::

    nextflow node -bg -cluster.join ip:10.0.2.1,10.0.2.2,10.0.2.4

AWS S3 bucket
===============

Creates an Amazon AWS S3 bucket that will hold the cluster members IP addresses. For example::

   nextflow node -bg -cluster.join s3:cluster_bucket




Advanced options
=====================

The following cluster node configuration options can be used.

=========================== ================
Name                        Description
=========================== ================
join                        IP address(es) of one more more cluster nodes to which the daemon make part. By using the value ``multicast``
group                       Cluster name of which this node makes part. It allows to create separate clusters. Default: ``nextflow``
slots                       Number of slots this damon node provides i.e. of process that can execute in parallel. By it is equal to the number of CPU cores.
interface                   Network interfaces that Ignite has to use. It can be the interface IP address or name.
config.file                 The file path of the Ignite configuration file (optional)
config.url                  the url of the Ignite configuration file (optional)
tcp.localAddress            Sets local host IP address.
tcp.localPort               Sets local port to listen to.
tcp.localPortRange          Range for local ports.
tcp.heartbeatFrequency      Gets delay between heartbeat messages sent by coordinator.
tcp.maxMissedHeartbeats     Gets max heartbeats count node can miss without initiating status check.
tcp.reconnectCount          Number of times node tries to (re)establish connection to another node.
tcp.networkTimeout          Gets network timeout.
tcp.socketTimeout           Sets socket operations timeout. This timeout is used to limit connection time and write-to-socket time. Note that when running Ignite on Amazon EC2, socket timeout must be set to a value significantly greater than the default (e.g. to 30000).
tcp.ackTimeout              Sets timeout for receiving acknowledgement for sent message.
tcp.maxAckTimeout           Sets maximum timeout for receiving acknowledgement for sent message.
tcp.joinTimeout             Sets join timeout.
=========================== ================

These options can be specified as command line parameters by pre-pending them the prefix ``-cluster.``, as shown below::

    nextflow node -bg -cluster.slots 4 -cluster.interface eth0

The same options can be entered in the nextflow configuration file named ``nextflow.config``, as shown below::


  cluster {
    join = 'ip:192.168.1.104'
    interface = 'eth0'
  }

Finally daemon options can be provided also as environment variables having the name in upper-case and by pre-pending
them with the prefix ``NXF_CLUSTER_``, for example::

    export NXF_CLUSTER_JOIN='ip:192.168.1.104'
    export NXF_CLUSTER_INTERFACE='eth0'


Pipeline execution
-----------------------

The pipeline should be launched in a `head` node i.e. a cluster node where the Nextflow node daemon is **not** running.
In order to execute your pipeline in the Ignite cluster you will need to specify to use the Ignite executor,
as shown below::

   nextflow run <your pipeline> -process.executor Ignite


If your network do no support multicast discovery, you will need to specify the `joining` strategy as you did for the
cluster daemons. For example, using a shared path::

    nextflow run <your pipeline> -process.executor ignite -cluster.join path:/net/shared/cluster



Execution with MPI
-------------------

Nextflow is able to deploy and self-configure a Ignite cluster on-demand, taking advantage the Open `MPI <https://en.wikipedia.org/wiki/Message_Passing_Interface>`_
standard that is commonly available in grid and supercomputer facilities.

In this scenario a Nextflow workflow needs to be executed as a MPI job. Under the hood Nextflow will launch a `driver`
process in the first of the nodes allocated by your job request and a Ignite daemon in the remaining nodes.

In practice you will need a wrapper script to submit a MPI job request to your batch scheduler/resource manager.
The batch scheduler must reserve the computing nodes in an exclusive manner to avoid to have multiple Ignite daemons
running on the same node. Nextflow must be launched using the ``mpirun`` utility, as if it were a MPI application,
specifying the ``--pernode`` option.

The following example shows a wrapper script for the `Platform LSF <https://en.wikipedia.org/wiki/Platform_LSF/>`_ resource manager::

    #!/bin/bash
    #BSUB -oo output_%J.out
    #BSUB -eo output_%J.err
    #BSUB -J <job name>
    #BSUB -q <queue name>
    #BSUB -W 02:00
    #BSUB -x
    #BSUB -n 80
    #BSUB -R "span[ptile=16]"
    export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)
    mpirun --pernode nextflow run <your-project-name> -with-mpi [pipeline parameters]

It requests 5 nodes (80 processes, with 16 cpus per node). The ``-x`` directive allocate the node in exclusive manner.
Nextflow needs to be executed using the ``-with-mpi`` command line option. It will automatically use ``ignite`` as executor.

The variable ``NXF_CLUSTER_SEED`` must contain an integer value (in the range 0-16777216) that will unequivocally identifies
your cluster instance. In the example it is randomly generated by using the ``shuf`` Linux command.

The following example shows a wrapper script for the Sun/`Univa grid engine <https://en.wikipedia.org/wiki/Univa_Grid_Engine>`_::

    #!/bin/bash
    #$ -cwd
    #$ -j y
    #$ -o <output file name>
    #$ -l virtual_free=120G
    #$ -q <queue name>
    #$ -N <job name>
    #$ -pe ompi 5
    export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)
    mpirun --pernode nextflow run <your-project-name> -with-mpi [pipeline parameters]

As in the previous example it allocates 5 processing nodes. SGE/UGE does not have an option to reserve a node in exclusive
manner. A common workaround is to request the maximum amount of memory or cpus available in the nodes of your cluster.
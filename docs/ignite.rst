.. _ignite-page:

*************
Apache Ignite
*************


Nextflow can be be deployed in a *cluster* mode by using `Apache Ignite <https://ignite.apache.org/>`_, an in-memory data-grid
and clustering platform.

Apache Ignite is packaged with Nextflow itself, so you won't need to install it separately or configure other third party
software.

.. _ignite-daemon:

Cluster daemon
--------------

In order to setup a cluster you will need to run a cluster daemon on each node within your cluster.
If you want to support the :ref:`Docker integration <docker-page>` feature provided by Nextflow, the Docker engine has
to be installed and must run in each node.

In its simplest form just launch the Nextflow daemon in each cluster node as shown below::

    nextflow node -bg

The command line option ``-bg`` launches the node daemon in the background. The output is stored in the log file ``.node-nextflow.log``. The daemon
process ``PID`` is saved in the file ``.nextflow.pid`` in the same folder.


Multicast discovery
===================

By default, the Ignite daemon tries to automatically discover all members in the cluster by using the network *multicast* discovery.
Note that NOT all networks support this feature (for example Amazon AWS does not).

.. tip::  To check if multicast is available in your network, use the `iperf <http://sourceforge.net/projects/iperf/>`_ tool.
  To test multicast, open a terminal on two machines within the network and run the following command in the first one
  ``iperf -s -u -B 228.1.2.4 -i 1`` and ``iperf -c 228.1.2.4 -u -T 32 -t 3 -i 1`` on the second.
  If data is being transferred then multicast is working.


Ignite uses the multicast group ``228.1.2.4`` and port ``47400`` by default. You can change these values, by using the
``cluster.join`` command line option, as shown below::

    nextflow node -cluster.join multicast:224.2.2.3:55701



In case multicast discovery is not available in your network, you can try one of the following alternative methods:

Shared file system
==================

Simply provide a path shared across the cluster by a network file system, as shown below::

    nextflow node -bg -cluster.join path:/net/shared/cluster


The cluster members will use that path to discover each other.


IP addresses
============

Provide a list of pre-configured IP addresses on the daemon launch command line, for example::

    nextflow node -cluster.join ip:10.0.2.1,10.0.2.2,10.0.2.4


Advanced options
=====================

The following cluster node configuration options can be used.

============================= ================
Name                          Description
============================= ================
join                          IP address(es) of one or more cluster nodes to which the daemon will join.
group                         The name of the cluster to which this node join. It allows the creation of separate cluster instances. Default: ``nextflow``
maxCpus                       Max number of CPUs that can be used by the daemon to run user tasks. By default it is equal to the number of CPU cores.
maxDisk                       Max amount of disk storage that can be used by user tasks eg. ``500 GB``.
maxMemory                     Max amount of memory that can be used by user tasks eg. ``64 GB``. Default total available memory.
interface                     Network interfaces that Ignite will use. It can be the interface IP address or name.
failureDetectionTimeout       Failure detection timeout is used to determine how long the communication or discovery SPIs should wait before considering a remote connection failed.
clientFailureDetectionTimeout Failure detection timeout is used to determine how long the communication or discovery SPIs should wait before considering a remote connection failed.
tcp.localAddress              Defines the local host IP address.
tcp.localPort                 Defines the local port to listen to.
tcp.localPortRange            Range for local ports.
tcp.reconnectCount            Number of times the node tries to (re)establish connection to another node.
tcp.networkTimeout            Defines the network timeout.
tcp.socketTimeout             Defines the socket operations timeout. This timeout is used to limit connection time and write-to-socket time. Note that when running Ignite on Amazon EC2, socket timeout must be set to a value significantly greater than the default (e.g. to 30000).
tcp.ackTimeout                Defines the timeout for receiving acknowledgement for sent message.
tcp.maxAckTimeout             Defines the maximum timeout for receiving acknowledgement for sent message.
tcp.joinTimeout               Defines the join timeout.
============================= ================



These options can be specified as command line parameters by adding the prefix ``-cluster.`` to them, as shown below::

    nextflow node -bg -cluster.maxCpus 4 -cluster.interface eth0

The same options can be entered into the Nextflow :ref:`configuration file<config-page>`, as shown below::

    cluster {
        join = 'ip:192.168.1.104'
        interface = 'eth0'
    }

Finally daemon options can be provided also as environment variables having the name in upper-case and by adding
the prefix ``NXF_CLUSTER_`` to them, for example::

    export NXF_CLUSTER_JOIN='ip:192.168.1.104'
    export NXF_CLUSTER_INTERFACE='eth0'


Pipeline execution
------------------

The pipeline execution needs to be launched in a `head` node i.e. a cluster node where the Nextflow node daemon
is **not** running. In order to execute your pipeline in the Ignite cluster you will need to use the Ignite executor,
as shown below::

   nextflow run <your pipeline> -process.executor ignite


If your network does no support multicast discovery, you will need to specify the `joining` strategy as you did for the
cluster daemons. For example, using a shared path::

    nextflow run <your pipeline> -process.executor ignite -cluster.join path:/net/shared/path



Execution with MPI
------------------

Nextflow is able to deploy and self-configure an Ignite cluster on demand, taking advantage of the Open `MPI <https://en.wikipedia.org/wiki/Message_Passing_Interface>`_
standard that is commonly available in grid and supercomputer facilities.

In this scenario a Nextflow workflow needs to be executed as an MPI job. Under the hood Nextflow will launch a `driver`
process in the first of the nodes, allocated by your job request, and an Ignite daemon in the remaining nodes.

In practice you will need a launcher script to submit an MPI job request to your batch scheduler/resource manager.
The batch scheduler must reserve the computing nodes in an exclusive manner to avoid having multiple Ignite daemons
running on the same node. Nextflow must be launched using the ``mpirun`` utility, as if it were an MPI application,
specifying the ``--pernode`` option.

Platform LSF launcher
=====================

The following example shows a launcher script for the `Platform LSF <https://en.wikipedia.org/wiki/Platform_LSF/>`_ resource manager::

    #!/bin/bash
    #BSUB -oo output_%J.out
    #BSUB -eo output_%J.err
    #BSUB -J <job name>
    #BSUB -q <queue name>
    #BSUB -W 02:00
    #BSUB -x
    #BSUB -n 80
    #BSUB -M 10240
    #BSUB -R "span[ptile=16]"
    export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)
    mpirun --pernode nextflow run <your-project-name> -with-mpi [pipeline parameters]

It requests 5 nodes (80 processes, with 16 cpus per node). The ``-x`` directive allocates the node in an exclusive manner.
Nextflow needs to be executed using the ``-with-mpi`` command line option. It will automatically use ``ignite`` as the executor.

The variable ``NXF_CLUSTER_SEED`` must contain an integer value (in the range 0-16777216) that will unequivocally identify
your cluster instance. In the above example it is randomly generated by using the `shuf <http://linux.die.net/man/1/shuf>`_ Linux tool.

Univa Grid Engine launcher
==========================

The following example shows a launcher script for the `Univa Grid Engine <https://en.wikipedia.org/wiki/Univa_Grid_Engine>`_ (aka SGE)::

    #!/bin/bash
    #$ -cwd
    #$ -j y
    #$ -o <output file name>
    #$ -l virtual_free=10G
    #$ -q <queue name>
    #$ -N <job name>
    #$ -pe ompi 5
    export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)
    mpirun --pernode nextflow run <your-project-name> -with-mpi [pipeline parameters]

As in the previous script it allocates 5 processing nodes. UGE/SGE does not have an option to reserve a node in an exclusive
manner. A common workaround is to request the maximum amount of memory or cpus available in the nodes of your cluster.


Linux SLURM launcher
====================

When using Linux SLURM you will need to use ``srun`` instead ``mpirun`` in your launcher script. For example::

    #!/bin/bash
    #SBATCH --job-name=<job name>
    #SBATCH --output=<log file %j>
    #SBATCH --ntasks=5
    #SBATCH --cpus-per-task=16
    #SBATCH --tasks-per-node=1
    export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)
    srun nextflow run hello.nf -with-mpi

As before, this allocates 5 processing nodes (``--ntasks=5``) and each node will be able to use up to 16 cpus
(``--cpus-per-task=16``). When using SLURM it's not necessary to allocate computing nodes in an exclusive manner.
It's even possible to launch more than one Nextflow daemon instance per node, though not suggested.

To submit the pipeline execution create a file like the above, then use the following command::

    sbatch <launcher script name>


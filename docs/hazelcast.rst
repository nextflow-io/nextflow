.. _hazelcast-page:

*******************
Hazelcast cluster
*******************


Nextflow can be be deployed in a *cluster* mode by setting up an `Hazelcast <http://www.hazelcast.com>`_ cluster,
a light-weight and efficient data-grid technology.

Hazelcast is packaged with Nextflow itself so you won't need to install it separately or configure other third party
software.

.. warning:: This an experimental feature.

.. _hazelcast-daemon:

Cluster daemon
---------------------

In order to setup a cluster you will need to run a cluster daemon on each mode that made up your cluster. It only
requires Java 7 to be installed.

In the case your network supports *multicast* discovery, simply launch the Nextflow daemon in each cluster node
as shown below::

    nextflow -daemon.name hazelcast -bg

The process will run in the background. The daemon output is stored in the log file ``.nxf-daemon.log``. The daemon
process ``PID`` is saved in the file ``.nextflow.pid`` in the same folder.

If you don't see members joining, then it is likely because multicast is not available.

.. tip::  To check if multicast is available, `iperf <http://sourceforge.net/projects/iperf/>`_ is a useful tool which is available for Windows/\*NIX/OSX.
  To test multicast open a terminal one 2 machines within the network and run the following in the first terminal
  ``iperf -s -u -B 224.2.2.3 -i 1`` and ``iperf -c 224.2.2.3 -u -T 32 -t 3 -i 1`` in the other terminal.
  If data is being transferred then multicast is working.


In the case the multicast discovery is not available in your network, you will need to provide at least one (or more)
IP address of well known members to connect to. For example::

    nextflow -daemon.name hazelcast -daemon.join 192.168.1.104 -bg

Configuration options
^^^^^^^^^^^^^^^^^^^^^^^

The following daemon configuration options can be used.

=========================== ================
Name                        Description
=========================== ================
port                        TCP port to which the daemon will connect. Default: ``5701``
join                        IP address(es) of one more more cluster nodes to which the daemon make part. By using the value ``multicast``
group                       Cluster name of which this node makes part. It allows to create separate clusters. Default: ``nextflow``
slots                       TNumber of slots this damon node provides i.e. of process that can execute in parallel. By it is equal to the number of CPU cores.
interface                   Network interfaces that Hazelcast has to use. It can be the interface IP address or name.
connectionMonitorInterval   Minimum interval to consider a connection error as critical in milliseconds (default: 300).
connectionMonitorMaxFaults  Maximum IO error count before disconnecting from a node (default: 5).

=========================== ================

These options can be specified as command line parameters by pre-pending them the prefix ``-daemon.``, as shown below::

    nextflow -daemon.name hazelcast -daemon.port 5705 -daemon.interface eth0 -bg

The same options can be entered in the nextflow configuration file named ``nextflow.config``, as shown below::


  daemon {
    name = 'hazelcast'
    port = '5705'
    join = '192.168.1.104'
    interface = 'eth0'
  }

Finally daemon options can be provided also as environment variables having the name in upper-case and by pre-pending
them the prefix ``NXF_DAEMON_``, for example::

    export NXF_DAEMON_PORT=5705
    export NXF_DAEMON_JOIN=192.168.1.104
    export NXF_DAEMON_INTERFACE=eth0


Pipeline execution
-----------------------

The pipeline should be launched in a `head` node i.e. a cluster node where the Nextflow daemon is NOT running.
In order to join and execute your pipeline in the Hazelcast cluster, you need to specify the Hazelcast executor
and an IP address of the cluster to join, as shown below::

   nextflow -executor.name hazelcast -executor.join 192.168.1.104 <your pipeline script>


In place of the command line parameters, the following environment variables can be used in a
semantically equivalent manner::

    export NXF_EXECUTOR_NAME=hazelcast
    export NXF_EXECUTOR_JOIN=192.168.1.104



Read more about the executor configuration options in :ref:`Hazelcast executor <hazelcast-executor>` page.



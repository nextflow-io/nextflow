.. _gridgain-page:

*******************
GridGain cluster
*******************


Nextflow can be be deployed in a *cluster* mode by using `GridGain <http://www.gridgain.com>`_, an In-Memory data-grid
and clustering platform.

GridGain is packaged with Nextflow itself, so you won't need to install it separately or configure other third party
software.

.. warning:: GridGain cluster integration is an incubating feature.

.. _gridgain-daemon:

Cluster daemon
---------------------

In order to setup a cluster you will need to run a cluster daemon on each mode that made up your cluster. It only
requires Java 7 to be installed.

In its simplest form just launch the Nextflow daemon in each cluster node as shown below::

    nextflow -d -bg

The process will run in the background. The daemon output is stored in the log file ``.nxf-daemon.log``. The daemon
process ``PID`` is saved in the file ``.nextflow.pid`` in the same folder.

If you don't see members joining, then it is likely because multicast is not available.

Multicast discovery
====================

By default, it tries to discover the cluster members by using the network *multicast* discovery. Note that NOT all
networks support this feature.

.. tip::  To check if multicast is available, `iperf <http://sourceforge.net/projects/iperf/>`_ is a useful tool which is available for Windows/\*NIX/OSX.
  To test multicast open a terminal one 2 machines within the network and run the following in the first terminal
  ``iperf -s -u -B 228.1.2.4 -i 1`` and ``iperf -c 228.1.2.4 -u -T 32 -t 3 -i 1`` in the other terminal.
  If data is being transferred then multicast is working.


GridGain uses the multicast group ``228.1.2.4`` and port ``47400`` by default. You can change these values, by using the
``daemon.join`` command line option, as shown below::

    nextflow -bg -daemon.join multicast:224.2.2.3:55701



In the case that multicast discovery is not available in your network, you can try one of the following alternative methods:

Shared file system
====================

Simply provide a path shared across the cluster by a network file system, as shown below::

    nextflow -bg -daemon.join path:/net/shared/cluster


The cluster members will use that path to discover each other.


IP address
============

Provide the list of pre-configured IP addresses on the daemon launch command line, for example::

    nextflow -bg -daemon.join ip:10.0.2.1,10.0.2.2,10.0.2.4

AWS S3 bucket
===============

Creates an Amazon AWS S3 bucket that will hold the cluster members IP addresses. For example::

   nextflow -bg -daemon.join s3:cluster_bucket




Advanced options
=====================

The following daemon configuration options can be used.

=========================== ================
Name                        Description
=========================== ================
join                        IP address(es) of one more more cluster nodes to which the daemon make part. By using the value ``multicast``
group                       Cluster name of which this node makes part. It allows to create separate clusters. Default: ``nextflow``
slots                       Number of slots this damon node provides i.e. of process that can execute in parallel. By it is equal to the number of CPU cores.
interface                   Network interfaces that GridGain has to use. It can be the interface IP address or name.
config.file                 The file path of the GridGain configuration file (optional)
config.url                  the url of the GridGain configuration file (optional)
tcp.localAddress            Sets local host IP address.
tcp.localPort               Sets local port to listen to.
tcp.localPortRange          Range for local ports.
tcp.heartbeatFrequency      Gets delay between heartbeat messages sent by coordinator.
tcp.maxMissedHeartbeats     Gets max heartbeats count node can miss without initiating status check.
tcp.reconnectCount          Number of times node tries to (re)establish connection to another node.
tcp.networkTimeout          Gets network timeout.
tcp.socketTimeout           Sets socket operations timeout. This timeout is used to limit connection time and write-to-socket time. Note that when running GridGain on Amazon EC2, socket timeout must be set to a value significantly greater than the default (e.g. to 30000).
tcp.ackTimeout              Sets timeout for receiving acknowledgement for sent message.
tcp.maxAckTimeout           Sets maximum timeout for receiving acknowledgement for sent message.
tcp.joinTimeout             Sets join timeout.
=========================== ================

These options can be specified as command line parameters by pre-pending them the prefix ``-daemon.``, as shown below::

    nextflow -bg -daemon.slots 4 -daemon.interface eth0

The same options can be entered in the nextflow configuration file named ``nextflow.config``, as shown below::


  daemon {
    join = 'ip:192.168.1.104'
    interface = 'eth0'
  }

Finally daemon options can be provided also as environment variables having the name in upper-case and by pre-pending
them the prefix ``NXF_DAEMON_``, for example::

    export NXF_DAEMON_JOIN='ip:192.168.1.104'
    export NXF_DAEMON_INTERFACE='eth0'


Pipeline execution
-----------------------

The pipeline should be launched in a `head` node i.e. a cluster node where the Nextflow daemon is not running.
In order to execute your pipeline in the GridGain cluster you will need to specify the GridGain executor,
as shown below::

   nextflow -executor.name gridgain <your pipeline script>


If your network do no support multicast discovery, you will need to specify the `joining` strategy as you did for the
cluster daemons. For example, using a shared path::

    nextflow -executor.name gridgain -executor.join path:/net/shared/cluster <your pipeline script>


In place of the command line parameters, the following environment variables can be used in a
semantically equivalent manner::

    export NXF_EXECUTOR_NAME='gridgain'
    export NXF_EXECUTOR_JOIN='path:/net/shared/cluster'






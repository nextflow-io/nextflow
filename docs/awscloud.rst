.. _amazonscloud-page:

************
Amazon Cloud
************

Nextflow provides out of the box support for the Amazon AWS cloud allowing you to setup a computing cluster,
deploy it and run your pipeline in the AWS infrastructure in a few commands.


Configuration
=============

Cloud configuration attributes are provided in the ``nextflow.config`` file as shown in the example below::

    cloud {
        imageId = 'ami-43f49030'
        instanceType = 'm4.xlarge'
        subnetId = 'subnet-05222a43'
    }

The above attributes define the virtual machine ID and type to be used and the VPC subnet ID to be applied
in you cluster. Replace these values with the ones of your choice.

Nextflow only requires a Linux image that provides support for `Cloud-init <http://cloudinit.readthedocs.io/>`_
bootstrapping mechanism and includes a Java runtime (version 7 or 8) and a Docker engine (version 1.11 or higher).

For your convenience the following pre-configured Amazon Linux AMI is available in the *EU Ireland* region:
``ami-43f49030``.


AWS credentials
---------------

Nextflow will use the AWS credentials defined in your environment, using the standard AWS variables shown below:

    * ``AWS_ACCESS_KEY_ID``
    * ``AWS_SECRET_ACCESS_KEY``
    * ``AWS_DEFAULT_REGION``


Alternatively AWS credentials can be specified in the Nextflow configuration file.
See :ref:`AWS configuration<config-aws>` for more details.


User & SSH key
--------------

By default Nextflow creates in each EC2 instance a user with the same name as the one in your local computer and install
the SSH public key available at the path ``$HOME/.ssh/id_rsa.pub``. A different user/key can be specified as shown below::

    cloud {
        userName = 'the-user-name'
        keyFile = '/path/to/ssh/key.pub'
    }

If you want to use a *key-pair* defined in your AWS account and the default user configured in the AMI, specify the
attribute ``keyName`` in place of ``keyFile`` and the name of the existing user specifying the ``userName`` attribute.


Storage
-------

The following storage types can be defined in the cloud instances: *boot*, *instance* and *shared*.

Boot storage
------------

You can set the size of the `root device volume <http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/RootDeviceStorage.html>`_
by specifying the attribute ``bootStorageSize``. For example::

    cloud {
        imageId = 'ami-xxx'
        bootStorageSize = '10 GB'
    }


Instance storage
----------------

Amazon instances can provide one or more `ephemeral volume storage <http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/InstanceStorage.html>`_,
depending the instance type chosen. This storage can be made available by using the ``instanceStorageMount``
and ``instanceStorageDevice`` configuration attributes, as shown below::

    cloud {
            imageId = 'xxx'
            instanceStorageMount = '/mnt/scratch'
            instanceStorageDevice = '/dev/xvdc'
    }


The mount path can be any of your choice. The storage device name should follow a pattern that may change depending
the AMI and machine virtualization type. See the `Amazon documentation for details <http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/device_naming.html>`_.

.. note:: Currently Nextflow can setup only one instance storage volume.


Shared file system
------------------

A shared file system can easily made available in your cloud cluster by using the `Amazon EFS <https://aws.amazon.com/efs/>`_
storage. You will only need to specify in the cloud configuration the EFS file system ID and optionally the
mount path. For example::

    cloud {
        imageId = 'ami-xxx'
        sharedStorageId = 'fs-1803efd1'
        sharedStorageMount = '/mnt/shared'
    }

.. note:: When the attribute ``sharedStorageMount`` is omitted the path ``/mnf/efs`` is used by default.


Cluster deployment
==================

Once defined the configuration settings in the ``nextflow.config`` file you can create the cloud cluster
by using the following command::

    nextflow cloud create my-cluster -c <num-of-nodes>

The string ``my-cluster`` identifies the cluster instance. Replace it with a name of your choice.

Finally replace ``num-of-nodes`` with the actual number of instances that will made-up the cluster.
One node is created as *master*, the remaining as *workers*. If the option ``-c`` is omitted only the *master* node
is created.

.. warning:: You will be charged accordingly the type and the number of instances chosen.


Pipeline execution
==================

Once the cluster initialization is complete, connect to the *master* node by using the SSH command that will be showed by
Nextflow.

You can run your Nextflow pipeline as usual, the environment is automatically configured to use the :ref:`Ignite<ignite-page>`
executor. If the Amazon EFS storage is specified in the cloud configuration the Nextflow work directory will
automatically be set in a shared folder in that file system.

The suggested approach is to run your pipeline downloading it from a public repository such
GitHub and to pack the binaries dependencies in a Docker container as described in the
:ref:`Pipeline sharing <sharing-page>` section.

Cluster shutdown
================

When completed shutdown the cluster instances by using the following command::

    nextflow cloud shutdown my-cluster


Cluster auto-scaling
====================

Nextflow integration for AWS cloud provides a native support auto-scaling that allows the computing cluster
to scale-out or scale-down i.e. add or remove computing nodes dynamically at runtime.

This is a critical feature, especially for pipelines crunching not homogeneous dataset, because it allows the
cluster to adapt dynamically to the actual workload computing resources need as they change over the time.

Cluster auto-scaling is enabled by adding the ``autoscale`` option group in the configuration file as shown below::

    cloud {
        imageId = 'xxx'
        autoscale {
            enabled = true
            maxInstances = 10
        }
    }


The above example enables automatic cluster scale-out i.e. new instances are automatically launched and added to the
cluster when tasks remain too long in wait status because there aren't enough computing resources available. The
``maxInstances`` attribute defines the upper limit to which the cluster can grow.

By default unused instances are not removed when are not utilised. If you want to enable automatic cluster scale-down
specify the ``terminateWhenIdle`` attribute in the ``autoscale`` configuration group.

It is also possible to define a different AMI image ID, type and spot price for instances launched by the Nextflow autoscaler.
For example::

    cloud {
        imageId = 'ami-xxx'
        instanceType = 'm4.large'

        autoscale {
            enabled = true
            spotPrice = 0.15
            minInstances = 5
            maxInstances = 10
            imageId = 'ami-yyy'
            instanceType = 'm4.4xlarge'
            terminateWhenIdle = true
        }
    }

By doing that it's is possible to create a cluster with a single node i.e. the master node. Then the autoscaler will
automatically add the missing instances, up to the number defined by the ``minInstances`` attributes. These will have a
different image and type from the master node and will be launched a *spot instances* because the ``spotPrice``
attribute has been specified.


Spot prices
===========

Nextflow includes an handy command to list the current price of EC2 spot instances. Simply type the following
command in your shell terminal::

    nextflow cloud spot-prices

It will print the current spot price for all available instances type, similar to the example below::

    TYPE        PRICE  PRICE/CPU ZONE       DESCRIPTION             CPUS   MEMORY DISK
    t1.micro    0.0044    0.0044 eu-west-1c Linux/UNIX                 1 627.7 MB -
    m4.4xlarge  0.1153    0.0072 eu-west-1a Linux/UNIX (Amazon VPC)   16    64 GB -
    m4.10xlarge 0.2952    0.0074 eu-west-1b Linux/UNIX (Amazon VPC)   40   160 GB -
    m4.large    0.0155    0.0077 eu-west-1b Linux/UNIX (Amazon VPC)    2     8 GB -
    m4.2xlarge  0.0612    0.0077 eu-west-1a Linux/UNIX (Amazon VPC)    8    32 GB -
    m4.xlarge   0.0312    0.0078 eu-west-1a Linux/UNIX (Amazon VPC)    4    16 GB -
    c4.8xlarge  0.3406    0.0095 eu-west-1c Linux/UNIX (Amazon VPC)   36    60 GB -
    m1.xlarge   0.0402    0.0100 eu-west-1b Linux/UNIX                 4    15 GB 4 x 420 GB
    c4.4xlarge  0.1652    0.0103 eu-west-1b Linux/UNIX (Amazon VPC)   16    30 GB -
    c1.xlarge   0.0825    0.0103 eu-west-1a Linux/UNIX                 8     7 GB 4 x 420 GB
    m1.medium   0.0104    0.0104 eu-west-1b Linux/UNIX (Amazon VPC)    1   3.8 GB 1 x 410 GB
    c3.8xlarge  0.3370    0.0105 eu-west-1a Linux/UNIX                32    60 GB 2 x 320 GB
    c3.2xlarge  0.0860    0.0108 eu-west-1c Linux/UNIX                 8    15 GB 2 x 80 GB
    c3.4xlarge  0.1751    0.0109 eu-west-1c Linux/UNIX (Amazon VPC)   16    30 GB 2 x 160 GB
    m3.2xlarge  0.0869    0.0109 eu-west-1c Linux/UNIX (Amazon VPC)    8    30 GB 2 x 80 GB
    r3.large    0.0218    0.0109 eu-west-1c Linux/UNIX                 2  15.2 GB 1 x 32 GB
    :


It's even possible to refine the showed data by specifying a filtering and ordering criteria. For example::

    nextflow cloud spot-prices -sort pricecpu -filter "cpus==4"


It will only print instance types having 4 cpus and sorting them by the best price per cpu.





Advanced configuration
======================

Read :ref:`Cloud configuration<config-cloud>` section to learn more about advanced cloud configuration options.








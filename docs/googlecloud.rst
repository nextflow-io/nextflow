.. _googlecloud-page:

************
Google Cloud
************

Nextflow provides out of the box support for the Google Cloud Platform allowing you to setup a computing cluster, deploy it, and run your pipeline in the Google Cloud infrastructure with a few commands.


Configuration
=============

Cloud configuration attributes are provided in the ``nextflow.config`` file as shown in the example below::

    cloud {
        imageId = 'centos-cloud/global/images/centos-7-v20181011'
        instanceType = 'n1-highcpu-8'
    }

    gce {
        project = 'your-project-id'
        zone = 'us-central1-f'
    }

The above attributes define the image ID and instance type to be used along with the project and zone to be used. Replace these values with the ones of your choice.

Nextflow only requires a Linux image that provides support for `Cloud-init <http://cloudinit.readthedocs.io/>`_ bootstrapping mechanism and includes a Java runtime (version 8) and a Docker engine (version 1.11 or higher).


Credentials
-----------

Download the Google Cloud CLI and initialise using the following command::

    gcloud init

`Enable the Cloud Genomics, Compute Engine, and Cloud Storage APIs <https://console.cloud.google.com/flows/enableapi?apiid=genomics,compute,storage_api>`_

Create and download credentials for google cloud in the console using the following procedure:

* Go to APIs & Services → Credentials →  Create Credentials
* Select Service account key.
* Choose the following Roles: Genomics > Genomics Admin, Storage > Storage Object Admin, & Service Account > Service Account User
* Download the json file and save as "creds.json".

Export as env variable::

    export GOOGLE_APPLICATION_CREDENTIALS=$PWD/creds.json

User & SSH key management
-------------------------

By default Nextflow creates in each GCE instance a user with the same name as the one in your local computer and install the SSH public key available at the path $HOME/.ssh/id_rsa.pub. A different user/key can be specified as shown below::

    cloud {
        userName = 'the-user-name'
        keyFile = '/path/to/ssh/key.pub'
        driver = 'gce'
    }

.. note:: The configuration of your instance on GCP requires a driver parameter as shown in the example above.


Cluster deployment
------------------

Once you have defined the configuration settings in the ``nextflow.config`` file you can create the cloud cluster by using the following command::

    ./nextflow cloud create my-cluster -c <num-of-nodes>

The string ``my-cluster`` identifies the cluster instance. Replace it with a name of your choice.

Finally replace ``<num-of-nodes>`` with the actual number of instances that will make-up the cluster. One node is created as master, the remaining as workers. If the option ``-c`` is omitted only the **master** node
is created.

.. warning:: You will be charged accordingly the type and the number of instances chosen.

The console will display the configuration that you have defined and ask you to confirm creation of the cluster with the configuration displayed. It will take some time for the cluster to deploy. You should be able to track the status of the deployed in the administration console for the Google Cloud Platform under VM instances in the Compute Engine section.


Pipeline execution
==================

Once the cluster initialisation is complete, Nextflow will display the SSH command to connect to the master node. Use that command to connect to the cluster.

.. note:: On MacOS, use the following command to avoid being asked for a pass-phrase even
  you haven't defined one::

    ssh-add -K [private key file]

The suggested approach is to run your pipeline downloading it from a public repository such as GitHub and to pack the binaries dependencies in a Docker container as described in the :ref:`Pipeline sharing <sharing-page>` section.

Cluster shutdown
================

When completed shutdown the cluster instances by using the following command::

    nextflow cloud shutdown my-cluster

Preemptible instances 
=====================

An optional parameter allows you to set the instance to be preemptible. Both master and worker instances can be set to be preemptible. The following example shows a cluster configuration with a preemptible setting::

    cloud {
        imageId = 'centos-cloud/global/images/centos-7-v20181011'
        instanceType = 'n1-highcpu-8'
        preemptible = true
    }

Setting an instance to preemptible allows the administrator to kill the VM at will and may affect the pricing of the instance.

Cluster auto-scaling
====================

Nextflow integration for GCP provides a native support auto-scaling that allows the computing cluster to scale out or scale down i.e., add or remove computing nodes dynamically at runtime.

This is a critical feature, especially for pipelines crunching non-homogeneous datasets, because it allows the cluster to adapt dynamically to the actual workload computing resources need as they change over the time.

Cluster auto-scaling is enabled by adding the autoscale option group in the configuration file as shown below::

    cloud {
        imageId = 'centos-cloud/global/images/centos-7-v20181011'
        autoscale {
            enabled = true
            maxInstances = 10
        }
    }


The above example enables automatic cluster scale-out i.e. new instances are automatically launched and added to the
cluster when tasks remain too long in wait status because there aren't enough computing resources available. The
``maxInstances`` attribute defines the upper limit to which the cluster can grow.

By default unused instances are not removed when they are not utilised. If you want to enable automatic cluster scale-down
specify the ``terminateWhenIdle`` attribute in the ``autoscale`` configuration group.

It is also possible to define a different AMI image ID, type and spot price for instances launched by the Nextflow autoscaler.
For example::

    cloud {
        imageId = 'instance-xxx'
        instanceType = 'n1-highcpu-8'
        preemptible = false
        autoscale {
            enable = true
            preemptible = true
            minInstances = 5
            maxInstances = 10
            imageId = 'instance-yyy'
            instanceType = 'n1-highcpu-8'
            terminateWhenIdle = true
        }
    }

By doing that it's is possible to create a cluster with a single node i.e. the master node. Then the autoscaler will
automatically add the missing instances, up to the number defined by the ``minInstances`` attributes. 


Advanced configuration
======================

Read :ref:`Cloud configuration<config-cloud>` section to learn more about advanced cloud configuration options.


.. _google-pipelines-API:

Google Pipelines API
====================

Google Pipelines API is a managed computing service that allows the execution of containerized workloads in the Google Cloud Platform infrastructure.

Nextflow provides built-in support for Google Pipelines API which allows the seamless deployment of a Nextflow pipeline in the cloud, offloading the process executions as pipelines.

Download Google Cloud CLI and initialize using the following command::

    gcloud init

`Enable the Cloud Genomics, Compute Engine, and Cloud Storage APIs <https://console.cloud.google.com/flows/enableapi?apiid=genomics,compute,storage_api>`_

Create and download credentials for google cloud in the console:

* Go to APIs & Services → Credentials → Create Credentials
* Select Service account key.
* Choose the following Roles: Genomics > Genomics Admin, Storage > Storage Object Admin, & Service Account > Service Account User
* Download the json file and save as "creds.json".

Export as env variable::

    export GOOGLE_APPLICATION_CREDENTIALS=$PWD/creds.json

Create a nextflow.config file in the project root directory. The config must specify the following parameters:

* workDir - **Must be a GS bucket**
* process.executor - **googlepipelines**
* gce.project - **GCP project to run in**
* gce.zone *or* gce.region - You need to specify either one, **not** both. Multiple regions og zones can be specified by separating them with a comma (,).

Example::

    workDir = 'gs://<your bucket>/<directory>/'
     
    process {
        executor = 'googlepipelines'
    }
    
    cloud {
        instanceType = 'n1-standard-1'
    }
     
    gce {
        project = 'your-project-id'
        zone = 'us-central1-f,us-central-1-b'
    }

Note that all processes defined in your Nextflow scripts must define the following directives (or you can configure it globally for all processes in your ``nextflow.config``:

============ ===================================== ==================================================
Directive    Global Value                          Description
============ ===================================== ==================================================
container    process.container                      The docker container image to run the process in.
instanceType process.instanceType/gce.instanceType  The vm instanceType to run the process in.
============ ===================================== ==================================================

Finally, run your Nextflow script.



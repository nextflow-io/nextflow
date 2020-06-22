.. _google-page:

************
Google Cloud
************

Requirements
============

Nextflow
--------
The support for Google Cloud requires Nextflow version ``20.01.0``. To install it define the following variables
in your system environment::

    export NXF_VER=20.01.0
    export NXF_MODE=google


Credentials
-----------

To allow the deployment in the Google Cloud you need to configure the security credentials using
a *Security account key* JSON file.

Nextflow looks for this file using the ``GOOGLE_APPLICATION_CREDENTIALS`` variable that
has to be defined in the launching environment.

If you don't have it, download the credentials file from the Google Cloud Console following these steps:

* Open the `Google Cloud Console <https://console.cloud.google.com>`_
* Go to APIs & Services → Credentials
* Click on the *Create credentials* (blue) drop-down and choose *Service account key*, in the following page
* Select an existing *Service account* or create a new one if needed
* Select JSON as *Key type*
* Click the *Create* button and download the JSON file giving a name of your choice e.g. ``creds.json``.

Finally define the following variable replacing the path in the example with the one of your
credentials file just downloaded::

    export GOOGLE_APPLICATION_CREDENTIALS=/path/your/file/creds.json


.. _google-lifesciences:

Cloud Life Sciences
===================

`Cloud Life Sciences <https://cloud.google.com/life-sciences/>`_ is a managed computing service that allows the execution of
containerized workloads in the Google Cloud Platform infrastructure.

Nextflow provides built-in support for Cloud Life Sciences API which allows the seamless deployment of a Nextflow pipeline
in the cloud, offloading the process executions through the Google Cloud service.

.. note::
  This features requires Nextflow ``20.01.0-edge`` or later.

.. warning::
  This API works well for coarse-grained workloads i.e. long running jobs. It's not suggested the use
  this feature for pipelines spawning many short lived tasks.

.. _google-lifesciences-config:

Configuration
-------------

Make sure to have defined in your environment the ``GOOGLE_APPLICATION_CREDENTIALS`` variable.
See the section `Requirements`_ for details.

.. tip:: Make sure to have enabled Cloud Life Sciences API to use this feature. To learn how to enable it
  follow `this link <https://cloud.google.com/life-sciences/docs/quickstart>`_.

Create a ``nextflow.config`` file in the project root directory. The config must specify the following parameters:

* Google Life Sciences as Nextflow executor i.e. ``process.executor = 'google-lifesciences'``.
* The Docker container images to be used to run pipeline tasks e.g. ``process.container = 'biocontainers/salmon:0.8.2--1'``.
* The Google Cloud `project` ID to run in e.g. ``google.project = 'rare-lattice-222412'``.
* The Google Cloud `region` or `zone`. This is where the Compute Engine VMs will be started.
  You need to specify either one, **not** both. Multiple regions or zones can be specified by
  separating them with a comma e.g. ``google.zone = 'us-central1-f,us-central-1-b'``.

Example::

    process {
        executor = 'google-lifesciences'
        container = 'your/container:latest'
    }

    google {
        project = 'your-project-id'
        zone = 'europe-west1-b'
    }


.. warning:: Make sure to specify in the above setting the project ID not the project name.

.. Note:: A container image must be specified to deploy the process execution. You can use a different Docker image for
  each process using one or more :ref:`config-process-selectors`.

The following configuration options are available:

============================================== =================
Name                                           Description
============================================== =================
google.project                                 The Google Project Id to use for the pipeline execution.
google.region                                  The Google *region* where the computation is executed in Compute Engine VMs. Multiple regions can be provided separating them by a comma. Do not specify if a zone is provided. See  `available Compute Engine regions and zones <https://cloud.google.com/compute/docs/regions-zones/>`_
google.zone                                    The Google *zone* where the computation is executed in Compute Engine VMs. Multiple zones can be provided separating them by a comma. Do not specify if a region is provided. See  `available Compute Engine regions and zones <https://cloud.google.com/compute/docs/regions-zones/>`_
google.location                                The Google *location* where the job executions are deployed to Cloud Life Sciences API. See  `available Cloud Life Sciences API locations <https://cloud.google.com/life-sciences/docs/concepts/locations>`_ (default: the same as the region or the zone specified).
google.enableRequesterPaysBuckets              When ``true`` uses the configured Google project id as the billing project for storage access. This is required when accessing data from *reqester pays enabled* buckets. See `Requester Pays on Google Cloud Storage documentation  <https://cloud.google.com/storage/docs/requester-pays>`_ (default: ``false``)
google.lifeSciences.cpuPlatform                Set the minimum CPU Platform e.g `'Intel Skylake'`  See `Specifying a minimum CPU Platform for VM instances  <https://cloud.google.com/compute/docs/instances/specify-min-cpu-platform#specifications>`_ (default: none).
google.lifeSciences.bootDiskSize               Set the size of the virtual machine boot disk e.g `50.GB` (default: none).
google.lifeSciences.copyImage                  The container image run to copy input and output files. It must include the ``gsutil`` tool (default: ``google/cloud-sdk:alpine``).
google.lifeSciences.debug                      When ``true`` copies the `/google` debug directory in that task bucket directory (default: ``false``)
google.lifeSciences.preemptible                When ``true`` enables the usage of *preemptible* virtual machines or ``false`` otherwise (default: ``true``)
google.lifeSciences.usePrivateAddress          When ``true`` the VM will NOT be provided with a public IP address, and only contain an internal IP. If this option is enabled, the associated job can only load docker images from Google Container Registry, and the job executable cannot use external services other than Google APIs (default: ``false``). Requires version `20.03.0-edge` or later.
google.lifeSciences.sshDaemon                  When ``true`` runs SSH daemon in the VM carrying out the job to which it's possible to connect for debugging purposes (default: ``false``).
google.lifeSciences.sshImage                   The container image used to run the SSH daemon (default: ``gcr.io/cloud-genomics-pipelines/tools``).
============================================== =================


Process definition
------------------
Processes can be defined as usual and by default the ``cpus`` and ``memory`` directives are used to instantiate a custom
machine type with the specified compute resources.  If ``memory`` is not specified, 1GB of memory is allocated per cpu.
A persistent disk will be created with size corresponding to the ``disk`` directive.  If ``disk`` is not specified, the
instance default is chosen to ensure reasonable I/O performance.

The process ``machineType`` directive may optionally be used to specify a predefined Google Compute Platform `machine type <https://cloud.google.com/compute/docs/machine-types>`_
If specified, this value overrides the ``cpus`` and ``memory`` directives.
If the ``cpus`` and ``memory`` directives are used, the values must comply with the allowed custom machine type `specifications <https://cloud.google.com/compute/docs/instances/creating-instance-with-custom-machine-type#specifications>`_ .  Extended memory is not directly supported, however high memory or cpu predefined
instances may be utilized using the ``machineType`` directive

Examples::

    process custom_resources_task {
        cpus 8
        memory '40 GB'
        disk '200 GB'

        """
        <Your script here>
        """
    }

    process predefined_resources_task {
        machineType 'n1-highmem-8'

        """
        <Your script here>
        """
    }

.. note:: This feature requires Nextflow 19.07.0 or later.

Pipeline execution
------------------

The pipeline can be launched either in a local computer or a cloud instance. Pipeline input data can be stored either
locally or in a Google Storage bucket.

The pipeline execution must specify a Google Storage bucket where the workflow's intermediate results are stored using
the ``-work-dir`` command line options. For example::

    nextflow run <script or project name> -work-dir gs://my-bucket/some/path


.. tip:: Any input data **not** stored in a Google Storage bucket will automatically be transferred to the
  pipeline work bucket. Use this feature with caution being careful to avoid unnecessary data transfers.

Preemptible instances
---------------------

Preemptible instances are supported adding the following setting in the Nextflow config file::

    google {
        lifeSciences.preemptible = true
    }

Since this type of virtual machines can be retired by the provider before the job completion, it is advisable
to add the following retry strategy to your config file to instruct Nextflow to automatically re-execute a job
if the virtual machine was terminated preemptively::

    process {
      errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
      maxRetries = 5
    }


Hybrid execution
----------------

Nextflow allows the use of multiple executors in the same workflow application. This feature enables the deployment
of hybrid workloads in which some jobs are executed in the local computer or local computing cluster and
some other jobs are offloaded to Google Pipelines service.

To enable this feature use one or more :ref:`config-process-selectors` in your Nextflow configuration file to apply
the Google Pipelines *executor* only to a subset of processes in your workflow.
For example::


    process {
        withLabel: bigTask {
            executor = 'google-lifesciences'
            container = 'my/image:tag'
        }
    }

    google {
        project = 'your-project-id'
        zone = 'europe-west1-b'
    }


Then deploy the workflow execution using the ``-bucket-dir`` to specify a Google Storage path
for the jobs computed by the Google Pipeline service and, optionally, the ``-work-dir`` to
specify the local storage for the jobs computed locally::

    nextflow run <script or project name> -bucket-dir gs://my-bucket/some/path

.. warning:: The Google Storage path needs to contain at least sub-directory. Don't use only the
  bucket name e.g. ``gs://my-bucket``.

Limitation
----------

* Currently it's not possible to specify a disk type different from the default one assigned
  by the service depending the chosen instance type.



Troubleshooting
---------------

* Make sure to have enabled Compute Engine API, Life Sciences API and Cloud Storage Service in the
  `APIs & Services Dashboard <https://console.cloud.google.com/apis/dashboard>`_ page.

* Make sure to have enough compute resources to run your pipeline in your project
  `Quotas <https://console.cloud.google.com/iam-admin/quotas>`_ (i.e. Compute Engine CPUs,
  Compute Engine Persistent Disk, Compute Engine In-use IP addresses, etc).

* Make sure your security credentials allows you to access any Google Storage bucket
  where input data and temporary files are stored.

* Check the directory ``google/`` created in the task work directory (in the bucket storage) created
  when on job failure and containing useful information of the job execution. The creation
  can be enabled as default setting the option ``google.lifeSciences.debug = true`` in the
  Nextflow config file

* Enable the optional SSH daemon in the job VM using the option ``google.lifeSciences.sshDaemon = true``

* Make sure you are choosing a `location` where  `Cloud Life Sciences API is available <https://cloud.google.com/life-sciences/docs/concepts/locations>`_,
  and a `region` or `zone` where `Compute Engine is available <https://cloud.google.com/compute/docs/regions-zones/>`_.


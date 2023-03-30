.. _google-page:

************
Google Cloud
************

Credentials
===========

Credentials for submitting requests to the Google Cloud Batch and Cloud LifeSciences API are picked up from your
environment using `Application Default Credentials <https://github.com/googleapis/google-auth-library-java#google-auth-library-oauth2-http>`_.
Application Default Credentials are designed to use the credentials most natural to the
environment in which a tool runs.

The most common case will be to pick up your end-user Google credentials from your
workstation. You can create these by running the command::

    gcloud auth application-default login

and running through the authentication flow. This will write a credential file to your gcloud
configuration directory that will be used for any tool you run on your workstation that
picks up default credentials.

The next most common case would be when running on a Compute Engine VM. In this case,
Application Default Credentials will pick up the Compute Engine Service Account
credentials for that VM.

See the `Application Default Credentials <https://github.com/googleapis/google-auth-library-java#google-auth-library-oauth2-http>`_ documentation for how to enable other use cases.


Finally, the ``GOOGLE_APPLICATION_CREDENTIALS`` environment variable can be used to specify location
of the Google credentials file.

If you don't have it, the credentials file can be download from the Google Cloud Console following these steps:

* Open the `Google Cloud Console <https://console.cloud.google.com>`_
* Go to APIs & Services → Credentials
* Click on the *Create credentials* (blue) drop-down and choose *Service account key*, in the following page
* Select an existing *Service account* or create a new one if needed
* Select JSON as *Key type*
* Click the *Create* button and download the JSON file giving a name of your choice e.g. ``creds.json``.

Then, define the following variable replacing the path in the example with the one of your
credentials file just downloaded::

    export GOOGLE_APPLICATION_CREDENTIALS=/path/your/file/creds.json

.. _google-batch:

Cloud Batch
============

`Google Cloud Batch <https://cloud.google.com/batch>`_ is a managed computing service that allows the execution of containerized workloads in the
Google Cloud Platform infrastructure.

Nextflow provides built-in support for Cloud Batch which allows the seamless deployment of a Nextflow pipeline
in the cloud, offloading the process executions through the Google Cloud service.


Requirements
------------

The support for Google Batch requires Nextflow version ``22.07.1-edge`` or later. If you have already Nextflow
installed make sure to update to the latest `edge` release using these commands::

    export NXF_EDGE=1
    nextflow -self-update

If you don't have Nextflow, install it with command below::

    curl get.nextflow.io | bash

when done, make sure to use the latest `edge` release running the snippet in the previous paragraph.

.. _google-batch-config:

Configuration
-------------

Make sure to have defined in your environment the ``GOOGLE_APPLICATION_CREDENTIALS`` variable.
See the section `Credentials`_ for details.

.. note::
    Make sure your Google account is allowed to access the Google Cloud Batch service by checking
    the `APIs & Services <https://console.cloud.google.com/apis/dashboard>`_ dashboard.

Create or edit the file ``nextflow.config`` in your project root directory. The config must specify the following parameters:

* Google Cloud Batch as Nextflow executor i.e. ``process.executor = 'google-batch'``.
* The Docker container image to be used to run pipeline tasks e.g. ``process.container = 'biocontainers/salmon:0.8.2--1'``.
* The Google Cloud `project` ID to run in e.g. ``google.project = 'rare-lattice-222412'``.
* The Google location e.g. ``google.location = 'us-central1'``.

Example::

    process {
        executor = 'google-batch'
        container = 'your/container:latest'
    }

    google {
        project = 'your-project-id'
        location = 'us-central1'
    }

.. note::
  Make sure to specify the project ID, not the project name.

.. note::
  Make sure to specify a location where Google Batch is available. Refer to the
  `Google Batch documentation <https://cloud.google.com/batch/docs/get-started#locations>`_
  for region availability.

.. Note::
  A container image must be specified to deploy the process execution. You can use a different Docker image for
  each process using one or more :ref:`config-process-selectors`.

The following configuration options are available:

============================================== =================
Name                                           Description
============================================== =================
google.project                                 The Google Project Id to use for the pipeline execution.
google.location                                The Google *location* where the job executions are deployed (default: ``us-central1``).
google.enableRequesterPaysBuckets              When ``true`` uses the configured Google project id as the billing project for storage access. This is required when accessing data from *requester pays enabled* buckets. See `Requester Pays on Google Cloud Storage documentation  <https://cloud.google.com/storage/docs/requester-pays>`_ (default: ``false``).
google.batch.allowedLocations                  Define the set of allowed locations for VMs to be provisioned. See `Google documentation <https://cloud.google.com/batch/docs/reference/rest/v1/projects.locations.jobs#locationpolicy>`_ for details (default: no restriction. Requires version ``22.12.0-edge`` or later).
google.batch.bootDiskSize                      Set the size of the virtual machine boot disk, e.g ``50.GB`` (default: none).
google.batch.cpuPlatform                       Set the minimum CPU Platform, e.g. ``'Intel Skylake'``. See `Specifying a minimum CPU Platform for VM instances <https://cloud.google.com/compute/docs/instances/specify-min-cpu-platform#specifications>`_ (default: none).
google.batch.spot                              When ``true`` enables the usage of *spot* virtual machines or ``false`` otherwise (default: ``false``).
google.batch.usePrivateAddress                 When ``true`` the VM will NOT be provided with a public IP address, and only contain an internal IP. If this option is enabled, the associated job can only load docker images from Google Container Registry, and the job executable cannot use external services other than Google APIs (default: ``false``).
google.batch.network                           Set network name to attach the VM's network interface to. The value will be prefixed with global/networks/ unless it contains a /, in which case it is assumed to be a fully specified network resource URL. If unspecified, the global default network is used.
google.batch.serviceAccountEmail               Define the Google service account email to use for the pipeline execution. If not specified, the default Compute Engine service account for the project will be used.
google.batch.subnetwork                        Define the name of the subnetwork to attach the instance to must be specified here, when the specified network is configured for custom subnet creation. The value is prefixed with `regions/subnetworks/` unless it contains a `/`, in which case it is assumed to be a fully specified subnetwork resource URL.
============================================== =================


Process definition
------------------
Processes can be defined as usual and by default the ``cpus`` and ``memory`` directives are used to find the cheapest machine
type available at current location that fits the requested resources. If ``memory`` is not specified, 1GB of memory is allocated per cpu.

As of version ``23.02.0-edge``, the process ``machineType`` directive can be a list of patterns separated by comma. The pattern can contain a ``*`` to match
any number of characters and ``?`` to match any single character. Examples of valid patterns: ``c2-*``, ```m?-standard*``, ```n*``.

Alternatively it can also be used to define a specific predefined Google Compute Platform `machine type <https://cloud.google.com/compute/docs/machine-types>`_
or a custom machine type.

Examples::

    process automatic_resources_task {
        cpus 8
        memory '40 GB'

        """
        <Your script here>
        """
    }

    process allowing_some_series {
        cpus 8
        memory '20 GB'
        machineType 'n2-*,c2-*,m3-*'

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

Pipeline execution
------------------

The pipeline can be launched either in a local computer or a cloud instance. Pipeline input data can be stored either
locally or in a Google Storage bucket.

The pipeline execution must specify a Google Storage bucket where the workflow's intermediate results are stored using
the ``-work-dir`` command line options. For example::

    nextflow run <script or project name> -work-dir gs://my-bucket/some/path

.. tip::
  Any input data **not** stored in a Google Storage bucket will automatically be transferred to the
  pipeline work bucket. Use this feature with caution being careful to avoid unnecessary data transfers.


.. warning::
  The Google Storage path needs to contain at least sub-directory. Don't use only the bucket name e.g. ``gs://my-bucket``.


Spot instances
---------------

Spot instances are supported adding the following setting in the Nextflow config file::

    google {
        batch.spot = true
    }

Since this type of virtual machines can be retired by the provider before the job completion, it is advisable
to add the following retry strategy to your config file to instruct Nextflow to automatically re-execute a job
if the virtual machine was terminated preemptively::

    process {
        errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
        maxRetries = 5
    }

Fusion file system
------------------

As of version ``23.02.0-edge``, Google Batch executor supports the use of :ref:`fusion-page`.

Fusion allows the use of Google Storage as a virtual distributed file system, optimising the data transfer
and speeding up most job I/O operations.

To enable the use of Fusion file system in your pipeline, add the following snippet in your Nextflow configuration file::

    fusion.enabled = true
    wave.enabled = true
    process.scratch = false
    tower.accessToken = '<YOUR ACCESS TOKEN>'

The `Tower <https://cloud.tower.nf>`_ access token is optional, but it enables higher API rate limits for the
:ref:`wave-page` service required by Fusion.

.. tip::
  When Fusion is enabled, by default, only machine types that allow to attach local SSD disks will be used. If you specify your own
  machine type or machine series they should allow to attach local SSD disks, otherwise the job scheduling will fail.



Supported directives
--------------------

The integration with Google Batch is a developer preview feature. Currently the following Nextflow directives are
supported:

* :ref:`process-accelerator`
* :ref:`process-container`
* :ref:`process-containeroptions`
* :ref:`process-cpus`
* :ref:`process-disk`
* :ref:`process-executor`
* :ref:`process-machinetype`
* :ref:`process-memory`
* :ref:`process-time`


.. _google-lifesciences:

Cloud Life Sciences
===================

Requirements
------------
The support for Google Cloud requires Nextflow version ``20.01.0`` or later. To install it define the following variables
in your system environment::

    export NXF_VER=20.01.0
    export NXF_MODE=google

.. note:: As of version ``21.04.0`` or later the above variables are not required anymore and therefore should not be used.


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
See the section `Credentials`_ for details.

.. tip:: Make sure to have enabled Cloud Life Sciences API to use this feature. To learn how to enable it
  follow `this link <https://cloud.google.com/life-sciences/docs/quickstart>`_.

Create a ``nextflow.config`` file in the project root directory. The config must specify the following parameters:

* Google Life Sciences as Nextflow executor i.e. ``process.executor = 'google-lifesciences'``.
* The Docker container image to be used to run pipeline tasks e.g. ``process.container = 'biocontainers/salmon:0.8.2--1'``.
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
google.enableRequesterPaysBuckets              When ``true`` uses the configured Google project id as the billing project for storage access. This is required when accessing data from *requester pays enabled* buckets. See `Requester Pays on Google Cloud Storage documentation  <https://cloud.google.com/storage/docs/requester-pays>`_ (default: ``false``)
google.lifeSciences.bootDiskSize               Set the size of the virtual machine boot disk e.g ``50.GB`` (default: none).
google.lifeSciences.copyImage                  The container image run to copy input and output files. It must include the ``gsutil`` tool (default: ``google/cloud-sdk:alpine``).
google.lifeSciences.cpuPlatform                Set the minimum CPU Platform e.g. ``'Intel Skylake'``. See `Specifying a minimum CPU Platform for VM instances <https://cloud.google.com/compute/docs/instances/specify-min-cpu-platform#specifications>`_ (default: none).
google.lifeSciences.debug                      When ``true`` copies the ``/google`` debug directory in that task bucket directory (default: ``false``)
google.lifeSciences.preemptible                When ``true`` enables the usage of *preemptible* virtual machines or ``false`` otherwise (default: ``true``)
google.lifeSciences.usePrivateAddress          When ``true`` the VM will NOT be provided with a public IP address, and only contain an internal IP. If this option is enabled, the associated job can only load docker images from Google Container Registry, and the job executable cannot use external services other than Google APIs (default: ``false``). Requires version ``20.03.0-edge`` or later.
google.lifeSciences.network                    Set network name to attach the VM's network interface to. The value will be prefixed with global/networks/ unless it contains a /, in which case it is assumed to be a fully specified network resource URL. If unspecified, the global default network is used. Requires version ``21.03.0-edge`` or later.
google.lifeSciences.serviceAccountEmail        Define the Google service account email to use for the pipeline execution. If not specified, the default Compute Engine service account for the project will be used. Requires version ``20.05.0-edge`` or later.
google.lifeSciences.subnetwork                 Define the name of the subnetwork to attach the instance to must be specified here, when the specified network is configured for custom subnet creation. The value is prefixed with `regions/subnetworks/` unless it contains a `/`, in which case it is assumed to be a fully specified subnetwork resource URL. Requires version ``21.03.0-edge`` or later.
google.lifeSciences.sshDaemon                  When ``true`` runs SSH daemon in the VM carrying out the job to which it's possible to connect for debugging purposes (default: ``false``).
google.lifeSciences.sshImage                   The container image used to run the SSH daemon (default: ``gcr.io/cloud-genomics-pipelines/tools``).
google.lifeSciences.keepAliveOnFailure         When ``true`` and a task complete with an unexpected exit status the associated computing node is kept up for 1 hour. This options implies ``sshDaemon=true`` (default: ``false``, requires Nextflow version ``21.06.0-edge`` or later).
google.storage.delayBetweenAttempts            Delay between download attempts from Google Storage (default `10 sec`, requires version ``21.06.0-edge`` or later).
google.storage.maxParallelTransfers            Max parallel upload/download transfer operations *per job* (default: ``4``, requires version ``21.06.0-edge`` or later).
google.storage.maxTransferAttempts             Max number of downloads attempts from Google Storage (default: `1`, requires version ``21.06.0-edge`` or later).
google.storage.parallelThreadCount             Defines the value for the option ``GSUtil:parallel_thread_count`` used by ``gsutil`` for transfer input and output data (default: ``1``, requires version ``21.06.0-edge`` or later).
google.storage.downloadMaxComponents           Defines the value for the option ``GSUtil:sliced_object_download_max_components`` used by ``gsutil`` for transfer input and output data (default: ``8``, requires version ``21.06.0-edge`` or later).
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

.. note:: Preemptible instances have a `runtime limit <https://cloud.google.com/compute/docs/instances/preemptible>`_ of 24 hours.

.. tip:: For an exhaustive list of all possible error codes, please refer to the official Google LifeSciences `documentation <https://cloud.google.com/life-sciences/docs/troubleshooting#error_codes>`_.

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

Quotas
------

Compute resources in Google Cloud are subject to `resource quotas <https://cloud.google.com/compute/quotas>`_ which may affect your ability to run pipelines at scale. You can request quota increases, and your quotas may automatically increase over time as you use the platform. In particular, GPU quotas are initially set to 0, so you must explicitly request a quota increase in order to use GPUs. Initially you can request an increase to 1 GPU at a time, and after one billing cycle you may be able to increase it further.

Limitations
-----------

* Currently it's not possible to specify a disk type different from the default one assigned
  by the service depending on the chosen instance type.



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


.. _charliecloud-page:

************************
Charliecloud containers
************************

`Charliecloud <https://hpc.github.io/charliecloud>`_ is an alternative Container engine to Docker, that is better
suited for use in HPC environments. Its main advantage is that it can be used with unpriviledged permissions,
making use of user namespaces in the Linux kernel. Charliecloud is able to pull from Docker registries.

By using this feature any process in a Nextflow script can be transparently executed into a Charliecloud container. This may
be extremely useful to package the binary dependencies of a script into a standard and portable format that can be 
executed on any platform supporting the Charliecloud engine.

.. note::
    This feature requires Nextflow version ``21.03.0-edge`` or later and Charliecloud ``v0.22` or later.

.. warning::
    This is an incubating feature. The use in production environment is not recommended.

Prerequisites
==============

You will need Charliecloud version ``0.22`` or later installed on your execution environment e.g. your computer or a
distributed cluster, depending on where you want to run your pipeline.

How it works
=============

You won't need to modify your Nextflow script in order to run it with Charliecloud. Simply specify the docker image from
where the containers are started by using the ``-with-charliecloud`` command line option. For example::

  nextflow run <your script> -with-charliecloud [container]

Every time your script launches a process execution, Nextflow will run it into a charliecloud container created by using the
specified container image. In practice Nextflow will automatically wrap your processes and run them by executing the ``ch-run``
command with the container you have provided.

.. note:: A container image can contain any tool or piece of software you may need to carry out a process execution. Moreover the
  container is run in such a way that the process result files are created in the hosting file system, thus
  it behaves in a completely transparent manner without requiring extra steps or affecting the flow in your pipeline.

If you want to avoid entering the Container image as a command line parameter, you can define it in the Nextflow configuration
file. For example you can add the following lines in the ``nextflow.config`` file::

    process.container = '/path/to/container'
    charliecloud.enabled = true

.. warning::
    If an absolute is provided, the container needs to be in the Charliecloud flat directory format.
    See section below for compatibility with Docker registries.

Read the :ref:`config-page` page to learn more about the ``nextflow.config`` file and how to use it to configure
your pipeline execution.

.. warning::
    Nextflow automatically manages the file system mounts each time a container is launched depending on the process
    input files. Note, however, that when a process input is a *symbolic link* file, the linked file **must** be stored
    in the same folder where the symlink is located, or any its sub-folder. Otherwise the process execution will fail because the
    launched container won't be able to access the linked file.

Charliecloud & Docker Hub
=========================

Nextflow is able to transparently pull remote container images stored in any Docker compatible registry and converts
them to the Charliecloud compatible format.

By default when a container name is specified, Nextflow checks if a container with that name exists in the local file
system. If it exists, it's used to execute the container. If a matching file does not exist,
Nextflow automatically tries to pull an image with the specified name from the Docker Hub.

To pull images from a third party Docker registry simply provide the URL to the image. If no URL is provided,
Docker Hub is assumed. For example this can be used to pull an image from quay.io and convert it automatically
to the Charliecloud container format::

    process.container = 'https://quay.io/biocontainers/multiqc:1.3--py35_2'
    charliecloud.enabled = true
 
Whereas this would pull from Docker Hub::

    process.container = 'nextflow/examples:latest'
    charliecloud.enabled = true


Multiple containers
===================

It is possible to specify a different container for each process definition in your pipeline script. For example,
let's suppose you have two processes named ``foo`` and ``bar``. You can specify two different containers
in the ``nextflow.config`` file as shown below::

    process {
        withName:foo {
            container = 'image_name_1'
        }
        withName:bar {
            container = 'image_name_2'
        }
    }
    charliecloud {
        enabled = true
    }

Read the :ref:`Process scope <config-process>` section to learn more about processes configuration.


Advanced settings 
=================

Charliecloud advanced configuration settings are described in :ref:`config-charliecloud` section in the Nextflow
configuration page.

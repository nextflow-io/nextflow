.. _singularity-page:

**********************
Singularity containers
**********************

`Singularity <http://singularity.lbl.gov/>`_ is a container engine alternative to Docker. The main advantages
of Singularity is that it can be used with unprivileged permissions and doesn't require a separate daemon process.
These, along other features, like for example the support for autofs mounts, makes Singularity a container engine
better suited the requirements of HPC workloads.

Nextflow provides built-in support for Singularity. This allows you to precisely control the execution environment
of the processes in your pipeline by running them in isolated containers along all their dependencies.

Moreover the support provided by Nextflow for different container technologies, allows the same pipeline to be
transparently executed both with Docker or Singularity containers, depending the available engine in the target
execution platforms.


Prerequisites
=============

You will need Singularity installed on your execution environment e.g. your computer or a distributed cluster, depending
on where you want to run your pipeline.

Images
======

Singularity makes use of a container image file, which physically contains the container. Refer to the `Singularity
documentation <http://singularity.lbl.gov/user-guide>`_ to learn how create Singularity images.

Docker images can automatically be converted to Singularity images by using the
`docker2singularity <https://github.com/singularityware/docker2singularity>`_ tool.


Singularity allows paths that do not currently exist within the container to be created
and mounted dynamically by specifying them on the command line. However this feature is only supported on hosts
that support the `Overlay file system <https://en.wikipedia.org/wiki/OverlayFS>`_ and is not enabled by default.

.. note::
    Nextflow expects that data paths are defined system wide, and your Singularity images need to be created having the
    mount paths defined in the container file system.

If your Singularity installation support the `user bind control` feature,
enable the Nextflow support for that by defining the ``singularity.autoMounts = true`` setting in the Nextflow
configuration file.


How it works
============

The integration for Singularity follows the same execution model implemented for Docker. You won't need to modify your
Nextflow script in order to run it with Singularity. Simply specify the Singularity image
file from where the containers are started by using the ``-with-singularity`` command line option. For example::

  nextflow run <your script> -with-singularity [singularity image file]

Every time your script launches a process execution, Nextflow will run it into a Singularity container created by using the
specified image. In practice Nextflow will automatically wrap your processes and launch them by running the
``singularity exec`` command with the image you have provided.

.. note:: A Singularity image can contain any tool or piece of software you may need to carry out a process execution.
  Moreover the container is run in such a way that the process result files are created in the hosting file system, thus
  it behaves in a completely transparent manner without requiring extra steps or affecting the flow in your pipeline.

If you want to avoid entering the Singularity image as a command line parameter, you can define it in the Nextflow
configuration file. For example you can add the following lines in the ``nextflow.config`` file::

    process.container = '/path/to/singularity.img'
    singularity.enabled = true

In the above example replace ``/path/to/singularity.img`` with any Singularity image of your choice.

Read the :ref:`config-page` page to learn more about the ``nextflow.config`` file and how to use it to configure
your pipeline execution.

.. note::
   Unlike Docker, Nextflow does not mount automatically host paths in the container when using Singularity.
   It expects they are configure and mounted system wide by the Singularity runtime. If your Singularity installation
   allows `user defined bind points` read the :ref:`Singularity configuration <config-singularity>` section to learn
   how to enable Nextflow auto mounts.

.. warning::
    When a process input is a *symbolic link* file, make sure the linked file **must** is stored in a host folder
    that is accessible from a bind path defined in your Singularity installation. Otherwise the process execution
    will fail because the launched container won't be able to access the linked file.


Multiple containers
===================

It is possible to specify a different Singularity image for each process definition in your pipeline script. For example,
let's suppose you have two processes named ``foo`` and ``bar``. You can specify two different Singularity images
specifing them in the ``nextflow.config`` file as shown below::

    process.$foo.container = 'image_name_1'
    process.$bar.container = 'image_name_2'
    singularity.enabled = true


Read the :ref:`Process scope <config-process>` section to learn more about processes configuration.


Singularity & Docker Hub
========================

Nextflow is able to transparently pull remote container images stored in the `Singularity-Hub <https://singularity-hub.org/>`_
or any Docker compatible registry.

.. note:: This feature requires you have installed Singularity 2.3.x or higher

Remote container images can specified in your Nextflow script or configuration file by simply prefixing the image name
with the ``shub://`` or ``docker://`` pseudo-protocol as required by the Singularity tool. For example::

    process.container = 'docker://debian:wheezy'
    singularity.enabled = true


Nextflow caches those images in the ``singularity`` directory in the pipeline work directory by default. However it is
suggest to provide a centralised caching directory by using either the ``NXF_SINGULARITY_CACHEDIR`` environment variable
or the ``singularity.cacheDir`` setting in the Nextflow config file.

.. warning:: When using a computing cluster the Singularity cache directory must be a shared folder accessible
  to all computing nodes.

.. error::  When pulling Docker images Singularity may be unable to determine the container size if the image was
  stored by using an old Docker format, resulting in a pipeline execution error. See the Singularity documentation for details.

Advanced settings
=================

Singularity advanced configuration settings are described in :ref:`config-singularity` section in the Nextflow
configuration page.














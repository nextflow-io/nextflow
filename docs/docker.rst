.. _docker-page:

*******************
Docker containers
*******************

Nextflow integration with `Docker containers <http://www.docker.io>`_ technology allows you to write self-contained
and truly reproducible computational pipelines.

By using this feature any process in a Nextflow script can be transparently executed into a Docker container. This may
be extremely useful to package the binaries dependencies of a script into a standard and portable format that can be 
executed on any platform supporting the Docker engine.

Prerequisites
==============

You will need Docker installed in your execution environment e.g. your computer or a distributed cluster, depending
where you want to run your pipeline.

If you are running Docker on Mac OSX make sure you are mounting your local ``/Users`` directory into the Docker VM as
explained in this excellent tutorial: `How to use Docker on OSX <http://viget.com/extend/how-to-use-docker-on-os-x-the-missing-guide>`_.


How it works
=============

In order to execute your Nextflow script with Docker you won't need to modify it. You simply have to specify the Docker
image from which containers have to be started by using the ``-with-docker`` command line option. For example::

  nextflow run <your script> -with-docker <docker image>

Every time your script launches a process execution Nextflow will run it into a Docker container created by using the image
specified. In practice Nextflow will automatically wrap your processes and run them by executing the ``docker run``
command with the image you have provided.

.. note:: A Docker image can contain any tool or piece of software you may need to carry out a process execution. Moreover the
  container is run is a such a way that the process result files are created in the hosting file system, in other words
  it behaves in a complete transparent manner without requiring extra steps or affecting the flow in your pipeline.

If you want to avoid to enter the Docker image as a command line parameter, you can define it in the Nextflow configuration
file. For example you can add the following lines in the ``nextflow.config`` file::

    process.container = 'nextflow/examples:latest'
    docker.enabled = true

In about example replace ``nextflow/examples:latest`` with any Docker image of your choice.

Read the :ref:`config-page` page to learn more about the ``nextflow.config`` file and how use it to configure
your pipeline execution.

Multiple containers
=====================

It is possible to specify a different Docker image for each process definition in your pipeline script. Let's
suppose you have two processes named ``foo`` and ``bar``. You can specify two different Docker images for them
in Nextflow script as shown below::

    process foo {
      container 'image_name_1'

      '''
      do this
      '''
    }

    process bar {
      container 'image_name_2'

      '''
      do that
      '''
    }


Alternatively, the same containers definition can be provided by using the ``nextflow.config`` as shown below::

    process.$foo.container = 'image_name_1'
    process.$bar.container = 'image_name_2'
    docker.enabled = true


Read the :ref:`Process scope <config-process>` section to learn more about processes configuration.


Advanced settings 
==================

Docker advanced configuration settings are described in :ref:`config-docker` section in the Nextflow configuration page.














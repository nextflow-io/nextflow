.. _docker-page:

******************
Docker containers
******************

Nextflow integration with `Docker containers <http://www.docker.io>`_ technology allows you to write self-contained
and truly reproducible computational pipelines.

By using this feature any process in a Nextflow script can be transparently executed into a Docker container. This may
be extremely useful to package the binary dependencies of a script into a standard and portable format that can be 
executed on any platform supporting the Docker engine.

Prerequisites
==============

You will need Docker installed on your execution environment e.g. your computer or a distributed cluster, depending
on where you want to run your pipeline.

If you are running Docker on Mac OSX make sure you are mounting your local ``/Users`` directory into the Docker VM as
explained in this excellent tutorial: `How to use Docker on OSX <http://viget.com/extend/how-to-use-docker-on-os-x-the-missing-guide>`_.


How it works
=============

You won't need to modify your Nextflow script in order to run it with Docker. Simply specify the Docker image from
where the containers are started by using the ``-with-docker`` command line option. For example::

  nextflow run <your script> -with-docker [docker image]

Every time your script launches a process execution, Nextflow will run it into a Docker container created by using the
specified image. In practice Nextflow will automatically wrap your processes and run them by executing the ``docker run``
command with the image you have provided.

.. note:: A Docker image can contain any tool or piece of software you may need to carry out a process execution. Moreover the
  container is run in such a way that the process result files are created in the hosting file system, thus
  it behaves in a completely transparent manner without requiring extra steps or affecting the flow in your pipeline.

If you want to avoid entering the Docker image as a command line parameter, you can define it in the Nextflow configuration
file. For example you can add the following lines in the ``nextflow.config`` file::

    process.container = 'nextflow/examples:latest'
    docker.enabled = true

In the above example replace ``nextflow/examples:latest`` with any Docker image of your choice.

Read the :ref:`config-page` page to learn more about the ``nextflow.config`` file and how to use it to configure
your pipeline execution.

.. warning::
    Nextflow automatically manages the file system mounts each time a container is launched depending on the process
    input files. Note, however, that when a process input is a *symbolic link* file, the linked file **must** be stored
    in the same folder where the symlink is located, or any its sub-folder. Otherwise the process execution will fail because the
    launched container won't be able to access the linked file.


Multiple containers
=====================

It is possible to specify a different Docker image for each process definition in your pipeline script. Let's
suppose you have two processes named ``foo`` and ``bar``. You can specify two different Docker images for them
in the Nextflow script as shown below::

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


Alternatively, the same containers definitions can be provided by using the ``nextflow.config`` file as shown below::

    process {
        withName:foo {
            container = 'image_name_1'
        }
        withName:bar {
            container = 'image_name_2'
        }
    }
    docker {
        enabled = true
    }


Read the :ref:`Process scope <config-process>` section to learn more about processes configuration.


Advanced settings 
==================

Docker advanced configuration settings are described in :ref:`config-docker` section in the Nextflow configuration page.














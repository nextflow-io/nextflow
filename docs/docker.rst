.. _docker-page:

*******************
Docker integration
*******************

Nextflow integration with `Docker container <http://www.docker.io>`_ technology allows you to write truly self-contained
and reproducible computational pipelines.

By using this feature any process in a pipeline script can be transparently executed into a Docker container. This may
be extremely useful to package the binaries dependencies of a script into a standard and portable for format
that can be executed on any platform supporting Docker.


Prerequisites
==============

You will need Docker installed in your target execution environment e.g. your computer or a distributed cluster depending
where you want to run your scripts.

If you are running Docker on Mac OSX make sure you are mounting your local ``/Users`` directory into the Docker VM as
explained in this tutorial `How to use Docker on OSX: The Missing Guide <http://viget.com/extend/how-to-use-docker-on-os-x-the-missing-guide>`_.


How it works
=============

In order to execute your pipeline script with Docker you won't need to modify it. You simply need to specify the Docker
image from which container has to be started by using the ``-with-docker`` command line option. For example::

  nextflow run <your script> -with-docker <docker image>

Each time a process execution is spawned, Nextflow will run it into a Docker contained created by using the specified image.
In practice Nextflow will automatically wrap your processes and run them by executing a ``docker run`` by the image you
have specified.

The Docker image can contain any tool or piece of software you may need to carry out the process execution. Moreover the
container is run is a such a way that the process result files are created in the hosting execution environment, in other words
in behaves in complete transparent manner without requiring extra steps or affecting the flow in your pipeline.

If you want to avoid to enter the Docker image as a command line parameter, you can define it in the Nextflow configuration
file. For example you can add the following lines in the ``nextflow.config`` file::

    process.container = 'nextflow/examples:latest'
    docker.enabled = true

In about example replace ``nextflow/examples:latest`` with any Docker image of your choice.

Multiple containers
=====================

It is also possible to specify a different Docker image for each process defined in your pipeline script. Let's
suppose you have have defined two processes named ``foo`` and ``bar``. You can specify two different Docker images
in the ``nextflow.config`` file as shown below::

    process.$boo.container = 'image_name_1'
    process.$bar.container = 'image_name_2'
    docker.enabled = true


Read :ref:`Process scope <config-process>` section to learn more about processes configuration.















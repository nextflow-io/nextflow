.. _shifter-page:

******************
Shifter Containers
******************

`Shifter <https://docs.nersc.gov/programming/shifter/overview/>`__ is container engine alternative to
`Docker <https://www.docker.com>`__. Shifter works by converting Docker images to a common format that can then be
distributed and launched on HPC systems. The user interface to Shifter enables a user to select an image
from the `dockerhub registry <https://hub.docker.com/>`__ and then submit jobs which run entirely within the container.

Nextflow provides built-in support for Shifter. This allows you to control the execution environment of the processes
in your pipeline by running them in isolated containers along-side their dependencies.

Moreover the support provided by Nextflow for different container technologies, allows the same pipeline to be
transparently executed both with :ref:`Docker <_docker-page>`, :ref:`Singularity <_singularity-page>` or
:ref:`Shifter <_shifter-page>` containers, depending on the available engine in target execution platforms.

Prerequisites
=============

You need Shifter and Shifter image gateway installed in your execution environment, i.e: your personal computed or the
entry node of a distributed cluster. In the case of the distributed cluster case, you should have Shifter installed on
all of the compute nodes and the ``shifterimg`` command should also be available and Shifter properly setup to access the
Image gateway, for more information see the
`official documentation <https://github.com/NERSC/shifter/tree/master/doc>`__.

.. note:: This feature requires Shifter version 18.03 (or later) and Nextflow 19.10.0 (or later).

Images
======

Shifter converts a docker image to squashfs layers which are distributed and launched in the cluster. For more info on
how to Build Shifter images see the
`official documentation <https://docs.nersc.gov/programming/shifter/how-to-use/#building-shifter-images>`__.

How it works
============

The integration for Shifter, at this time, requires you to set up the following parameters in your config file::

  process.container = "dockerhub_user/image_name:image_tag"
  shifter.enabled = true

and it will always try to search the Docker Hub registry for the images.

.. note:: if you do not specify an image tag it will fetch the 'latest' tag by default.

Multiple Containers
===================

It is possible to specify a different Shifter image for each process definition in your pipeline script. For example,
let's suppose you have two processes named ``foo`` and ``bar``. You can specify two different Shifter images
specifying them in the ``nextflow.config`` file as shown below::

    process {
        withName:foo {
            container = 'image_name_1'
        }
        withName:bar {
            container = 'image_name_2'
        }
    }
    shifter {
        enabled = true
    }

Read the :ref:`Process scope <config-process>` section to learn more about processes configuration.

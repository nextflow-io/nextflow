.. _podman-page:

******************
Podman containers
******************

Nextflow integration with `Podman containers <http://www.podman.io>`_ technology allows you to write self-contained
and truly reproducible computational pipelines.

By using this feature any process in a Nextflow script can be transparently executed into a Podman container. This may
be extremely useful to package the binary dependencies of a script into a standard and portable format that can be 
executed on any platform supporting the Podman engine.

.. note::
    This feature requires Nextflow version ``20.01.0`` or later.

.. warning::
    This is an incubating feature. The use in production environment is not recommended.

Prerequisites
==============

You will need Podman installed on your execution environment e.g. your computer or a distributed cluster, depending
on where you want to run your pipeline. Running in rootless mode requires appropriate OS configuration. Due to current
Podman limits using cpuset for cpus and memory such is only possible using sudo.


How it works
=============

You won't need to modify your Nextflow script in order to run it with Podman. Simply specify the Podman image from
where the containers are started by using the ``-with-podman`` command line option. For example::

  nextflow run <your script> -with-podman [OCI container image]

Every time your script launches a process execution, Nextflow will run it into a Podman container created by using the
specified image. In practice Nextflow will automatically wrap your processes and run them by executing the ``podman run``
command with the image you have provided.

.. note:: A OCI container image can contain any tool or piece of software you may need to carry out a process execution. Moreover the
  container is run in such a way that the process result files are created in the hosting file system, thus
  it behaves in a completely transparent manner without requiring extra steps or affecting the flow in your pipeline.

If you want to avoid entering the Podman image as a command line parameter, you can define it in the Nextflow configuration
file. For example you can add the following lines in the ``nextflow.config`` file::

    process.container = 'nextflow/examples:latest'
    podman.enabled = true

In the above example replace ``nextflow/examples:latest`` with any Podman image of your choice.

Read the :ref:`config-page` page to learn more about the ``nextflow.config`` file and how to use it to configure
your pipeline execution.

.. warning::
    Nextflow automatically manages the file system mounts each time a container is launched depending on the process
    input files. Note, however, that when a process input is a *symbolic link* file, the linked file **must** be stored
    in the same folder where the symlink is located, or any its sub-folder. Otherwise the process execution will fail because the
    launched container won't be able to access the linked file.


Multiple containers
=====================

It is possible to specify a different container image for each process definition in your pipeline script. Let's
suppose you have two processes named ``foo`` and ``bar``. You can specify two different container images for them
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
    podman {
        enabled = true
    }


Read the :ref:`Process scope <config-process>` section to learn more about processes configuration.


Advanced settings 
==================

Podman advanced configuration settings are described in :ref:`config-podman` section in the Nextflow configuration page.



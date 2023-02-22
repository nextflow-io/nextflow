.. _container-page:

**********
Containers
**********

Nextflow supports a variety of container runtimes. Containerization allows you to write
self-contained and truly reproducible computational pipelines, by packaging the binary
dependencies of a script into a standard and portable format that can be executed on any
platform that supports a container runtime. Furthermore, the same pipeline can be transparently
executed with any of the supported container runtimes, depending on which runtimes are available
in the target compute environment.

.. note::
  When creating your container image to use with Nextflow make sure it includes Bash (3.x or later) and ``ps`` are installed in your
  image, along with other tools required for collecting metrics (See the list :ref:`here <trace-required-packages>`). Also, Bash should be available on the path `/bin/bash` and it should be the container entry-point. 

.. _container-charliecloud:

Charliecloud
============

`Charliecloud <https://hpc.github.io/charliecloud>`_ is an alternative container runtime to Docker, that is better
suited for use in HPC environments. Its main advantage is that it can be used without root privileges,
making use of user namespaces in the Linux kernel. Charliecloud is able to pull from Docker registries.

.. note::
    This feature requires at least Nextflow version ``21.03.0-edge`` and a Charliecloud version between ``v0.22`` and ``v0.27``.
    As of Nextflow version ``22.09.0-edge``, Charliecloud ``v0.28`` or later is required.

.. warning::
    This feature is experimental. Using it in a production environment is not recommended.

Prerequisites
-------------

You will need Charliecloud installed in your execution environment e.g. on your computer or a
distributed cluster, depending on where you want to run your pipeline.

How it works
------------

You won't need to modify your Nextflow script in order to run it with Charliecloud. Simply specify the docker image from
where the containers are started by using the ``-with-charliecloud`` command line option. For example::

  nextflow run <your script> -with-charliecloud [container]

Every time your script launches a process execution, Nextflow will run it into a charliecloud container created by using the
specified container image. In practice Nextflow will automatically wrap your processes and run them by executing the ``ch-run``
command with the container you have provided.

.. note::
    A container image can contain any tool or piece of software you may need to carry out a process execution. Moreover, the
    container is run in such a way that the process result files are created in the host file system, thus
    it behaves in a completely transparent manner without requiring extra steps or affecting the flow in your pipeline.

If you want to avoid entering the Container image as a command line parameter, you can define it in the Nextflow configuration
file. For example you can add the following lines in the ``nextflow.config`` file::

    process.container = '/path/to/container'
    charliecloud.enabled = true

.. warning::
    If an absolute path is provided, the container needs to be in the Charliecloud flat directory format.
    See the section below on compatibility with Docker registries.

Read the :ref:`config-page` page to learn more about the ``nextflow.config`` file and how to use it to configure
your pipeline execution.

.. warning::
    Nextflow automatically manages the file system mounts whenever a container is launched depending on the process
    input files. However, when a process input is a *symbolic link*, the linked file **must** be stored
    in the same folder where the symlink is located, or a sub-folder of it. Otherwise the process execution will fail because the
    launched container won't be able to access the linked file.

Charliecloud & Docker Hub
-------------------------

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
-------------------

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
-----------------

Charliecloud advanced configuration settings are described in :ref:`config-charliecloud` section in the Nextflow
configuration page.


.. _container-docker:

Docker
======

`Docker <http://www.docker.io>`_ is the industry standard container runtime.

Prerequisites
-------------

You will need Docker installed on your execution environment e.g. your computer or a distributed cluster, depending
on where you want to run your pipeline.

If you are running Docker on Mac OSX make sure you are mounting your local ``/Users`` directory into the Docker VM as
explained in this excellent tutorial: `How to use Docker on OSX <http://viget.com/extend/how-to-use-docker-on-os-x-the-missing-guide>`_.

How it works
------------

You won't need to modify your Nextflow script in order to run it with Docker. Simply specify the Docker image from
where the containers are started by using the ``-with-docker`` command line option. For example::

  nextflow run <your script> -with-docker [docker image]

Every time your script launches a process execution, Nextflow will run it into a Docker container created by using the
specified image. In practice Nextflow will automatically wrap your processes and run them by executing the ``docker run``
command with the image you have provided.

.. note::
    A Docker image can contain any tool or piece of software you may need to carry out a process execution. Moreover, the
    container is run in such a way that the process result files are created in the host file system, thus
    it behaves in a completely transparent manner without requiring extra steps or affecting the flow in your pipeline.

If you want to avoid entering the Docker image as a command line parameter, you can define it in the Nextflow configuration
file. For example you can add the following lines in the ``nextflow.config`` file::

    process.container = 'nextflow/examples:latest'
    docker.enabled = true

In the above example replace ``nextflow/examples:latest`` with any Docker image of your choice.

Read the :ref:`config-page` page to learn more about the ``nextflow.config`` file and how to use it to configure
your pipeline execution.

.. warning::
    Nextflow automatically manages the file system mounts whenever a container is launched depending on the process
    input files. However, when a process input is a *symbolic link*, the linked file **must** be stored
    in the same folder where the symlink is located, or a sub-folder of it. Otherwise the process execution will fail because the
    launched container won't be able to access the linked file.

Multiple containers
-------------------

It is possible to specify a different Docker image for each process definition in your pipeline script. However, this
can't be done with ``-with-docker`` command line option, since it doesn't support process selectors in the ``nextflow.config`` file, and doesn't
check for the ``container`` process directive in the pipeline script file. Let's suppose you have two processes named ``foo``
and ``bar``. You can specify two different Docker images for them in the Nextflow script as shown below::

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

After running your pipeline, you can easily check what container image each process used with the following command::

    nextflow log last -f name,container

Advanced settings
-----------------

Docker advanced configuration settings are described in :ref:`config-docker` section in the Nextflow configuration page.


.. _container-podman:

Podman
======

`Podman <http://www.podman.io>`_ is a drop-in replacement for Docker that can run containers with or without root privileges.

.. note::
    This feature requires Nextflow version ``20.01.0`` or later.

.. warning::
    This feature is experimental. Using it in a production environment is not recommended.

Prerequisites
-------------

You will need Podman installed on your execution environment e.g. your computer or a distributed cluster, depending
on where you want to run your pipeline. Running in rootless mode requires appropriate OS configuration. Due to current
Podman limits using cpuset for cpus and memory such is only possible using sudo.

How it works
------------

You won't need to modify your Nextflow script in order to run it with Podman. Simply specify the Podman image from
where the containers are started by using the ``-with-podman`` command line option. For example::

  nextflow run <your script> -with-podman [OCI container image]

Every time your script launches a process execution, Nextflow will run it into a Podman container created by using the
specified image. In practice Nextflow will automatically wrap your processes and run them by executing the ``podman run``
command with the image you have provided.

.. note::
    An OCI container image can contain any tool or piece of software you may need to carry out a process execution. Moreover, the
    container is run in such a way that the process result files are created in the host file system, thus
    it behaves in a completely transparent manner without requiring extra steps or affecting the flow in your pipeline.

If you want to avoid entering the Podman image as a command line parameter, you can define it in the Nextflow configuration
file. For example you can add the following lines in the ``nextflow.config`` file::

    process.container = 'nextflow/examples:latest'
    podman.enabled = true

In the above example replace ``nextflow/examples:latest`` with any Podman image of your choice.

Read the :ref:`config-page` page to learn more about the ``nextflow.config`` file and how to use it to configure
your pipeline execution.

.. warning::
    Nextflow automatically manages the file system mounts whenever a container is launched depending on the process
    input files. However, when a process input is a *symbolic link*, the linked file **must** be stored
    in the same folder where the symlink is located, or a sub-folder of it. Otherwise the process execution will fail because the
    launched container won't be able to access the linked file.

Multiple containers
-------------------

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
-----------------

Podman advanced configuration settings are described in :ref:`config-podman` section in the Nextflow configuration page.


.. _container-sarus:

Sarus
=======

`Sarus <https://sarus.readthedocs.io>`_ is an alternative container runtime to
Docker. Sarus works by converting Docker images to a common format that can then be
distributed and launched on HPC systems. The user interface to Sarus enables a user to select an image
from `Docker Hub <https://hub.docker.com/>`_ and then submit jobs which run entirely within the container.

Prerequisites
-------------

You need Sarus installed in your execution environment,
i.e: your personal computer or a distributed cluster, depending
on where you want to run your pipeline.

.. note:: This feature requires Sarus version 1.5.1 (or later) and Nextflow 22.12.0-edge (or later).

Images
------

Sarus converts a docker image to squashfs layers which are distributed and launched in the cluster. For more information on
how to build Sarus images see the `official documentation <https://sarus.readthedocs.io/en/stable/user/user_guide.html#develop-the-docker-image>`_.

How it works
------------

The integration for Sarus, at this time, requires you to set up the following parameters in your config file::

  process.container = "dockerhub_user/image_name:image_tag"
  sarus.enabled = true

and it will always try to search the Docker Hub registry for the images.

.. note:: if you do not specify an image tag, the ``latest`` tag will be fetched by default.

Multiple containers
-------------------

It is possible to specify a different Sarus image for each process definition in your pipeline script. For example,
let's suppose you have two processes named ``foo`` and ``bar``. You can specify two different Sarus images
specifying them in the ``nextflow.config`` file as shown below::

    process {
        withName:foo {
            container = 'image_name_1'
        }
        withName:bar {
            container = 'image_name_2'
        }
    }
    sarus {
        enabled = true
    }

Read the :ref:`Process scope <config-process>` section to learn more about processes configuration.


.. _container-shifter:

Shifter
=======

`Shifter <https://docs.nersc.gov/programming/shifter/overview/>`_ is an alternative container runtime to
Docker. Shifter works by converting Docker images to a common format that can then be
distributed and launched on HPC systems. The user interface to Shifter enables a user to select an image
from `Docker Hub <https://hub.docker.com/>`_ and then submit jobs which run entirely within the container.

Prerequisites
-------------

You need Shifter and Shifter image gateway installed in your execution environment, i.e: your personal computer or the
entry node of a distributed cluster. In the case of the distributed cluster, you should have Shifter installed on
all of the compute nodes and the ``shifterimg`` command should also be available and Shifter properly setup to access the
Image gateway, for more information see the `official documentation <https://github.com/NERSC/shifter/tree/master/doc>`_.

.. note:: This feature requires Shifter version 18.03 (or later) and Nextflow 19.10.0 (or later).

Images
------

Shifter converts a docker image to squashfs layers which are distributed and launched in the cluster. For more information on
how to build Shifter images see the `official documentation <https://docs.nersc.gov/programming/shifter/how-to-use/#building-shifter-images>`_.

How it works
------------

The integration for Shifter, at this time, requires you to set up the following parameters in your config file::

  process.container = "dockerhub_user/image_name:image_tag"
  shifter.enabled = true

and it will always try to search the Docker Hub registry for the images.

.. note:: if you do not specify an image tag, the ``latest`` tag will be fetched by default.

Multiple containers
-------------------

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


.. _container-singularity:

Singularity
===========

`Singularity <http://singularity.lbl.gov/>`_ is an alternative container runtime to Docker. The main advantages
of Singularity are that it can be used without root privileges and it doesn't require a separate daemon process.
These, along with other features such as support for autofs mounts, makes Singularity
better suited to the requirements of HPC workloads. Singularity is able to use existing Docker images and can pull from Docker
registries.

Prerequisites
-------------

You will need Singularity installed on your execution environment e.g. your computer or a distributed cluster, depending
on where you want to run your pipeline.

Images
------

Singularity makes use of a container image file, which physically contains the container. Refer to the `Singularity
documentation <https://www.sylabs.io/docs/>`_ to learn how create Singularity images.

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
------------

The integration for Singularity follows the same execution model implemented for Docker. You won't need to modify your
Nextflow script in order to run it with Singularity. Simply specify the Singularity image
file from where the containers are started by using the ``-with-singularity`` command line option. For example::

  nextflow run <your script> -with-singularity [singularity image file]

Every time your script launches a process execution, Nextflow will run it into a Singularity container created by using the
specified image. In practice Nextflow will automatically wrap your processes and launch them by running the
``singularity exec`` command with the image you have provided.

.. note::
    A Singularity image can contain any tool or piece of software you may need to carry out a process execution.
    Moreover, the container is run in such a way that the process result files are created in the host file system, thus
    it behaves in a completely transparent manner without requiring extra steps or affecting the flow in your pipeline.

If you want to avoid entering the Singularity image as a command line parameter, you can define it in the Nextflow
configuration file. For example you can add the following lines in the ``nextflow.config`` file::

    process.container = '/path/to/singularity.img'
    singularity.enabled = true

In the above example replace ``/path/to/singularity.img`` with any Singularity image of your choice.

Read the :ref:`config-page` page to learn more about the ``nextflow.config`` file and how to use it to configure
your pipeline execution.

.. note::
    Unlike Docker, Nextflow does not automatically mount host paths in the container when using Singularity.
    It expects that the paths are configure and mounted system wide by the Singularity runtime. If your Singularity installation
    allows user defined bind points, read the :ref:`Singularity configuration <config-singularity>` section to learn
    how to enable Nextflow auto mounts.

.. warning::
    When a process input is a *symbolic link* file, make sure the linked file is stored in a host folder that is
    accessible from a bind path defined in your Singularity installation. Otherwise the process execution will fail
    because the launched container won't be able to access the linked file.

Multiple containers
-------------------

It is possible to specify a different Singularity image for each process definition in your pipeline script. For example,
let's suppose you have two processes named ``foo`` and ``bar``. You can specify two different Singularity images
specifying them in the ``nextflow.config`` file as shown below::

    process {
        withName:foo {
            container = 'image_name_1'
        }
        withName:bar {
            container = 'image_name_2'
        }
    }
    singularity {
        enabled = true
    }


Read the :ref:`Process scope <config-process>` section to learn more about processes configuration.

Singularity & Docker Hub
------------------------

Nextflow is able to transparently pull remote container images stored in the `Singularity-Hub <https://singularity-hub.org/>`_,
`Singularity Library <https://cloud.sylabs.io/library/>`_, or any Docker compatible registry.

.. note:: This feature requires Singularity 2.3.x or higher

By default when a container name is specified, Nextflow checks if an image file with that name exists in the local file
system. If that image file exists, it's used to execute the container. If a matching file does not exist,
Nextflow automatically tries to pull an image with the specified name from the Docker Hub.

If you want Nextflow to check only for local file images, prefix the container name with the ``file://`` pseudo-protocol.
For example::

    process.container = 'file:///path/to/singularity.img'
    singularity.enabled = true

.. warning::
    Use three ``/`` slashes to specify an **absolute** file path, otherwise the path will be interpreted as
    relative to the workflow launch directory.

To pull images from the Singularity Hub or a third party Docker registry simply prefix the image name
with the ``shub://``, ``docker://`` or ``docker-daemon://`` pseudo-protocol as required by the Singularity tool. For example::

    process.container = 'docker://quay.io/biocontainers/multiqc:1.3--py35_2'
    singularity.enabled = true

You do not need to specify ``docker://`` to pull from a Docker repository.
Nextflow will automatically prepend it to your image name when Singularity is enabled.
Additionally, the Docker engine will not work with containers specified as ``docker://``.

Nextflow version 18.10 introduced support for the `Singularity Library <https://cloud.sylabs.io/library/>`_ repository. This feature also requires Singularity 3.0::

    process.container = 'library://library/default/alpine:3.8'

The pseudo-protocol allows you to import Singularity using a local Docker installation instead of downloading
the container image from the Docker registry. It requires Nextflow 19.04.0 or later and Singularity 3.0.3 or later.

.. note::
    This feature requires the ``singularity`` tool to be installed
    where the workflow execution is launched (as opposed to the compute nodes).

Nextflow caches those images in the ``singularity`` directory in the pipeline work directory by default. However it is
suggested to provide a centralised cache directory by using either the ``NXF_SINGULARITY_CACHEDIR`` environment variable
or the ``singularity.cacheDir`` setting in the Nextflow config file.

As of version ``21.09.0-edge``, when looking for a Singularity image file, Nextflow first checks the *library* directory,
and if the image file is not found, the *cache* directory is used as usual. The library directory can be defined either using
the ``NXF_SINGULARITY_LIBRARYDIR`` environment variable or the ``singularity.libraryDir`` configuration setting (the
latter overrides the former).

.. warning::
    When using a compute cluster, the Singularity cache directory must reside in a shared filesystem accessible to all compute nodes.

.. danger::
    When pulling Docker images, Singularity may be unable to determine the container size if the image was
    stored using an old Docker format, resulting in a pipeline execution error. See the Singularity documentation for details.

Advanced settings
-----------------

Singularity advanced configuration settings are described in :ref:`config-singularity` section in the Nextflow
configuration page.

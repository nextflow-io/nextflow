.. _spack-page:

******************
Spack environments
******************


`Spack <https://spack.io/>`_ Spack is a package manager for supercomputers, Linux, and macOS.
It makes installing scientific software easy. Spack is not tied to a particular language;
you can build a software stack in Python or R, link to libraries written in C, C++, or Fortran,
and easily swap compilers or target specific microarchitectures.

Nextflow has built-in support for Spack that allows the configuration of workflow dependencies
using Spack recipes and environment files.

This allows Nextflow applications to build packages from source on the compute infrastructure in use,
whilst taking advantage of the configuration flexibility provided by Nextflow.
At shared compute facilities where Spack has been configured by the administrators,
this may result in optimised builds without user intervention. With appropriate options,
this also permits end users to customise binary optimisations by themselves.


Prerequisites
-------------

This feature requires the `Spack <https://spack.io>`_ package manager to be installed on your system.


How it works
------------

Nextflow  automatically creates and activates the Spack environment(s) given the dependencies
specified by each process.

Dependencies are specified by using the :ref:`process-spack` directive, providing either
the names of the required Spack packages, the path of a Spack environment yaml file or
the path of an existing Spack environment directory.

.. note:: Spack always installs the software packages in its own directories, regardless of the Nextflow specifications.
  The Spack environment created by Nextflow only contains symbolic links pointing to the appropriate package locations,
  and therefore it is relatively small in size.

You can specify the directory where the Spack environment is stored using the ``spack.cacheDir``
configuration property (see the :ref:`configuration page <config-spack>` for details).
When using a computing cluster, make sure to use a shared file system path
accessible from all compute nodes.

.. warning:: The Spack environment feature is not supported by executors that use
  remote object storage as a work directory e.g. AWS Batch.

Enabling Spack environment
==========================

The use of Spack recipes specified using the :ref:`process-spack`
directive needs to be enabled explicitly by setting the option shown below in the pipeline
configuration file (i.e. ``nextflow.config``)::

    spack.enabled = true


Alternatively it can be specified by setting the variable ``NXF_SPACK_ENABLED=true`` in your environment
or by using the ``-with-spack true`` command line option.


Use Spack package names
=======================

Spack package names can specified using the ``spack`` directive. Multiple package names can be specified
by separating them with a blank space.
For example::

  process foo {
    spack 'bwa samtools py-multiqc'

    '''
    your_command --here
    '''
  }

Using the above definition a Spack environment that includes BWA, Samtools and MultiQC tools is created and
activated when the process is executed.

The usual Spack package syntax and naming conventions can be used. The version of a package can be
specified after the package name as shown here ``bwa@0.7.15``.

Optimisation for the local CPU microarchitetcure can be requested by adding the option ``target=<LOCAL ARCH>``
after a package name. For instance, if the compute infrastructure uses AMD Zen3 microprocessors,
use the following to optimise all packages for it: ``bwa target=zen3 samtools target=zen3 py-multiqc target=zen3``.

Read the Spack documentation for more details about `package specifications <https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies>`_.


Use Spack environment files
===========================

Spack environments can also be defined using one or more Spack environment files. This is a file that
lists the required packages and channels structured using the YAML format. For example::

    spack:
      specs:
      - star@2.5.4a
      - bwa@0.7.15
    
      view: true
      concretizer:
        unify: true

Here, the ``view`` and ``concretizer`` options are sensible Spack defaults for environments.

There are concise ways to specify the target microarchitecture (and eventually other options) within a Spack environment file.
For instance, the following environment file specifies build optimisation for an AMD Zen3 target microprocessor::

    spack:
      packages:
        all:
          target: [zen3]
      specs:
      - star@2.5.4a
      - bwa@0.7.15
    
      view: true
      concretizer:
        unify: true

Read the Spack documentation for more details about how to create `environment files <https://spack.readthedocs.io/en/latest/environments.html>`_.

The path of an environment file can be specified using the ``spack`` directive::

  process foo {
    spack '/some/path/my-env.yaml'

    '''
    your_command --here
    '''
  }

.. warning:: The environment file name **must** have a ``.yaml`` extension or else it won't be properly recognised.


Use existing Spack environments
===============================

If you already have a local Spack environment, you can use it in your workflow specifying the
installation directory of such environment by using the ``spack`` directive::

  process foo {
    spack '/path/to/an/existing/env/directory'

    '''
    your_command --here
    '''
  }


Best practices
--------------

Building Spack packages for Nextflow pipelines
==============================================

Spack builds most software package from their source codes, and it does this for a request package
and for all its required dependencies. As a result, Spack builds can last for long, even several hours.
This can represent an inconvenience, in that it can significantly lengthen the duration of Nextflow processes.
Here we briefly discuss two strategies to mitigate this aspect, and render the usage of Spack more effective.

1.  Use a Spack yaml file, and pre-build the environment outside of Nextflow, prior to running the pipeline.
    Building packages outside of the Nextflow pipeline will work since Spack always installs packages in its own directories,
    and only creates symbolic links in the environment. This sequence of commands will do the trick in most cases::

    spack env create myenv /path/to/spack.yaml
    spack env activate myenv
    spack concretize -f
    spack install -y
    spack env deactivate

2.  Use the Nextflow stub functionality prior to running the pipeline for production.
    Nextflow will run the stub pipeline, skipping process executions but still setting up the required software packages.
    This option is useful if it is not possible to write a Spack yaml file for the environment.
    The stub functionality is described in the :ref:`Stub <process-stub>` section of the Processes page.


Configuration file
==================

When a ``spack`` directive is used in any ``process`` definition within the workflow script, Spack tool is required for
the workflow execution.

Specifying the Spack environments in a separate configuration :ref:`profile <config-profiles>` is therefore
recommended to allow the execution via a command line option and to enhance the workflow portability. For example::

  profiles {
    spack {
      process.spack = 'samtools'
    }

    docker {
      process.container = 'biocontainers/samtools'
      docker.enabled = true
    }
  }

The above configuration snippet allows the execution either with Spack or Docker specifying ``-profile spack`` or
``-profile docker`` when running the workflow script.


Advanced settings
-----------------

Spack advanced configuration settings are described in the :ref:`Spack <config-spack>` section on the Nextflow configuration page.

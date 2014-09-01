.. _docker-page:

*******************
Docker integration
*******************

Nextflow integration with `Docker container <http://www.docker.io>`_ technology allows you to write truly self-contained
and reproducible computational pipelines.

By using this feature any process in a pipeline script can be seamlessly executed through a Docker container. This may
be extremely useful to package the binaries dependencies in a pipeline script into a standard and portable for format
that can be executed on any platform supporting Docker.


Prerequisites
==============

You will need Docker installed in your target execution environment e.g. your computer or a distributed cluster depending
where you want to run your scripts.



How it works
=============

In order your to execute your pipeline script with Docker you won't need to modify it. In it's simplest form add
the ``-with-docker`` command line option by specifying the name of a Docker image you want to use to run the process
in your pipeline. For example::

  nextflow run <your script> -with-docker <docker image>


Each time a process execution is spawned, Nextflow will run it into a Docker contained created by using the specified image.
The image can contain any tool or piece of software you may need to carry out the process execution. The container
is run is a such a way that the process result files are created in the hosting execution environment, in other words
in behaves in complete transparent manner.



Configuration
=============

Add in the `nextflow.config` file for your pipeline a

Available Docker configuration options

================== ================
Name                Description
================== ================
image               Name of a Docker image used to start
enabled             By default execution through is disabled. You need to enable it by turning this flag to ``true`` or by using the `-with-docker` command line options  (default: ``false``)
sudo                (default: ``false``)
temp                It allows to mount
remove              Clean-up the container after the execution (default: ``true``). See also http://docs.docker.com/reference/run/#clean-up-rm
runOptions          Provide Docker run o
================== ================












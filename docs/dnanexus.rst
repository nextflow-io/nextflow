.. _dnanexus-page:

****************
DNAnexus cloud
****************

DNAnexus is a cloud platform specifically designed for genomics applications. Nextflow provides experimental support
for the DNAnexus platform allowing you to execute your pipeline in the cloud.


Prerequisites
==============

You will need a DNAnexus account and the `DNAnexus SDK toolkit <https://wiki.dnanexus.com/Downloads#DNAnexus-Platform-SDK>`_
installed in your computer.

Moreover, since extra runtime packages are required by Nextflow to enable it to run in the cloud, you will need to download
the full Nextflow source tree in your computer. You can do this by using the following `git` command::

   git clone https://github.com/nextflow-io/nextflow.git


App packing
============

In order to execute a Nextflow pipeline in the DNAnexus cloud you need to package it as a DNAnexus app (or applet).
Thus you will need to use the ``dx-app-wizard`` command to create the application structure and copy the pipeline script,
along with the Nextflow runtime, in the app project folder..

.. tip:: Read more about DNAnexus applications at `this link <https://wiki.dnanexus.com/Developer-Tutorials/Intro-to-Building-Apps>`_.

Nextflow provides a shortcut that helps you create the DNAnexus application skeleton pre-configured to work with it.

Simply change to the Nextflow source project folder, that you have cloned with `git`, and type the command::

   ./gradlew dnanexus


This will:

#. Build the Nextflow runtime with the required dependencies
#. Create the DNAnexus app skeleton in the folder named ``dx-project``
#. Create the ``dxapp.json`` application manifest
#. Create the ``dxapp.sh`` pipeline launcher
#. Copy the Nextflow runtime to the ``dx-project/resources/usr/bin`` folder
#. Copy some example scripts to the ``dx-project/resources`` folder


Then, you will need to add your pipeline script into the ``resources`` folder (adding any other eventually needed
binary dependencies) and modify the ``dxapp.json`` file providing your application name.

When you have finished, package and upload the app by using the DNAnexus `build` command, as shown below::

    dx build -f


Execution
==========

When ready, you can launch it like any other DNAnexus app, by using the `run` command, for example::

    dx run nextflow

.. tip:: Replace ``nextflow`` with your app name in the above example.

By running this example and confirming the default options, a simple `Hello world` script will run in the cloud, producing
the output shown below. ::


    Job Log
    -------
    Watching job job-BFY5Pz00kF9B7FY891602V7g. Press Ctrl+C to stop.
    * Nextflow pipelines framework (nextflow:main) (running) job-BFY5Pz00kF9B7FY891602V7g
      pditommaso 2014-01-27 13:58:36 (running for 0:00:05)
    2014-01-27 14:02:18 Nextflow pipelines framework INFO Logging initialized (priority)
    2014-01-27 14:02:22 Nextflow pipelines framework STDOUT Installing apt packages openjdk-7-jre-headless, dx-toolkit
    2014-01-27 14:02:47 Nextflow pipelines framework STDOUT >>> Unpacking resources.tar.gz to /
    2014-01-27 14:02:47 Nextflow pipelines framework STDOUT *** Downloading bundled file resources.tar.gz
    2014-01-27 14:02:48 Nextflow pipelines framework STDOUT bash running (job ID job-BFY5Pz00kF9B7FY891602V7g)
    2014-01-27 14:02:51 Nextflow pipelines framework STDOUT N E X T F L O W  ~  version 0.6.0
    2014-01-27 14:02:54 Nextflow pipelines framework STDOUT [63/9a33ce] Running process > sayhello (1)
    2014-01-27 14:04:09 Nextflow pipelines framework STDOUT Hello world!
    2014-01-27 14:04:23 Nextflow pipelines framework STDOUT nextflow exitstatus > 0
    2014-01-27 14:04:24 Nextflow pipelines framework STDOUT hello.log > file-BFY5Xf00kF97JFXyzJ3Xg637
    * Nextflow pipelines framework (nextflow:main) (done) job-BFY5Pz00kF9B7FY891602V7g
      pditommaso 2014-01-27 13:58:36 (runtime 0:02:11)
      Output: -


How it works
=============

The DNAnexus app is executed by the BASH script called ``dxapp.sh`` which defines two functions: `main` and `process`.

The `main` function is the one that launches the Nextflow pipeline script execution in the cloud, specifying the required
configuration parameters.


Nextflow will be able to manage the processes scheduling and execution in the cloud because it is launched by specifying
the ``dnanexus`` executor.

When a new process needs to be executed, Nextflow will use the DNAnexus API to run a new job, which in turn calls
the function `process`, defined in the ``dxapp.sh`` script, specifying the target command to run.

Read more about Nextflow's `dnanexus` executor configuration properties :ref:`here <dnanexus-executor>`.








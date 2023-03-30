.. _flux-page:

**************
Flux Framework
**************

`Flux Framework <https://flux-framework.org/>`_ "flux" is modern resource manager that can span the space between
cloud and HPC. If your center does not provide Flux for you, you can `build Flux on your own <https://flux-framework.readthedocs.io/en/latest/quickstart.html#building-the-code>`_
and launch it as a job under your resource manager of choice (e.g., SLURM or on the cloud.)


Tutorial
========

In the `docker/flux <https://github.com/nextflow-io/nextflow/tree/master/docker/flux>`_ directory we provide 
a `Dockerfile for interacting with Flux <https://github.com/nextflow-io/nextflow/tree/master/docker/flux/.devcontainer/Dockerfile>`_

along with a `VSCode Developer Container <https://code.visualstudio.com/docs/devcontainers/containers>`_
environment that you can put at the root of the project to be provided with a Flux agent and the dependencies
needed to build Nextflow. There are two ways to use this:

 - Build a container from scratch and bind your code to it (e.g., for development or testing)
 - Use VSCode and DevContainers to create a more seamless environment.

Both strategies are described below. For this tutorial, you will generally want 
to prepare a workflow to use the ``flux`` executor, create an environment with Flux, start a Flux instance, and
interact with it.

Prepare Nextflow Workflow 
-------------------------

For your workflow, to use Flux you'll want to specify it in your config. Here is an example 
``nextflow.config``


.. code-block:: console

    // This is an example config for running with flux.

    manifest {
      mainScript = 'demo.nf'
      homePage = 'https://github.com/nextflow-io/nextflow/tree/master/docker/flux'
      description = 'Demo using Nextflow with Flux'
    } 

    process {
      executor = 'flux'
    }

Note that the executor is ``flux``. For additional settings for Flux, see :ref:`flux-executor`.
Here is an a demo workflow ``demo.nf`` of a job we want to run!

.. code-block:: console

    breakfast = Channel.from '🥞️', '🥑️', '🥧️', '🍵️', '🍞️'
 
    workflow {
      haveMeal(breakfast)
    }
    
    process haveMeal {
      debug true
      input:
      val food
      """
      printf '$food for breakfast!'
      """
    }

We will be using these files to run our test workflow. Next, assuming you don't have one handy,
let's set up an environment with Flux.

Container Environment
---------------------

You can choose just to build the Docker image from the root of the repository:

.. code-block:: console
    
    $ docker build -f docker/flux/.devcontainer/Dockerfile --platform linux/amd64 -o type=docker -t nextflow-flux .


And then shell into the container for a development environment. You'd need to bind
the present working directory to ``/code`` to see your local changes in the container:

.. code-block:: console

    $ docker run -it -v $PWD:/code nextflow-flux 


You can also move the .devcontainer directory to the root of your repository, and just
open it in VSCode.

.. code-block:: console

    $ cp -R docker/flux/.devcontainer .devcontainer


Then open in vscode, and select to "Re-open in container"

.. code-block:: console

    $ code .

Then you should be able to open a terminal (Terminal -> New Terminal)
to interact with the command line. Try running `make`` again!
Whichever of these two approaches you take, you should be in a container 
environment with ``flux`` on the path.

Start a Flux Instance 
---------------------

Once in your container, you can start an interactive flux instance (from which you can submit jobs on
the command line to test with Nextflow) as follows:

.. code-block:: console

    $ flux start --test-size=4


Getting Familiar with Flux
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

    This step is optional!

Here is an example of submitting a job and getting the log for it.

.. code-block:: console

    $ flux mini submit echo "HELLO MOTO"
    ƒEzWqspb


And then getting the log for it:


.. code-block:: console

    $ flux job attach ƒEzWqspb
    HELLO MOTO


Try submitting a longer job:

.. code-block:: console

    $ flux mini submit sleep 60


And then seeing it in the jobs listing.

.. code-block:: console

    $ flux jobs
           JOBID USER     NAME       ST NTASKS NNODES     TIME INFO
       ƒ4tkMUAAT root     sleep       R      1      1   2.546s ab6634a491bb


We are going to be issuing similar submit commands via Nextflow!

Submitting with Nextflow
------------------------

Prepare your ``nextflow.config`` and ``demo.nf`` in the same directory.

.. code-block:: console

    $ ls .
    demo.nf    nextflow.config 


If you've installed Nextflow already with the flux executor, you are good to go! If you are working 
with development code and need to build nextflow: 

.. code-block:: console

    $ make assemble

Make sure ``nextflow`` is on your PATH (here we are in the root of Nextflow):

.. code-block:: console

    $ export PATH=$PWD:$PATH
    $ which nextflow
    /workspaces/nextflow/nextflow

Then cd into the directory with your config and demo file:

.. code-block:: console

    $ cd docker/flux


And then run the workflow with flux!

.. code-block:: console

    $ nextflow -c nextflow.config run demo.nf

    # nextflow -c nextflow.config run demo.nf 
    N E X T F L O W  ~  version 22.10.0
    Launching `demo.nf` [clever_blackwell] DSL2 - revision: f8cda838cb
    executor >  flux (5)
    [4c/f162db] process > haveMeal (3) [100%] 5 of 5 ✔
    🥞️ for breakfast!
    🍞️ for breakfast!
    🍵️ for breakfast!
    🥑️ for breakfast!
    🥧️ for breakfast!

And that's it! You've just run a workflow using nextflow and Flux. 

(flux-page)=

# Flux Framework

:::{versionadded} 22.11.0-edge
:::

The [Flux Framework](https://flux-framework.org/) is a modern resource manager that can span the space between cloud and HPC. If your center does not provide Flux for you, you can [build Flux on your own](https://flux-framework.readthedocs.io/en/latest/quickstart.html#building-the-code) and launch it as a job with your resource manager of choice (e.g. SLURM or a cloud provider).

## Tutorial

In the [`docker/flux`](https://github.com/nextflow-io/nextflow/tree/master/docker/flux) directory we provide a [Dockerfile for interacting with Flux](https://github.com/nextflow-io/nextflow/tree/master/docker/flux/.devcontainer/Dockerfile) along with a [VSCode Developer Container](https://code.visualstudio.com/docs/devcontainers/containers) environment that you can put at the root of the project to be provided with a Flux agent and the dependencies needed to build Nextflow. There are two ways to use this:

- Build a container from scratch and bind your code to it (e.g. for development or testing)
- Use VSCode and DevContainers to create a more seamless environment

Both strategies are described below. For this tutorial, you will generally want to prepare a pipeline to use the `flux` executor, create an environment with Flux, start a Flux instance, and interact with it.

### Prepare your pipeline

To run your pipeline with Flux, you'll want to specify it in your config. Here is an example `nextflow.config`:

```groovy
manifest {
    mainScript = 'demo.nf'
    homePage = 'https://github.com/nextflow-io/nextflow/tree/master/docker/flux'
    description = 'Demo using Nextflow with Flux'
}

process {
    executor = 'flux'
}
```

For additional Flux settings, see the {ref}`flux-executor` section.

Here is an example pipeline that we will use:

```groovy
workflow {
    breakfast = Channel.from '🥞️', '🥑️', '🥧️', '🍵️', '🍞️'
    haveMeal(breakfast)
}

process haveMeal {
    debug true
    input:
    val food
    script:
    """
    printf '$food for breakfast!'
    """
}
```

### Container Environment

You can either build the Docker image from the root of the Nextflow repository:

```console
$ docker build -f docker/flux/.devcontainer/Dockerfile --platform linux/amd64 -o type=docker -t nextflow-flux .
```

And then shell into the container for a development environment. You'll need to bind the present working directory to `/code` to see your local changes in the container:

```console
$ docker run -it -v $PWD:/code nextflow-flux
```

You can also move the `.devcontainer` directory to the root of your repository, and open it in VSCode:

```console
$ cp -R docker/flux/.devcontainer .devcontainer
```

Then open in VSCode, and select **Re-open in container**:

```console
$ code .
```

Then you should be able to open a terminal (**Terminal** -> **New Terminal**) to interact with the command line. Try running `make` again! Whichever of these two approaches you take, you should be in a container environment with the `flux` command available.

### Start a Flux Instance

Once in your container, you can start an interactive Flux instance (from which you can submit jobs on the command line to test with Nextflow) as follows:

```console
$ flux start --test-size=4
```

#### Getting Familiar with Flux

:::{note}
This step is optional!
:::

Here is an example of submitting a job and getting the log for it.

First submit the job:

```console
$ flux mini submit echo "HELLO MOTO"
ƒEzWqspb
```

Then get the log for it:

```console
$ flux job attach ƒEzWqspb
HELLO MOTO
```

Try submitting a longer job:

```console
$ flux mini submit sleep 60
```

And then seeing it in the jobs listing.

```console
$ flux jobs
       JOBID USER     NAME       ST NTASKS NNODES     TIME INFO
   ƒ4tkMUAAT root     sleep       R      1      1   2.546s ab6634a491bb
```

### Submitting with Nextflow

Prepare your `nextflow.config` and `demo.nf` in the same directory.

```console
$ ls .
demo.nf    nextflow.config
```

If you've installed Nextflow already, you are good to go! If you are working with development code and need to build Nextflow:

```console
$ make assemble
```

Make sure `nextflow` is on your PATH (here we are in the root of the Nextflow repository):

```console
$ export PATH=$PWD:$PATH
$ which nextflow
/workspaces/nextflow/nextflow
```

Then change to the directory with your config and demo file:

```console
$ cd docker/flux
```

And then run the pipeline with Flux!

```console
$ nextflow -c nextflow.config run demo.nf

N E X T F L O W  ~  version 22.10.0
Launching `demo.nf` [clever_blackwell] DSL2 - revision: f8cda838cb
executor >  flux (5)
[4c/f162db] process > haveMeal (3) [100%] 5 of 5 ✔
🥞️ for breakfast!
🍞️ for breakfast!
🍵️ for breakfast!
🥑️ for breakfast!
🥧️ for breakfast!
```

And that's it! You've just run a pipeline using nextflow and Flux.

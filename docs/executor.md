(executor-page)=

# Executors

In the Nextflow framework architecture, the *executor* is the component that determines the system where a pipeline process is run and supervises its execution.

The executor provides an abstraction between the pipeline processes and the underlying execution system. This allows you to write the pipeline functional logic independently from the actual processing platform.

In other words, you can write your pipeline script once and have it running on your computer, a cluster resource manager, or the cloud — simply change the executor definition in the Nextflow configuration file.

(awsbatch-executor)=

## AWS Batch

Nextflow supports the [AWS Batch](https://aws.amazon.com/batch/) service that allows job submission in the cloud without having to spin out and manage a cluster of virtual machines. AWS Batch uses Docker containers to run tasks, which greatly simplifies pipeline deployment.

The pipeline processes must specify the Docker image to use by defining the `container` directive, either in the pipeline script or the `nextflow.config` file.

To enable this executor, set `process.executor = 'awsbatch'` in the `nextflow.config` file.

The pipeline can be launched either in a local computer, or an EC2 instance. EC2 is suggested for heavy or long-running workloads. Additionally, an S3 bucket must be used as the pipeline work directory.

Resource requests and other job characteristics can be controlled via the following process directives:

- {ref}`process-accelerator`
- {ref}`process-container`
- {ref}`process-containerOptions`
- {ref}`process-cpus`
- {ref}`process-memory`
- {ref}`process-queue`
- {ref}`process-resourcelabels`
- {ref}`process-time`

See the {ref}`AWS Batch<aws-batch>` page for further configuration details.

(azurebatch-executor)=

## Azure Batch

:::{versionadded} 21.04.0
:::

Nextflow supports the [Azure Batch](https://azure.microsoft.com/en-us/services/batch/) service that allows job submission in the cloud without having to spin out and manage a cluster of virtual machines. Azure Batch uses Docker containers to run tasks, which greatly simplifies pipeline deployment.

The pipeline processes must specify the Docker image to use by defining the `container` directive, either in the pipeline script or the `nextflow.config` file.

To enable this executor, set `process.executor = 'azurebatch'` in the `nextflow.config` file.

The pipeline can be launched either in a local computer, or a cloud virtual machine. The cloud VM is suggested for heavy or long-running workloads. Additionally, an Azure Blob storage container must be used as the pipeline work directory.

Resource requests and other job characteristics can be controlled via the following process directives:

- {ref}`process-container`
- {ref}`process-containerOptions`
- {ref}`process-cpus`
- {ref}`process-machineType`
- {ref}`process-memory`
- {ref}`process-queue`
- {ref}`process-resourcelabels`
- {ref}`process-time`

See the {ref}`Azure Batch <azure-batch>` page for further configuration details.

(bridge-executor)=

## Bridge

:::{versionadded} 22.09.1-edge
:::

[Bridge](https://github.com/cea-hpc/bridge) is an abstraction layer to ease batch system and resource manager usage in heterogeneous HPC environments.

It is open source software that can be installed on top of existing classical job schedulers such as Slurm, LSF, or other schedulers. Bridge allows you to submit jobs, get information on running jobs, stop jobs, get information on the cluster system, etc.

For more details on how to install the Bridge system, see the [documentation](https://github.com/cea-hpc/bridge).

To enable the Bridge executor, set `process.executor = 'bridge'` in the `nextflow.config` file.

Resource requests and other job characteristics can be controlled via the following process directives:

- {ref}`process-clusterOptions`
- {ref}`process-cpus`
- {ref}`process-memory`
- {ref}`process-queue`
- {ref}`process-time`

(flux-executor)=

## Flux Executor

:::{versionadded} 22.11.0-edge
:::

The `flux` executor allows you to run your pipeline script using the [Flux Framework](https://flux-framework.org).

Nextflow manages each process as a separate job that is submitted to the cluster by using the `flux mini submit` command.

To enable the Flux executor, set `process.executor = 'flux'` in the `nextflow.config` file.

Resource requests and other job characteristics can be controlled via the following process directives:

- {ref}`process-clusterOptions`
- {ref}`process-cpus`
- {ref}`process-queue`
- {ref}`process-time`

:::{note}
Flux does not support the `memory` directive.
:::

:::{note}
By default, Flux will send all output to the `.command.log` file. To send this output to stdout and stderr instead, set `flux.terminalOutput = true` in your config file.
:::

(ga4ghtes-executor)=

## GA4GH TES

:::{warning} *Experimental: may change in a future release.*
:::

The [Task Execution Schema](https://github.com/ga4gh/task-execution-schemas) (TES) project by the [GA4GH](https://www.ga4gh.org) standardization initiative is an effort to define a standardized schema and API for describing batch execution tasks in a portable manner.

Nextflow supports the TES API via the `tes` executor, which allows the submission of workflow tasks to a remote execution backend exposing a TES API endpoint.

To use this feature, define the following variables in the workflow launching environment:

```bash
export NXF_MODE=ga4gh
export NXF_EXECUTOR=tes
export NXF_EXECUTOR_TES_ENDPOINT='http://back.end.com'
```

It is important that the endpoint is specified without the trailing slash; otherwise, the resulting URLs will not be normalized and the requests to TES will fail.

You will then be able to run your workflow over TES using the usual Nextflow command line. Be sure to specify the Docker image to use, i.e.:

```bash
nextflow run rnaseq-nf -with-docker alpine
```

:::{note}
If the variable `NXF_EXECUTOR_TES_ENDPOINT` is omitted, the default endpoint is `http://localhost:8000`.
:::

:::{tip}
You can use a local [Funnel](https://ohsu-comp-bio.github.io/funnel/) server using the following launch command line:

```bash
./funnel server --Server.HTTPPort 8000 --LocalStorage.AllowedDirs $HOME run
```

(tested with version 0.8.0 on macOS)
:::

:::{warning}
Make sure the TES backend can access the Nextflow work directory when data is exchanged using a local or shared file system.
:::

### Known Limitations

- Automatic deployment of workflow scripts in the `bin` folder is not supported.

  :::{versionchanged} 23.07.0-edge
  Automatic upload of the `bin` directory is now supported.
  :::

- Process output directories are not supported. For details see [#76](https://github.com/ga4gh/task-execution-schemas/issues/76).

- Glob patterns in process output declarations are not supported. For details see [#77](https://github.com/ga4gh/task-execution-schemas/issues/77).

(google-batch-executor)=

## Google Cloud Batch

:::{versionadded} 22.07.1-edge
:::

[Google Cloud Batch](https://cloud.google.com/batch) is a managed computing service that allows the execution of containerized workloads in the Google Cloud Platform infrastructure.

Nextflow provides built-in support for the Cloud Batch API, which allows the seamless deployment of a Nextflow pipeline in the cloud, offloading the process executions as pipelines.

The pipeline processes must specify the Docker image to use by defining the `container` directive, either in the pipeline script or the `nextflow.config` file. Additionally, the pipeline work directory must be located in a Google Storage bucket.

To enable this executor, set `process.executor = 'google-batch'` in the `nextflow.config` file.

Resource requests and other job characteristics can be controlled via the following process directives:

- {ref}`process-accelerator`
- {ref}`process-container`
- {ref}`process-containerOptions`
- {ref}`process-cpus`
- {ref}`process-disk`
- {ref}`process-machineType`
- {ref}`process-memory`
- {ref}`process-resourcelabels`
- {ref}`process-time`

See the {ref}`Google Cloud Batch <google-batch>` page for further configuration details.

(google-lifesciences-executor)=

## Google Life Sciences

:::{versionadded} 20.01.0
:::

[Google Cloud Life Sciences](https://cloud.google.com/life-sciences) is a managed computing service that allows the execution of containerized workloads in the Google Cloud Platform infrastructure.

Nextflow provides built-in support for the Life Sciences API, which allows the seamless deployment of a Nextflow pipeline in the cloud, offloading the process executions as pipelines.

The pipeline processes must specify the Docker image to use by defining the `container` directive, either in the pipeline script or the `nextflow.config` file. Additionally, the pipeline work directory must be located in a Google Storage bucket.

To enable this executor, set `process.executor = 'google-lifesciences'` in the `nextflow.config` file.

Resource requests and other job characteristics can be controlled via the following process directives:

- {ref}`process-accelerator`
- {ref}`process-cpus`
- {ref}`process-disk`
- {ref}`process-machineType`
- {ref}`process-memory`
- {ref}`process-resourcelabels`
- {ref}`process-time`

See the {ref}`Google Life Sciences <google-lifesciences>` page for further configuration details.

(htcondor-executor)=

## HTCondor

:::{warning} *Experimental: may change in a future release.*
:::

The `condor` executor allows you to run your pipeline script by using the [HTCondor](https://research.cs.wisc.edu/htcondor/) resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster using the `condor_submit` command.

The pipeline must be launched from a node where the `condor_submit` command is available, which is typically the cluster login node.

:::{note}
The HTCondor executor for Nextflow does not currently support HTCondor's ability to transfer input/output data to the corresponding job's compute node. Therefore, the data must be made accessible to the compute nodes through a shared file system directory from where the Nextflow workflow is executed (or specified via the `-w` option).
:::

To enable the HTCondor executor, set `process.executor = 'condor'` in the `nextflow.config` file.

Resource requests and other job characteristics can be controlled via the following process directives:

- {ref}`process-clusterOptions`
- {ref}`process-cpus`
- {ref}`process-disk`
- {ref}`process-memory`
- {ref}`process-time`

(hyperqueue-executor)=

## HyperQueue

:::{versionadded} 22.05.0-edge
:::

:::{warning} *Experimental: may change in a future release.*
:::

The `hyperqueue` executor allows you to run your pipeline script by using the [HyperQueue](https://github.com/It4innovations/hyperqueue) job scheduler.

Nextflow manages each process as a separate job that is submitted to the cluster using the `hq` command line tool.

The pipeline must be launched from a node where the `hq` command is available, which is typically the cluster login node.

To enable the HTCondor executor, set `process.executor = 'hyperqueue'` in the `nextflow.config` file.

Resource requests and other job characteristics can be controlled via the following process directives:

- {ref}`process-accelerator`
- {ref}`process-clusterOptions`
- {ref}`process-cpus`
- {ref}`process-memory`
- {ref}`process-time`

(ignite-executor)=

## Ignite

:::{warning}
This feature is no longer maintained.
:::

:::{versionchanged} 22.01.0-edge
The `ignite` executor must be enabled via the `nf-ignite` plugin.
:::

The `ignite` executor allows you to run a pipeline on an [Apache Ignite](https://ignite.apache.org/) cluster.

To enable this executor, set `process.executor = 'ignite'` in the `nextflow.config` file.

Resource requests and other job characteristics can be controlled via the following process directives:

- {ref}`process-cpus`
- {ref}`process-disk`
- {ref}`process-memory`

See the {ref}`ignite-page` page to learn how to configure Nextflow to deploy and run an Ignite cluster in your infrastructure.

(k8s-executor)=

## Kubernetes

The `k8s` executor allows you to run a pipeline on a [Kubernetes](http://kubernetes.io/) cluster.

Resource requests and other job characteristics can be controlled via the following process directives:

- {ref}`process-accelerator`
- {ref}`process-cpus`
- {ref}`process-disk`
- {ref}`process-memory`
- {ref}`process-pod`
- {ref}`process-resourcelabels`
- {ref}`process-time`

See the {ref}`Kubernetes <k8s-page>` page to learn how to set up a Kubernetes cluster to run Nextflow pipelines.

(local-executor)=

## Local

The `local` executor is used by default. It runs the pipeline processes on the computer where Nextflow is launched. The processes are parallelised by spawning multiple threads, taking advantage of the multi-core architecture of the CPU.

The `local` executor is useful for developing and testing a pipeline script on your computer, before switching to a cluster or cloud environment with production data.

(lsf-executor)=

## LSF

The `lsf` executor allows you to run your pipeline script using a [Platform LSF](http://en.wikipedia.org/wiki/Platform_LSF) cluster.

Nextflow manages each process as a separate job that is submitted to the cluster using the `bsub` command.

The pipeline must be launched from a node where the `bsub` command is available, which is typically the cluster login node.

To enable the LSF executor, set `process.executor = 'lsf'` in the `nextflow.config` file.

Resource requests and other job characteristics can be controlled via the following process directives:

- {ref}`process-clusterOptions`
- {ref}`process-cpus`
- {ref}`process-memory`
- {ref}`process-queue`
- {ref}`process-time`

:::{note}
LSF supports both *per-core* and *per-job* memory limits. Nextflow assumes that LSF works in the *per-core* mode, thus it divides the requested {ref}`process-memory` by the number of requested {ref}`process-cpus`.

When LSF is configured to work in the *per-job* memory limit mode, you must specify this limit with the `perJobMemLimit` option in the {ref}`config-executor` scope of your Nextflow config file.

See also the [Platform LSF documentation](https://www.ibm.com/support/knowledgecenter/SSETD4_9.1.3/lsf_config_ref/lsf.conf.lsb_job_memlimit.5.dita).
:::

(moab-executor)=

## Moab

:::{versionadded} 19.07.0
:::

:::{warning} *Experimental: may change in a future release.*
:::

The `moab` executor allows you to run your pipeline script using the [Moab](https://en.wikipedia.org/wiki/Moab_Cluster_Suite) resource manager by [Adaptive Computing](http://www.adaptivecomputing.com/).

Nextflow manages each process as a separate job that is submitted to the cluster using the `msub` command provided by the resource manager.

The pipeline must be launched from a node where the `msub` command is available, which is typically the cluster login node.

To enable the `Moab` executor, set `process.executor = 'moab'` in the `nextflow.config` file.

Resource requests and other job characteristics can be controlled via the following process directives:

- {ref}`process-clusterOptions`
- {ref}`process-cpus`
- {ref}`process-memory`
- {ref}`process-queue`
- {ref}`process-time`

(nqsii-executor)=

## NQSII

The `nsqii` executor allows you to run your pipeline script using the [NQSII](https://www.rz.uni-kiel.de/en/our-portfolio/hiperf/nec-linux-cluster) resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster using the `qsub` command provided by the scheduler.

The pipeline must be launched from a node where the `qsub` command is available, which is typically the cluster login node.

To enable the NQSII executor, set `process.executor = 'nqsii'` in the `nextflow.config` file.

Resource requests and other job characteristics can be controlled via the following process directives:

- {ref}`process-clusterOptions`
- {ref}`process-cpus`
- {ref}`process-memory`
- {ref}`process-queue`
- {ref}`process-time`

(oar-executor)=

## OAR

:::{versionadded} 19.11.0-edge
:::

The `oar` executor allows you to run your pipeline script using the [OAR](https://oar.imag.fr) resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster using the `oarsub` command.

The pipeline must be launched from a node where the `oarsub` command is available, which is typically the cluster login node.

To enable the OAR executor set `process.executor = 'oar'` in the `nextflow.config` file.

Resource requests and other job characteristics can be controlled via the following process directives:

- {ref}`process-clusterOptions`
- {ref}`process-cpus`
- {ref}`process-memory`
- {ref}`process-queue`
- {ref}`process-time`

### Known Limitations

- Multiple `clusterOptions` should be semicolon-separated to ensure that the OAR job script is accurately formatted:
  ```groovy
  clusterOptions = '-t besteffort;--project myproject'
  ```

(pbs-executor)=

## PBS/Torque

The `pbs` executor allows you to run your pipeline script using a resource manager from the [PBS/Torque](http://en.wikipedia.org/wiki/Portable_Batch_System) family of batch schedulers.

Nextflow manages each process as a separate job that is submitted to the cluster using the `qsub` command provided by the scheduler.

The pipeline must be launched from a node where the `qsub` command is available, which is typically the cluster login node.

To enable the PBS executor, set `process.executor = 'pbs'` in the `nextflow.config` file.

Resource requests and other job characteristics can be controlled via the following process directives:

- {ref}`process-clusterOptions`
- {ref}`process-cpus`
- {ref}`process-memory`
- {ref}`process-queue`
- {ref}`process-time`

(pbspro-executor)=

## PBS Pro

The `pbspro` executor allows you to run your pipeline script using the [PBS Pro](https://www.pbspro.org/) resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster using the `qsub` command provided by the scheduler.

The pipeline must be launched from a node where the `qsub` command is available, which is typically the cluster login node.

To enable the PBS Pro executor, set `process.executor = 'pbspro'` in the `nextflow.config` file.

Resource requests and other job characteristics can be controlled via the following process directives:

- {ref}`process-clusterOptions`
- {ref}`process-cpus`
- {ref}`process-memory`
- {ref}`process-queue`
- {ref}`process-time`

(sge-executor)=

## SGE

The `sge` executor allows you to run your pipeline script using a [Sun Grid Engine](http://en.wikipedia.org/wiki/Oracle_Grid_Engine) cluster or a compatible platform ([Open Grid Engine](http://gridscheduler.sourceforge.net/), [Univa Grid Engine](http://www.univa.com/products/grid-engine.php), etc).

Nextflow manages each process as a separate grid job that is submitted to the cluster using the `qsub` command.

The pipeline must be launched from a node where the `qsub` command is available, which is typically the cluster login node.

To enable the SGE executor, set `process.executor = 'sge'` in the `nextflow.config` file.

Resource requests and other job characteristics can be controlled via the following process directives:

- {ref}`process-clusterOptions`
- {ref}`process-cpus`
- {ref}`process-memory`
- {ref}`process-penv`
- {ref}`process-queue`
- {ref}`process-time`

(slurm-executor)=

## SLURM

The `slurm` executor allows you to run your pipeline script using the [SLURM](https://slurm.schedmd.com/documentation.html) resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster using the `sbatch` command.

The pipeline must be launched from a node where the `sbatch` command is available, which is typically the cluster login node.

To enable the SLURM executor, set `process.executor = 'slurm'` in the `nextflow.config` file.

Resource requests and other job characteristics can be controlled via the following process directives:

- {ref}`process-clusterOptions`
- {ref}`process-cpus`
- {ref}`process-memory`
- {ref}`process-queue`
- {ref}`process-time`

:::{note}
SLURM partitions can be specified with the `queue` directive.
:::

:::{note}
Nextflow does not provide direct support for SLURM multi-clusters. If you need to submit workflow executions to a cluster other than the current one, specify it with the `SLURM_CLUSTERS` variable in the launch environment.
:::

:::{versionadded} 23.07.0-edge
Some SLURM clusters require memory allocations to be specified with `--mem-per-cpu` instead of `--mem`. You can specify `executor.perCpuMemAllocation = true` in the Nextflow configuration to enable this behavior. Nextflow will automatically compute the memory per CPU for each task (by default 1 CPU is used).
:::

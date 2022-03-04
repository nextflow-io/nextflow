.. _executor-page:

*********
Executors
*********

In the Nextflow framework architecture, the `executor` is the component that determines the system where a pipeline
process is run and supervises its execution.

The `executor` provides an abstraction between the pipeline processes and the underlying execution system. This
allows you to write the pipeline functional logic independently from the actual processing platform.

In other words you can write your pipeline script once and have it running on your computer, a cluster resource manager
or the cloud by simply changing the executor definition in the Nextflow configuration file.


.. toctree::
   :maxdepth: 1

   aws
   azure
   ga4gh
   google
   htcondor
   ignite
   kubernetes
   local
   lsf
   moab
   nsqii
   oar
   pbs
   pbspro
   sge
   slurm

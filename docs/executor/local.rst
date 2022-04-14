.. _local-executor:

Local
=====

The ``local`` executor is used by default. It runs the pipeline processes in the computer where Nextflow
is launched. The processes are parallelised by spawning multiple `threads` and by taking advantage of multi-cores
architecture provided by the CPU.

In a common usage scenario, the ``local`` executor can be useful to develop and test your pipeline script in your computer,
switching to a cluster facility when you need to run it on production data.

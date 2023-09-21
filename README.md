![Nextflow logo](https://github.com/nextflow-io/trademark/blob/master/nextflow2014_no-bg.png)

*"Dataflow variables are spectacularly expressive in concurrent programming"*
<br>[Henri E. Bal , Jennifer G. Steiner , Andrew S. Tanenbaum](https://dl.acm.org/doi/abs/10.1145/72551.72552)

![Nextflow CI](https://github.com/nextflow-io/nextflow/workflows/Nextflow%20CI/badge.svg)
[![Nextflow version](https://img.shields.io/github/release/nextflow-io/nextflow.svg?colorB=26af64&style=popout)](https://github.com/nextflow-io/nextflow/releases/latest)
[![Nextflow Twitter](https://img.shields.io/twitter/url/https/nextflowio.svg?colorB=26af64&&label=%40nextflow&style=popout)](https://twitter.com/nextflowio)
[![Nextflow Publication](https://img.shields.io/badge/Published-Nature%20Biotechnology-26af64.svg?colorB=26af64&style=popout)](https://www.nature.com/articles/nbt.3820)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?colorB=26af64&style=popout)](http://bioconda.github.io/recipes/nextflow/README.html)
[![Nextflow license](https://img.shields.io/github/license/nextflow-io/nextflow.svg?colorB=26af64&style=popout)](https://github.com/nextflow-io/nextflow/blob/master/COPYING)

Nextflow is a workflow manager for developing scalable, portable, and reproducible workflows. It can deploy workflows on a variety of execution platforms including local, HPC schedulers, AWS Batch, Azure Batch, Google Cloud Batch, and Kubernetes. Additionally, it supports many ways to manage your software dependencies, including Conda, Spack, Docker, Podman, Singularity, and more.

Quick start
===========

Install Nextflow with a single command:

```bash
curl -fsSL https://get.nextflow.io | bash
```

It creates the `nextflow` executable file in the current directory. You can then move it to a directory in your `$PATH` to run it from anywhere.

Nextflow can also be installed from Bioconda:

```bash
conda install -c bioconda nextflow
```

Documentation
=============

Nextflow's documentation is available at https://nextflow.io/docs/latest/.

Community
=========

You can post questions, or report problems by using the Nextflow [discussions](https://github.com/nextflow-io/nextflow/discussions)
or the Nextflow [Slack community chat](https://www.nextflow.io/slack-invite.html).

Nextflow also hosts a yearly workshop showcasing researcher's workflows and advancements in the language. Talks from the past workshops are available on the [Nextflow YouTube Channel](https://www.youtube.com/@Nextflow)

The [nf-core](https://nf-co.re/) project is a community effort aggregating high quality Nextflow workflows which can be used by the community.

Contributing
============

Contributions are more than welcome. See the [CONTRIBUTING](CONTRIBUTING.md) file for details.

License
=======

Nextflow is released under the Apache 2.0 license.

Citations
=========

If you use Nextflow in your work, please cite:

P. Di Tommaso, et al. Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319 (2017) doi:[10.1038/nbt.3820](http://www.nature.com/nbt/journal/v35/n4/full/nbt.3820.html)

Credits
=======

Nextflow is built on two excellent open source software projects: <a href='http://groovy-lang.org' target='_blank'>Groovy</a>
and <a href='http://www.gpars.org/' target='_blank'>Gpars</a>.

<a href='http://www.yourkit.com' target='_blank'>YourKit</a> is kindly supporting Nextflow with its fully-featured Java Profiler.


# Nextflow's documentation!

[![Nextflow CI](https://github.com/nextflow-io/nextflow/workflows/Nextflow%20CI/badge.svg)](https://github.com/nextflow-io/nextflow/actions/workflows/build.yml?query=branch%3Amaster+event%3Apush)
[![Nextflow version](https://img.shields.io/github/release/nextflow-io/nextflow.svg?colorB=58bd9f&style=popout)](https://github.com/nextflow-io/nextflow/releases/latest)
[![Nextflow Twitter](https://img.shields.io/twitter/url/https/nextflowio.svg?colorB=58bd9f&&label=%40nextflow&style=popout)](https://twitter.com/nextflowio)
[![Nextflow Publication](https://img.shields.io/badge/Published-Nature%20Biotechnology-26af64.svg?colorB=58bd9f&style=popout)](https://www.nature.com/articles/nbt.3820)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?colorB=58bd9f&style=popout)](http://bioconda.github.io/recipes/nextflow/README.html)
[![Nextflow license](https://img.shields.io/github/license/nextflow-io/nextflow.svg?colorB=58bd9f&style=popout)](https://github.com/nextflow-io/nextflow/blob/master/COPYING)

Nextflow is a workflow system for creating scalable, portable, and reproducible workflows.

## Rationale

The rise of big data has made it increasingly necessary to be able to analyze and perform experiments on large datasets in a portable and reproducible manner.

Parallelization and distributed computing are the best ways to tackle this challenge, but the tools commonly available to computational scientists often lack good support for these techniques, or they provide a model that fits poorly with the needs of computational scientists and often require knowledge of complex tools and APIs.

The Nextflow language is inspired by [the Unix philosophy](https://en.wikipedia.org/wiki/Unix_philosophy), in which many simple command line tools can be chained together into increasingly complex tasks. Similarly, a Nextflow script consists of composing many simple processes into increasingly complex pipelines. Each process executes a given tool or scripting language, and by specifying the process inputs and outputs, Nextflow coordinates the execution of tasks for you.

The Nextflow runtime integrates with many popular execution platforms (HPC schedulers, cloud providers) and software tools (Git, Docker, Conda), allowing you to fully describe a computational pipeline with all of its dependencies and run it in nearly any environment -- write once, run anywhere.

```{toctree}
:hidden:
:caption: User Guide
:maxdepth: 1

userguide/running-pipelines/index
userguide/developing-pipelines/index
userguide/sharing-pipelines
userguide/software-deps/index
userguide/compute/index
userguide/storage/index
```

```{toctree}
:hidden:
:caption: Language
:maxdepth: 1

reference/syntax
reference/stdlib
reference/process
reference/channel
reference/operator
reference/workflow
reference/module
reference/config
reference/dsl1
```

```{toctree}
:hidden:
:caption: Execution
:maxdepth: 1

reference/cli
reference/executor
reference/container
reference/git
reference/plugins
```

```{toctree}
:hidden:
:caption: Tutorials
:maxdepth: 1

tutorials/flux
tutorials/metrics
```

```{toctree}
:hidden:
:caption: Contributing
:maxdepth: 1

reference/developer/index
reference/developer/packages
reference/developer/plugins
```

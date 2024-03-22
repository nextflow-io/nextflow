
# Nextflow reference documentation

[![Nextflow CI](https://github.com/nextflow-io/nextflow/workflows/Nextflow%20CI/badge.svg)](https://github.com/nextflow-io/nextflow/actions/workflows/build.yml?query=branch%3Amaster+event%3Apush)
[![Nextflow version](https://img.shields.io/github/release/nextflow-io/nextflow.svg?colorB=58bd9f&style=popout)](https://github.com/nextflow-io/nextflow/releases/latest)
[![Nextflow Twitter](https://img.shields.io/twitter/url/https/nextflowio.svg?colorB=58bd9f&&label=%40nextflow&style=popout)](https://twitter.com/nextflowio)
[![Nextflow Publication](https://img.shields.io/badge/Published-Nature%20Biotechnology-26af64.svg?colorB=58bd9f&style=popout)](https://www.nature.com/articles/nbt.3820)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?colorB=58bd9f&style=popout)](http://bioconda.github.io/recipes/nextflow/README.html)
[![Nextflow license](https://img.shields.io/github/license/nextflow-io/nextflow.svg?colorB=58bd9f&style=popout)](https://github.com/nextflow-io/nextflow/blob/master/COPYING)

Nextflow is a workflow system for creating scalable, portable, and reproducible workflows.

- Get an {ref}`overview <overview-page>` of Nextflow and its key concepts.
- Get started with Nextflow by {ref}`installing <install-page>` it and running {ref}`your first script <your-first-script>`.
- Check out [this blog post](https://www.nextflow.io/blog/2023/learn-nextflow-in-2023.html) for even more resources on how to learn Nextflow.

```{toctree}
:hidden:
:caption: Introduction
:maxdepth: 1

overview
install
your-first-script
```

```{toctree}
:hidden:
:caption: Language
:maxdepth: 1

script
process
channel
operator
workflow
module
config
dsl1
```

```{toctree}
:hidden:
:caption: Execution
:maxdepth: 1

cli
executor
cache-and-resume
tracing
metrics
sharing
metadata
mail
plugins
secrets
```

```{toctree}
:hidden:
:caption: Software dependencies
:maxdepth: 1

container
conda
spack
wave
```

```{toctree}
:hidden:
:caption: Compute & storage platforms
:maxdepth: 1

aws
amazons3
azure
fusion
google
kubernetes
```

```{toctree}
:hidden:
:caption: Additional integrations
:maxdepth: 1

flux
ignite
```

```{toctree}
:hidden:
:caption: Contributing
:maxdepth: 1

developer/index
developer/packages
developer/plugins
```

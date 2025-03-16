
# Nextflow

*"Dataflow variables are spectacularly expressive in concurrent programming"*
<br>[Henri E. Bal , Jennifer G. Steiner , Andrew S. Tanenbaum](https://dl.acm.org/doi/abs/10.1145/72551.72552)

[![Nextflow CI](https://github.com/nextflow-io/nextflow/workflows/Nextflow%20CI/badge.svg)](https://github.com/nextflow-io/nextflow/actions/workflows/build.yml?query=branch%3Amaster+event%3Apush)
[![Nextflow version](https://img.shields.io/github/release/nextflow-io/nextflow.svg?colorB=58bd9f&style=popout)](https://github.com/nextflow-io/nextflow/releases/latest)
[![Nextflow Twitter](https://img.shields.io/twitter/url/https/nextflowio.svg?colorB=58bd9f&&label=%40nextflow&style=popout)](https://twitter.com/nextflowio)
[![Nextflow Publication](https://img.shields.io/badge/Published-Nature%20Biotechnology-26af64.svg?colorB=58bd9f&style=popout)](https://www.nature.com/articles/nbt.3820)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?colorB=58bd9f&style=popout)](http://bioconda.github.io/recipes/nextflow/README.html)
[![Nextflow license](https://img.shields.io/github/license/nextflow-io/nextflow.svg?colorB=58bd9f&style=popout)](https://github.com/nextflow-io/nextflow/blob/master/COPYING)

Nextflow is a workflow system for creating scalable, portable, and reproducible workflows. It uses a dataflow programming model that simplifies writing parallel and distributed pipelines by allowing you to focus on data flow and computation. Nextflow can deploy workflows on a variety of execution platforms, including your local machine, HPC schedulers, and cloud. Additionally, Nextflow supports a range of compute environments, software container runtimes, and package managers, allowing workflows to be executed in reproducible and isolated environments.

## Get started

To get started with Nextflow:

1. See the Nextflow {ref}`overview <overview-page>` to learn key concepts.
2. Download and {ref}`install <install-page>` Nextflow.
3. Set up an {ref}`environment <devenv-page>` with the {ref}`Nextflow VS Code extension <devenv-nextflow>`.
4. Run {ref}`your first script <your-first-script>`.

To continue learning about Nextflow, visit the [Nextflow community training portal](https://training.nextflow.io/latest/) and find a training course that is right for you. Seqera, the company that develops Nextflow, also runs a variety of training events. See [Seqera Events](https://seqera.io/events/) for more information.

## Community

You can post questions in the [Nextflow community forum](https://community.seqera.io) or the [Nextflow Slack](https://www.nextflow.io/slack-invite.html). Bugs and feature requests should be reported as [GitHub issues](https://github.com/nextflow-io/nextflow/issues/new/choose).

The Nextflow community is highly active with regular community meetings, events, a podcast, and more. You can view this material on the [Nextflow](https://www.youtube.com/@Nextflow) YouTube channel.

The [nf-core](https://nf-co.re/) project is a community effort aggregating high-quality Nextflow workflows that can be used by everyone.

## Contributing

Contributions to Nextflow are welcome. See {ref}`Contributing <contributing-page>` for more details.

## License

Nextflow is released under the Apache 2.0 license. Nextflow is a [registered trademark](https://github.com/nextflow-io/trademark).

## Citations

If you use Nextflow in your work, please cite:

P. Di Tommaso, et al. Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319 (2017) doi:[10.1038/nbt.3820](http://www.nature.com/nbt/journal/v35/n4/full/nbt.3820.html)

```{toctree}
:hidden:
:caption: Get started
:maxdepth: 1

overview
install
developer-env
your-first-script
```

```{toctree}
:hidden:
:caption: Running pipelines
:maxdepth: 1

cli
config
executor
cache-and-resume
reports
plugins
```

```{toctree}
:hidden:
:caption: Developing pipelines
:maxdepth: 1

script
working-with-files
process
channel
workflow
module
notifications
secrets
sharing
vscode
```

```{toctree}
:hidden:
:caption: Software dependencies
:maxdepth: 1

git
container
conda
spack
wave
```

```{toctree}
:hidden:
:caption: Compute & storage
:maxdepth: 1

aws
amazons3
azure
google
kubernetes
fusion
```

```{toctree}
:hidden:
:caption: Reference
:maxdepth: 1

reference/syntax
reference/cli
reference/config
reference/env-vars
reference/feature-flags
reference/stdlib
reference/process
reference/channel
reference/operator
```

```{toctree}
:hidden:
:caption: Updates
:maxdepth: 1

updating-nextflow
updating-syntax
dsl1
```

```{toctree}
:hidden:
:caption: Contributing
:maxdepth: 1

developer/index
developer/diagram
developer/packages
developer/plugins
```

```{toctree}
:hidden:
:caption: Guides
:maxdepth: 1

updating-spot-retries
metrics
flux
```


# Nextflow

*"Dataflow variables are spectacularly expressive in concurrent programming"*
<br />[Henri E. Bal , Jennifer G. Steiner , Andrew S. Tanenbaum](https://dl.acm.org/doi/abs/10.1145/72551.72552)

[![Nextflow CI](https://github.com/nextflow-io/nextflow/workflows/Nextflow%20CI/badge.svg)](https://github.com/nextflow-io/nextflow/actions/workflows/build.yml?query=branch%3Amaster+event%3Apush)
[![Nextflow version](https://img.shields.io/github/release/nextflow-io/nextflow.svg?colorB=58bd9f&style=popout)](https://github.com/nextflow-io/nextflow/releases/latest)
[![Nextflow Twitter](https://img.shields.io/twitter/url/https/nextflowio.svg?colorB=58bd9f&&label=%40nextflow&style=popout)](https://twitter.com/nextflowio)
[![Nextflow Publication](https://img.shields.io/badge/Published-Nature%20Biotechnology-26af64.svg?colorB=58bd9f&style=popout)](https://www.nature.com/articles/nbt.3820)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?colorB=58bd9f&style=popout)](http://bioconda.github.io/recipes/nextflow/README.html)
[![Nextflow license](https://img.shields.io/github/license/nextflow-io/nextflow.svg?colorB=58bd9f&style=popout)](https://github.com/nextflow-io/nextflow/blob/master/COPYING)

Nextflow is a workflow system for creating scalable, portable, and reproducible workflows. It uses a dataflow programming model that simplifies writing parallel and distributed pipelines by allowing you to focus on data flow and computation. Nextflow can deploy workflows on a variety of execution platforms, including your local machine, HPC schedulers, and cloud. Additionally, Nextflow supports a range of compute environments, software container runtimes, and package managers, allowing workflows to be executed in reproducible and isolated environments.

## Get started

To get started with Nextflow:

1. See the Nextflow [overview][overview-page] to learn key concepts.
2. Download and [install][install-page>] Nextflow.
3. Set up an [environment][devenv-page] with the [Nextflow VS Code extension][devenv-nextflow].
4. Run [your first script][your-first-script].

To continue learning about Nextflow, visit the [Nextflow community training portal](https://training.nextflow.io/latest/) and find a training course that is right for you. Seqera, the company that develops Nextflow, also runs a variety of training events. See [Seqera Events](https://seqera.io/events/) for more information.

## Community

You can post questions in the [Nextflow community forum](https://community.seqera.io) or the [Nextflow Slack](https://www.nextflow.io/slack-invite.html). Bugs and feature requests should be reported as [GitHub issues](https://github.com/nextflow-io/nextflow/issues/new/choose).

The Nextflow community is highly active with regular community meetings, events, a podcast, and more. You can view this material on the [Nextflow](https://www.youtube.com/@Nextflow) YouTube channel.

The [nf-core](https://nf-co.re/) project is a community effort aggregating high-quality Nextflow workflows that can be used by everyone.

## Contributing

Contributions to Nextflow are welcome. See [Contributing][contributing-page] for more details.

## License

Nextflow is released under the Apache 2.0 license. Nextflow is a [registered trademark](https://github.com/nextflow-io/trademark).

## Citations

If you use Nextflow in your work, please cite:

P. Di Tommaso, et al. Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319 (2017) doi:[10.1038/nbt.3820](http://www.nature.com/nbt/journal/v35/n4/full/nbt.3820.html)

[contributing-page]: /nextflow_docs/nextflow_repo/docs/developer/index
[overview-page]: /nextflow_docs/nextflow_repo/docs/overview
[install-page]: /nextflow_docs/nextflow_repo/docs/install.md
[devenv-page]: /nextflow_docs/nextflow_repo/docs/developer-env.mdx
[devenv-nextflow]: /nextflow_docs/nextflow_repo/docs/developer-env#nextflow
[your-first-script]: /nextflow_docs/nextflow_repo/docs/your-first-script
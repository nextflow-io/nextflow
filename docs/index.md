
# Nextflow's documentation!

The rise of big data has made it increasingly necessary to be able to analyze and perform experiments on large datasets in a portable and reproducible manner.

Parallelization and distributed computing are the best ways to tackle this challenge, but the tools commonly available to computational scientists often lack good support for these techniques, or they provide a model that fits poorly with the needs of computational scientists and often require knowledge of complex tools and APIs.

Nextflow is a workflow system for creating scalable, portable, and reproducible workflows. It is based on the dataflow programming model, which greatly simplifies writing parallel and distributed pipelines, allowing you to focus on the actual computation and flow of data.

The Nextflow language is inspired by [the Unix philosophy](https://en.wikipedia.org/wiki/Unix_philosophy), in which many simple command line tools can be chained together into increasingly complex tasks. Similarly, a Nextflow script consists of composing many simple processes into increasingly complex pipelines. Each process executes a given tool or scripting language, and by specifying the process inputs and outputs, Nextflow coordinates the execution of tasks for you.

The Nextflow runtime integrates with many popular execution platforms (HPC schedulers, cloud providers) and software tools (Git, Docker, Conda), allowing you to fully describe a computational pipeline with all of its dependencies and run it in nearly any environment -- write once, run anywhere.

```{toctree}
:hidden:

getstarted
basic
script
process
channel
operator
executor
config
dsl2
cli
container
conda
spack
wave
fusion
aws
amazons3
azure
flux
google
ignite
kubernetes
tracing
metrics
sharing
metadata
mail
plugins
secrets
```

(install-page)=

# Installation

(install-requirements)=

## Requirements

Nextflow can be used on any POSIX compatible system (Linux, macOS, etc). It requires Bash 3.2 (or later) and [Java 11 (or later, up to 21)](http://www.oracle.com/technetwork/java/javase/downloads/index.html) to be installed.

For the execution in a cluster of computers, the use of a shared file system is required to allow the sharing of tasks input/output files.

Nextflow can also be run on Windows through [WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux).

:::{tip}
We recommend that you install Java through [SDKMAN!](https://sdkman.io/), and that you use the latest LTS version of Corretto or Temurin. See [this website](https://whichjdk.com/) for more information. While other Java distros may work at first or even most of the time, many users have experienced issues that are difficult to debug and are usually resolved by using one of the recommended distros.

To install Corretto 17:

```bash
sdk install java 17.0.6-amzn
```

To install Temurin 17:

```bash
sdk install java 17.0.6-tem
```
:::

(install-nextflow)=

## Installation

Nextflow is distributed as a self-installing package, which means that it does not require any special installation procedure.

It only needs two easy steps:

1. Download the executable package by copying and pasting either one of the following commands in your terminal window: `wget -qO- https://get.nextflow.io | bash`

   Or, if you prefer `curl`: `curl -s https://get.nextflow.io | bash`

   This will create the `nextflow` main executable file in the current directory.

2. Make the binary executable on your system by running `chmod +x nextflow`.

3. Optionally, move the `nextflow` file to a directory accessible by your `$PATH` variable (this is only required to avoid remembering and typing the full path to `nextflow` each time you need to run it).

:::{tip}
Set `export CAPSULE_LOG=none` to make the dependency installation logs less verbose.
:::

:::{tip}
If you don't have `curl` or `wget`, you can also download the Nextflow launcher script from the [project releases page](https://github.com/nextflow-io/nextflow/releases/latest) on GitHub, in lieu of step 1.
:::

:::{tip}
To avoid downloading the dependencies, you can also use the `nextflow-VERSION-all` distribution available for every Nextflow release on Github.

1. Go to the [Github releases page](https://github.com/nextflow-io/nextflow/releases) and expand the `Assets` section for a specific release.
2. Copy the URL of the `nextflow-VERSION-all` asset and enter the download command in your terminal, e.g. `wget -qO- ASSET-URL`. It will create the completely self-contained `nextflow-VERSION-all` executable file in the current directory.
:::

## Updates

Having Nextflow installed in your computer you can update to the latest version using the following command:

```bash
nextflow self-update
```

:::{tip}
You can temporarily switch to a specific version of Nextflow by prefixing the `nextflow` command with the `NXF_VER` environment variable. For example:

```bash
NXF_VER=20.04.0 nextflow run hello
```
:::

## Stable and Edge releases

A *stable* version of Nextflow is released on a six-months basic schedule, in the 1st and 3rd quarter of every year.

Along with the stable release, an *edge* version is released on a monthly basis. This version is useful to test and use most recent updates and experimental features.

To use the latest edge release run the following snippet in your shell terminal:

```bash
export NXF_EDGE=1
nextflow self-update
```

(updating-nextflow-page)=

# Updating Nextflow

This page describes Nextflow release cadence, how to self-update Nextflow, and how select your version of Nextflow.

## Releases

A stable version of Nextflow is released in the 4th and 10th month of each year. An edge version of Nextflow is released on a monthly basis. The edge version can be used to access the latest updates and experimental features.

Nextflow uses [Calendar Versioning](https://calver.org). Versions are numbered as `<year>.<month>.<patch>`. For example, `23.10.1` corresponds to the first patch of the October 2023 stable release.

You can find an exhaustive list of releases and updates in the [Nextflow changelog](https://github.com/nextflow-io/nextflow/blob/master/changelog.txt).

## Self-update

To update to the latest stable release of Nextflow, run the `self-update` command:

```{code-block} bash
:class: copyable
nextflow self-update
```

To use the latest edge release, set `NXF_EDGE=1` when you self-update Nextflow:

```{code-block} bash
:class: copyable
NXF_EDGE=1 nextflow self-update
```

:::{warning}
Nextflow will update its executable during the self-update process. The update can fail if the Nextflow executable is in a directory with restricted permissions.
:::

## Version selection

The `NXF_VER` environment variable can be used to define which version of Nextflow to use. To switch to a specific version of Nextflow for a single run, set the `NXF_VER` environment variable in your execution command. For example:

```{code-block} bash
:class: copyable
NXF_VER=23.10.0 nextflow info
```

To set a specific version of Nextflow for a terminal session, export the `NXF_VER` environment variable. For example:

```{code-block} bash
:class: copyable
export NXF_VER=23.10.0
```

To set a specific version of Nextflow for your user profile, add the above `NXF_VER` export command to your shell configuration file, such as `~/.bashrc` or `~/.zshrc`, and restart your session.

:::{tip}
You can use `NXF_VER` to switch to an edge release. For example:

```{code-block} bash
:class: copyable
NXF_VER=24.06.0-edge nextflow info
```
:::

:::{warning}
Nextflow will update its executable during the self-update process. The update can fail if the Nextflow executable is in a directory with restricted permissions.
:::

(install-page)=

# Installation

(install-requirements)=

## Requirements

Nextflow can be used on any POSIX-compatible system (Linux, macOS, etc), and on Windows through [WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux). It requires Bash 3.2 (or later) and [Java 17 (or later, up to 23)](http://www.oracle.com/technetwork/java/javase/downloads/index.html) to be installed. You can see which version you have using the following command:

```{code-block} bash
:class: copyable
java -version
```

:::{versionchanged} 24.11.0-edge
Support for Java versions prior to 17 was dropped.
:::

If you don't have a compatible version of Java installed in your computer, it is recommended that you install it through [SDKMAN!](https://sdkman.io/), and that you use the latest LTS version of Temurin. See [this website](https://whichjdk.com/) for more information.

To install Java with SDKMAN:

1. [Install SDKMAN](https://sdkman.io/install):

    ```{code-block} bash
    :class: copyable
    curl -s https://get.sdkman.io | bash
    ```

2. Open a new terminal.

3. Install Java:

    ```{code-block} bash
    :class: copyable
    sdk install java 17.0.10-tem
    ```

4. Confirm that Java is installed correctly:

    ```{code-block} bash
    :class: copyable
    java -version
    ```

(install-nextflow)=

## Install Nextflow

Nextflow is distributed as a self-installing package, in order to make the installation process as simple as possible:

1. Install Nextflow:

    ```{code-block} bash
    :class: copyable
    curl -s https://get.nextflow.io | bash
    ```

    This will create the `nextflow` executable in the current directory.

    :::{tip}
    Set `export CAPSULE_LOG=none` to make the installation logs less verbose.
    :::

2. Make Nextflow executable:

    ```{code-block} bash
    :class: copyable
    chmod +x nextflow
    ```

3. Move Nextflow into an executable path. For example:

    ```{code-block} bash
    :class: copyable
    mkdir -p $HOME/.local/bin/
    mv nextflow $HOME/.local/bin/
    ```

    :::{tip}
    Ensure the directory `$HOME/.local/bin/` is included in your `PATH` variable. Temporarily add this directory to `PATH` by setting `export PATH="$PATH:$HOME/.local/bin"`. Add the directory to `PATH` permanently by adding the export command to your shell configuration file, such as `~/.bashrc` or `~/.zshrc`. Alternatively, move the `nextflow` executable to a directory already in your `PATH`.
    :::

    :::{warning}
    Nextflow will update its executable during the self update process, therefore the update can fail if the executable is placed in a directory with restricted permissions.
    :::

4. Confirm that Nextflow is installed correctly:

    ```{code-block} bash
    :class: copyable
    nextflow info
    ```

## Updates

With Nextflow installed in your environment, you can update to the latest version using the following command:

```{code-block} bash
:class: copyable
nextflow self-update
```

You can also temporarily switch to a specific version of Nextflow with the `NXF_VER` environment variable. For example:

```{code-block} bash
:class: copyable
NXF_VER=23.10.0 nextflow info
```

## Seqera Platform

You can launch workflows directly from [Seqera Platform](https://seqera.io/platform/) without installing Nextflow locally.

Launching from Seqera Platform provides you with:
- User-friendly launch interfaces.
- Automated cloud infrastructure creation.
- Organizational user management.
- Advanced analytics with resource optimization.

Seqera Cloud Basic is free for small teams. Researchers at qualifying academic institutions can apply for free access to Seqera Cloud Pro.
See the [Seqera Platform documentation](https://docs.seqera.io/platform) for set-up information and tutorials to get started.

## Stable and edge releases

A *stable* version of Nextflow is released every six months, in the 4th and 10th month of each year.

Additionally, an *edge* version is released on a monthly basis. The edge releases can be used to access the latest updates and experimental features.

To use the latest edge release, set `NXF_EDGE=1` when updating:

```{code-block} bash
:class: copyable
NXF_EDGE=1 nextflow self-update
```

You can also use `NXF_VER` to temporarily switch to any edge release. For example:

```{code-block} bash
:class: copyable
NXF_VER=24.06.0-edge nextflow info
```

## Standalone distribution

The Nextflow standalone distribution (i.e. the `dist` distribution) consists of self-contained `nextflow` executable file
that includes all the application dependencies for core functionalities, and it can run without downloading third parties
libraries. This distribution is mainly useful for offline environments.

Note however the support for cloud services e.g. AWS, Seqera Platform, Wave, etc. still require the download
of the corresponding Nextflow plugins.

To use the standalone distribution:

1. Download the standalone distribution from Assets section of the [GitHub releases page](https://github.com/nextflow-io/nextflow/releases).

2. Grant execution permissions to the downloaded file. For example:

    ```{code-block} bash
    :class: copyable
    chmod -x nextflow-24.10.1-dist
    ```

3. Use it as a drop-in replacement for `nextflow` command. For example:

    ```{code-block} bash
    :class: copyable
    ./nextflow-24.10.1-dist run hello
    ```

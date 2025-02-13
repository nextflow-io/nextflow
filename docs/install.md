(install-page)=

# Installation

Nextflow can be used on any POSIX-compatible system (Linux, macOS, etc), and on Windows through [WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux). This page describes how to install Nextflow.

:::{note}
New versions of Nextflow are released regularly. See {ref}`updating-nextflow-page` for more information about Nextflow release cadence, how to update Nextflow, and how select your version of Nextflow.
:::

(install-requirements)=

## Requirements

Nextflow requires Bash 3.2 (or later) and [Java 17 (or later, up to 23)](http://www.oracle.com/technetwork/java/javase/downloads/index.html) to be installed. To see which version of Java you have, run the following command:

```{code-block} bash
:class: copyable
java -version
```

:::{versionchanged} 24.11.0-edge
Support for Java versions prior to 17 was dropped.
:::

If you don't have a compatible version of Java installed, it is recommended that you install it through [SDKMAN!](https://sdkman.io/), and that you use the latest Long-Term-Support (LTS) version of Temurin. See [Which version of JDK should I use?](https://whichjdk.com/) for more information about different versions of Java.

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

To install Nextflow:

1. Download Nextflow:

    ```{code-block} bash
    :class: copyable
    curl -s https://get.nextflow.io | bash
    ```

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

## Seqera Platform

You can launch workflows directly from [Seqera Platform](https://seqera.io/platform/) without installing Nextflow locally.

Launching from Seqera Platform provides you with:

- User-friendly launch interfaces.
- Automated cloud infrastructure creation.
- Organizational user management.
- Advanced analytics with resource optimization.

Seqera Cloud Basic is free for small teams. Researchers at qualifying academic institutions can apply for free access to Seqera Cloud Pro.
See the [Seqera Platform documentation](https://docs.seqera.io/platform) for set-up information and tutorials to get started.

## Standalone distribution

The Nextflow standalone distribution (i.e. the `dist` release) is a self-contained `nextflow` executable that can run without needing to download core dependencies at runtime. This distribution is useful for offline environments, as well as building and testing Nextflow locally.

The standalone distribution will still download core and third-party plugins as needed at runtime.

To use the standalone distribution:

1. Download the standalone distribution from Assets section of the [GitHub releases page](https://github.com/nextflow-io/nextflow/releases).

2. Grant execution permissions to the downloaded file. For example:

    ```{code-block} bash
    :class: copyable
    chmod +x nextflow-24.10.1-dist
    ```

3. Use it as a drop-in replacement for `nextflow` command. For example:

    ```{code-block} bash
    :class: copyable
    ./nextflow-24.10.1-dist run hello
    ```

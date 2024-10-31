(install-page)=

# Installation

(install-requirements)=

## Requirements

Nextflow can be used on any POSIX-compatible system (Linux, macOS, etc), and on Windows through [WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux). It requires Bash 3.2 (or later) and [Java 17 (or later, up to 23)](http://www.oracle.com/technetwork/java/javase/downloads/index.html) to be installed. You can see which version you have using the following command:

```bash
java -version
```

:::{versionchanged} 24.11.0-edge
Support for Java versions prior to 17 was dropped.
:::

If you don't have a compatible version of Java installed in your computer, it is recommended that you install it through [SDKMAN!](https://sdkman.io/), and that you use the latest LTS version of Temurin. See [this website](https://whichjdk.com/) for more information.

To install Java with SDKMAN:

1. [Install SDKMAN](https://sdkman.io/install):

    ```bash
    curl -s https://get.sdkman.io | bash
    ```

2. Open a new terminal.

3. Install Java:

    ```bash
    sdk install java 17.0.10-tem
    ```

4. Confirm that Java is installed correctly:

    ```bash
    java -version
    ```

(install-nextflow)=

## Install Nextflow

Nextflow is distributed as a self-installing package, in order to make the installation process as simple as possible:

1. Install Nextflow:

    ```bash
    curl -s https://get.nextflow.io | bash
    ```

    This will create the `nextflow` executable in the current directory.

    :::{tip}
    You can set `export CAPSULE_LOG=none` to make the installation logs less verbose.
    :::

2. Make Nextflow executable:

    ```bash
    chmod +x nextflow
    ```

3. Move Nextflow into an executable path. For example:

    ```bash
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

    ```bash
    nextflow info
    ```

## Updates

With Nextflow installed in your environment, you can update to the latest version using the following command:

```bash
nextflow self-update
```

You can also temporarily switch to a specific version of Nextflow with the `NXF_VER` environment variable. For example:

```bash
NXF_VER=23.10.0 nextflow info
```

## Stable and edge releases

A *stable* version of Nextflow is released every six months, in the 4th and 10th month of each year.

Additionally, an *edge* version is released on a monthly basis. The edge releases can be used to access the latest updates and experimental features.

To use the latest edge release, set `NXF_EDGE=1` when updating:

```bash
NXF_EDGE=1 nextflow self-update
```

You can also use `NXF_VER` to temporarily switch to any edge release. For example:

```bash
NXF_VER=24.06.0-edge nextflow info
```

## Standalone distribution

Nextflow has a set of {ref}`core plugins <plugins-core>` which are downloaded at runtime by default. There is also a standalone distribution (i.e. the `all` distribution) which comes pre-packaged with all core plugins. This distribution is mainly useful for offline environments.

The installer for the `all` distribution can be found on the [GitHub releases page](https://github.com/nextflow-io/nextflow/releases), under the "Assets" section for a specific release. The installation procedure is the same as for the standard distribution, only using this URL instead of `https://get.nextflow.io`:

```bash
export NXF_VER=23.10.0
curl -s https://github.com/nextflow-io/nextflow/releases/download/v$NXF_VER/nextflow-$NXF_VER-all
```

:::{warning}
The `all` distribution does not support third-party plugins. Only the {ref}`core plugins <plugins-core>` are supported.
:::

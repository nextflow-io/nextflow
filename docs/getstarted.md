(getstarted-page)=

# Get started
 
 This guide covers the installation pThis guide will help you get started quickly. 

- Requirements
- Installation
- Your first script

Also take a look at:

- {ref}`config-page`
- {ref}`workflow-page`
- {ref}`operator-page`

If you're just starting out with Nextflow, [take a look at the tutorials on this page](https://www.nextflow.io/blog/2023/learn-nextflow-in-2023.html#nf-core). Then come back here. 

(getstarted-requirement)=

## Requirements

Nextflow can be used on any POSIX compatible system (Linux, macOS, etc). It requires Bash 3.2 (or later) and [Java 11 (or later, up to 21)](http://www.oracle.com/technetwork/java/javase/downloads/index.html) to be installed. You can check what version you have using the following command:

```bash
java -version
```

We recommend that you install Java through [SDKMAN!](https://sdkman.io/), and that you use the latest LTS version of Temurin. See [this website](https://whichjdk.com/) for more information.

To install Java with SDKMAN:

1. Install SDKMAN:

    ```bash
    curl -s https://get.sdkman.io | bash
    ```

2. Open a new terminal and install Java:

    ```bash
    sdk install java 17.0.10-tem
    ```

3. Confirm that Java is installed correctly:

    ```bash
    java -version
    ```

4. To install Temurin 17:

   ```bash
   sdk install java 17.0.6-tem
   ```

## Updates

With Nextflow installed in your environment, you can update to the latest version using the following command:

```bash
nextflow self-update
```

You can also temporarily switch to a specific version of Nextflow with the `NXF_VER` environment variable. For example:

```bash
NXF_VER=22.10.0 nextflow run hello
```

## Stable and Edge releases

A *stable* version of Nextflow is released every six months, in the 4th and 10th month of each year.

Additionally, an *edge* version is released on a monthly basis. The edge releases can be used to access the latest updates and experimental features.

To use the latest edge release, set `NXF_EDGE=1` when updating:

```bash
NXF_EDGE=1 nextflow self-update
```

You can also use `NXF_VER` to switch to any edge release:

```bash
$ nextflow info
```

For the execution in a cluster of computers, the use of a shared file system is required to allow the sharing of tasks input/output files. Nextflow can also be run on Windows through [WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux).

(getstarted-install)=

## Installation

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

3. Move Nextflow into an executable path:

    ```bash
    sudo mv nextflow /usr/local/bin
    ```

4. Confirm that Nextflow is installed correctly:

    ```bash
    nextflow info
    ```

### Standalone distribution

Nextflow has a set of {ref}`core plugins <plugins-core>` which are downloaded at runtime by default. There is also a standalone distribution (i.e. the `all` distribution) which comes pre-packaged with all core plugins. This distribution is mainly useful for offline environments.

The installer for the `all` distribution can be found on the [GitHub releases page](https://github.com/nextflow-io/nextflow/releases), under the "Assets" section for a specific release. The installation procedure is the same as for the standard distribution, only using this URL instead of `https://get.nextflow.io`:

```bash
export NXF_VER=23.10.0
curl -s https://github.com/nextflow-io/nextflow/releases/download/v$NXF_VER/nextflow-$NXF_VER-all
```

:::{warning}
The `all` distribution does not support third-party plugins. Only the {ref}`core plugins <plugins-core>` are supported.
:::

## Updates

With Nextflow installed in your environment, you can update to the latest version using the following command:

```bash
nextflow self-update
```

You can also temporarily switch to a specific version of Nextflow with the `NXF_VER` environment variable. For example:

```bash
NXF_VER=22.10.0 nextflow run hello
```

(getstarted-first)=

## Your first script

:::{note}
For versions of Nextflow prior to `22.10.0`, you must explicitly enable DSL2 by adding `nextflow.enable.dsl=2` to the top of the script or by using the `-dsl2` command-line option.
:::

Copy the following example into your favorite text editor and save it to a file named `tutorial.nf`:

```{literalinclude} snippets/your-first-script.nf
:language: groovy
```

This script defines two processes. The first splits a string into 6-character chunks, writing each one to a file with the prefix `chunk_`, and the second receives these files and transforms their contents to uppercase letters. The resulting strings are emitted on the `result` channel and the final output is printed by the `view` operator.

Execute the script by entering the following command in your terminal:

```console
$ nextflow run tutorial.nf

N E X T F L O W  ~  version 22.10.0
executor >  local (3)
[69/c8ea4a] process > splitLetters   [100%] 1 of 1 ✔
[84/c8b7f1] process > convertToUpper [100%] 2 of 2 ✔
HELLO
WORLD!
```

You can see that the first process is executed once, and the second twice. Finally the result string is printed.

It's worth noting that the process `convertToUpper` is executed in parallel, so there's no guarantee that the instance processing the first split (the chunk `Hello`) will be executed before the one processing the second split (the chunk `world!`).

Thus, it is perfectly possible that you will get the final result printed out in a different order:

```
WORLD!
HELLO
```

:::{tip}
The hexadecimal string, e.g. `22/7548fa`, is the unique hash of a task, and the prefix of the directory where the task is executed. You can inspect a task's files by changing to the directory `$PWD/work` and using this string to find the specific task directory.
:::

(getstarted-resume)=

### Modify and resume

Nextflow keeps track of all the processes executed in your pipeline. If you modify some parts of your script, only the processes that are actually changed will be re-executed. The execution of the processes that are not changed will be skipped and the cached result used instead.

This helps a lot when testing or modifying part of your pipeline without having to re-execute it from scratch.

For the sake of this tutorial, modify the `convertToUpper` process in the previous example, replacing the process script with the string `rev $x`, so that the process looks like this:

```groovy
process convertToUpper {
  input:
    path x
  output:
    stdout

  """
  rev $x
  """
}
```

Then save the file with the same name, and execute it by adding the `-resume` option to the command line:

```bash
nextflow run tutorial.nf -resume
```

It will print output similar to this:

```
N E X T F L O W  ~  version 22.10.0
executor >  local (2)
[69/c8ea4a] process > splitLetters   [100%] 1 of 1, cached: 1 ✔
[d0/e94f07] process > convertToUpper [100%] 2 of 2 ✔
olleH
!dlrow
```

You will see that the execution of the process `splitLetters` is actually skipped (the process ID is the same), and its results are retrieved from the cache. The second process is executed as expected, printing the reversed strings.

:::{tip}
The pipeline results are cached by default in the directory `$PWD/work`. Depending on your script, this folder can take up a lot of disk space. It's a good idea to clean this folder periodically, as long as you know you won't need to resume any pipeline runs.
:::

For more information, see the {ref}`cache-resume-page` page.

(getstarted-params)=

### Pipeline parameters

Pipeline parameters are simply declared by prepending to a variable name the prefix `params`, separated by dot character. Their value can be specified on the command line by prefixing the parameter name with a double dash character, i.e. `--paramName`

For the sake of this tutorial, you can try to execute the previous example specifying a different input string parameter, as shown below:

```bash
nextflow run tutorial.nf --str 'Bonjour le monde'
```

The string specified on the command line will override the default value of the parameter. The output will look like this:

```
N E X T F L O W  ~  version 22.10.0
executor >  local (4)
[8b/16e7d7] process > splitLetters   [100%] 1 of 1 ✔
[eb/729772] process > convertToUpper [100%] 3 of 3 ✔
m el r
edno
uojnoB
```

:::{versionchanged} 20.11.0-edge
Any `.` (dot) character in a parameter name is interpreted as the delimiter of a nested scope. For example, `--foo.bar Hello` will be interpreted as `params.foo.bar`. If you want to have a parameter name that contains a `.` (dot) character, escape it using the back-slash character, e.g. `--foo\.bar Hello`.
:::


## Stable and Edge releases

A *stable* version of Nextflow is released on a six-months basic schedule, in the 1st and 3rd quarter of every year.

Along with the stable release, an *edge* version is released on a monthly basis. This version is useful to test and use most recent updates and experimental features.

To use the latest edge release run the following snippet in your shell terminal:

```bash
export NXF_EDGE=1
nextflow self-update
```


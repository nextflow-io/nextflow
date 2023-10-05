# Overview

This section provides a high-level overview of the Nextflow source code for users who want to understand or contribute to it. Rather than a comprehensive API documentation, these docs simply provide a conceptual map to help you understand the key concepts of the Nextflow implementation, and to quickly find code sections of interest for further investigation.

Before you dive into code, be sure to check out the [CONTRIBUTING.md](https://github.com/nextflow-io/nextflow/blob/master/CONTRIBUTING.md) for Nextflow to learn about the many ways to contribute to the project.

## IntelliJ IDEA

The suggested development environment is [IntelliJ IDEA](https://www.jetbrains.com/idea/download/). Nextflow development with IntelliJ IDEA requires a recent version of the IDE (2019.1.2 or later).

After installing IntelliJ IDEA, use the following steps to use it with Nextflow:

1. Clone the Nextflow repository to a directory in your computer.

2. Open IntelliJ IDEA and go to **File > New > Project from Existing Sources...**.

3. Select the Nextflow project root directory in your computer and click **OK**.

4. Select **Import project from external model > Gradle** and click **Finish**.

5. After the import process completes, select **File > Project Structure...**.

6. Select **Project**, and make sure that the **SDK** field contains Java 11 (or later).

7. Go to **File > Settings > Editor > Code Style > Groovy > Imports** and apply the following settings:

   * Use single class import
   * Class count to use import with '*': `99`
   * Names count to use static import with '*': `99`
   * Imports layout:
      * `import java.*`
      * `import javax.*`
      * *blank line*
      * all other imports
      * all other static imports

New files must include the appropriate license header boilerplate and the author name(s) and contact email(s) ([see for example](https://github.com/nextflow-io/nextflow/blob/e8945e8b6fc355d3f2eec793d8f288515db2f409/modules/nextflow/src/main/groovy/nextflow/Const.groovy#L1-L15)).

## Groovy

Nextflow is written in [Groovy](http://groovy-lang.org/), which is itself a programming language based on [Java](https://www.java.com/). Groovy is designed to be highly interoperable with Java -- Groovy programs compile to Java bytecode, and nearly any Java program is also a valid Groovy program. However, Groovy adds several language features (e.g. closures, list and map literals, optional typing, optional semicolons, meta-programming) and standard libraries (e.g. JSON and XML parsing) that greatly improve the overall experience of developing for the Java virtual machine.

Recommended resources for Groovy, from most reference-complete to most user-friendly, are listed below:

- [Groovy documentation](http://groovy-lang.org/documentation.html)
- [Groovy in Action](https://www.manning.com/books/groovy-in-action-second-edition)
- [Groovy: The Awesome Parts](https://www.slideshare.net/paulk_asert/awesome-groovy)
- [Groovy cheat sheet](http://www.cheat-sheets.org/saved-copy/rc015-groovy_online.pdf)

## Software Dependencies

Nextflow depends on a variety of libraries and frameworks, the most prominent of which are listed below:

- [AWS SDK for Java 1.x](https://aws.amazon.com/sdk-for-java/): AWS integration
- [Azure SDK for Java](https://learn.microsoft.com/en-us/azure/developer/java/sdk/): Azure integration
- [Google Cloud Client Libraries for Java](https://cloud.google.com/java/docs/reference): Google Cloud integration
- [GPars](http://gpars.org/1.2.1/guide/guide/dataflow.html): dataflow concurrency
- [Gradle](https://gradle.org/): build automation
- [JCommander](https://jcommander.org/): command line interface
- [JGit](https://www.eclipse.org/jgit/): Git integration
- [Kryo](https://github.com/EsotericSoftware/kryo): serialization
- [LevelDB](https://mvnrepository.com/artifact/org.iq80.leveldb/leveldb): key-value store for the cache database
- [Logback](https://logback.qos.ch/): application logging
- [PF4J](https://pf4j.org/): plugin extensions
- [Spock](https://spockframework.org/): unit testing framework

Any other integrations are likely implemented using a CLI (e.g. Conda, Docker, HPC schedulers) or REST API (e.g. Kubernetes).

## Class Diagrams

Each package has a class diagram, abridged and annotated for relevance and ease of use.

Each node is a class. Fields are selectively documented in order to show only the core data structures and the classes that "own" them. Methods are not explicitly documented, but they are mentioned in certain links where appropriate. Links are selectively documented in order to show only the most important classes and relationships.

Links between classes denote one of the following relationships:

- Inheritance (`A <|-- B`): `B` is a subclass of `A`
- Composition (`A --* B`): `A` contains `B`
- Instantiation (`A --> B : f`): `A` creates instance(s) of `B` at runtime via `A::f()`

See {ref}`packages-page` for the list of Nextflow packages.

```{warning}
Class diagrams are manually curated, so they might not always reflect the latest version of the source code.
```

## Building from source

If you are interested in modifying the source code, you only need Java 11 or later to build Nextflow from source. Nextflow uses the [Gradle](http://www.gradle.org/) build automation system, but you do not need to install Gradle to build Nextflow. In other words, if you can run Nextflow, then you can probably build it too!

To build locally from a branch (useful for testing PRs):

```bash
git clone -b <branch> git@github.com:nextflow-io/nextflow.git
cd nextflow
make compile
```

The build system will automatically download all of the necessary dependencies on the first run, which may take several minutes.

Once complete, you can run your local build of Nextflow using the `launch.sh` script in place of the `nextflow` command:

```bash
./launch.sh run <script> ...
```

A self-contained executable Nextflow package can be created with the following command:

```bash
make pack
```

Again, use `launch.sh` in place of the `nextflow` command to use your local build.

## Testing

To run the unit tests:

```bash
# run all tests
make test

# run individual test
make test module=<nextflow|plugins:nf-amazon|...> class=<package>.<class>.<method>

# refer to the Makefile for all build rules
```

When a test fails, it will give you a report that you can open in your browser to view the reason for each failed test. The **Standard output** tab is particularly useful as it shows the console output of each test.

Refer to the [build.yml](https://github.com/nextflow-io/nextflow/tree/master/.github/workflows/build.yml) configuration to see how to run integration tests locally, if you are interested.

## Installing from source

The `nextflow` command is just a Bash script that downloads and executes the Nextflow JAR. When you install Nextflow using `get.nextflow.io`, it only downloads this launcher script, while the Nextflow JAR is downloaded on the first Nextflow run.

You can run `make install` to install a local build of the Nextflow JAR to `$NXF_HOME`. Note that this command will overwrite any existing Nextflow packages with the same version. This approach is useful for testing non-core plugins with a local build of Nextflow.

If you need to test changes to the `nextflow` launcher script, you can run it directly as `./nextflow`, or you can install it using `cp nextflow $(which nextflow)` and then run it as `nextflow`.

## Debugging

### Groovy REPL

The `groovysh` command provides a command-line REPL that you can use to play around with Groovy code independently of Nextflow. The `groovyConsole` command provides a graphical REPL similar to `nextflow console`. These commands require a standalone Groovy distribution, which can be installed as described for Java in {ref}`Getting started <getstarted-requirement>`.

:::{note}
If you are using WSL, you must also install an X server for Windows, such as [VcXsrv](https://sourceforge.net/projects/vcxsrv/) or [Xming](http://www.straightrunning.com/XmingNotes/), in order to use these commands.
:::

### IntelliJ IDEA

:::{versionadded} 23.09.0-edge
:::

You can perform limited breakpoint debugging on a Nextflow script using IntelliJ IDEA.

1. Set a breakpoint in your Nextflow script by clicking on a line number.

2. Run `nextflow -remote-debug run <script>`

3. Select the **Run / Debug Configurations** dropdown, select **Edit Configurations...**, and create a new configuration of type **Remote JVM Debug**. Set the port that appeared in the terminal when you launched your Nextflow script. Click **OK**.

4. Select the green bug icon to begin the remote debug session. The Debug window will appear and allow you to step through and inspect your script as it runs.

Note that this approach can only be used to debug the *script* execution, which does not include the *pipeline* execution.

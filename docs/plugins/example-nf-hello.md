(nf-hello-page)=

# Example: nf-hello

[nf-hello](https://github.com/nextflow-io/nf-hello) is a simple Nextflow plugin that uses the Nextflow Gradle plugin and provides examples for various plugin extension points.

The `nf-hello` plugin has the following structure:

```
nf-hello
├── COPYING
├── Makefile
├── README.md
├── build.gradle
├── gradle
│   └── wrapper
│       ├── gradle-wrapper.jar
│       └── gradle-wrapper.properties
├── gradlew
├── settings.gradle
└── src
    ├── main
    │   └── groovy
    │       └── nextflow
    │           └── hello
    │               ├── HelloConfig.groovy
    │               ├── HelloExtension.groovy
    │               ├── HelloFactory.groovy
    │               ├── HelloObserver.groovy
    │               └── HelloPlugin.groovy
    └── test
        └── groovy
            └── nextflow
                └── hello
                    ├── HelloDslTest.groovy
                    └── HelloFactoryTest.groovy
```

It includes the following extensions:

- A custom `hello` config scope (see `HelloConfig`).
- A custom trace observer that prints a message when the workflow starts and when the workflow completes (see `HelloObserver`).
- A custom channel factory called `reverse` (see `HelloExtension`).
- A custom operator called `goodbye` (see `HelloExtension`).
- A custom function called `randomString` (see `HelloExtension`).

The `nf-hello` plugin can be configured via nextflow configuration files or at runtime. For example:

```bash
nextflow run hello -plugins nf-hello@0.5.0
```

See the [nf-hello plugin repository](https://github.com/nextflow-io/nf-hello) for the plugin source code.

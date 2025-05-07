(nf-hello-page)=

# Example: nf-hello

[`nf-hello`](https://github.com/nextflow-io/nf-hello/tree/gradle-plugin-example) is a simple Nextflow plugin that uses the Gradle plugin for Nextflow plugins and is commonly used as a starting point for third-party plugin development.

The [`nf-hello` plugin](https://github.com/nextflow-io/nf-hello/tree/gradle-plugin-example) has the following structure:

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

It includes examples of different plugin extensions:

- A custom trace observer that prints a message when the workflow starts and when the workflow completes.
- A custom channel factory called reverse.
- A custom operator called goodbye.
- A custom function called randomString.

It also includes several classes that demonstrate different plugin functionality:

- `HelloConfig`: An example of how to handle options from the Nextflow configuration.
- `HelloExtension`: An example of how to create custom channel factories, operators, and functions that can be included in pipeline scripts.
- `HelloFactory`: An example of a workflow event with custom behavior.
- `HelloObserver`: An example of a workflow event with custom behavior.
- `HelloPlugin`: An example of a plugin entry point.

The `nf-hello` plugin can be configured via nextflow configuration files or at runtime. For example:

```bash
nextflow run hello -plugins nf-hello@0.5.0
```

See the [nf-hello plugin repository](https://github.com/nextflow-io/nf-hello/tree/gradle-plugin-example) for the plugin source code.

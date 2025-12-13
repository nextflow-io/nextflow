# Groovy console plugin for Nextflow

This plugin provides an interactive console and graphical user interface for Nextflow development, testing, and workflow exploration.

## Get Started

To use this plugin, add it to your `nextflow.config`:

```groovy
plugins {
    id 'nf-console'
}
```

Launch the console:

```bash
nextflow console
```

This opens a Groovy console with the Nextflow runtime pre-loaded, allowing you to interactively develop and test workflow components.

## Examples

### Launching the Console

```bash
nextflow console
```

### Loading a Script in the Console

```bash
nextflow console my-script.nf
```

### Interactive Testing

In the console, you can test Nextflow constructs:

```groovy
// Test a channel operation
Channel.of(1, 2, 3)
    .map { it * 2 }
    .view()

// Test file operations
file('data.txt').text = 'Hello World'

// Test process inputs/outputs
def reads = Channel.fromPath('*.fastq')
reads.view()
```

### Exploring the Runtime

```groovy
// Access workflow metadata
println workflow.projectDir
println workflow.launchDir

// Check available executors
println executorFactory.getExecutorNames()
```

## Resources

- [Nextflow Console Documentation](https://nextflow.io/docs/latest/cli.html#console)
- [Groovy Console Guide](https://groovy-lang.org/groovyconsole.html)

## License

[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0)

(structure-page)=

# Structure

## The `templates` directory

The `templates` directory in the Nextflow project root can be used to store scripts.

```
├── templates
│   └── sayhello.py
└── main.nf
```

It allows custom scripts to be invoked like regular scripts from any process in your pipeline using the `template` function:

```
process sayHello {
    
    input:
    val x

    output:
    stdout

    script:
    template 'sayhello.py'
}

workflow {
    Channel.of("Foo") | sayHello | view
}
```

Variables prefixed with the dollar character (`$`) are interpreted as Nextflow variables when the template script is executed by Nextflow:

```
#!/usr/bin/env python

print("Hello ${x}!")
```

The pipeline will fail if a template variable is missing, regardless of where it occurs in the template.

Templates can be tested independently of pipeline execution by providing each input as an environment variable. For example:

```bash
STR='foo' bash templates/my_script.sh
```

Template scripts are only recommended for Bash scripts. Languages that do not prefix variables with `$` (e.g. Python and R) can't be executed directly as a template script from the command line as variables prefixed with `$` are interpreted as Bash variables. Similarly, template variables escaped with `\$` will be interpreted as Bash variables when executed by Nextflow but not the command line.

:::{warning}
Template variables are evaluated even if they are commented out in the template script.
:::

:::{tip}
The best practice for using a custom script is to first embed it in the process definition and transfer it to a separate file with its own command line interface once the code matures.
:::

(bundling-executables)=

## The `bin` directory

The `bin` directory in the Nextflow project root can be used to store executable scripts.

```
├── bin
│   └── sayhello.py
└── main.nf
```

It allows custom scripts to be invoked like regular commands from any process in your pipeline without modifying the `PATH` environment variable or using an absolute path. Each script should include a shebang to specify the interpreter. Inputs should be supplied as arguments.

```python
#!/usr/bin/env python

import argparse

def main():
    parser = argparse.ArgumentParser(description="A simple argparse example.")
    parser.add_argument("name", type=str, help="Person to greet.")
    
    args = parser.parse_args()
    print(f"Hello {args.name}!")

if __name__ == "__main__":
    main()
```

:::{tip}
Use `env` to resolve the interpreter's location instead of hard-coding the interpreter path.
:::

Scripts placed in the `bin` directory must have executable permissions. Use `chmod` to grant the required permissions. For example:

```
chmod a+x bin/sayhello.py
```

Like modifying a process script, changing the executable script will cause the task to be re-executed on a resumed run.

:::{warning}
When using containers and the Wave service, Nextflow will send the project-level `bin` directory to the Wave service for inclusion as a layer in the container. Any changes to scripts in the `bin` directory will change the layer md5sum and the hash for the final container. The container identity is a component of the task hash calculation and will force re-calculation of all tasks in the workflow.

When using the Wave service, use module-specific bin directories instead. See {ref}`module-binaries` for more information.
:::

## The `lib` directory

The `lib` directory can be used to add utility code or external libraries without cluttering the pipeline scripts. The `lib` directory in the Nextflow project root is added to the classpath by default.

```
├── lib
│   └── DNASequence.groovy
└── main.nf
```

Classes or packages defined in the `lib` directory will be available in the execution context. Scripts or functions defined outside of classes will not be available in the execution context.

For example, `lib/DNASequence.groovy` defines the `DNASequence` class:

```groovy
// lib/DNASequence.groovy
class DNASequence {
    String sequence

    // Constructor
    DNASequence(String sequence) {
        this.sequence = sequence.toUpperCase() // Ensure sequence is in uppercase for consistency
    }

    // Method to calculate melting temperature using the Wallace rule
    double getMeltingTemperature() {
        int g_count = sequence.count('G')
        int c_count = sequence.count('C')
        int a_count = sequence.count('A')
        int t_count = sequence.count('T')

        // Wallace rule calculation
        double tm = 4 * (g_count + c_count) + 2 * (a_count + t_count)
        return tm
    }

    String toString() {
        return "DNA[$sequence]"
    }
}
```

The `DNASequence` class is available in the execution context:

```nextflow
// main.nf
workflow {
    Channel.of('ACGTTGCAATGCCGTA', 'GCGTACGGTACGTTAC')
    .map { seq -> new DNASequence(seq) }
    .view { dna -> 
        def meltTemp = dna.getMeltingTemperature()
        "Found sequence '$dna' with melting temperature ${meltTemp}°C" 
    }
}
```

It returns:

```
Found sequence 'DNA[ACGTTGCAATGCCGTA]' with melting temperaure 48.0°C
Found sequence 'DNA[GCGTACGGTACGTTAC]' with melting temperaure 50.0°C
```

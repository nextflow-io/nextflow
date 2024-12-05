(structure-page)=

# Project structure

(structure-templates)=

## The `templates` directory

The `templates` directory in the Nextflow project root can be used to store template files.

```
├── templates
│   └── sayhello.sh
└── main.nf
```

Template files can be invoked like regular scripts from any process in your pipeline using the `template` function. Variables prefixed with the dollar character (`$`) are interpreted as Nextflow variables when the template file is executed by Nextflow.

See {ref}`process-template` for more information about utilizing template files.

(structure-bin)=

## The `bin` directory

The `bin` directory in the Nextflow project root can be used to store executable scripts.

```
├── bin
│   └── sayhello.py
└── main.nf
```

The `bin` directory allows binary scripts to be invoked like regular commands from any process in your pipeline without using an absolute path of modifying the `PATH` environment variable. Each script should include a shebang to specify the interpreter and inputs should be supplied as arguments to the executable. For example:

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

Binary scripts placed in the `bin` directory must have executable permissions. Use `chmod` to grant the required permissions. For example:

```
chmod a+x bin/sayhello.py
```

Binary scripts in the `bin` directory can then be invoked like regular commands.

```
process sayHello {
    
    input:
    val x

    output:
    stdout

    script:
    """
    sayhello.py --name $x
    """
}

workflow {
    Channel.of("Foo") | sayHello | view
}
```

Like modifying a process script, modifying the binary script will cause the task to be re-executed on a resumed run.

:::{note}
Binary scripts require a local or shared file system for the pipeline work directory or {ref}`wave-page` when using cloud-based executors.
:::

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

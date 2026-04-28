(process-typed-page)=

# Processes (typed)

:::{versionadded} 25.10.0
:::

:::{warning}
Typed processes are a preview feature. The syntax and behavior may change in future releases.
:::

Typed processes use a new syntax for inputs and outputs that supports static typing.

```nextflow
nextflow.preview.types = true

process hello {
    input:
    message: String

    output:
    file('hello.txt')

    script:
    """
    echo '${message}' > hello.txt
    """
}
```

To use this feature:

1. Enable the {ref}`strict syntax <strict-syntax-page>` by setting the `NXF_SYNTAX_PARSER` environment variable to `v2`:

    ```bash
    export NXF_SYNTAX_PARSER=v2
    ```

2. Set `nextflow.preview.types = true` in every script that uses typed processes.

See {ref}`syntax-process-typed` for the complete syntax reference and {ref}`migrating-static-types` to migrate existing code to static typing.

## Inputs

The `input:` section declares process inputs. In typed processes, each input declaration consists of a name and type:

```nextflow
process fastqc {
    input:
    meta: Map
    fastq: Path
    extra_args: String

    script:
    """
    echo 'meta: ${meta}'
    echo 'fastq: ${fastq}'
    echo 'extra_args: ${extra_args}'
    """
}
```

All {ref}`standard types <stdlib-types>` except for the dataflow types (`Channel` and `Value`) can be used as type annotations in processes.

### File inputs

Nextflow automatically stages `Path` inputs and `Path` collections (such as `Set<Path>`) into the task directory.

### Nullable inputs

By default, tasks fail if any input receives a `null` value. To allow `null` values, add `?` to the type annotation:

```nextflow
process cat_opt {
    input:
    input: Path?

    stage:
    stageAs input, 'input.txt'

    output:
    stdout()

    script:
    '''
    [[ -f input.txt ]] && cat input.txt || echo 'empty input'
    '''
}
```

### Record inputs

Record inputs can be declared using a record type:

```nextflow
process fastqc {
    input:
    sample: Sample

    script:
    """
    echo 'id: ${sample.id}'
    echo 'fastq: ${sample.fastq}'
    """
}

record Sample {
    id: String
    fastq: Path
}
```

In this example, the record input is staged as `sample`, and `sample.fastq` is staged as an input file since it is declared with type `Path` in the `Sample` record type. Each field in the record type is staged into the task the same way as an individual input.

When the process is invoked, the incoming record should contain the specified fields, or else the run will fail. If the incoming record has additional fields not declared by the process input, they are ignored.

Record inputs can also be declared as a *destructured* input:

```nextflow
process fastqc {
    input:
    record(
        id: String,
        fastq: Path
    )

    script:
    """
    echo 'id: ${id}'
    echo 'fastq: ${fastq}'
    """
}
```

This pattern mirrors the standard `record()` function used to construct records. In this example, `fastq` is staged as an input file since the `fastq` field is declared with type `Path`.

:::{tip}
Record inputs are a useful way to select a subset of fields from a larger record. This way, the process stages only what it needs, keeping related data together in your workflow logic.
:::

### Tuple inputs

Tuple inputs can be declared as a *destructured* input:

```nextflow
process fastqc {
    input:
    tuple(id: String, fastq: Path)

    script:
    """
    echo 'id: ${id}'
    echo 'fastq: ${fastq}'
    """
}
```

This pattern mirrors the standard `tuple()` function used to construct tuples. Each tuple component is staged into the task the same way as an individual input.

## Stage directives

The `stage:` section defines custom staging behavior using *stage directives*.  It should be specified after the `input:` section. These directives serve the same purpose as input qualifiers such as `env` and `stdin` in the legacy syntax. 

### Environment variables

The `env` directive declares an environment variable in terms of task inputs:

```nextflow
process echo_env {
    input:
    hello: String

    stage:
    env 'HELLO', hello

    script:
    '''
    echo "$HELLO world!"
    '''
}
```

### Standard input (stdin)

The `stdin` directive defines the standard input of the task script:

```nextflow
process cat {
    input:
    message: String

    stage:
    stdin message

    script:
    """
    cat -
    """
}
```

### Custom file staging

The `stageAs` directive stages an input file (or files) under a custom file pattern:

```nextflow
process blast {
    input:
    fasta: Path

    stage:
    stageAs fasta, 'query.fa'

    script:
    """
    blastp -query query.fa -db nr
    """
}
```

The file pattern can also reference task inputs:

```nextflow
process grep {
    input:
    id: String
    fasta: Path

    stage:
    stageAs fasta, "${id}.fa"

    script:
    """
    cat ${id}.fa | grep '>'
    """
}
```

:::{versionadded} 26.04.0
:::

Input files can be staged using a *staging closure* instead of a file pattern:

```nextflow
process ls {
    input:
    slice: Set<Path>

    stage:
    stageAs(slice) { file -> "${file.parent.name}/${file.name}" }

    script:
    """
    ls -1 */*.txt | sort
    """
}
```

The staging closure should define the stage name for a given input file.

See {ref}`process-reference-typed` for available stage directives.

## Outputs

The `output:` section declares the outputs of a typed process. Each output declaration consists of a name, an optional type, and an output value:

```nextflow
process echo {
    input:
    message: String

    output:
    out_env: String = env('MESSAGE')
    out_file: Path = file('message.txt')
    out_std: String = stdout()

    script:
    """
    export MESSAGE='${message}'

    echo \$MESSAGE > message.txt

    cat message.txt
    """
}
```

When there is only one output, the name can be omitted:

```nextflow
process echo {
    input:
    message: String

    output:
    stdout()

    script:
    """
    echo '${message}'
    """
}
```

See {ref}`process-reference-typed` for available output functions.

### File outputs

You can use the `file()` and `files()` functions in the `output:` section to get a single file or collection of files from the task directory.

By default, the `file()` function fails if the specified file is not present in the task directory. You can specify `optional: true` to allow missing files. The `file()` function returns `null` for missing files. For example:

```nextflow
process foo {
    output:
    file('output.txt', optional: true)

    script:
    """
    exit 0
    """
}
```

### Structured outputs

Whereas legacy process outputs could only be structured using specific qualifiers like `val` and `tuple`, typed process outputs are regular values.

The `record()` standard library function can be used to create a record:

```nextflow
process fastqc {
    input:
    record(
        id: String,
        fastq: Path
    )

    output:
    record(
        id: id,
        fastqc: file('fastqc_logs')
    )

    script:
    // ...
}
```

The `tuple()` standard library function can be used to create a tuple:

```nextflow
process fastqc {
    input:
    tuple(id: String, fastq: Path)

    output:
    tuple(id, file('fastqc_logs'))

    script:
    // ...
}
```

(process-typed-topics)=

## Topics

The `topic:` section emits values to {ref}`topic channels <channel-topic>`. A topic emission consists of an output value and a topic name:

```nextflow
process cat {
    input:
    message: Path

    output:
    stdout()

    topic:
    tuple('bash', eval('bash --version')) >> 'versions'
    tuple('cat', eval('cat --version')) >> 'versions'

    script:
    """
    cat ${message}
    """
}
```

Topic emissions can use the same {ref}`output functions <process-reference-typed>` as the `output:` section.

## Script

The `script:` and `exec:` sections behave the same way as {ref}`legacy processes <process-script>`.

## Stub

The `stub:` section behaves the same way as {ref}`legacy processes <process-stub>`.

## Directives

Directives behave the same way as {ref}`legacy processes <process-directives>`.

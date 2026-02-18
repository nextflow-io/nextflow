(process-typed-page)=

# Processes (typed)

:::{versionadded} 25.10.0
:::

:::{warning}
Typed processes are a preview feature. The syntax and behavior may change in future releases.
:::

Typed processes use a new syntax for inputs and outputs based on static types.

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

See {ref}`syntax-process-typed` for the complete syntax reference and {ref}`migrating-static-types` to migrate existing code to static types.

## Inputs

The `input:` section declares process inputs. In typed processes, each input declaration consists of a name and type:

```nextflow
process fastqc {
    input:
    (meta, fastq): Tuple<Map,Path>
    extra_args: String

    script:
    """
    echo 'meta: ${meta}`
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
    stageAs 'input.txt', input

    output:
    stdout()

    script:
    '''
    [[ -f input.txt ]] && cat input.txt || echo 'empty input'
    '''
}
```

### Record inputs

Inputs with type `Record` can declare the name and type of each record field:

```nextflow
process fastqc {
    input:
    sample: Record {
        id: String
        fastq: Path
    }

    script:
    """
    echo 'id: ${sample.id}`
    echo 'fastq: ${sample.fastq}'
    """
}
```

In this example, the record is staged into the task as `sample`, and `sample.fastq` is staged as an input file since the `fastq` field is declared with type `Path`.

When the process is invoked, the incoming record should contain the specified fields, or else the run will fail. If the record has additional fields not declared by the process input, they are ignored.

:::{tip}
Record inputs are a useful way to select a subset of fields from a larger record. This way, the process only stages what it needs, allowing you to keep related data together in your workflow logic.
:::

You can achieve the same behavior using an external record type:

```nextflow
process fastqc {
    input:
    sample: Sample

    script:
    """
    echo 'id: ${sample.id}`
    echo 'fastq: ${sample.fastq}'
    """
}

record Sample {
    id: String
    fastq: Path
}
```

This approach is useful when the record type can be re-used elsewhere in the pipeline.

### Tuple inputs

Inputs with type `Tuple` can declare the name of each tuple component:

```nextflow
process fastqc {
    input:
    (id, fastq): Tuple<String,Path>

    script:
    """
    echo 'id: ${id}`
    echo 'fastq: ${fastq}'
    """
}
```

This pattern is called *tuple destructuring*. Each tuple component is staged into the task the same way as an individual input.

The generic types inside the `Tuple<...>` annotation specify the type of each tuple compomnent and should match the component names. In the above example, `id` has type `String` and `fastq` has type `Path`.

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
    stageAs 'query.fa', fasta

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
    stageAs "${id}.fa", fasta

    script:
    """
    cat ${id}.fa | grep '>'
    """
}
```

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
    sample: Record {
        id: String
        fastq: Path
    }

    output:
    record(id: sample.id, fastqc: file('fastqc_logs'))

    script:
    // ...
}
```

The `tuple()` standard library function can be used to create a tuple:

```nextflow
process fastqc {
    input:
    (id, fastq): Tuple<String,Path>

    output:
    tuple(id, file('fastqc_logs'))

    script:
    // ...
}
```

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

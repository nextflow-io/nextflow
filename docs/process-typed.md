(process-typed-page)=

# Processes (typed)

:::{versionadded} 25.10.0
:::

:::{note}
This feature requires the {ref}`strict syntax <strict-syntax-page>` to be enabled (`NXF_SYNTAX_PARSER=v2`).
:::

:::{note}
Typed processes require the `nextflow.preview.types` feature flag to be enabled in every script that uses them.
:::

Typed processes use a new syntax for inputs and outputs that is based on static types:

```nextflow
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

See {ref}`syntax-process-typed` for a full description of the process syntax. See {ref}`migrating-static-types` for more information on migrating existing code to static types.

## Inputs

The `input:` section is used to declare the inputs of a process. An input declaration in a typed process consists of a name and a type:

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

Any of the {ref}`standard types <stdlib-types>` can be used as type annotations, except for the dataflow types (`Channel` and `Value`) which can only be used in workflows.

### File inputs

Inputs of type `Path` or a collection of `Path` (e.g. `Set<Path>`) are automatically staged into the task directory.

By default, the task will fail if any input receives a `null` value. You can mark an input as nullable by appending `?` to the type annotation:

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

### Stage directives

The `stage:` section can be specified after the `input:` section. You can use it to specify custom staging behavior using *stage directives*. These directives serve the same purpose as input qualifiers such as `env` and `stdin` in the legacy syntax.

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

See {ref}`process-reference-typed` for the set of available stage directives.

## Outputs

The `output:` section is used to declare the outputs of a typed process. An output declaration in a typed process consists of a name, an optional type, and an output value:

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

See {ref}`process-reference-typed` for the set of available output functions.

### File outputs

You can use the `file()` and `files()` functions in the `output:` section to get a single file or collection of files from the task directory.

By default, the `file()` function will fail if the specified file is not present in the task directory. You can specify `optional: true` to allow the file to be missing, in which case the `file()` function will return `null`. For example:

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

## Topics

The `topic:` section is used to emit values to a {ref}`topic channel <channel-topic>`. A topic emission consists of an output value and a topic name:

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

Topic emissions can use the same {ref}`output functions <process-reference-typed>` that are available in the `output:` section.

## Script

The `script:` and `exec:` sections behave the same way as {ref}`legacy processes <process-script>`.

## Stub

The `stub:` section behaves the same way as {ref}`legacy processes <process-stub>`.

## Directives

Directives behave the same way as {ref}`legacy processes <process-directives>`.

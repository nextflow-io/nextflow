# Typed processes

- Authors: Ben Sherman
- Status: accepted
- Deciders: Ben Sherman, Paolo Di Tommaso
- Date: 2025-10-17
- Tags: lang, static-types, processes

## Updates

### Version 1.1 (2026-03-23)

- Changed the method signature for `stageAs` from `(filePattern, value)` to `(value, filePattern)` to mirror commands like `cp`, `mv`, etc.

- Replaced annotation-based tuple syntax (`(...): Tuple<...>`) with destructuring syntax (`record(...)`) for better continuity with legacy syntax and record input syntax.

## Summary

Support static typing for process inputs and outputs.

## Problem Statement

The legacy process syntax uses *qualifiers* to describe both the type and staging behavior of each input and output:

```groovy
process FASTQC {
    input:
    tuple val(id), path(fastq_1), path(fastq_2)
    path index

    output:
    tuple val(id), path("fastqc_${id}_logs")

    script:
    // ...
}
```

This syntax has several limimtations:

- **No static typing**: The `val` qualifier can not specify a type, so there is no way to validate input values. The `path` qualifier can not distinguish between a file and a file collection. The `arity` option was introduced to address this ambiguity, but it is cumbersome and rarely used.

- **Type and staging behavior are coupled**: Qualifiers like `path` describe both the type *and* the staging behavior (link into task directory). There is no way to specify staging behavior separately, such as staging a tuple element or record field as an environment variable.

- **No nullability**: There is no way to declare that an input may be `null`. The `path` qualifier raises a runtime error if a null value is received. Outputs can be marked optional, but optional outputs are handled by emitting nothing rather than emitting `null`. A tuple output can be optional, but a tuple element can not.

- **Limited output expressiveness**: Outputs must be expressed in terms of qualifiers that mirror the input qualifiers. It is difficult to express many kinds of output values, and it is unclear to the user whether a given expression is valid or not.

## Goals

- Provide a way to model process inputs and outputs with types from the Nextflow standard library.

- Separate the *type* of an input from its *staging behavior*.

- Provide first-class support for nullable inputs and outputs.

- Allow outputs to be arbitrary expressions, ensuring consistency with the rest of the language.

- Enable compile-time type checking for processes.

## Non-goals

- Removing the legacy qualifier syntax -- legacy processes must continue to work without modification.

- Enforcing type checking -- static type checking will be introduced progressively as an opt-in feature.

## Decision

Introduce **typed processes**, which use a new syntax for inputs and outputs based on type annotations instead of qualifiers.

All other process sections (directives, script, stub, etc) are supported by typed processes without changes. Only the `input:` and `output:` sections are changed.

## Core Capabilities

### Typed inputs

Each input is declared as `name: Type`:

```groovy
process FASTQC {
    input:
    id: String
    fastq_1: Path
    fastq_2: Path

    script:
    """
    fastqc -o fastqc_${id}_logs ${fastq_1} ${fastq_2}
    """
}
```

All standard library types except `Channel` and `Value` are valid input types. Inputs of type `Path` (or `Path` collections such as `Set<Path>`) are automatically staged into the task directory.

### Nullable inputs

Appending `?` to a type annotation allows the input to be `null`:

```groovy
process CAT_OPT {
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

By default, a task fails if any input receives `null`.

### Stage directives

Staging behavior is moved to a dedicated `stage:` section that appears after `input:`. This replaces the staging aspects of legacy qualifiers:

| Legacy qualifier  | Stage directive    |
|-------------------|--------------------|
| `env('NAME')`     | `env 'NAME', value` |
| `stdin`           | `stdin value`       |
| `path('name.fa')` | `stageAs file, 'name.fa'` |

For example:

```groovy
process BLAST {
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

Separating staging from type declaration keeps the inputs clean and makes it easier to specify staging behavior independently of the input type.

### Tuple inputs

Tuples are declared inline using `tuple(name: Type, ...)`:

```groovy
process FASTQC {
    input:
    tuple(id: String, fastq_1: Path, fastq_2: Path)

    script:
    // ...
}
```

Each component is destructured into a local variable. This mirrors the `tuple()` constructor used in the output section and in workflow logic, making the syntax consistent.

### Typed outputs

Each output declaration consists of an optional name and type, and a value expression:

```groovy
process ECHO {
    input:
    message: String

    output:
    out_file: Path = file('message.txt')
    out_std: String = stdout()

    script:
    """
    echo '${message}' | tee message.txt
    """
}
```

When there is only one output, the name and type can be omitted:

```groovy
process ECHO {
    input:
    message: String

    output:
    file('message.txt')

    script:
    """
    echo '${message}' > message.txt
    """
}
```

Outputs can be arbitrary expressions, rather that being restricted to specific qualifiers such as `tuple` and `val`. Special functions such as `file()`, `files()`, `env()`, and `stdout()` can be composed into the desired output structure.

### Nullable outputs

By default, the `file()` and `files()` function raise an error if the given file is missing. These functions can be called with `optionel: true` to allow missing files. This way, it is possible to declare a tuple output that contains nullable values:

```groovy
process MAYBE {
    input:
    id: String

    output:
    tuple(id, file('result.txt'))

    script:
    """
    [[ '$id' == 42 ]] && touch result.txt
    """
}
```

### Topic emissions

A `topic:` section emits values to topic channels using the `>>` operator:

```groovy
process CAT {
    input:
    message: Path

    output:
    stdout()

    topic:
    tuple('cat', eval('cat --version')) >> 'versions'

    script:
    """
    cat ${message}
    """
}
```

Moving topic emissions to a dedicated section allows them to be defined without having to include them in the process outputs.

## Distinguishing between typed and legacy processes

Typed processes are gated behind the `nextflow.enable.types` feature flag, in order to distinguish between typed and legacy processes in the language.

When a script enables this feature flag, its processes are treated as typed processes; otherwise, its processes are treated as legacy processes. This way, typed and legacy processes cannot be mixed in the same script, but they can be used together as long as they are declared in different scripts.

While typed and legacy processes are syntactically distinct and could theoretically be allowed in the same script, the feature flag helps distinguish typed vs legacy to the reader (whether human or agent).

## Alternatives

### Implicit tuple input

The syntax for typed process inputs aims to be consistent with typed inputs throughout the rest of the language, such as the `params` block and workflow inputs, which use the pattern of `<name>: <type>`. The `tuple` input qualifier does not fit neatly into this pattern, since it specifies multiple tuple *components*:

```groovy
process QUANT {
    input:
    tuple(id: String, fastq_1: Path, fastq_2: Path)
    index: Path

    // ...
}

workflow {
    ch_samples = channel.of( tuple('1', file('1_1.fq'), file('1_2.fq')) )
    index = file('index.fa')
    QUANT(ch_samples, index)
}
```

One alternative is to remove tuple inputs altogether and treat the entire `input:` section as an implicit tuple input:

```groovy
process QUANT {
    input:
    id: String
    fastq_1: Path
    fastq_2: Path
    index: Path

    // ...
}

workflow {
    ch_samples = channel.of( tuple('1', file('1_1.fq'), file('1_2.fq')) )
    index = file('index.fa')
    QUANT( ch_samples.combine(index) )
}
```

With this approach, a process would always be called with a single input, and multiple sources (e.g. `ch_samples` and `index`) would need to be combined into a single input. This could be done explicitly with the `combine` operator or implicitly by the runtime.

However, this approach would be a significant change to process call semantics, even if only applied to typed processes. It would likely be difficult to validate for processes with many inputs.

The tuple destructuring syntax makes it possible to migrate legacy processes to typed processes without changing workflow logic or call semantics. While the `tuple(...)` syntax is a deviation from the typed input syntax used by the rest of the language, such deviations can be appropriate and even advantageous when used judiciously in a custom language.

### Type annotation syntax for tuple inputs

Another alternative for tuple inputs is to use a type annotation:

```groovy
process FASTQC {
    input:
    (id, fastq_1, fastq_2): Tuple<String,Path,Path>
    index: Path

    // ...
}
```

This approach attempts to bring the syntax closer to the `<name>: <type>` pattern while maintaining support for tuple destructuring. This syntax was used in the first preview of typed processes in Nextflow 25.10.

However, this syntax needlessly separates the component name from its corresponding type, making it harder to read and validate. Although it is semantically equivalent to the legacy syntax, it looks and feels very different, which can be jarring for users.

With the introduction of records, the `tuple(...)` destructuring syntax emerged as a clear pattern to follow for both records and tuples:

**Legacy process:**
```groovy
process FASTQC {
    input:
    tuple val(id), path(fastq_1), path(fastq_2)

    output:
    tuple val(id), path("fastqc_${id}_logs")

    // ...
}
```

**Typed process (tuple):**
```groovy
process FASTQC {
    input:
    tuple(id: String, fastq_1: Path, fastq_2: Path)

    output:
    tuple(id, file("fastqc_${id}_logs"))

    // ...
}
```

**Typed process (record):**
```groovy
process FASTQC {
    input:
    record(
        id: String,
        fastq_1: Path,
        fastq_2: Path
    )

    output:
    record(
        id: id,
        fastqc: file("fastqc_${id}_logs")
    )

    // ...
}
```

This pattern provides the best balance of continuity with the old way and consistency with static typing:

- A legacy process can be migrated to a typed process by replacing the `tuple` input/output qualifier with the `tuple` destructor/constructor.
- A typed process can be migrated from tuples to records by replacing `tuple` with `record` and adding fields to the record output.
- The `tuple` and `record` destructors use the same `<name>: <type>` pattern used by the rest of the language.
- At each stage, the inputs and outputs mirror each other without creating syntactic confusion.

## Consequences

**Positive:**

- Type annotations make processes self-documenting and provide the information needed to perform static type checking.

- Separating type from staging behavior (the `stage:` section) makes each concern independently clear.

- Nullable types (`?`) provide first-class support for nullable input files.

- Outputs can be structured arbitrarily and can contain nullable files.

**Negative:**

- The `each` qualifier is not supported; pipelines using it must be refactored to use the `combine` operator before migrating to typed processes.

- The typed syntax must be maintained alongside the legacy syntax, which makes the codebase more complex and may cause confusion as the community transitions to the new syntax.

**Neutral:**

- Typed processes use the same standard types as the rest of the language, so no additional type vocabulary is introduced.

- Typed processes are enabled by a feature flag, which introduces new functionality without breaking existing code and helps distinguish between typed and legacy code.

## Links

- [Nextflow standard types](https://nextflow.io/docs/latest/reference/stdlib-types.html)
- Community issues: #1694, #2678
- Related nf-core discussion: https://github.com/nf-core/modules/issues/4311
- Original implementation: #4553

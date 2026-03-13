# Record types

- Authors: Ben Sherman
- Status: accepted
- Deciders: Ben Sherman, Paolo Di Tommaso
- Date: 2026-03-06
- Tags: lang, static-types

## Summary

Provide a way to model composite data types in the Nextflow language.

## Problem Statement

Nextflow pipelines typically need a way to model *aggregates*, or "data that travels together", such as a paired-end read consisting of a sample ID and two FASTQ files.

Primary use cases:

- Model complex pipeline inputs and outputs (e.g. samplesheets as collections of records)

- Generate JSON schemas for pipeline inputs and outputs from source code (e.g. to facilitate pipeline chaining with external validation)

- Model directory outputs as records, enabling more fine-grained validation

## Goals

- Introduce a alternative data structure to tuples that allows users to model their data domain more precisely.

- Introduce a way to validate custom data types at compile-time while keeping pipeline code concise and readable.

## Non-goals

- Removing support for tuples -- tuples should continue to work as before

- Introducing new dataflow operators or changing existing ones -- any changes to operators will be addressed in future efforts as needed

- Introducing type inheritance -- the Nextflow type system avoids type inheritance as much as possible in order to not introduce unnecessary complexity

- Introducing object methods -- this can be handled by standalone functions for now, and may be improved in a future effort (e.g. namespaces)

## Decision

Add support for **records**, an alternative data structure to tuples, and **record types**, a way to declare custom types that can be applied to records a la carte.

## Core Capabilities

### Records

Records are an attempt to take the best qualities of maps and custom classes while avoiding the downsides.

A record can be created using the `record()` function:

```groovy
sample = record(
    id: '1',
    fastq_1: file('1_1.fastq'),
    fastq_2: file('1_2.fastq')
)

sample.id = '2' // error: record cannot be modified
sample += record(id: '2') // ok

println sample.id
```

This function effectively creates an immutable map (`Map<String,?>`):

- The keys are just field names
- The values can have any type
- The record can’t be modified -- use the `+` operator instead

Records can have arbitrary fields, unlike custom classes, which makes them easy to use with dataflow operators.

For example, a future version of the `join` operator could join records as follows:

```groovy
ch_bam = channel.of( record(id: '1', bam: file('1.bam')) )
ch_bai = channel.of( record(id: '1', bai: file('1.bai')) )

ch_bam.join(ch_bai, by: 'id').view()

// -> record(id: '1', bam: file('1.bam'), bai: file('1.bai'))
```

Whereas tuples are joined by a matching index, records would be joined by a matching key.

### Record types

A record type is a user-defined type that consists of a name and a set of fields:

```groovy
record Sample {
    id: String
    fastq_1: Path
    fastq_2: Path?
}
```

Fields in a record type are declared the same way as [typed parameters](https://nextflow.io/docs/latest/workflow.html#typed-parameters). All [standard types](https://nextflow.io/docs/latest/reference/stdlib-types.html) can be used. Fields can be marked as optional by appending a `?` to the field type.

The purpose of a record type is to specify a *minimum set of requirements* for a record *in a particular context*. A record created with the `record()` function simply has the type `Record` -- it makes no guarantees about which fields it provides. A *record type* can be used (e.g. in a workflow input) to make a stronger guarantee.

For example:

```groovy
workflow RNASEQ {
    take:
    samples: Channel<Sample>

    main:
    // ...
}

workflow {
    ch_samples = channel.of( record(id: '1', fastq_2: file('1_2.fastq')) )
    RNASEQ(ch_samples) // error: `ch_samples` is missing `fastq_1` field required by Sample
}
```

This workflow definition specifies that the `samples` input should be a channel of records, where each record has at minimum the fields specified by the `Sample` record type. The records can still have additional fields, but only the fields in `Sample` are guaranteed to be present.

In other words, records are *duck-typed*. Duck-typing semantics are used whenever a record is validated against a record type:

- Supplying a record argument to a workflow, process, or function with a record type input (as shown above)

- Casting a record to a record type (e.g. `record(...) as Sample`)

Record types can be included across modules:

```groovy
include { Sample } from './module'
```

Because of duck-typing, two record types with the same fields and field types are effectively equivalent, even if they have different names:

```groovy
record Sample {
    id: String
    fastq_1: Path
    fastq_2: Path?
}

record FastqPair {
    id: String
    fastq_1: Path
    fastq_2: Path?
}
```

This makes it easier to compose modules and workflows that use their own record types.

### Process inputs

When a record is supplied as input to a process, the process needs to know how to stage input files from the record, like it does with the `path` qualifier.

Typed processes can stage inputs using the `stage:` section, but ideally the files in a record should be automatically detected and staged.

A typed process can declare a record using an *inline record type*:

```groovy
process FASTQC {
    input:
    sample: Record {
        id: String
        fastq_1: Path
        fastq_2: Path
    }

    // ...
}
```

All record fields that are a `Path` or `Path` collection (e.g. `Set<Path>`) are automatically staged. The record itself is declared in the process body as `sample`, like any other input, and record fields are accessed as `sample.id`, `sample.fastq_1`, and so on.

A typed process can also use an explicit record type to achieve the same behavior:

```groovy
process FASTQC {
    input:
    sample: FastqPair

    // ...
}

record FastqPair {
    id: String
    fastq_1: Path
    fastq_2: Path
}
```

The only difference between these two aprooaches is that the `FastqPair` type can be used elsewhere in pipeline code because it is declared externally.

### Process outputs

Typed processes can declare outputs with arbitrary expressions, so no new syntax is required to support record outputs. Simply use the `record()` function to create a record:

```groovy
process FASTQC {
  // ...

  output:
  record(id: id, fastqc: file('fastqc_logs'))

  // ...
}
```

The type of this output is an *implicit* record type that is inferred from the code: `Record { id: String ; fastqc: Path }`.

## Alternatives

### Custom classes

Define Groovy-style classes (see #2085) and use them to model composite data:

```groovy
@nextflow.io.ValueObject
class Sample {
    String id
    Path fastq_1
    Path fastq_2
}

workflow {
    sample = new Sample('1', file('1_1.fastq'), file('1_2.fastq'))

    println sample.id
}
```

This approach can be used, but in practice it requires a lot of extra dataflow logic around process calls to convert between custom types and tuples, because processes don’t know how to stage input files from custom types.

Custom classes are not very flexible. For example, joining two channels of custom classes would be more complicated than joining two tuples by a matching key, because you would need to define an additional class for the “joined” type and explicitly construct it from the two joining classes.

### Maps

Use maps to create composite data structures dynamically (see #2127):

```groovy
sample = [
    id: '1',
    fastq_1: file('1_1.fastq'),
    fastq_2: file('1_2.fastq')
]

println sample.id
```

Maps are flexible because you can store arbitrary fields rather than being restricted to a fixed set of fields. However, maps are meant to be used for a single value type (e.g. `Map<String,Integer>`).

Unlike tuples, maps are mutable (i.e. they can be modified). Modifying maps can lead to race conditions if done improperly. As a best practice, maps should be modified by adding another map, which creates a copy:

```groovy
sample2 = sample + [id: '2']

println sample.id // -> '1'
println sample2.id // -> '2'
```

### Implicit process record output

A process record output can be defined using the `record()` function as shown above:

```groovy
process FASTQC {
    // ...

    output:
    record(
        id: id,
        fastqc: file('fastqc_logs')
    )

    // ...
}
```

One alternative is to re-interpret the existing typed output syntax as an implicit record, treating each line as a record field:

```groovy
process FASTQC {
    // ...

    output:
    id: String = id
    fastqc: Path = file('fastqc_logs')

    // ...
}
```

This approach is syntactically more concise, and it re-uses the typed output syntax that was introduced in Nextflow 25.10.

However, with this approach, the same syntax can have different meanings depending on the surrounding context (e.g. presence/absence of the `nextflow.preview.types` feature flag), which can be confusing for both users and agents.

The `record()` approach works "out of the box", and it isn't much more verbose, so we decided that it is sufficient for now.

## Consequences

**Positive:**

- Records replace positional elements with named fields, making data structures self-documenting and less error-prone (e.g. no more accidentally swapping `fastq_1` and `fastq_2` by position).

- Duck-typing makes module composition easier: downstream processes declare only the fields they need, and record types from different modules are interchangeable if they share the same fields.

- Immutability by default eliminates the race conditions that can occur when mutable maps are improperly used in workflow logic.

- A single record output replaces multiple per-file tuple output channels (as shown in the prokka example), reducing the total number of channels in a workflow.

- Record types provide a foundation for generating JSON schemas for pipeline inputs and outputs directly from source code, enabling external validation and pipeline chaining.

- Backward compatibility is fully preserved: existing tuple-based pipelines continue to work without modification.

**Negative:**

- Since records must match based on field name rather than element input, users must be careful to use consistent naming conventions or write additional adaptor logic in their workflow. Type checking for records will be essential to streamline the developer experience as much as possible.

- Tuples and records coexist as parallel data model options, which may cause confusion about which to use for a given situation. Guidelines will be needed to help users make the right choice.

- Dataflow operators such as `cross`, `groupTuple`, and `join` need to be updated to support records natively.

**Neutral:**

- Record types use structural (duck) typing rather than nominal typing. Two record types with identical fields are interchangeable regardless of their names. This is intentional and enables flexible module composition, but it differs from the nominal typing that most users encounter in other languages, so it may be surprising at first.

- Records have no methods; behavior must be expressed via standalone functions. This is consistent with the functional style of Nextflow pipelines, and may be improved in the future with namespaces.

## Links

- Community issues: #2085, #2127
- Related nf-core discussion: https://github.com/nf-core/modules/issues/4311
- Original implementation: #4553
- nf-core/fetchngs POC: https://github.com/nf-core/fetchngs/pull/309
- Inspired by: [Simple Made Easy](https://github.com/matthiasn/talk-transcripts/blob/master/Hickey_Rich/SimpleMadeEasy.md)
- Type systems: [Nominal typing](https://en.wikipedia.org/wiki/Nominal_type_system) vs [Structural typing](https://en.wikipedia.org/wiki/Structural_type_system) vs [Duck typing](https://en.wikipedia.org/wiki/Duck_typing)

## Appendix

### Example: nf-core/prokka

The [nf-core/prokka](https://github.com/nf-core/modules/blob/master/modules/nf-core/prokka/main.nf) module produces several output files, emitting a tuple channel for each output file. If a downstream process requires multiple outputs, the individual output channels must be joined as needed to match the process input tuple.

Here is how the process might look using records:

```groovy
process PROKKA {
    // ...

    input:
    sample: Record {
        meta: Map
        fasta: Path
    }
    proteins: Path
    prodigal_tf: Path

    output:
    record(
        meta: sample.meta,
        gff: file("${prefix}/*.gff"),
        gbk: file("${prefix}/*.gbk"),
        fna: file("${prefix}/*.fna"),
        faa: file("${prefix}/*.faa"),
        ffn: file("${prefix}/*.ffn"),
        sqn: file("${prefix}/*.sqn"),
        fsa: file("${prefix}/*.fsa"),
        tbl: file("${prefix}/*.tbl"),
        err: file("${prefix}/*.err"),
        log: file("${prefix}/*.log"),
        txt: file("${prefix}/*.txt"),
        tsv: file("${prefix}/*.tsv")
    )

    topic:
    file("versions.yml")  >> 'versions'

    script:
    prefix = sample.meta.id
    // ...
}
```

The tuple input is refactored as a record input with an inline record type. The tuple outputs are combined into a single record output. No external record types are needed, although they could be used if desired.

*NOTE:* The `meta` map has not been changed in this example for brevity. However, it could be modeled with a record type instead of the generic `Map` type, or it could even be replaced with explicit fields such as `id: String`.

Now, suppose there are two downstream processes that want to use the outputs of `PROKKA`:

1. Process `FOO` only needs the `gff` file
2. Process `BAR` only needs the `fna`, `faa`, and `tbl` files

These processes would be defined as follows:

```groovy
process FOO {

    input:
    sample: Record {
        meta: Map
        gff: Path
    }

    // ...
}

process BAR {

    input:
    sample: Record {
        meta: Map
        fna: Path
        faa: Path
        tbl: Path
    }

    // ...
}
```

And the calling workflow would be written as follows:

```groovy
workflow {
    ch_inputs = channel.of( /* ... */ )
    proteins = // ...
    prodigal_tf = // ...
    cn_prokka = PROKKA( ch_inputs, proteins, prodigal_tf )

    FOO(ch_prokka)
    BAR(ch_prokka)
}
```

Each process declares a record input containing only the fields that it needs. When the output of `PROKKA` is passed to `FOO` and `BAR`, each process stages only the files that it declared in the record input.

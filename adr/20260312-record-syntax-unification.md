# Unified record syntax for process inputs and outputs

- Authors: Paolo Di Tommaso
- Status: proposed
- Deciders: Paolo Di Tommaso, Ben Sherman
- Date: 2026-03-12
- Tags: lang, records, syntax

Technical Story: Follow-up to [Record types ADR](20260306-record-types.md)

## Summary

The current record types implementation uses two different syntactic forms for records in process inputs (block syntax) vs outputs (function-call syntax). This RFC proposes using the `record()` function-call notation uniformly for both inputs and outputs, combined with standard assignment and type annotations.

## Problem Statement

The accepted record types ADR ([20260306-record-types](20260306-record-types.md)) introduces two distinct syntactic forms for records within process definitions:

**Input** — a `Record { ... }` block syntax unique to inputs:
```nextflow
process FASTQC {
    input:
    sample: Record {
        id: String
        fastq_1: Path
        fastq_2: Path
    }
    ...
}
```

**Output** — a `record()` function call:
```nextflow
process FASTQC {
    ...
    output:
    record(
        id: sample.id,
        html: file('*.html'),
        zip: file('*.zip')
    )
}
```

This asymmetry means the same concept (a record) is expressed with two different syntactic forms depending on context. The block syntax `Record { ... }` exists only in process input declarations and has no counterpart elsewhere in the language. Meanwhile, the `record()` function call used in outputs is already a general-purpose construct usable in any expression context.

## Goals

- **Syntactic consistency** — use a single notation for records across inputs and outputs.
- **Alignment with existing syntax** — reuse assignment (`=`) and type annotation (`: Type`) patterns already present in process I/O, rather than introducing new block syntax.
- **Standard type semantics** — record assignments should follow the same type compatibility rules as any other typed assignment in the language.

## Non-goals

- Changing the top-level `record` type definition syntax — the `record Name { field: Type }` declaration form is a type-level construct and is not affected by this proposal.
- Changing the `record()` function runtime behavior or the `RecordMap` implementation.
- Removing support for external type references (e.g. `sample: Sample`).

## Considered Options

### Option 1: Current syntax (status quo)

Input uses a dedicated block syntax, output uses the `record()` function call:

```nextflow
process FASTQC {
    input:
    sample: Record {
        id: String
        fastq_1: Path
        fastq_2: Path
    }

    output:
    record(
        id: sample.id,
        html: file('*.html'),
        zip: file('*.zip')
    )
}
```

The output can optionally be an assignment with a type annotation:

```nextflow
    output:
    result: FastqcResult = record(
        id: sample.id,
        html: file('*.html'),
        zip: file('*.zip')
    )
```

- Good, because input block syntax mirrors the top-level `record` definition.
- Bad, because two different notations for the same concept in the same process definition.
- Bad, because `Record { ... }` block syntax only exists in input declarations — it is not a general-purpose construct.

### Option 2: Block syntax for both inputs and outputs

Use `record { ... }` blocks in both input and output:

```nextflow
process FASTQC {
    input:
    record sample {
        id: String
        fastq_1: Path
        fastq_2: Path
    }

    output:
    record {
        id: String = sample.id
        html: Path = file('*.html')
        zip: Path = file('*.zip')
    }
}
```

- Good, because symmetric — same block form on both sides.
- Bad, because the output block mixes type declarations with value assignments (`Path = file(...)`).
- Bad, because block syntax in process I/O diverges from the function-call style already established for `record()`.

### Option 3: Unified `record()` function notation with assignment

Use the `record()` function-call syntax for both inputs and outputs, with standard assignment:

```nextflow
process FASTQC {
    input:
    sample = record(id: String, fastq_1: Path, fastq_2: Path)

    output:
    result = record(id: sample.id, html: file('*.html'), zip: file('*.zip'))
}
```

With optional explicit type annotation:

```nextflow
process FASTQC {
    input:
    sample: Sample = record(id: String, fastq_1: Path, fastq_2: Path)

    output:
    result: QcResult = record(id: sample.id, html: file('*.html'), zip: file('*.zip'))
}
```

- Good, because same notation on both sides — `name = record(...)`.
- Good, because reuses existing assignment and type annotation patterns.
- Good, because `record()` is already a general-purpose function, no new syntax needed.
- Good, because type annotations follow standard rules — `sample: Sample = record(...)` works like any typed assignment.
- Bad, because input `record()` arguments are types rather than values, which is a different usage of the function.

## Solution or decision outcome

**Option 3**: Use the `record()` function-call notation uniformly for both process inputs and outputs, combined with standard assignment (`=`) and optional type annotation (`: Type`).

## Rationale & discussion

The key insight is that the `record()` call is just a constructor, and everything else is standard Nextflow assignment and type annotation. This eliminates the need for a dedicated `Record { ... }` block syntax in process inputs.

### Syntax pattern

The unified pattern is `name: Type = record(...)` for both inputs and outputs:

- **Input**: `sample = record(id: String, fastq_1: Path, fastq_2: Path)` — declares the fields and their types being received.
- **Output**: `result = record(id: sample.id, html: file('*.html'))` — declares the fields and their values being produced.

The only difference is what goes inside the `record()` call — types on input (declaring structure), expressions on output (producing values). This parallels how assignment works elsewhere: the left side declares, the right side provides.

### Type annotations

Type annotations are optional and follow standard semantics:

```nextflow
// Inferred type from record fields
sample = record(id: String, fastq_1: Path, fastq_2: Path)

// Explicit type — compiler checks compatibility with Sample
sample: Sample = record(id: String, fastq_1: Path, fastq_2: Path)
```

This is the same as writing `x: Integer = 42` vs `x = 42` — nothing record-specific about the assignment semantics.

### Alignment with existing process syntax

The proposed syntax reuses patterns that already exist in Nextflow process definitions:

| Existing pattern | Example | Record equivalent |
|-----------------|---------|-------------------|
| Assignment in output | `id = sample.id` | `result = record(...)` |
| Typed assignment in output | `id: String = sample.id` | `result: QcResult = record(...)` |
| Type annotation in input | `id: String` | `sample: Sample = record(...)` |

### External type reference

When using a pre-defined record type, the syntax naturally simplifies:

```nextflow
// With inline fields
sample: Sample = record(id: String, fastq_1: Path, fastq_2: Path)

// With external type only (no inline fields needed)
sample: Sample
```

The `sample: Sample` shorthand remains valid — the `record()` call is only needed when defining fields inline.

### Full example

```nextflow
nextflow.preview.types = true

record Sample {
    id: String
    fastq_1: Path
    fastq_2: Path
}

process TOUCH {
    input:
    id: String

    output:
    result = record(id: id, fastq_1: file('*_1.fastq'), fastq_2: file('*_2.fastq'))

    script:
    """
    touch ${id}_1.fastq
    touch ${id}_2.fastq
    """
}

process FASTQC {
    input:
    sample: Sample = record(id: String, fastq_1: Path, fastq_2: Path)

    output:
    result = record(id: sample.id, html: file('*.html'), zip: file('*.zip'))

    script:
    """
    touch ${sample.id}.html
    touch ${sample.id}.zip
    """
}

workflow {
    ch_samples = TOUCH(channel.of('a', 'b', 'c'))
    ch_fastqc = FASTQC(ch_samples)
    ch_fastqc.view()
}
```

## Links

- Supersedes input syntax in [Record types ADR](20260306-record-types.md)
- Related: [Record types syntax summary](../plans/record-types-syntax-new.md)

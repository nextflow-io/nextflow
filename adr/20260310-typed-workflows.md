# Typed workflows

- Authors: Ben Sherman
- Status: accepted
- Date: 2026-03-10
- Tags: lang, static-types, workflows

## Summary

Extend workflows and dataflow logic to provide first-class support for static typing and records.

## Problem Statement

Workflow logic in Nextflow consists of composing processes, channels, and *dataflow operators* (or just *operators*). Operators are essential for transforming, filtering, and combining channels to control the flow of data through a pipeline.

However, workflows were not originally designed with static typing in mind. The introduction of static typing throughout the rest of the language has revealed several gaps in the design of workflow logic.

### Operators

Many operators cannot be statically type-checked because they do not have well-defined argument types and return types. Operators such as `combine` and `join` are designed to work with tuples, but do not support records.

A broader issue is that the operator library is very large, which makes it difficult to find the right operator for a given situation. Several operators deal with additional concerns such as reading/writing specific data formats, which blurs the distinction between dataflow logic and the domain-specific aspects of a workflow. Several operators rely on the ordering of values in a channel, which can cause non-deterministic behavior and hinder reproducibility.

The need for static typing is also an opportunity to address these issues by encouraging the use of a core subset of operators that provide all necessary functionality and support static typing.

### Dataflow syntax

There are many different ways to express the same dataflow logic. Consider the following example:

```groovy
ch_input = channel.of('Hello', 'Hola', 'Ciao')

// alt 1
ch_input
    | GREET
    | map { v -> v.toUpperCase() }
    | view
    | set { ch_upper }

// alt 2
GREET(ch_input)
GREET.out
    .map { v -> v.toUpperCase() }
    .tap { ch_upper }
    .view()

// alt 3
ch_greet = GREET(ch_input)
ch_upper = ch_greet
    .map { v -> v.toUpperCase() }
    .view()
```

Here we see several syntax variants:

- Processes and operators can be composed with pipes (alt 1) or with method calls (alt 2, alt 3).

- Channels can be assigned using `set`, `tap`, or a regular assignment.

- Process outputs can be accessed using the `.out` property on the process name (alt 2) or by assignment (alt 3). The `.out` property can refer to a single output or a record of outputs, depending on the process definition.

Every syntax variant has a cost -- it make code look less familiar to new users, it can cause counterproductive debates over which variant is "better", and it makes Nextflow code less consistent overall. Even if you stick to your preferred syntax, you still have to learn the other variants because you might encounter them when reading someone else's code.

Therefore, syntax sugar should be used judiciously -- it should provide some value that makes adding it worth the aforementioned cost. The variants shown in alt 1 and alt 2 do not add much value relative to their cost.

Even the pipe (`|`), which is loved by many users, can rarely be used in its ideal form because processes usually have additional arguments that can’t be specified in a pipe chain.

## Goals

- Introduce first-class support for static typing and records with dataflow operators

- Encourage the use of a core set of operators (`map`, `filter`, `join`, etc)

- Discourage the use of non-deterministic operators (`buffer`, `distinct`, `first`, etc)

- Discourage the use of syntax variants that do not provide sufficient value to the language

## Non-goals

- Remove support for existing workflow syntax and semantics -- static typing should be opt-in

- Change the way that processes are called -- processes are still called directly with channels, preserving the common mental model of "processes connected by channels"

## Solution

Introduce **typed workflows**, which provide a streamlined syntax for workflows that supports static typing.

Typed workflows can be used with the `nextflow.enable.types` feature flag:

```groovy
// typed workflow
nextflow.enable.types = true

workflow HELLO {
    take:
    ch_names: Channel<String>

    main:
    ch_names.subscribe { name ->
        println "Hello, $name!"
    }
}
```

```groovy
// legacy workflow
workflow HELLO {
    take:
    ch_names

    main:
    ch_names.subscribe { name ->
        println "Hello, $name!"
    }
}
```

This flag behaves the same way for typed processes and typed workflows:

- The flag must be specified in every script that uses typed processes/workflows
- Typed processes/workflows cannot be mixed with legacy processes/workflows in the same script
- Typed and non-typed scripts can be used in the same pipeline

### Operators

The operator library is extended to support static typing and records:

- The `combine` and `join` operators are extended to support both tuples and records.

- The `groupBy` operator is introduced as a statically-typed replacement for `groupTuple`

All operators can be used with or without static typing, with some caveats:

- Some operators have stricter semantics when static typing is enabled via `nextflow.enable.types`. These changes are necessary in order to support static typing effectively. They should not affect the majority of existing code.

- Some operators are discouraged from use with static typing. While they can still be used, the type checker will not be able to validate them. Users should be encouraged to migrate away from them in favor of the *core operators* that are statically typed.

The accompanying reference documentation and best practices guide explain these updates in detail. Here we highlight the most important changes.

The *core operators* are:

- `collect`: collect the channel values into a collection (dataflow value)
- `combine`: emit the combinations of two channels
- `filter`: emit only the channel values that satisfy a condition
- `flatMap`: emit multiple values for each channel value with a closure
- `groupBy`: group channel values by a grouping key
- `join`: relational join of two channels based on a matching key
- `map`: transform each channel value with a closure
- `mix`: concatenate two channels
- `reduce`: reduce channel values into a single value with an accumulator
- `subscribe`: perform an action for each channel value
- `unique`: emit unique values
- `until`: emit each channel value until a stopping condition is satisfied
- `view`: print each channel value

The core operators provide a minimal subset (13 out of ~50) that covers practically all use cases and supports static typing. Encouraging this subset as best practice makes it easier to find the right operator for a given situation, while preserving existing code that uses other operators.

The *legacy operators* are:

| Operator | Problem | Migration strategy |
|---|---|---|
| `branch`                        | Redundant | Use `filter` and `map` for each branch instead |
| `buffer`, `collate`             | Non-deterministic | Use `List::collate()` instead |
| `collectFile`                   | Not statically typed | Use `collect`, `groupBy`, and `Iterable::toSorted()` instead |
| `concat`                        | Redundant | Use `mix` instead |
| `count`, `max`, `min`, `sum`    | Redundant, rarely used | Use `collect` and the corresponding `Iterable` method instead |
| `cross`                         | Redundant | Use `join` instead |
| `distinct`                      | Non-deterministic | Use `unique` instead |
| `dump`                          | Redundant | Use `view` with `tag` option instead |
| `first`, `last`, `take`         | Non-deterministic | Use a list instead |
| `flatten`                       | Not statically typed | Use `flatMap` instead |
| `groupTuple`                    | Not statically typed | Use `groupBy` instead |
| `ifEmpty`                       | Not statically typed | Use `map` with `?:` instead |
| `merge`                         | Non-deterministic | Use `join` instead |
| `multiMap`                      | Redundant | Use `map` instead |
| `randomSample`                  | Non-deterministic | - |
| `set`, `tap`                    | Redundant | Use a regular assignment instead |
| `splitCsv`, `splitFasta`, `splitFastq`, `splitJson`, `splitText` | Not statically typed | Use `flatMap` with the equivalent `Path` method instead |
| `countCsv`, `countFasta`, `countFastq`, `countJson`, `countLines` | Not statically typed | Use `flatMap` with the equivalent `Path` method instead |
| `toDouble`, `toFloat`, `toInteger`, `toLong` | Redundant, rarely used | Use `map` and the corresponding `String` method instead |
| `toList`                        | Redundant | Use `collect` instead |
| `toSortedList`                  | Redundant | Use `collect` and `Iterable::toSorted()` instead |
| `transpose`                     | Not statically typed | Use `flatMap` instead |

In most cases, a legacy operator can be rewritten in terms of core operators and standard library functions. The accompanying best practices guide provides detailed examples for each operator. Since legacy operators can still be used in typed workflows, users can migrate away from legacy operators at their own pace.

### Fewer syntax variants

Typed workflows do not support the following syntax variants:

- Implicit `it` closure parameter → declare an explicit parameter instead
  - `it` can still be used as a variable name as long as it is explicitly declared

- Using `Channel` to access channel factories → use `channel` instead
  - `Channel` should be used only in type annotations

- Using `set` or `tap` to assign channels → use assignments instead

- Special dataflow operators `|` and `&` → use assignments and method calls instead
  - The equivalent bitwise operators are still allowed

- Using the `.out` property to access process and workflow outputs → use assignments instead  

These restrictions are designed to make Nextflow code more consistent across the board and more familiar to users from other programming languages. Things like variable assignments and method calls in Nextflow look and feel the same as most other languages, whereas things like `set` assignments and the `.out` property make Nextflow code feel more unfamiliar without adding much value.

This aspect of the language is becoming more salient as code is increasingly read and written by AI agents. Agents need many examples of a programming language in order to use it effectively, so when a niche language has many syntax variants or syntax that deviates heavily from the common patterns used by other languages, it hurts the agent's ability to read and write code in that language.

## Distinguishing between typed and legacy workflows

Static typing has been introduced as multiple independent features:

- Typed parameters (`params` block)
- Typed outputs (`output` block)
- Typed processes
- Record types
- Typed workflows (this proposal)

This incremental approach was done in contrast to DSL2, which was a monolithic change that required an entire pipeline to be updated at once. With static typing, each new feature can be adopted independently of the others, rather than requiring all new features to be adopted at once (e.g. "DSL3").

Most of the features for static typing are new concepts that can be used alongside existing code. However, typed processes and typed workflows modify existing concepts (`process` and `workflow` definitions), so they require a feature flag.

The `nextflow.enable.types` feature flag will be used to distinguish between typed and legacy code, indefinitely. It would only be removed if the support for legacy syntax was removed, which is unlikely since DSL2 has been the standard Nextflow syntax for many years.

To help distinguish between typed and legacy workflows, the use of type annotations should be allowed only for typed workflows:

```groovy
// legacy workflow
workflow greet {
    take:
    greetings

    main:
    messages = greetings.map { v -> "$v world!" }

    emit:
    messages
}
```

```groovy
// typed workflow
nextflow.enable.types = true

workflow greet {
    take:
    greetings: Channel<String>

    main:
    messages = greetings.map { v -> "$v world!" }

    emit:
    messages: Channel<String>
}
```

## Interoperability between typed and legacy workflows

Typed and legacy workflows use different underlying dataflow types:

- **Legacy workflows (v1)** use raw GPars types: `DataflowBroadcast` (queue channel) and `DataflowVariable` (value channel).

- **Typed workflows (v2)** use wrapper types: `ChannelImpl` (wraps a `DataflowBroadcast`) and `ValueImpl` (wraps a `DataflowVariable`). These wrappers implement the new operators and integrate with the type system.

While a given script must be entirely typed or entirely legacy (controlled by the `nextflow.enable.types` flag), **typed and legacy workflows can call each other across different scripts**. This interoperability enables incremental migration -- individual scripts can be migrated to static typing without having to update the entire pipeline at once.

### Normalization at call sites

When a workflow calls another workflow, the Nextflow runtime automatically converts dataflow arguments and return values to the appropriate type for each side of the call site.

Normalization can occur in either direction:

- **v2 → v1 (unwrap)**: when passing typed channels to a legacy component, `ChannelImpl` / `ValueImpl` are unwrapped to the underlying `DataflowBroadcast` / `DataflowVariable`.

- **v1 → v2 (wrap)**: when passing legacy channels to a typed component, `DataflowBroadcast` / `DataflowVariable` are wrapped as `ChannelImpl` / `ValueImpl`.

The normalization is applied twice per call: once to the arguments (converted to match the *callee's* semantics), and once to the return value (converted to match the *caller's* semantics).

### Example: typed workflow calling a legacy workflow

**`legacy.nf`**
```groovy
workflow LEGACY_ALIGN {
    take:
    reads // DataflowBroadcast

    main:
    ALIGN(reads)

    emit:
    bam = ALIGN.out // DataflowBroadcast
}
```

**`typed.nf`**
```groovy
nextflow.enable.types = true

include { LEGACY_ALIGN } from './legacy'

workflow {
    reads = channel.fromPath('*.fastq') // ChannelImpl

    // `reads` is unwrapped to DataflowBroadcast when passed to LEGACY_ALIGN
    // The return value (DataflowBroadcast) is wrapped to ChannelImpl
    bam = LEGACY_ALIGN(reads)

    bam.view() // ChannelImpl
}
```

### Example: legacy workflow calling a typed workflow

**`typed.nf`**
```groovy
nextflow.enable.types = true

workflow TYPED_TRIM {
    take:
    reads: Channel<Record>

    main:
    ch_trimmed = TRIM(reads)

    emit:
    trimmed = ch_trimmed
}
```

**`legacy.nf`**
```groovy
include { TYPED_TRIM } from './typed'

workflow {
    reads = Channel.fromPath('*.fastq') // DataflowBroadcast

    // `reads` is wrapped as ChannelImpl when passed to TYPED_TRIM
    // The return value (ChannelImpl) is unwrapped to DataflowBroadcast
    trimmed = TYPED_TRIM(reads)

    trimmed.view() // DataflowBroadcast
}
```

### Process and workflow outputs (`ChannelOut`)

Processes and workflows -- regardless of whether they are legacy or typed -- always return a `ChannelOut`, a specialized class that can contain one or more named outputs (`DataflowBroadcast` / `DataflowVariable`).

When a `ChannelOut` is returned to a typed workflow, it is normalized as follows:

- If the `ChannelOut` contains only one output, it is unwrapped to the underlying `DataflowBroadcast` / `DataflowVariable` and then wrapped as a `ChannelImpl` / `ValueImpl`.

- If the `ChannelOut` contains multiple outputs, it is converted to a record (`RecordMap`), where each named output is wrapped as a `ChannelImpl` / `ValueImpl`.

For example:

**`legacy.nf`**
```groovy
workflow LEGACY_QC {
    take:
    reads

    main:
    FASTQC(reads)
    MULTIQC(FASTQC.out)

    emit:
    fastqc = FASTQC.out     // DataflowBroadcast
    multiqc = MULTIQC.out   // DataflowBroadcast
}
```

**`typed.nf`**
```groovy
nextflow.enable.types = true

include { LEGACY_QC } from './legacy'

workflow {
    reads = channel.fromPath('*.fastq')

    // LEGACY_QC returns a ChannelOut with two outputs
    // which is converted to a record:
    //   Record { fastqc: ChannelImpl ; multiqc: ChannelImpl }
    qc = LEGACY_QC(reads)

    // RecordMap provides same semantics as ChannelOut
    qc.fastqc.view()
    qc.multiqc.view()
}
```

## Alternatives

### Processes in operator closures

A process call is essentially a task function wrapped in a `map` operation. But processes are called directly on channels, which has a few implications:

- The true structure of process calls are somewhat obscured
- Process calls can have a different return type (`Channel` or `Value`) depending on how they are called
- Processes can only be called as a `map` operation, not with other operators like `reduce`
- Processes can not be chained like operator calls (without additional syntax like `|`)

These limitations could be addressed by calling processes in operator closures instead of calling them directly with channels:

```groovy
ch_samples = channel.of(...)
fasta = file(...)

// before
SALMON(ch_samples, fasta)

// after
ch_samples.map { sample -> SALMON(sample, fasta) }
```

Where `SALMON` is defined as follows:

```groovy
process SALMON {
    input:
    record(
        id: String,
        fastq: Path
    )
    fasta: Path

    // ...
}
```

This syntax brings a number of benefits:

- The `map` operation is explicit, making process calls consistent with other operator logic
- The process call matches the process definition -- it accepts and returns regular values, not channels or dataflow values
- Processes can be chained without needing a pipe syntax (e.g. `ch.map(FOO).map(BAR).map(BAZ) ...`)
- Processes could theoretically be composed with other operators (e.g. an iterative process with the `reduce` operator)
- The closure around the process call can be used to handle process inputs and outputs without additional operator calls

Decoupling the process lifecycle from an implicit `map` operation, however, breaks a key assumption of the Nextflow runtime:

- Handling a process call in an arbitrary closure instead of as an implicit `map` operation is significantly more complex, and would likely require new language semantics and compiler transformations to implement.
- Alternatively, the compiler could restrict such closures to specific patterns (e.g. a single process call with some statements before and after), but this would add complexity for developer experience (i.e. having to remember which patterns are allowed in which cases).

Additionally, some of the problems that motivated this approach have been addressed by type checking and records:

- The type checker can infer the return type of a direct process call from the call arguments (e.g. `Channel` vs `Value`)
- Records and record types provide additional flexibility that eliminates much of the adaptor logic that was required between tuple channels and processes

Ultimately, this change would mostly be a cosmetic syntax improvement that would do little to improve the developer experience, but would introduce a great deal of complexity to the compiler and runtime. It would also be a significant break from the way that Nextflow workflows have been written since the introduction of DSL2.

### Processes as operator closures

A moderated version of calling processes in operator closures is to call them *as* operator closures:

```groovy
ch_samples = channel.of(...)
fasta = file(...)

// before
SALMON(ch_samples, fasta)

// after
ch_samples.map(SALMON, index: fasta)
```

The process name takes the place of the `map` closure. The channel calling `map` is supplied as the first process input, and any additional inputs are supplied as named arguments to `map`.

This approach avoids much of the aforementioned complexity risk while retaining many of the benefits.

For example, process calls can be chained with other operator calls:

```groovy
ch_input
    .map(GREET, greeting: "Hello")
    .map { v -> v.toUpperCase() }
    .view()
```

And processes can be called with other operators such as `reduce`:

```groovy
process ACCUMULATE {
    input:
    result: Path
    input: Path

    script:
    """
    cat ${input} >> ${result}
    """

    output:
    file('result.txt')
}

workflow {
    channel.fromPath("*.txt").reduce(ACCUMULATE).view()
}
```

This particular pattern was proposed as a cleaner alternative to the experimental [recursion](https://nextflow.io/docs/latest/workflow.html#process-and-workflow-recursion) feature. As long as the process matches the signature of the accumulator closure (two inputs and one output), the process can be executed iteratively.

While this approach avoids most of the potential complexity that would be required to call processes in operator closures, it is still a significant syntax change with dubious relative benefit.

Investigating these approaches revealed an important trade-off -- Nextflow sacrifices a small amount of syntactic precision in order to make process calls prominent in the workflow logic. While calling processes in an operator would be more correct and provide some additional flexibility (e.g. using processes with other operators), it would make workflows feel much more like "operators that call processes in closures" instead of "processes connected by channels".

The reality is that most Nextflow users think of their pipelines as "processes connected by channels", and operator logic is a minor detail at best and a confusing distraction at worst. While we can and should make channel operators as simple and pleasant to use as possible, it should be in service of making them less prominent in the language, not more.

### Implicit dataflow values

Dataflow values (a.k.a. *value channels*) are analogous to Futures or Promises in other languages. For example, given a `CompletableFuture` in Java, you can either call `get()` to await the value or `thenAccept()` / `thenApply()` to invoke a callback when the value is ready.

Dataflow values can similarly call `subscribe` or `map`, but it is not possible to "await" a dataflow value directly. For example, it is not possible to use a dataflow value in an `if` statement:

```groovy
vals = channel.of(1..10).collect()
if( vals.size() > 2 )
    println 'More than two!'
```

Instead, you must use `subscribe` to act on the value asynchronously:

```groovy
vals = channel.of(1..10).collect()
vals.subscribe { _vals ->
    if( _vals.size() > 2 )
        println 'More than two!'
}
```

This is a common frustration for many users, that dataflow values don't quite work like regular values, even though it seems like they should.

A solution could be to make dataflow values *implicit* -- users would use them like regular values (i.e. the first example above) and the compiler would translate the user's code into explicit dataflow logic (i.e. the second example).

To do this, the compiler would need to:

1. distinguish implicit dataflow values from regular values via type inference (e.g. the result of a `collect` operator),

2. wrap downstream code in `map` and/or `subscribe` operators as needed to produce the desired dataflow logic.

In the end, however, this change does not seem worthwhile:

- It makes type inference an essential part of the compilation process rather than an optional enhancement.

- The above example is simple to understand, but it is easy to construct more complicated examples that quickly cast doubt on whether the compiler could solve this problem in general.

- Even if there is a general solution, any mistake or edge case would likely lead to unexpected behavior that would be extremely difficult to debug (e.g. a low-level compiler error, compiled code that is silently incorrect).

Additionally, most of the problems that motivated this idea have been effectively solved by type checking:

- There is now an explicit `Value` type which allows both developers and the type checker to distinguish between channels and dataflow values.

- While users still can't use a dataflow value in an `if` statement, they can get clear and early feedback on whether their code is valid, which is what ultimately matters.

- Being transparent about regular values vs dataflow values in the language may be for the best anyway -- it provides a clear picture of how things are working "under the hood", and it is still far simpler than the async programming models employed by most languages.

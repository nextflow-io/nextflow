# Typed workflows

- Authors: Ben Sherman
- Status: accepted
- Date: 2026-03-10
- Tags: lang, static-types

## Summary

Provide a way to write workflow logic in a way that supports static typing and records as first-class features.

## Problem Statement

Workflow logic in Nextflow consists of composing processes, channels, and *dataflow operators* (or just *operators*). Operators are essential for transforming, filtering, and combining channels to control the flow of data through a pipeline.

### Operators

Over the years, the dataflow operators have accrued a number of problems:

- The operator library is bloated. There are over [50 operators](https://nextflow.io/docs/latest/reference/operator.html), and there are several cases of similar operators such as `collect` / `toList` / `toSortedList` and `combine` / `cross` / `join`, which makes it difficult to determine which operator to use for a given situation.

- There are several operators that mix dataflow with other concerns, such as `collectFile` and `splitCsv`. These operators tend to be complicated, with many options and behaviors, making them hard to understand. They also blur the scope and purpose of dataflow operators, making it harder to separate dataflow logic conceptually from other aspects of a pipeline.

- There are several operators that rely on the ordering of values in a channel, even though channels are unordered by design, such as `buffer`, `distinct`, `first`, and `last`. These operators can cause non-deterministic behavior when used improperly. It is usually more appropriate to use an equivalent function, such as `List::first()` instead of the `first` operator.

The introduction of static typing and records has brought more problems:

- Some operators such as `flatten` cannot be statically type-checked, because they do not have a well-defined return type.

- Some operators such as `groupTuple` and `join` are designed to work with tuples, but do not have first-class support for records.

Simply updating the operators, even to only add support for records, would require either a breaking change or a lot of added complexity to maintain backwards compatibility. Instead, the need for static typing and records is a good opportunity to introduce a new operator library that addresses all of the above issues.

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

- Simplify the operators to a core set of stream operations (`map`, `filter`, `join`, etc)

- Discourage the use of non-deterministic operators (`buffer`, `distinct`, `first`, etc)

- Discourage the use of syntax variants that do not provide sufficient value to the language

## Non-goals

- Remove support for existing workflow semantics -- typed workflows should be opt-in

- Change the way that processes are called -- processes are still called directly with channels, preserving the common mental model of "processes connected by channels"

## Solution

Introduce new workflow semantics with the `nextflow.preview.types` feature flag:

```groovy
// hello.nf
nextflow.preview.types = true

workflow HELLO {
    // typed workflow semantics
}
```

```groovy
// goodbye.nf
workflow GOODBYE {
    // legacy workflow semantics
}
```

This flag applies the same rules for both typed processes and typed workflows:

- The flag must be specified in every script that uses typed processes/workflows
- Typed processes/workflows cannot be mixed with legacy processes/workflows in the same script
- Typed and non-typed scripts can be used in the same pipeline

### New operators

Typed workflows use a new set of operators for the `Channel` type. When `nextflow.preview.types` is enabled, channels will use these operators instead of the legacy operators.

The accompanying reference documentation and migration guide explain these operators in detail and how to update legacy code to use them. Here we will highlight the most important changes.

The new operators are:

- `collect`: collect the channel values into a collection (dataflow value)

- `cross`: cross-product of two channels

- `filter`: emit only the channel values that satisfy a condition

- `flatMap`: emit multiple values for each channel value with a closure (scatter)

- `groupBy`: group channel values by a grouping key (gather)

- `join`: relational join of two channels based on a matching key

- `map`: transform each channel value with a closure

- `mix`: concatenate two channels

- `reduce`: reduce channel values into a single value with an accumulator

- `subscribe`: perform an action for each channel value

- `unique`: emit unique values

- `until`: emit each channel value until a stopping condition is satisfied

- `view`: print each channel value

These operators are essentially the core subset of operators that are needed to implement dataflow logic. They are much fewer (13 instead of ~50) and it is much clearer which operator should be used for a given situation.

Some of the new operators have slightly different semantics:

- `collect` does not flatten lists (use `flatMap` instead)

- `cross` does not join on a matching key (use `join` instead)

- `flatMap` does not flatten maps (use `Map::entrySet()` with `flatMap` instead)

- `groupBy` accepts a 3-tuple of `(key, size, value)` instead of wrapping the grouping key with `groupKey()`

- `join` operates on records instead of tuples, matching on record fields instead of tuple indices

- `join` emits cross products for duplicate matches (consistent with relational join semantics)

These changes are intended to make these operators easy to use while enabling first-class support for static typing and records.

Legacy operators that are not supported in typed workflows can be migrated as follows:

| Operator | Migration strategy |
|---|---|
| `branch`                        | Use `filter` and `map` for each branch instead |
| `buffer`, `collate`             | Use `List::collate()` instead |
| `collectFile`                   | Use `collect`, `groupBy`, and `Iterable::toSorted()` instead |
| `combine`                       | Use `cross` or `join` instead |
| `concat`                        | Use `mix` instead |
| `count`, `max`, `min`, `sum`    | Use `collect` and the corresponding `Iterable` method instead |
| `distinct` [^nondeterministic]  | Use `unique` instead |
| `dump`                          | Use `view` instead (it will have the `tag` option) |
| `first`, `last`, `randomSample`, `take` [^nondeterministic] | Use a list instead |
| `flatten`                       | Use `flatMap` instead |
| `ifEmpty`                       | Use `map` with `?:` instead |
| `merge` [^nondeterministic]     | Use `join` instead |
| `multiMap`                      | Use `map` instead |
| `set`, `tap`                    | Use a regular assignment instead |
| `splitCsv`, `splitFasta`, `splitFastq`, `splitJson`, `splitText` | Use `flatMap` with the equivalent `Path` method instead |
| `countCsv`, `countFasta`, `countFastq`, `countJson`, `countLines` | Use `flatMap` with the equivalent `Path` method instead |
| `toDouble`, `toFloat`, `toInteger`, `toLong` | Use `map` and the corresponding `String` method instead |
| `toList`                        | Use `collect` instead |
| `toSortedList`                  | Use `collect` and `Iterable::toSorted()` instead |
| `transpose`                     | Use `flatMap` instead |

[^nondeterministic]: This operator is non-deterministic -- you probably shouldn't be using it anyway.

In most cases, a legacy operator can be rewritten in terms of existing operators and standard library functions. The accompanying migration guide provides detailed examples for each operator.

### Fewer syntax variants

Typed workflows do not support the following syntax variants:

- Implicit `it` closure parameter → declare an explicit parameter instead
  - `it` can still be used as a variable name as long as it is explicitly declared

- Using `Channel` to access channel factories → use `channel` instead
  - `Channel` should be used only in type annotations

- Special dataflow operators `|` and `&` → use assignments and method calls instead
  - The equivalent bitwise operators are still allowed

- The `.out` property for processes and workflows → use assignments instead  

THese restrictions are designed to make Nextflow code more consistent across the board and more familiar to users from other programming languages. Things like variable assignments and method calls in Nextflow look and feel the same as most other languages, whereas things like `set` assignments and the `.out` property make Nextflow code feel more unfamiliar without adding much value.

This aspect of the language is becoming more salient as code is increasingly read and written by AI agents. Agents need many examples of a programming language in order to use it effectively, so when a niche language has many syntax variants or syntax that deviates heavily from the common patterns used by other languages, it hurts the agent's ability to read and write code in that language.

## How to distinguish between typed and legacy workflows?

Static typing has been introduced as multiple independent features:

- Type annotations
- Typed parameters (`params` block)
- Typed outputs (`output` block)
- Typed processes
- Record types
- Typed workflows (this proposal)

This approach was done in contrast to DSL2, which was a monolithic change that required an entire pipeline to be updated at once. With static typing, each new feature can be adopted independently of the others, rather than requiring all new features to be adopted at once (e.g. "DSL3").

However, the challenge with this approach is to make sure that it is easy for users (and agents) to distinguish between new and old syntax.

Several alternative approaches are considered below:

### Option 1: Use `nextflow.enable.types` to enable typed processes and typed workflows

Most of the features for static typing are *purely additive* -- they are new concepts that can be used alongside existing code. However, typed processes and typed workflows modify existing concepts (`process` and `workflow` definitions), so they require the `nextflow.preview.types` feature flag.

The `nextflow.enable.types` feature flag will replace the preview flag once the feature set is stable, and it will likely be used indefinitely to distinguish between typed and legacy code. It would only be removed if the support for legacy syntax was removed, which is unlikely since DSL2 has been the standard Nextflow syntax for many years.

However, while typed processes look significantly different from legacy processes, typed workflows do not. Typed workflows look very similar but have slightly different semantics. A feature flag may not be enough to signal the difference to users and agents, even if it is sufficient for the compiler and language server.

### Option 2: Use `nextflow.enable.types` to enable all static typing

Now that the entire language has been updated to support static typing, it could make sense to provide it as a single feature controlled by a single feature flag:

```groovy
// "dynamically typed" code
// nextflow.enable.types = true

// "statically typed" code
nextflow.enable.types = true
```

Even though these features can be adopted independently in principle, they are designed to work together, and in practice it is difficult to adopt one feature without the others:

- Migrating a large pipeline to workflow outputs is very difficult without also migrating to typed processes and record types.
- Adopting type annotations (e.g. for workflow takes and emits) can provide some basic documentation and validation, but most workflow logic still cannot be effectively validated by the type checker without typed processes.

Enabling all static typing features via `nextflow.enable.types` would establish a clear boundary between *statically typed* code and *dynamically typed* code. This way, the poor distinction between typed workflows and legacy workflows is made up by the clear distinction of type annotations, record types, etc in the same context.

Since type annotations and typed parameters were introduced in Nextflow 25.10 as stable, requiring a feature flag for them in 26.04 would be a breaking change. However, this break might be acceptable for now since these features are still new and many users are waiting for full static typing anyway. These features could be allowed with a warning in 26.04 to ease the transition.

See also: static compilation in Groovy via `@CompileStatic`

### Option 3: Enable new operators via `include` declaration

Since operators are methods of the `Channel` type, new operators can be understood as a new `Channel` type / `channel` namespace. Therefore, the new operators could be introduced simply by including a different version of `Channel` or `channel`:

```groovy
// legacy operators (default)
// include { channel } from 'dataflow/v1'

// typed operators
include { channel } from 'dataflow/v2'
```

This approach is similar to using a feature flag, but it more clearly expresses the intent of using the new operators. The other aspects of typed workflows -- removal of certain syntax patterns, process named arguments -- could also be enabled by this include or by the `nextflow.enable.types` feature flag.

Either way, the feature flag will still be needed to enable typed processes, so users will end up using both the feature flag and include across their scripts. This might be more complicated than just using a feature flag.

### Option 4: Use new operator names in typed workflows

Typed workflows could simply rename all operators that were changed. This would clearly distinguish typed workflows from legacy workflows.

The problem of semantic changes essentially comes down to `cross` and `join`:

- `groupBy` was renamed from `groupTuple`
- all other operators are effectively identical, with minor differences that amount to bug fixes

Even `join` will be distinct in typed workflows because it will join on record fields instead of tuple indices:

```groovy
// legacy workflow
left.join(right, by: 0)

// typed workflow
left.join(right, by: 'id')
```

Some possible names:

- `cross` -> `crossV2`, `crossProduct`, `crossJoin`, `combine`
- `join` -> `joinV2`, `joinBy`, `joinInner` (remainder: false), `joinOuter` (remainder: true), `combineBy`

Ironically, the new operators are more true to their names than the old ones:

- `cross` now performs a true cross product (the legacy `cross` implicitly joined on matching keys)
- `join` now performs a true relational join (the legacy `join` did not handle duplicates correctly)

Ultimately, new operator names alone might not be enough to signal the other aspects of typed workflows, such as the removal of many other operators and syntax patterns.

### Option 5: Replace `process` and `workflow` with `task` and `flow`

The key issue is that typed processes and typed workflows modify existing concepts (processes and workflows) rather than adding to them. Instead, we could make these features purely additive by introducing them as new top-level concepts:

```groovy
// legacy semantics
process FASTQC { /* ... */ }
workflow RNASEQ { /* ... */ }

// typed semantics
task FASTQC { /* ... */ }
flow RNASEQ { /* ... */ }
```

This approach would make the distinction very clear. However, this change would be a significant break, since *processes* and *workflows* are long-standing and fundamental ideas in Nextflow. Users think about Nextflow pipelines in terms of *processes* and *workflows*, so introducing new terminology would be confusing.

See also: [Prefect](https://www.prefect.io/prefect/open-source) (uses `@task` and `@flow` in their DSL)

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
    sample: Record {
        id: String
        fastq: Path
    }
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

For example, process calls can be chained with other opartor calls:

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

While this approach avoids most of the potential complexity that would be required to call processes in operator closures, it is still a significant syntax change with dubious relatively benefit.

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

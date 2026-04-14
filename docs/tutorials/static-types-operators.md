(migrating-static-types-operators)=

# Using operators with static typing

Nextflow 26.04 brings updates to the operator library in order to support static tying and records. This page provides best practices for using operators with static typing.

See {ref}`migrating-static-types` for more information about migrating pipelines to static typing.

## Overview

All operators can be used with or without static typing (i.e. {ref}`typed workflows <workflow-typed-page>`). However, only a core subset of operators are recommended for use with static typing, while the rest are discouraged. They are distinguished here as *core operators* and *legacy operators*.

## Core operators

The {ref}`core operators <operator-typed-page>` are recommended for use with static typing. When static typing is enabled (via `nextflow.preview.types`), some of these operators have stricter semantics which may require minor changes to pipeline code. These cases are described below.

### collect

When using `collect` with static typing, it has the same semantics as `toList`. Collected values are not flattened, and when the source channel is empty, an empty list is emitted.

### combine

When using `combine` with static typing, the right operand should be a channel, dataflow value, or named arguments corresponding to record fields.

For uses of `combine` with the `by` option, use `join` instead:

```nextflow
// before
left = channel.of( [1, 'alpha'], [2, 'beta'] )
right = channel.of( [1, 'x'], [1, 'y'], [2, 'p'] )
left.combine(right, by: 0).view()

// [1, alpha, x]
// [1, alpha, y]
// [2, beta, p]

// after (static typing enabled)
left = channel.of(
    record(id: 1, name: 'alpha'),
    record(id: 2, name: 'beta')
)
right = channel.of(
    record(id: 1, code: 'x'),
    record(id: 1, code: 'y'),
    record(id: 2, code: 'p')
)
left.join(right, by: 'id').view()

// [id:1, name:alpha, code:x]
// [id:1, name:alpha, code:y]
// [id:2, name:beta, code:p]
```

### filter

When using `filter` with static typing, the predicate should be a closure.

```nextflow
// before
ch.filter( ~/^a.*/ )

// after (static typing enabled)
ch.filter { v -> v == ~/^a.*/ }
```

### flatMap

When using `flatMap` with static typing, the mapping closure should always return a collection. Maps and tuples are not automatically flattened because they are not collection types.

```nextflow
// before
channel.of( 1, 2, 3 )
    .flatMap { n -> [ number: n, square: n*n, cube: n*n*n ] }
    .view { entry -> "${entry.key}: ${entry.value}" }

// after (static typing enabled)
channel.of( 1, 2, 3 )
    .flatMap { n -> [ tuple('number', n), tuple('square', n*n), tuple('cube', n*n*n) ] }
    .view { key, value -> "${key}: ${value}" }
```

### groupBy

The `groupBy` operator is a replacement for `groupTuple` that is statically typed.

While `groupTuple` accepts tuples of arbitrary length, `groupBy` accepts either a 2-tuple of `(<key>, <value>)` or a 3-tuple of `(<key>, <size>, <value>)`. Specifying the group size with each input tuple provides the same behavior as using the `size` option (or wrapping each key with `groupKey()`) does with `groupTuple`.

While `groupTuple` can group multiple lists in a group, `groupBy` always emits 2-tuples of the form `(<key>, <values>)`, where `<values>` is an unordered collection (`Bag`). This approach avoids a pitfall with `groupTuple` where the grouped lists can be ordered inconsistently.

### join

When using `join` with static typing, the `by` option is required. It should be either an integer (for joining tuples by index) or a string (for joining records by field name).

When using `join` with records, the `failOnDuplicate` and `failOnMismatch` options are not supported. Duplicate matches are handled by emitting each matching combination (like a relational join). Unmatched records are either emitted or discarded depending on whether the `remainder` option is set. To fail on mismatches, use the `remainder` option and check for unmatched records in downstream logic.

```nextflow
// tuples
left  = channel.of( ['X', 1], ['Y', 2], ['Z', 3], ['P', 7] )
right = channel.of( ['Z', 6], ['Y', 5], ['X', 4] )
left.join(right).view()

// [X, 1, 4]
// [Y, 2, 5]
// [Z, 3, 6]

// records
left  = channel.of(
    record(id: 'X', a: 1),
    record(id: 'Y', a: 2),
    record(id: 'Z', a: 3),
    record(id: 'P', a: 7)
)
right = channel.of(
    record(id: 'Z', b: 6),
    record(id: 'Y', b: 5),
    record(id: 'X', b: 4)
)
left.join(right, by: 'id').view()

// [id: X, a: 1, b: 4]
// [id: Y, a: 2, b: 5]
// [id: Z, a: 3, b: 6]
```

### map

When using `map` with static typing, `null` values are not automatically discarded. Use `filter` to discard `null` values explicitly.

```nextflow
// before
ch.map { r -> r.id }

// after (static typing enabled)
ch.map { r -> r.id }.filter { id -> id != null }
```

### mix

When using `mix` with static typing, only one argument should be supplied for each `mix` call.

```nextflow
// before
ch1.mix(ch2, ch3)

// after (static typing enabled)
ch1.mix(ch2).mix(ch3)
```

## Legacy operators

The {ref}`legacy operators <operator-page>` are discouraged from use with static typing. They can still be used, but the type checker will not be able to validate them.

This section describes how to rewrite each legacy operator with core operators.

### branch

Use `filter` and `map` for each branch instead. Using records instead of tuples can eliminate much of the need for `branch`.

Example requiring only `filter`:

```nextflow
// before
ch_gvcf_branch = ch_gvcf.branch { meta, gvcf, tbi ->
    no_tbi: !tbi
        return tuple(meta, gvcf)
    tbi: tbi
        return tuple(meta, gvcf, tbi)
}

ch_gvcf_branch.no_tbi.view()
ch_gvcf_branch.tbi.view()

// after
ch_gvcf_no_tbi = ch_gvcf.filter { s -> !s.tbi }
ch_gvcf_tbi = ch_gvcf.filter { s -> s.tbi }
```

Example requiring `filter` and `map`:

```nextflow
// before
ch_input_by_type = ch_input.branch { meta, platform, fastq_1, fastq_2 ->
    fastq: meta.single_end || fastq_2
        return tuple(meta + [type: "short"], fastq_2 ? [fastq_1, fastq_2] : [fastq_1])

    nanopore: platform == 'OXFORD_NANOPORE'
        meta.single_end = true
        return tuple(meta + [type: "long"], [fastq_1])

    pacbio: platform == 'PACBIO_SMRT'
        meta.single_end = true
        return tuple(meta + [type: "long"], [fastq_1])
}

ch_input_by_type.fastq.view()
ch_input_by_type.nanopore.view()
ch_input_by_type.pacbio.view()

// after -- no more fastq_1/fastq_2 wrangling
ch_input_fastq = ch_input
    .filter { s -> s.single_end || s.fastq_2 }
    .map { s -> s + record(type: 'short') }

ch_input_nanopore = ch_input
    .filter { s -> s.platform == 'OXFORD_NANOPORE'}
    .map { s -> s + record(single_end: true, type: 'long') }

ch_input_pacbio = ch_input
    .filter { s -> s.platform == 'PACBIO_SMRT' }
    .map { s -> s + record(single_end: true, type: 'long') }

ch_input_fastq.view()
ch_input_nanopore.view()
ch_input_pacbio.view()
```

### buffer, collate

These operators are {ref}`non-deterministic <cache-nondeterministic-inputs>`. Use `List::collate()` instead.

### collectFile

The `collectFile` operator is useful for collecting intermediate results into a final output file, or writing a samplesheet. In many cases, `collectFile` can be replaced by a {ref}`workflow output <workflow-output-def>`, which can generate an index file for a published channel.

For other cases, consider the following alternatives:

- Use an `exec` process to write text files (see {ref}`working-with-files`)
- Use the `groupBy` operator to group  
- Use `Iterable::toSorted` to sort

You can compose these functions and operators as needed to achieve the desired functionality.

### concat

Use `mix` instead.

### count, max, min, sum

Use `collect` and the corresponding `Iterable <stdlib-types-iterable>` methods instead.

### cross

Use `join` with records instead.

```nextflow
// before
left = channel.of( [1, 'alpha'], [2, 'beta'] )
right = channel.of( [1, 'x'], [1, 'y'], [2, 'p'] )
left.cross(right).view()

// [[1, alpha], [1, x]]
// [[1, alpha], [1, y]]
// [[2, beta], [2, p]]

// after
left = channel.of(
    record(id: 1, name: 'alpha'),
    record(id: 2, name: 'beta')
)
right = channel.of(
    record(id: 1, code: 'x'),
    record(id: 1, code: 'y'),
    record(id: 2, code: 'p')
)
left.join(right, by: 'id').view()

// [id:1, name:alpha, code:x]
// [id:1, name:alpha, code:y]
// [id:2, name:beta, code:p]
```

### distinct

This operator is {ref}`non-deterministic <cache-nondeterministic-inputs>`. Use `unique` instead.

### dump

Use `view` instead. The `view` operator now supports the `tag` option, allowing it to be used like `dump`.

### first, last, take

These operators are {ref}`non-deterministic <cache-nondeterministic-inputs>`. Use the corresponding {ref}`List <stdlib-types-list>` methods instead.

### flatten

Use `flatMap` instead.

### ifEmpty

The `ifEmpty` operator is typically used to either (1) raise an error if a channel is empty or (2) provide a fallback for a null dataflow value.

With static typing, both cases can be implemented without `ifEmpty`:

```nextflow
// (1) fail if channel is empty
files_ch = channel.fromPath('*.txt')
files_ch.collect().subscribe { files ->
    if( files.isEmpty() )
        error 'no input files were found'
}
files_ch.view()

// (2) provide a fallback for dataflow value
index_file = FETCH_INDEX().map { index ->
    index ?: file('index_default.txt')
}
index_file.view()
```

The example for (2) assumes that `FETCH_INDEX` is a typed process. Typed processes emit `null` when an optional output is missing, whereas legacy processes emit nothing.

### merge

This operator is {ref}`non-deterministic <cache-nondeterministic-inputs>`. Use `join` instead.

### multiMap

Use `map` for each branch instead. Using records instead of tuples can eliminate much of the need for `branch`.

For example:

```nextflow
// before
ch_input_by_type = ch_input.multiMap { families, meta, cram, crai, gvcf, tbi, roi ->
    def new_meta = meta + [
        family_count: families[meta.family].size(),
        type: gvcf && cram ? "gvcf_cram" : gvcf ? "gvcf" : "cram"
    ]

    gvcf: tuple(new_meta, gvcf, tbi)
    cram: tuple(new_meta, cram, crai)
    roi:  tuple(new_meta, roi)
}

ch_input_by_type.gvcf.view()
ch_input_by_type.cram.view()
ch_input_by_type.roi.view()

// after -- just keep everything in a single record
ch_input = ch_input.map { s ->
    s + record(
        family_count: s.families[s.family].size(),
        type: s.gvcf && s.cram ? "gvcf_cram" : s.gvcf ? "gvcf" : "cram"
    )
}
```

### randomSample

This operator is {ref}`non-deterministic <cache-nondeterministic-inputs>`. It should not be used.

If needed, it is possible to implement a function that samples a collection (e.g., using `Math.random()` from the Java standard library).

### set

Use standard assignments instead:

```nextflow
// before
channel.of(10, 20, 30).set { my_channel }

// after
my_channel = channel.of(10, 20, 30)
```

### tap

Use standard assignments instead:

```nextflow
// before
channel.of(10, 20, 30)
    .tap { log1 }
    .map { v -> v * 2 }
    .tap { log2 }

// after
log1 = channel.of(10, 20, 30)
log2 = log1.map { v -> v * 2 }
```

### splitCsv, splitFasta, splitFastq, splitJson, splitText

Use the equivalent {ref}`stdlib-types-path` method with `flatMap` instead:

```nextflow
// before
channel.fromPath('samplesheet.csv')
    .splitCsv(sep: ',')
    .view()

// after
channel.fromPath('samplesheet.csv')
    .flatMap { csv -> csv.splitCsv(sep: ',') }
    .view()
```

### toList

Use `collect` instead.

### toSortedList

Use `collect` and `Iterable::toSorted` instead:

```nextflow
// before
channel.of(3, 2, 1, 4)
    .toSortedList()
    .view()

// after
channel.of(3, 2, 1, 4)
    .collect()
    .map { vals -> vals.toSorted() }
    .view()
```

### transpose

Use `flatMap` instead:

```nextflow
// before
channel.of(
        tuple(1, ['A', 'B', 'C']),
        tuple(2, ['C', 'A']),
        tuple(3, ['B', 'D']),
    )
    .transpose()
    .view()

// after
channel.of(
        tuple(1, ['A', 'B', 'C']),
        tuple(2, ['C', 'A']),
        tuple(3, ['B', 'D']),
    )
    .flatMap { key, values ->
        values.collect { value -> tuple(key, value) }
    }
    .view()
```

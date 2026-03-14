(migrating-typed-operators)=

# Migrating to typed operators

Nextflow 26.04 introduces {ref}`typed workflows <workflow-typed-page>` for scripts that enable the `nextflow.preview.types` feature flag. Among other changes, typed workflows introduce a streamlined operator library that has first-class support for static typing and records.

While some operators have been preserved with no changes or minor changes, a few operators have significant changes, and many operators are not supported at all. This tutorial describes the changes to each operator, and demonstrates how to migrate each operator from a legacy workflow to a typed workflow.

See {ref}`operator-typed-page}` for the full set of operators that are supported in typed workflows.

## Supported with no changes

The following operators are supported in typed workflows with no changes.

- `reduce`
- `subscribe`
- `unique`
- `until`
- `view`

## Supported with minor changes

The following operators are supported in typed workflows with minor changes. These operators can usually be carried over with no changes, but may require minor modifications in certain cases.

### collect

The `collect` operator is supported with the following changes:

- The optional mapping closure is no longer supported -- use `map` instead

- The `flat` option is no longer supported and collected values are no longer flattened by default -- use `flatMap` instead

- The `sort` option is no longer supported -- use `map` with `Iterable::toSorted` instead

### filter

The `filter` operator is supported with one change: filtering by literal value, regular expression, or type qualifier is no longer supported. Use a closure instead.

```nextflow
// legacy workflow
ch.filter( ~/^a.*/ )

// typed workflow
ch.filter { v -> v == ~/^a.*/ }
```

### flatMap

The `flatMap` operator is supported with one change: the mapping closure *must* return a collection. Maps and tuples are not automatically flattened, and individual values must be wrapped in a list.

```nextflow
// legacy workflow
channel.of( 1, 2, 3 )
    .flatMap { n -> [ number: n, square: n*n, cube: n*n*n ] }
    .view { entry -> "${entry.key}: ${entry.value}" }

// typed workflow
channel.of( 1, 2, 3 )
    .flatMap { n -> [ tuple('number', n), tuple('square', n*n), tuple('cube', n*n*n) ] }
    .view { key, value -> "${key}: ${value}" }
```

### map

The `map` operator is supported with one change: `null` values are no longer discarded. Use `filter` instead.

```nextflow
// legacy workflow
ch.map { r -> r.id }

// typed workflow
ch.map { r -> r.id }.filter { id -> id != null }
```

### mix

The `mix` operator is supported with one change: multiple arguments are no longer supported. Use a separate `mix` for each new operand.

```nextflow
// legacy workflow
ch1.mix(ch2, ch3)

// typed workflow
ch1.mix(ch2).mix(ch3)
```

## Supported with major changes

The following operators are supported in typed workflows with significant changes in behavior. Review the following changes carefully if you use any of these operators.

### cross

The `cross` operator now implements a full cross product. It does not filter combinations by matching key.

The mapping closure is no longer supported. Use `map` instead.

### groupBy (groupTuple)

The `groupTuple` operator is replaced by `groupBy` in typed workflows.

Whereas `groupTuple` accepts tuples of arbitrary length, `groupBy` accepts either a 2-tuple of `(<key>, <value>)` or a 3-tuple of `(<key>, <size>, <value>)`. Specifying the group size with each input tuple provides the same behavior as using the `size` option (or wrapping each key with `groupKey()`) does with `groupTuple`.

Whereas `groupTuple` can group multiple lists in a group, `groupBy` always emits 2-tuples of the form `(<key>, <values>)`, where `<values>` is an unordered collection (`Bag`). This approach avoids a pitfall with `groupTuple` where the grouped lists can be ordered inconsistently.

### join

The `join` operator now implements a relational join (i.e. *SQL join*) on channels of records. Instead of joining tuples on a matching tuple element (e.g. `by: 0`), it joins records on a matching record field (e.g. `by: 'id'`).

The `join` operator does not support channels of tuples in typed workflows. You must refactor incoming channels to use records instead of tuples. See {ref}`migrating-records` for more information on how to migrate to records.

The `failOnDuplicate` option is no longer supported. Duplicate matches are handled by emitting each matching combination (like a relational join).

The `failOnMismatch` option is no longer supported. Unmatched values are either emitted or discarded depending on whether the `remainder` option is set. To fail on mismatches, use the `remainder` option and check for unmatched values in downstream logic.

```nextflow
// legacy workflow
left  = channel.of( ['X', 1], ['Y', 2], ['Z', 3], ['P', 7] )
right = channel.of( ['Z', 6], ['Y', 5], ['X', 4] )
left.join(right).view()

// typed workflow
left  = channel.of( record(id: 'X', a: 1), record(id: 'Y', a: 2), record(id: 'Z', a: 3), record(id: 'P', a: 7) )
right = channel.of( record(id: 'Z', b: 6), record(id: 'Y', b: 5), record(id: 'X', b: 4) )
left.join(right, by: 'id').view()
```

## Unsupported operators

The following operators are not supported in typed workflows. This section describes how to replace each unsupported operator when migrating to typed workflows.

:::{note}
Some operators are noted as being *non-deterministic*. This means that the operator depends on the ordering of values in the source channel. These operators are generally discouraged in favor of equivalent standard library functions because they can lead to {ref}`non-deterministic behavior <cache-nondeterministic-inputs>` if used improperly.
:::

### branch

Use `filter` and `map` for each branch instead. Using records instead of tuples can eliminate much of the need for `branch`.

Examples requiring only `filter`:

```nextflow
// before
ch_samples_by_type = ch_samples.branch { meta, cram, crai ->
    cram: cram.extension == 'cram'
    bam: cram.extension == 'bam'
}

ch_samples_by_type.cram.view()
ch_samples_by_type.bam.view()

// after
ch_samples_cram = ch_samples.filter { s -> s.cram.extension == 'cram' }
ch_samples_bam = ch_samples.filter { s -> s.cram.extension == 'bam' }

ch_samples_cram.view()
ch_samples_bam.view()
```

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

*These operators are non-deterministic.*

Use `List::collate()` instead.

### collectFile

The `collectFile` operator is useful for collecting intermediate results into a final output file, or writing a samplesheet. In many cases, `collectFile` can be replaced by a {ref}`workflow output <workflow-output-def>`, which can generate an index file for a published channel.

For other cases, consider the following alternatives:

- Use an `exec` process to write text files (see {ref}`working-with-files`)
- Use the `groupBy` operator to group  
- Use `Iterable::toSorted` to sort

You can compose these functions and operators as needed to achieve the desired functionality.

### combine

For uses of `combine` with the `by` option, use `join` instead. It implements the equivalent relational join.

For uses of `combine` without the `by` option, use `cross` instead. It implements the equivalent cross product.

### concat

Use `mix` instead.

### count, max, min, sum

Use `collect` and the corresponding `Iterable <stdlib-types-iterable>` methods instead.

### distinct

*This operator is non-deterministic.*

Use `unique` instead.

### dump

Use `view` instead. The `view` operator in typed workflows supports the same `tag` option, allowing it to be used like `dump`.

### first, last, take

*These operators are non-deterministic.*

These operators should not be used downstream of a process because they would essentially throw away computations. Instead, the workflow should be designed such that the unnecessary computations aren’t performed in the first place.

In other cases such as the beginning of a pipeline, use the corresponding {ref}`List <stdlib-types-list>` methods instead.

### flatten

Use `flatMap` instead.

### ifEmpty

The `ifEmpty` operator is typically used to either (1) raise an error if a channel is empty or (2) provide a fallback for a null dataflow value.

Both cases can be implemented in typed workflows without `ifEmpty`:

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

*This operator is non-deterministic.*

Use `join` to combine channels in a way that doesn’t depend on the channel ordering.

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

*This operator is non-deterministic.*

The `randomSample` operator has no equivalent in typed workflows. In general, this operator should not be used because it leads to non-deterministic behavior.

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
log2 = log1..map { v -> v * 2 }
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

### countFasta, countFastq, countJson, countLines

Use the equivalent {ref}`stdlib-types-path` method with `flatMap` instead:

```nextflow
// before
channel.fromPath('sample.fastq')
  .countFastq( by: 10 )
  .view()

// after
channel.fromPath('sample.fastq')
  .flatMap { fastq -> fastq.countFastq(by: 10) }
  .view()
```

### toInteger, toLong, toFloat, toDouble

Use the equivalent {ref}`stdlib-types-string` method with `map` instead:

```nextflow
// before
channel.of('1', '7', '12')
    .toInteger()
    .view()

// after
channel.of('1', '7', '12')
    .map { v -> v.toInteger() }
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

(operator-page)=

# Operators

Nextflow **operators** are methods that allow you to manipulate channels. Every operator, with the exception of [set](#set) and [subscribe](#subscribe), produces one or more new channels, allowing you to chain operators to fit your needs.

This page is a comprehensive reference for all Nextflow operators. However, if you are new to Nextflow, here are some suggested operators to learn for common use cases:

- Filtering: [filter](#filter), [randomSample](#randomsample), [take](#take), [unique](#unique)
- Reduction: [collect](#collect), [groupTuple](#grouptuple), [reduce](#reduce)
- Text processing: [splitCsv](#splitcsv), [splitJson](#splitjson), [splitText](#splittext)
- Combining channels: [combine](#combine), [concat](#concat), [join](#join), [mix](#mix)
- Forking channels: [branch](#branch), [multiMap](#multimap)
- Maths: [count](#count), [max](#max), [min](#min), [sum](#sum)
- Other: [ifEmpty](#ifempty), [map](#map), [set](#set), [view](#view)

(operator-branch)=

## branch

:::{versionadded} 19.08.0-edge
:::

*Returns: map of queue channels*

The `branch` operator forwards each item from a source channel to one of multiple output channels, based on a selection criteria.

The selection criteria is a {ref}`closure <script-closure>` that defines, for each output channel, a unique label followed by a boolean expression. When an item is received, it is routed to the first output channel whose expression evaluates to `true`. For example:

```{literalinclude} snippets/branch.nf
:language: groovy
```

```{literalinclude} snippets/branch.out
:language: console
```

:::{note}
The above output may be printed in any order since the two `view` operations are executed asynchronously.
:::

A fallback condition can be specified using `true` as the last branch condition:

```{literalinclude} snippets/branch-with-fallback.nf
:language: groovy
```

```{literalinclude} snippets/branch-with-fallback.out
:language: console
```

The value emitted to each branch can be customized with an expression statement (or statements) after the branch condition:

```{literalinclude} snippets/branch-with-mapper.nf
:language: groovy
```

```{literalinclude} snippets/branch-with-mapper.out
:language: console
```

:::{tip}
When the `return` keyword is omitted, the value of the last expression statement is implicitly returned.
:::

The `branchCriteria()` method can be used to create a branch criteria as a variable that can be passed as an argument to any number of `branch` operations, as shown below:

```{literalinclude} snippets/branch-criteria.nf
:language: groovy
```

```{literalinclude} snippets/branch-criteria.out
:language: console
```

## buffer

*Returns: queue channel*

The `buffer` operator collects items from a source channel into subsets and emits each subset separately.

This operator has multiple variants:

`buffer( closingCondition )`

: Emits each subset when `closingCondition` is satisfied. The closing condition can be a literal value, a {ref}`regular expression <script-regexp>`, a type qualifier (i.e. Java class), or a boolean predicate. For example:

  ```{literalinclude} snippets/buffer-with-closing.nf
  :language: groovy
  ```

  ```{literalinclude} snippets/buffer-with-closing.out
  :language: console
  ```

`buffer( openingCondition, closingCondition )`

: Creates a new subset when `openingCondition` is satisfied and emits the subset when is `closingCondition` is satisfied. The opening and closing conditions can each be a literal value, a {ref}`regular expression <script-regexp>`, a type qualifier (i.e. Java class), or a boolean predicate. For example:

  ```{literalinclude} snippets/buffer-with-opening-closing.nf
  :language: groovy
  ```

  ```{literalinclude} snippets/buffer-with-opening-closing.out
  :language: console
  ```

`buffer( size: n )`

: Emits a new subset for every `n` items. Remaining items are discarded. For example:

  ```{literalinclude} snippets/buffer-with-size.nf
  :language: groovy
  ```

  ```{literalinclude} snippets/buffer-with-size.out
  :language: console
  ```

  The `remainder` option can be used to emit any remaining items as a partial subset:

  ```{literalinclude} snippets/buffer-with-size-remainder.nf
  :language: groovy
  ```

  ```{literalinclude} snippets/buffer-with-size-remainder.out
  :language: console
  ```

`buffer( size: n, skip: m )`

: Emits a new subset for every `n` items, skipping `m` items before collecting each subset. For example:

  ```{literalinclude} snippets/buffer-with-size-skip.nf
  :language: groovy
  ```

  ```{literalinclude} snippets/buffer-with-size-skip.out
  :language: console
  ```

  The `remainder` option can be used to emit any remaining items as a partial subset.

See also: [collate](#collate)

## collate

*Returns: queue channel*

The `collate` operator collects items from a source channel into groups of *N* items.

This operator has multiple variants:

`collate( size, remainder = true )`

: Collects items into groups of `size` items:

  ```{literalinclude} snippets/collate.nf
  :language: groovy
  ```

  ```{literalinclude} snippets/collate.out
  :language: console
  ```

  By default, any remaining items are emitted as a partial group. You can specify `false` as the second parameter to discard them instead:

  ```{literalinclude} snippets/collate-with-no-remainder.nf
  :language: groovy
  ```

  ```{literalinclude} snippets/collate-with-no-remainder.out
  :language: console
  ```

  :::{note}
  This version of `collate` is equivalent to `buffer( size: n, remainder: true | false )`.
  :::

`collate( size, step, remainder = true )`

: Collects items into groups of `size` items using a *sliding window* that moves by `step` items at a time:

  ```{literalinclude} snippets/collate-with-step.nf
  :language: groovy
  ```

  ```{literalinclude} snippets/collate-with-step.out
  :language: console
  ```

  You can specify `false` as the third parameter to discard any remaining items.

See also: [buffer](#buffer)

(operator-collect)=

## collect

*Returns: value channel*

The `collect` operator collects all items from a source channel into a list and emits it as a single item:

```{literalinclude} snippets/collect.nf
:language: groovy
```

```{literalinclude} snippets/collect.out
:language: console
```

An optional {ref}`closure <script-closure>` can be used to transform each item before it is collected:

```{literalinclude} snippets/collect-with-mapper.nf
:language: groovy
```

```{literalinclude} snippets/collect-with-mapper.out
:language: console
```

Available options:

`flat`
: When `true`, nested list structures are flattened and their items are collected individually (default: `true`).

`sort`
: When `true`, the collected items are sorted by their natural ordering (default: `false`). Can also be a {ref}`closure <script-closure>` or a [Comparator](https://docs.oracle.com/javase/8/docs/api/java/util/Comparator.html) which defines how items are compared during sorting.

See also: [toList](#tolist), [toSortedList](#tosortedlist)

## collectFile

*Returns: queue channel*

The `collectFile` operator collects the items from a source channel and saves them to one or more files, emitting the collected file(s).

This operator has multiple variants:

`collectFile( name: '...', options = [:] )`

: Collects the items and saves them to a single file specified by the `name` option:

  ```{literalinclude} snippets/collectfile.nf
  :language: groovy
  ```

`collectFile( closure, options = [:] )`

: Collects the items into groups and saves each group to a file, using a grouping criteria. The grouping criteria is a {ref}`closure <script-closure>` that maps each item to a pair, where the first element is the file name for the group and the second element is the content to be appended to that file. For example:

  ```{literalinclude} snippets/collectfile-closure.nf
  :language: groovy
  ```

  ```{literalinclude} snippets/collectfile-closure.out
  :language: console
  ```

  When the items from the source channel are files, the grouping criteria can be omitted. In this case, the items will be grouped by their source filename.

The following example shows how to use a closure to collect and sort all sequences in a FASTA file from shortest to longest:

```groovy
Channel
    .fromPath('/data/sequences.fa')
    .splitFasta( record: [id: true, sequence: true] )
    .collectFile( name: 'result.fa', sort: { it.size() } ) {
        it.sequence
    }
    .view { it.text }
```

:::{warning}
The `collectFile` operator needs to store files in a temporary directory that is automatically deleted on workflow completion. For performance reasons, this directory is located in the machine's local storage, and it should have as much free space as the data that is being collected. The `tempDir` option can be used to specify a different temporary directory.
:::

Available options:

`cache`
: Controls the caching ability of the `collectFile` operator when using the *resume* feature. It follows the same semantic of the {ref}`process-cache` directive (default: `true`).

`keepHeader`
: Prepend the resulting file with the header fetched in the first collected file. The header size (ie. lines) can be specified by using the `skip` option (default: `0`), to determine how many lines to remove from all collected files except for the first (where no lines will be removed).

`name`
: Name of the file where all received values are stored.

`newLine`
: Appends a `newline` character automatically after each entry (default: `false`).

`seed`
: A value or a map of values used to initialize the files content.

`skip`
: Skip the first `n` lines e.g. `skip: 1` (default: `0`).

`sort`
: Defines sorting criteria of content in resulting file(s). Can be one of the following values:

  - `false`: Disable content sorting. Entries are appended as they are produced.
  - `true`: Order the content by the entry's natural ordering i.e. numerical for number, lexicographic for string, etc. See the [Java documentation](http://docs.oracle.com/javase/tutorial/collections/interfaces/order.html) for more information.
  - `'index'`: Order the content by the incremental index number assigned to each entry while they are collected.
  - `'hash'`: (default) Order the content by the hash number associated to each entry
  - `'deep'`: Similar to the previous, but the hash number is created on actual entries content e.g. when the entry is a file the hash is created on the actual file content.
  - A custom sorting criteria can be specified with a {ref}`Closure <script-closure>` or a [Comparator](http://docs.oracle.com/javase/7/docs/api/java/util/Comparator.html) object.

  The file content is sorted in such a way that it does not depend on the order in which entries were added to it, which guarantees that it is consistent (i.e. does not change) across different executions with the same data.

`storeDir`
: Folder where the resulting file(s) are stored.

`tempDir`
: Folder where temporary files, used by the collecting process, are stored.

(operator-combine)=

## combine

*Returns: queue channel*

The `combine` operator produces the combinations (i.e. cross product, "Cartesian" product) of two source channels, or a channel and a list (as the right operand), emitting each combination separately.

For example:

```{literalinclude} snippets/combine.nf
:language: groovy
```

```{literalinclude} snippets/combine.out
:language: console
```

The `by` option can be used to combine items that share a matching key. The value should be the zero-based index of the tuple, or a list of indices. For example:

```{literalinclude} snippets/combine-by.nf
:language: groovy
```

```{literalinclude} snippets/combine-by.out
:language: console
```

:::{note}
The `combine` operator is similar to `cross` and `join`, making them easy to confuse. Their differences can be summarized as follows:

- `combine` and `cross` both produce an *outer product* or *cross product*, whereas `join` produces an *inner product*.

- `combine` filters pairs with a matching key only if the `by` option is used, whereas `cross` always filters pairs with a matching key.

- `combine` with the `by` option merges and flattens each pair, whereas `cross` does not. Compare the examples for `combine` and `cross` to see this difference.
:::

See also: [cross](#cross), [join](#join)

(operator-concat)=

## concat

*Returns: queue channel*

The `concat` operator emits the items from two or more source channels into a single output channel. Each source channel is emitted in the order in which it was specified.

In other words, given *N* channels, the items from the *i+1*-th channel are emitted only after all of the items from the *i*-th channel have been emitted.

For example:

```{literalinclude} snippets/concat.nf
:language: groovy
```

```{literalinclude} snippets/concat.out
:language: console
```

See also: [mix](#mix)

(operator-count)=

## count

*Returns: value channel*

The `count` operator computes the total number of items in a source channel and emits it:

```{literalinclude} snippets/count.nf
:language: groovy
```

```{literalinclude} snippets/count.out
:language: console
```

An optional filter can be provided to select which items to count. The selection criteria can be a literal value, a {ref}`regular expression <script-regexp>`, a type qualifier (i.e. Java class), or a boolean predicate. For example:

```{literalinclude} snippets/count-with-filter-number.nf
:language: groovy
```

```{literalinclude} snippets/count-with-filter-number.out
:language: console
```

```{literalinclude} snippets/count-with-filter-regex.nf
:language: groovy
```

```{literalinclude} snippets/count-with-filter-regex.out
:language: console
```

```{literalinclude} snippets/count-with-filter-closure.nf
:language: groovy
```

```{literalinclude} snippets/count-with-filter-closure.out
:language: console
```

(operator-countfasta)=

## countFasta

*Returns: value channel*

Counts the total number of records in a channel of FASTA files, equivalent to `splitFasta | count`. See [splitFasta](#splitfasta) for the list of available options.

(operator-countfastq)=

## countFastq

*Returns: value channel*

Counts the total number of records in a channel of FASTQ files, equivalent to `splitFastq | count`. See [splitFastq](#splitfastq) for the list of available options.

(operator-countjson)=

## countJson

*Returns: value channel*

Counts the total number of records in a channel of JSON files, equivalent to `splitJson | count`. See [splitJson](#splitjson) for the list of available options.

(operator-countlines)=

## countLines

*Returns: value channel*

Counts the total number of lines in a channel of text files, equivalent to `splitText | count`. See [splitLines](#splittext) for the list of available options.

(operator-cross)=

## cross

*Returns: queue channel*

The `cross` operator emits every pairwise combination of two channels for which the pair has a matching key.

By default, the key is defined as the first entry in a list or map, or the value itself for any other data type. For example:

```{literalinclude} snippets/cross.nf
:language: groovy
```

```{literalinclude} snippets/cross.out
:language: console
```

An optional closure can be used to define the matching key for each item:

```{literalinclude} snippets/cross-with-mapper.nf
:language: groovy
```

```{literalinclude} snippets/cross-with-mapper.out
:language: console
```

There are two important caveats when using the `cross` operator:

1. The operator is not *commutative*, i.e. `a.cross(b)` is not the same as `b.cross(a)`
2. Each source channel should not emit any items with duplicate keys, i.e. each item should have a unique key.

See also: [combine](#combine)

## distinct

*Returns: queue channel*

The `distinct` operator forwards a source channel with *consecutively* repeated items removed, such that each emitted item is different from the preceding one:

```{literalinclude} snippets/distinct.nf
:language: groovy
```

```{literalinclude} snippets/distinct.out
:language: console
```

An optional {ref}`closure <script-closure>` can be used to transform each value before it is evaluated for distinct-ness:

```{literalinclude} snippets/distinct-with-mapper.nf
:language: groovy
```

```{literalinclude} snippets/distinct-with-mapper.out
:language: console
```

See also: [unique](#unique)

(operator-dump)=

## dump

*Returns: queue channel or value channel, depending on the input*

The `dump` operator prints each item in a source channel when the pipeline is executed with the `-dump-channels` command-line option, otherwise it does nothing. It is a useful way to inspect and debug channels quickly without having to modify the pipeline script.

The `tag` option can be used to select which channels to dump:

```{literalinclude} snippets/dump.nf
:language: groovy
```

Then, you can run your pipeline with `-dump-channels foo` or `-dump-channels bar` to dump the content of either channel. Multiple tag names can be specified as a comma-separated list.

Available options:

`pretty`
: :::{versionadded} 22.10.0
  :::
: When `true`, format the output as pretty-printed JSON (default: `false`).

`tag`
: Associate the channel with a tag that can be specified with the `-dump-channels` option to select which channels to dump.

## filter

*Returns: queue channel*

The `filter` operator emits the items from a source channel that satisfy a condition, discarding all other items. The filter condition can be a literal value, a {ref}`regular expression <script-regexp>`, a type qualifier (i.e. Java class), or a boolean predicate.

The following example filters a channel with a regular expression that only matches strings beginning with `a`:

```{literalinclude} snippets/filter-regex.nf
:language: groovy
```

```{literalinclude} snippets/filter-regex.out
:language: console
```

The following example filters a channel with the `Number` type qualifier so that only numbers are emitted:

```{literalinclude} snippets/filter-type.nf
:language: groovy
```

```{literalinclude} snippets/filter-type.out
:language: console
```

The following example filters a channel using a boolean predicate, which is a {ref}`closure <script-closure>` that returns a boolean value. In this case, the predicate is used to select only odd numbers:

```{literalinclude} snippets/filter-closure.nf
:language: groovy
```

```{literalinclude} snippets/filter-closure.out
:language: console
```

(operator-first)=

## first

*Returns: value channel*

The `first` operator emits the first item in a source channel, or the first item that matches a condition. The condition can be a {ref}`regular expression<script-regexp>`, a type qualifier (i.e. Java class), or a boolean predicate. For example:

```{literalinclude} snippets/first.nf
:language: groovy
```

(operator-flatmap)=

## flatMap

*Returns: queue channel*

The `flatMap` operator applies a *mapping function* to each item from a source channel.

When the mapping function returns a list, each element in the list is emitted separately:

```{literalinclude} snippets/flatmap-list.nf
:language: groovy
```

```{literalinclude} snippets/flatmap-list.out
:language: console
```

When the mapping function returns a map, each key-value pair in the map is emitted separately:

```{literalinclude} snippets/flatmap-map.nf
:language: groovy
```

```{literalinclude} snippets/flatmap-map.out
:language: console
```

(operator-flatten)=

## flatten

*Returns: queue channel*

The `flatten` operator flattens each item from a source channel that is a list or other collection, such that each element in each collection is emitted separately:

```{literalinclude} snippets/flatten.nf
:language: groovy
```

```{literalinclude} snippets/flatten.out
:language: console
```

As shown in the above example, deeply nested collections are also flattened.

See also: [flatMap](#flatmap)

(operator-grouptuple)=

## groupTuple

*Returns: queue channel*

The `groupTuple` operator collects lists (i.e. *tuples*) from a source channel into groups based on a grouping key. A new tuple is emitted for each distinct key.

To be more precise, the operator transforms a sequence of tuples like *(K, V, W, ..)* into a sequence of tuples like *(K, list(V), list(W), ..)*.

For example:

```{literalinclude} snippets/grouptuple.nf
:language: groovy
```

```{literalinclude} snippets/grouptuple.out
:language: console
```

By default, the first element of each tuple is used as the grouping key. The `by` option can be used to specify a different index, or list of indices. For example, to group by the second element of each tuple:

```{literalinclude} snippets/grouptuple-by.nf
:language: groovy
```

```{literalinclude} snippets/grouptuple-by.out
:language: console
```

By default, if you don't specify a size, the `groupTuple` operator will not emit any groups until *all* inputs have been received. If possible, you should always try to specify the number of expected elements in each group using the `size` option, so that each group can be emitted as soon as it's ready. In cases where the size of each group varies based on the grouping key, you can use the built-in `groupKey()` function, which allows you to define a different expected size for each group:

```{literalinclude} snippets/grouptuple-groupkey.nf
:language: groovy
```

```{literalinclude} snippets/grouptuple-groupkey.out
:language: console
```

Available options:

`by`
: The zero-based index of the element to use as the grouping key. Can also be a list of indices, e.g. `by: [0,2]` (default: `[0]`).

`remainder`
: When `true`, incomplete tuples (i.e. groups with less than `size` items) are emitted as partial groups, otherwise they are discarded (default: `false`). This option can only be used with `size`.

`size`
: The required number of items for each group. When a group reaches the required size, it is emitted.

`sort`
: Defines the sorting criteria for the grouped items. Can be one of the following values:

  - `false`: No sorting is applied (default).
  - `true`: Order the grouped items by the item's natural ordering i.e. numerical for number, lexicographic for string, etc. See the [Java documentation](http://docs.oracle.com/javase/tutorial/collections/interfaces/order.html) for more information.
  - `'hash'`: Order the grouped items by the hash number associated to each entry.
  - `'deep'`: Similar to the previous, but the hash number is created on actual entries content e.g. when the item is a file, the hash is created on the actual file content.
  - A custom sorting criteria used to order the nested list elements of each tuple. It can be a {ref}`Closure <script-closure>` or a [Comparator](http://docs.oracle.com/javase/7/docs/api/java/util/Comparator.html) object.

(operator-ifempty)=

## ifEmpty

*Returns: value channel*

The `ifEmpty` operator emits a source channel, or a default value if the source channel is *empty* (doesn't emit any value):

```{literalinclude} snippets/ifempty-1.nf
:language: groovy
```

```{literalinclude} snippets/ifempty-1.out
:language: console
```

```{literalinclude} snippets/ifempty-2.nf
:language: groovy
```

```{literalinclude} snippets/ifempty-2.out
:language: console
```

The default value can also be a {ref}`closure <script-closure>`, in which case the closure is evaluated and the result is emitted when the source channel is empty.

See also: {ref}`channel-empty` channel factory

(operator-join)=

## join

*Returns: queue channel*

The `join` operator emits the inner product of two source channels using a matching key.

To be more precise, the operator transforms a sequence of tuples like *(K, V1, V2, ..)* and *(K, W1, W1, ..)* into a sequence of tuples like *(K, V1, V2, .., W1, W2, ..)*. It is equivalent to an *inner join* in SQL, or an *outer join* when `remainder` is `true`.

For example:

```{literalinclude} snippets/join.nf
:language: groovy
```

```{literalinclude} snippets/join.out
:language: console
```

By default, the first element of each item is used as the key. The `by` option can be used to specify a different index, or list of indices.

By default, unmatched items are discarded. The `remainder` option can be used to emit them at the end:

```{literalinclude} snippets/join-with-remainder.nf
:language: groovy
```

```{literalinclude} snippets/join-with-remainder.out
:language: console
```

Available options:

`by`
: The zero-based index of each item to use as the matching key. Can also be a list of indices, e.g. `by: [0, 2]` (default: `[0]`).

`failOnDuplicate`
: When `true`, an error is reported when the operator receives multiple items from the same channel with the same key (default: `true` if {ref}`strict mode <config-feature-flags>` is enabled, `false` otherwise).

`failOnMismatch`
: When `true`, an error is reported when the operator receives an item from one channel for which there no matching item from the other channel (default: `true` if {ref}`strict mode <config-feature-flags>` is enabled, `false` otherwise). This option cannot be used with `remainder`.

`remainder`
: When `true`, unmatched items are emitted at the end, otherwise they are discarded (default: `false`). 

See also: [combine](#combine), [cross](#cross)

(operator-last)=

## last

*Returns: value channel*

The `last` operator emits the last item from a source channel:

```{literalinclude} snippets/last.nf
:language: groovy
```

```{literalinclude} snippets/last.out
:language: console
```

(operator-map)=

## map

*Returns: queue channel*

The `map` operator applies a *mapping function* to each item from a source channel:

```{literalinclude} snippets/map.nf
:language: groovy
```

```{literalinclude} snippets/map.out
:language: console
```

(operator-max)=

## max

*Returns: value channel*

The `max` operator emits the item with the greatest value from a source channel:

```{literalinclude} snippets/max.nf
:language: groovy
```

```{literalinclude} snippets/max.out
:language: console
```

An optional {ref}`closure <script-closure>` can be used to control how the items are compared. The closure can be a *mapping function*, which transforms each item before it is compared, or a *comparator function*, which defines how to compare two items more generally.

The following examples show how to find the longest string in a channel:

```{literalinclude} snippets/max-with-mapper.nf
:language: groovy
```

```{literalinclude} snippets/max-with-mapper.out
:language: console
```

```{literalinclude} snippets/max-with-comparator.nf
:language: groovy
```

```{literalinclude} snippets/max-with-comparator.out
:language: console
```

(operator-merge)=

## merge

*Returns: queue channel*

The `merge` operator joins the items from two or more channels into a new channel:

```{literalinclude} snippets/merge.nf
:language: groovy
```

```{literalinclude} snippets/merge.out
:language: console
```

An optional closure can be used to control how two items are merged:

```{literalinclude} snippets/merge-with-mapper.nf
:language: groovy
```

```{literalinclude} snippets/merge-with-mapper.out
:language: console
```

:::{danger}
In general, the use of the `merge` operator is discouraged. Processes and channel operators are not guaranteed to emit items in the order that they were received, as they are executed concurrently. Therefore, if you try to merge output channels from different processes, the resulting channel may be different on each run, which will cause resumed runs to {ref}`not work properly <cache-nondeterministic-inputs>`.

You should always use a matching key (e.g. sample ID) to merge multiple channels, so that they are combined in a deterministic way. For this purpose, you can use the [join](#join) operator.
:::

(operator-min)=

## min

*Returns: value channel*

The `min` operator emits the item with the lowest value from a source channel:

```{literalinclude} snippets/min.nf
:language: groovy
```

```{literalinclude} snippets/min.out
:language: console
```

An optional {ref}`closure <script-closure>` can be used to control how the items are compared. The closure can be a *mapping function*, which transforms each item before it is compared, or a *comparator function*, which defines how to compare two items more generally.

The following examples show how to find the shortest string in a channel:

```{literalinclude} snippets/min-with-mapper.nf
:language: groovy
```

```{literalinclude} snippets/min-with-mapper.out
:language: console
```

```{literalinclude} snippets/min-with-comparator.nf
:language: groovy
```

```{literalinclude} snippets/min-with-comparator.out
:language: console
```

(operator-mix)=

## mix

*Returns: queue channel*

The `mix` operator emits the items from two or more source channels into a single output channel:

```{literalinclude} snippets/mix.nf
:language: groovy
```

```{literalinclude} snippets/mix.out
:language: console
```

The items in the mixed output channel may appear in any order, regardless of which source channel they came from. Thus, the previous example could also output the following:

```console
z
1
a
2
b
3
```

See also: [concat](#concat)

(operator-multimap)=

## multiMap

:::{versionadded} 19.11.0-edge
:::

*Returns: map of queue channels*

The `multiMap` operator applies a set of mapping functions to a source channel, producing a separate output channel for each mapping function.

The multi-map criteria is a {ref}`closure <script-closure>` that defines, for each output channel, a label followed by a mapping expression.

For example:

```{literalinclude} snippets/multimap.nf
:language: groovy
```

```{literalinclude} snippets/multimap.out
:language: console
```

Multiple labels can share the same mapping expression using the following shorthand:

```{literalinclude} snippets/multimap-shared.nf
:language: groovy
```

```{literalinclude} snippets/multimap-shared.out
:language: console
```

The above example creates two channels as before, but now they both receive the same items.

You can use the `multiMapCriteria()` method to create a multi-map criteria as a variable that can be passed as an argument to any number of `multiMap` operations, as shown below:

```{literalinclude} snippets/multimap-criteria.nf
:language: groovy
```

:::{note}
If you use `multiMap` to split a tuple or map into multiple channels, it is recommended that you retain a matching key (e.g. sample ID) with *each* new channel, so that you can re-combine these channels later on if needed. In general, you should not expect to be able to merge channels correctly without a matching key, due to the concurrent nature of Nextflow pipelines.
:::

(operator-randomsample)=

## randomSample

*Returns: queue channel*

The `randomSample` operator emits a randomly-selected subset of items from a source channel:

```{literalinclude} snippets/random-sample.nf
:language: groovy
```

The above snippet will print 10 randomly-selected numbers between 1 and 100 (without replacement).

An optional second parameter can be used to set the initial *seed* for the random number generator, which ensures that the `randomSample` operator produces the same pseudo-random sequence across runs:

```{literalinclude} snippets/random-sample-with-seed.nf
:language: groovy
```

The above example will print 10 randomly-selected numbers between 1 and 100 (without replacement). Each subsequent script execution will produce the same sequence.

(operator-reduce)=

## reduce

*Returns: value channel*

The `reduce` operator applies an *accumulator function* sequentially to each item in a source channel, and emits the final accumulated value. The accumulator function takes two parameters -- the accumulated value and the *i*-th emitted item -- and it should return the accumulated result, which is passed to the next invocation with the *i+1*-th item. This process is repeated for each item in the source channel.

For example:

```{literalinclude} snippets/reduce.nf
:language: groovy
```

```{literalinclude} snippets/reduce.out
:language: console
```

By default, the first item is used as the initial accumulated value. You can optionally specify a different initial value as shown below:

```{literalinclude} snippets/reduce-with-initial-value.nf
:language: groovy
```

```{literalinclude} snippets/reduce-with-initial-value.out
:language: console
```

(operator-set)=

## set

*Returns: nothing*

The `set` operator assigns a source channel to a variable, whose name is specified as a closure parameter:

```groovy
Channel.of(10, 20, 30).set { my_channel }
```

Using `set` is semantically equivalent to assigning a variable:

```groovy
my_channel = Channel.of(10, 20, 30)
```

See also: [tap](#tap)

(operator-splitcsv)=

## splitCsv

*Returns: queue channel*

The `splitCsv` operator parses and splits [CSV-formatted](http://en.wikipedia.org/wiki/Comma-separated_values) text from a source channel into records, or groups of records with a given size.

For example:

```{literalinclude} snippets/splitcsv.nf
:language: groovy
```

```{literalinclude} snippets/splitcsv.out
:language: console
```

The above example shows hows CSV text is parsed and split into individual rows, where each row is simply a list of columns.

When the CSV begins with a header line defining the column names, and the `header` option is `true`, each row is returned as a map instead:

```{literalinclude} snippets/splitcsv-with-header.nf
:language: groovy
```

```{literalinclude} snippets/splitcsv-with-header.out
:language: console
```

The `header` option can also just be a list of columns:

```{literalinclude} snippets/splitcsv-with-columns.nf
:language: groovy
```

```{literalinclude} snippets/splitcsv-with-columns.out
:language: console
```

Available options:

`by`
: When specified, group rows into *chunks* with the given size (default: none).

`charset`
: Parse the content with the specified charset, e.g. `UTF-8`. See the list of [standard charsets](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/nio/charset/StandardCharsets.html) for available options.

`decompress`
: When `true`, decompress the content using the GZIP format before processing it (default: `false`). Files with the `.gz` extension are decompressed automatically.

`elem`
: The index of the element to split when the source items are lists or tuples (default: first file object or first element).

`header`
: When `true`, the first line is used as the columns names (default: `false`). Can also be a list of columns names.

`limit`
: Limits the number of records to retrieve for each source item (default: no limit).

`quote`
: The character used to quote values (default: `''` or `""`).

`sep`
: The character used to separate values (default: `,`)

`skip`
: Number of lines to ignore from the beginning when parsing the CSV text (default: `0`).

`strip`
: When `true`, remove leading and trailing blanks from values (default: `false`).

(operator-splitfasta)=

## splitFasta

*Returns: queue channel*

The `splitFasta` operator splits [FASTA-formatted](http://en.wikipedia.org/wiki/FASTA_format) text from a source channel into individual sequences.

The `by` option can be used to group sequences into chunks of a given size. The following example shows how to read a FASTA file and split it into chunks of 10 sequences each:

```groovy
Channel
     .fromPath('misc/sample.fa')
     .splitFasta( by: 10 )
     .view()
```

:::{warning}
Chunks are stored in memory by default. When splitting large files, specify `file: true` to save the chunks into files in order to avoid running out of memory. See the list of options below for details.
:::

The `record` option can be used to split FASTA content into *records* instead of text chunks. Each record is a map that allows you to access the FASTA sequence data with ease. For example:

```groovy
Channel
     .fromPath('misc/sample.fa')
     .splitFasta( record: [id: true, seqString: true] )
     .filter { record -> record.id =~ /^ENST0.*/ }
     .view { record -> record.seqString }
```

The above example loads the `misc/sample.fa` file, splits it into records containing the `id` and `seqString` fields (i.e. the sequence id and the sequence data), filters records by their ID, and finally prints the sequence string of each record.

Available options:

`by`
: Defines the number of sequences in each chunk (default: `1`).

`charset`
: Parse the content with the specified charset, e.g. `UTF-8`. See the list of [standard charsets](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/nio/charset/StandardCharsets.html) for available options.

`compress`
: When `true`, resulting file chunks are GZIP compressed (default: `false`). The `.gz` suffix is automatically added to chunk file names.

`decompress`
: When `true`, decompress the content using the GZIP format before processing it (default: `false`). Files with the `.gz` extension are decompressed automatically.

`elem`
: The index of the element to split when the source items are lists or tuples (default: first file object or first element).

`file`
: When `true`, saves each split to a file. Use a string instead of `true` value to create split files with a specific name (split index number is automatically added). Finally, set this attribute to an existing directory, in order to save the split files into the specified directory.

`limit`
: Limits the number of sequences to retrieve for each source item (default: no limit).

`record`
: Parse each entry in the FASTA file into a record. The following fields are available:

  - `id`: The FASTA sequence identifier, i.e. the word following the `>` symbol up to the first blank or newline character
  - `header`: The first line in a FASTA sequence without the `>` character
  - `desc`: The text in the FASTA header following the ID value
  - `text`: The complete FASTA sequence including the header
  - `seqString`: The sequence data as a single-line string, i.e. containing no newline characters
  - `sequence`: The sequence data as a multi-line string, i.e. always ending with a newline character
  - `width`: Define the length of a single line when the `sequence` field is used, after which the sequence data continues on a new line.

`size`
: Defines the size of the expected chunks as a memory unit, e.g. `1.MB`.

See also: [countFasta](#countfasta)

(operator-splitfastq)=

## splitFastq

*Returns: queue channel*

The `splitFasta` operator splits [FASTQ formatted](http://en.wikipedia.org/wiki/FASTQ_format) text from a source channel into individual sequences.

The `by` option can be used to group sequences into chunks of a given size. The following example shows how to read a FASTQ file and split it into chunks of 10 sequences each:

```groovy
Channel
    .fromPath('misc/sample.fastq')
    .splitFastq( by: 10 )
    .view()
```

:::{warning}
Chunks are stored in memory by default. When splitting large files, specify `file: true` to save the chunks into files in order to avoid running out of memory. See the list of options below for details.
:::

The `record` option can be used to split FASTQ content into *records* instead of text chunks. Each record is a map that allows you to access the FASTQ sequence data with ease. For example:

```groovy
Channel
    .fromPath('misc/sample.fastq')
    .splitFastq( record: true )
    .view { record -> record.readHeader }
```

The `pe` option can be used to split paired-end FASTQ files. The source channel must emit tuples containing the file pairs. For example:

```groovy
Channel
    .fromFilePairs('/my/data/SRR*_{1,2}.fastq', flat: true)
    .splitFastq(by: 100_000, pe: true, file: true)
    .view()
```

:::{note}
`Channel.fromFilePairs()` requires the `flat: true` option in order to emit the file pairs as separate elements in the produced tuples.
:::

:::{note}
This operator assumes that the order of the paired-end reads correspond with each other and that both files contain the same number of reads.
:::

Available options:

`by`
: Defines the number of sequences in each chunk (default: `1`).

`charset`
: Parse the content with the specified charset, e.g. `UTF-8`. See the list of [standard charsets](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/nio/charset/StandardCharsets.html) for available options.

`compress`
: When `true`, resulting file chunks are GZIP compressed (default: `false`). The `.gz` suffix is automatically added to chunk file names.

`decompress`
: When `true`, decompress the content using the GZIP format before processing it (default: `false`). Files with the `.gz` extension are decompressed automatically.

`elem`
: The index of the element to split when the source items are lists or tuples (default: first file object or first element).

`file`
: When `true`, saves each split to a file. Use a string instead of `true` value to create split files with a specific name (split index number is automatically added). Finally, set this attribute to an existing directory, in order to save the split files into the specified directory.

`limit`
: Limits the number of sequences to retrieve for each source item (default: no limit).

`pe`
: When `true`, splits paired-end read files. Items emitted by the source channel must be tuples with the file pairs.

`record`
: Parse each entry in the FASTQ file into a record. The following fields are available:

  - `readHeader`: Sequence header (without the `@` prefix)
  - `readString`: The raw sequence data
  - `qualityHeader`: Base quality header (it may be empty)
  - `qualityString`: Quality values for the sequence

See also: [countFastq](#countfastq)

(operator-splitjson)=

## splitJson

*Returns: queue channel*

The `splitJson` operator splits [JSON formatted](https://en.wikipedia.org/wiki/JSON) text from a source channel into individual records.

If the source item is a JSON array, each element of the array will be emitted:

```{literalinclude} snippets/splitjson-array.nf
:language: groovy
```

```{literalinclude} snippets/splitjson-array.out
:language: console
```

If the source item is a JSON object, each key-value pair will be emitted as a map with the properties `key`  and `value`:

```{literalinclude} snippets/splitjson-object.nf
:language: groovy
```

```{literalinclude} snippets/splitjson-object.out
:language: console
```

The `path` option can be used to query a section of the JSON document to parse and split:

```{literalinclude} snippets/splitjson-with-path.nf
:language: groovy
```

```{literalinclude} snippets/splitjson-with-path.out
:language: console
```

Available options:

`limit`
: Limits the number of records to retrieve for each source item (default: no limit).

`path`
: Defines a query for a section of each source item to parse and split. The expression should be a path similar to [JSONPath](https://goessner.net/articles/JsonPath/). The empty string is the document root (default). An integer in brackets is a zero-based index in a JSON array. A string preceded by a dot `.` is a key in a JSON object.

See also: [countJson](#countjson)

(operator-splittext)=

## splitText

*Returns: queue channel*

The `splitText` operator splits multi-line text content from a source channel into chunks of *N* lines:

```groovy
Channel
    .fromPath('/some/path/*.txt')
    .splitText()
    .view()
```

The above example loads a collection of text files, splits the content of each file into individual lines, and prints each line.

The `by` option can be used to emit chunks of *N* lines:

```groovy
Channel
    .fromPath('/some/path/*.txt')
    .splitText( by: 10 )
    .subscribe {
        print it;
        print "--- end of the chunk ---\n"
    }
```

An optional {ref}`closure <script-closure>` can be used to transform each text chunk produced by the operator. The following example shows how to split text files into chunks of 10 lines and transform them to uppercase letters:

```groovy
Channel
    .fromPath('/some/path/*.txt')
    .splitText( by: 10 ) { it.toUpperCase() }
    .view()
```

:::{note}
Text chunks returned by the `splitText` operator are always terminated by a `\n` newline character.
:::

Available options:

`by`
: Defines the number of lines in each `chunk` (default: `1`).

`charset`
: Parse the content with the specified charset, e.g. `UTF-8`. See the list of [standard charsets](https://docs.oracle.com/en/java/javase/11/docs/api/java.base/java/nio/charset/StandardCharsets.html) for available options.

`compress`
: When `true`, resulting file chunks are GZIP compressed (default: `false`). The `.gz` suffix is automatically added to chunk file names.

`decompress`
: When `true`, decompresses the content using the GZIP format before processing it (default: `false`). Files with the `.gz` extension are decompressed automatically.

`elem`
: The index of the element to split when the source items are lists or tuples (default: first file object or first element).

`file`
: When `true`, saves each split to a file. Use a string instead of `true` value to create split files with a specific name (split index number is automatically added). Finally, set this attribute to an existing directory, in order to save the split files into the specified directory.

`keepHeader`
: Parses the first line as header and prepends it to each emitted chunk (default: `false`).

`limit`
: Limits the number of lines to retrieve for each source item (default: no limit).

See also: [countLines](#countlines)

(operator-subscribe)=

## subscribe

*Returns: nothing*

The `subscribe` operator invokes a custom function for each item from a source channel:

```{literalinclude} snippets/subscribe.nf
:language: groovy
```

```{literalinclude} snippets/subscribe.out
:language: console
```

The closure parameter can be defined explicitly if needed, using a name other than `it` and, optionally, the expected type:

```{literalinclude} snippets/subscribe-with-param.nf
:language: groovy
```

```{literalinclude} snippets/subscribe-with-param.out
:language: console
```

The `subscribe` operator supports multiple types of event handlers:

```{literalinclude} snippets/subscribe-with-on-complete.nf
:language: groovy
```

```{literalinclude} snippets/subscribe-with-on-complete.out
:language: console
```

Available options:

`onNext`
: Closure that is invoked when an item is emitted. Equivalent to providing a closure as the first argument.

`onComplete`
: Closure that is invoked after the last item is emitted by the channel.

`onError`
: Closure that is invoked when an exception is raised while handling the `onNext` event. It will not make further calls to `onNext` or `onComplete`. The `onError` method takes as its parameter the `Throwable` that caused the error.

(operator-sum)=

## sum

*Returns: value channel*

The `sum` operator emits the sum of all items in a source channel:

```{literalinclude} snippets/sum.nf
:language: groovy
```

```{literalinclude} snippets/sum.out
:language: console
```

An optional {ref}`closure <script-closure>` can be used to transform each item before it is added to the sum:

```{literalinclude} snippets/sum-with-mapper.nf
:language: groovy
```

```{literalinclude} snippets/sum-with-mapper.out
:language: console
```

## take

*Returns: queue channel*

The `take` operator takes the first *N* items from a source channel:

```{literalinclude} snippets/take.nf
:language: groovy
```

```{literalinclude} snippets/take.out
:language: console
```

:::{tip}
Specifying a size of `-1` causes the operator to take all values.
:::

See also: [until](#until)

## tap

*Returns: queue channel*

The `tap` operator assigns a source channel to a variable, and emits the source channel. It is a useful way to extract intermediate output channels from a chain of operators. For example:

```{literalinclude} snippets/tap.nf
:language: groovy
```

```{literalinclude} snippets/tap.out
:language: console
```

See also: [set](#set)

## toInteger

*Returns: queue channel*

The `toInteger` operator converts string values from a source channel to integer values:

```{literalinclude} snippets/tointeger.nf
:language: groovy
```

```{literalinclude} snippets/tointeger.out
:language: console
```

:::{note}
`toInteger` is equivalent to:

```groovy
map { it -> it as Integer }
```
:::

:::{note}
You can also use `toLong`, `toFloat`, and `toDouble` to convert to other numerical types.
:::

## toList

*Returns: value channel*

The `toList` operator collects all the items from a source channel into a list and emits the list as a single item:

```{literalinclude} snippets/tolist.nf
:language: groovy
```

```{literalinclude} snippets/tolist.out
:language: console
```

:::{note}
There are two main differences between `toList` and `collect`:

- When there is no input, `toList` emits an empty list whereas `collect` emits nothing.
- By default, `collect` flattens list items by one level.

In other words, `toList` is equivalent to:

```groovy
collect(flat: false).ifEmpty([])
```
:::

See also: [collect](#collect)

## toSortedList

*Returns: value channel*

The `toSortedList` operator collects all the items from a source channel into a sorted list and emits the list as a single item:

```{literalinclude} snippets/tosortedlist.nf
:language: groovy
```

```{literalinclude} snippets/tosortedlist.out
:language: console
```

An optional closure can be used to control how items are compared when sorting. For example, to sort tuples by their second element in descending order:

```{literalinclude} snippets/tosortedlist-with-comparator.nf
:language: groovy
```

```{literalinclude} snippets/tosortedlist-with-comparator.out
:language: console
```

:::{note}
`toSortedList` is equivalent to:

```groovy
collect(flat: false, sort: true).ifEmpty([])
```
:::

See also: [collect](#collect)

## transpose

*Returns: queue channel*

The `transpose` operator "transposes" each tuple from a source channel by flattening any nested list in each tuple, emitting each nested item separately.

To be more precise, the operator transforms a sequence of tuples like *(K, list(V), list(W), ..)* into a sequence of tuples like *(K, V, W, ..)*.

For example:

```{literalinclude} snippets/transpose-1.nf
:language: groovy
```

```{literalinclude} snippets/transpose-1.out
:language: console
```

If each source item has more than two elements, these will be flattened by the first element in the item, and a new item will be emitted only when it is complete:

```{literalinclude} snippets/transpose-2.nf
:language: groovy
```

```{literalinclude} snippets/transpose-2.out
:language: console
```

The `remainder` option can be used to emit any incomplete items:

```{literalinclude} snippets/transpose-2-with-remainder.nf
:language: groovy
```

```{literalinclude} snippets/transpose-2-with-remainder.out
:language: console
```

Available options:

`by`
: The zero-based index of the element to be transposed. Can also be a list of indices, e.g. `by: [0,2]`. By default, every list element is transposed.

`remainder`
: When `true`, incomplete tuples are emitted with `null` values for missing elements, otherwise they are discarded (default: `false`). 

See also: [groupTuple](#grouptuple)

## unique

*Returns: queue channel*

The `unique` operator emits the unique items from a source channel:

```{literalinclude} snippets/unique.nf
:language: groovy
```

```{literalinclude} snippets/unique.out
:language: console
```

An optional {ref}`closure <script-closure>` can be used to transform each item before it is evaluated for uniqueness:

```{literalinclude} snippets/unique-with-mapper.nf
:language: groovy
```

```{literalinclude} snippets/unique-with-mapper.out
:language: console
```

:::{note}
The difference between `unique` and `distinct` is that `unique` removes *all* duplicate values, whereas `distinct` removes only *consecutive* duplicate values. As a result, `unique` must process the entire source channel before it can emit anything, whereas `distinct` can emit each value immediately.
:::

See also: [distinct](#distinct)

## until

*Returns: queue channel*

The `until` operator emits each item from a source channel until a stopping condition is satisfied:

```{literalinclude} snippets/until.nf
:language: groovy
```

```{literalinclude} snippets/until.out
:language: console
```

See also: [take](#take)

(operator-view)=

## view

*Returns: queue channel*

The `view` operator prints each item from a source channel to standard output:

```{literalinclude} snippets/view.nf
:language: groovy
```

```{literalinclude} snippets/view.out
:language: console
```

An optional closure can be used to transform each item before it is printed:

```{literalinclude} snippets/view-with-mapper.nf
:language: groovy
```

```{literalinclude} snippets/view-with-mapper.out
:language: console
```

The `view` operator also emits every item that it receives, allowing it to be chained with other operators.

Available options:

`newLine`
: Print each item to a separate line (default: `true`).

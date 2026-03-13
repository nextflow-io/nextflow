(operator-typed-page)=

# Operators (typed)

This page describes the channel operators that can be used in {ref}`typed workflows <workflow-typed-page>`.

## collect

**`Channel<E> collect() -> Value<Bag<E>>`**

The `collect` operator collects all values from a source channel into a collection and emits it as a dataflow value:

```{literalinclude} ../snippets/collect.nf
:language: nextflow
```

```{literalinclude} ../snippets/collect.out
:language: console
```

## cross

**`Channel<L> cross( other: Channel<R> ) -> Channel<Tuple<L, R>>`**

**`Channel<L> cross( other: Value<R> ) -> Channel<Tuple<L, R>>`**

**`Channel<L> cross( other: R ) -> Channel<Tuple<L, R>>`**

The `cross` operator emits every pairwise combination of a source channel with another channel or value:

```{literalinclude} ../snippets/cross-v2.nf
:language: nextflow
```

```{literalinclude} ../snippets/cross-v2.out
:language: console
```

## filter

**`Channel<E> filter( condition: (E) -> Boolean ) -> Channel<E>`**

The `filter` operator emits the values from a source channel that satisfy a condition, discarding all other values:

```{literalinclude} ../snippets/filter-closure.nf
:language: nextflow
```

```{literalinclude} ../snippets/filter-closure.out
:language: console
```

## flatMap

**`Channel<E> flatMap( transform: (E) -> Iterable<R> ) -> Channel<R>`**

The `flatMap` operator applies a *mapping function* to each value from a source channel. The mapping function should return a collection, and each element in the collection is emmited separately.

For example:

```{literalinclude} ../snippets/flatmap-list.nf
:language: nextflow
```

```{literalinclude} ../snippets/flatmap-list.out
:language: console
```

## groupBy

**`Channel<Tuple<K, V>> groupBy() -> Channel<Tuple<K, Bag<V>>>`**

**`Channel<Tuple<K, Integer, V>> groupBy() -> Channel<Tuple<K, Bag<V>>>`**

The `groupBy` operator collects values from a source channel into groups based on a grouping key. A tuple is emitted for each group, containing the grouping key and collection of values.

The source channel should supply either 2-tuples of the form `(<key>, <value>)` or 3-tuples of the form `(<key>, <size>, <value>)`.

If the source tuples do not specify a size, `groupBy` will not emit any groups until *all* inputs have been received:

```{literalinclude} ../snippets/groupby-v2.nf
:language: nextflow
```

```{literalinclude} ../snippets/groupby-v2.out
:language: console
```

If the source tuples do specify a size, then `groupBy` will be able to emit each group as soon as it is ready:

```{literalinclude} ../snippets/groupby-size-v2.nf
:language: nextflow
```

```{literalinclude} ../snippets/groupby-size-v2.out
:language: console
```

:::{note}
When specifying the group size, make sure that the number of inputs for a given group matches the specified size for that group. Otherwise, the run will fail.
:::

## join

**`Channel<Record> join( other: Channel<Record>, [opts] ) -> Channel<Record>`**

The `join` operator emits the relational join of two channels of records, using a matching key given by the `by` option:

```{literalinclude} ../snippets/join-v2.nf
:language: nextflow
```

```{literalinclude} ../snippets/join-v2.out
:language: console
```

Duplicate matching keys are handled by emitting each matching combination (just like a relational join):

```{literalinclude} ../snippets/join-duplicates-v2.nf
:language: nextflow
```

```{literalinclude} ../snippets/join-duplicates-v2.out
:language: console
```

By default, unmatched values are discarded. The `remainder` option can be used to emit them at the end:

```{literalinclude} ../snippets/join-remainder-v2.nf
:language: nextflow
```

```{literalinclude} ../snippets/join-remainder-v2.out
:language: console
```

Available options:

`by: String`
: *This option is required.*
: The record field to use as the matching key.

`remainder: Boolean`
: When `true`, unmatched values are emitted at the end, otherwise they are discarded (default: `false`). 

## map

**`Channel<E> map( transform: (E) -> R ) -> Channel<R>`**

The `map` operator applies a *mapping function* to each value from a source channel:

```{literalinclude} ../snippets/map.nf
:language: nextflow
```

```{literalinclude} ../snippets/map.out
:language: console
```

## mix

**`Channel<E> mix( other: Channel<E> ) -> Channel<E>`**

**`Channel<E> mix( other: Value<E> ) -> Channel<E>`**

The `mix` operator emits the values from a channel and another channel or dataflow value into a single output channel:

```{literalinclude} ../snippets/mix.nf
:language: nextflow
```

```{literalinclude} ../snippets/mix.out
:language: console
```

The values in the output channel may be emitted in any order, for example:

```console
z
1
a
2
b
3
```

## reduce

**`Channel<E> reduce( seed: R, accumulator: (R, E) -> R ) -> Value<R>`**

**`Channel<E> reduce( accumulator: (E, E) -> E ) -> Value<E>`**

The `reduce` operator applies an *accumulator function* sequentially to each value in a source channel, and emits the final accumulated value.

The accumulator function takes two parameters -- the accumulated value and the *i*-th emitted value -- and it should return the accumulated result, which is passed to the next invocation with the *i+1*-th value. This process is repeated for each value in the source channel.

For example:

```{literalinclude} ../snippets/reduce.nf
:language: nextflow
```

```{literalinclude} ../snippets/reduce.out
:language: console
```

By default, the first value is used as the initial accumulated value (the *seed*). You can optionally specify a different initial value as shown below:

```{literalinclude} ../snippets/reduce-with-initial-value.nf
:language: nextflow
```

```{literalinclude} ../snippets/reduce-with-initial-value.out
:language: console
```

## subscribe

**`Channel<E> subscribe( action: (E) -> () )`**

**`Channel<E> subscribe( [opts] )`**

The `subscribe` operator invokes a custom function for each value from a source channel:

```{literalinclude} ../snippets/subscribe.nf
:language: nextflow
```

```{literalinclude} ../snippets/subscribe.out
:language: console
```

The `subscribe` operator supports multiple types of event handlers:

```{literalinclude} ../snippets/subscribe-with-on-complete.nf
:language: nextflow
```

```{literalinclude} ../snippets/subscribe-with-on-complete.out
:language: console
```

:::{note}
Unlike most operators, `subscribe` does not return anything. It should only be used for *side effects*, such as printing to the console, writing to a file, or making HTTP requests.
:::

Available options:

`onNext: (E) -> ()`
: Closure that is invoked when an value is emitted. Equivalent to providing a single closure argument.

`onComplete: () -> ()`
: Closure that is invoked after the last value is emitted by the channel.

`onError: (T) -> ()`
: Closure that is invoked when an exception is raised while handling the `onNext` event. It will not make further calls to `onNext` or `onComplete`. The `onError` method takes as its parameter the `Throwable` that caused the error.

## unique

**`Channel<E> unique( transform: (E) -> ? ) -> Channel<E>`**
**`Channel<E> unique() -> Channel<E>`**

The `unique` operator emits the unique values from a source channel:

```{literalinclude} ../snippets/unique.nf
:language: nextflow
```

```{literalinclude} ../snippets/unique.out
:language: console
```

An optional closure can be used to transform each value before it is evaluated for uniqueness:

```{literalinclude} ../snippets/unique-with-mapper.nf
:language: nextflow
```

```{literalinclude} ../snippets/unique-with-mapper.out
:language: console
```

## until

**`Channel<E> until( condition: (E) -> Boolean ) -> Channel<E>`**

The `until` operator emits each value from a source channel until a stopping condition is satisfied:

```{literalinclude} ../snippets/until.nf
:language: nextflow
```

```{literalinclude} ../snippets/until.out
:language: console
```

## view

**`Channel<E> view( transform: (E) -> String, [opts] ) -> Channel<E>`**

**`Channel<E> view( [opts] ) -> Channel<E>`**

The `view` operator prints each value from a source channel to standard output:

```{literalinclude} ../snippets/view.nf
:language: nextflow
```

```{literalinclude} ../snippets/view.out
:language: console
```

An optional closure can be used to transform each value before it is printed:

```{literalinclude} ../snippets/view-with-mapper.nf
:language: nextflow
```

```{literalinclude} ../snippets/view-with-mapper.out
:language: console
```

The `tag` option can be used to print the channel only when the `-dump-channels` command-line option is specified with the given tag:

```{literalinclude} ../snippets/view-tag-v2.nf
:language: nextflow
```

You can run this script with `-dump-channels plus1` or `-dump-channels exp2` to print either channel, or `-dump-channels plus1,exp2` to print both.

The `view` operator also emits every value that it receives, allowing it to be chained with other operators.

Available options:

`newLine: Boolean`
: Print each value to a separate line (default: `true`).

`tag: String`
: Print the channel values only when `-dump-channels` is specified on the command line with the given tag.

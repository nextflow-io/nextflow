(dataflow-page)=

# Dataflow

Nextflow uses a **dataflow** programming model to define workflows declaratively. In this model, {ref}`processes <process-page>` in a pipeline are connected to each other through *dataflow channels* and *dataflow values*.

(dataflow-type-channel)=

## Channels

A *dataflow channel* (or simply *channel*) is an asynchronous sequence of values.

The values in a channel cannot be accessed directly, but only through an operator or process. For example:

```nextflow
channel.of(1, 2, 3).view { v -> "channel emits ${v}" }
```

```console
channel emits 1
channel emits 2
channel emits 3
```

### Factories

A channel can be created by factories in the `channel` namespace. For example, the `channel.fromPath()` factory creates a channel from a file name or glob pattern, similar to the `files()` function:

```nextflow
channel.fromPath('input/*.txt').view()
```

See {ref}`channel-factory` for the full list of channel factories.

### Operators

Channel operators, or *operators* for short, are functions that consume and produce channels. Because channels are asynchronous, operators are necessary to manipulate the values in a channel. Operators are particularly useful for implementing glue logic between processes.

Commonly used operators include:

- {ref}`operator-combine`: emit the combinations of two channels

- {ref}`operator-collect`: collect the values from a channel into a list

- {ref}`operator-filter`: select the values in a channel that satisfy a condition

- {ref}`operator-flatMap`: transform each value from a channel into a list and emit each list element separately

- {ref}`operator-grouptuple`: group the values from a channel based on a grouping key

- {ref}`operator-join`: join the values from two channels based on a matching key

- {ref}`operator-map`: transform each value from a channel with a mapping function

- {ref}`operator-mix`: emit the values from multiple channels

- {ref}`operator-view`: print each value in a channel to standard output

See {ref}`operator-page` for the full list of operators.

(dataflow-type-value)=

## Values

A *dataflow value* is an asynchronous value.

Dataflow values can be created using the {ref}`channel.value <channel-value>` factory, and they are created by processes (under {ref}`certain conditions <process-out-singleton>`).

A dataflow value cannot be accessed directly, but only through an operator or process. For example:

```nextflow
channel.value(1).view { v -> "dataflow value is ${v}" }
```

```console
dataflow value is 1
```

See {ref}`stdlib-types-value` for the set of available methods for dataflow values.

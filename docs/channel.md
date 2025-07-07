(channel-page)=

# Channels

In Nextflow, **channels** are the key data structures that facilitate the dataflow dependencies between each step (i.e. {ref}`process <process-page>`) in a pipeline.

There are two kinds of channels, *queue channels* and *value channels*. Channels are created using *channel factories* and transformed using *channel operators*.

(channel-type-queue)=

## Queue channels

A *queue channel* is a channel that *emits* an asynchronous sequence of values.

A queue channel can be created by channel factories (e.g., {ref}`channel.of <channel-of>` and {ref}`channel.fromPath <channel-path>`), operators (e.g., {ref}`operator-map` and  {ref}`operator-filter`), and processes (see {ref}`Process outputs <process-output>`).

The values in a queue channel cannot be accessed directly -- they can only be accessed by passing the channel as input to an operator or process. For example:

```nextflow
channel.of(1, 2, 3).view { v -> "queue channel emits ${v}" }
```

```console
queue channel emits 1
queue channel emits 2
queue channel emits 3
```

(channel-type-value)=

## Value channels

A *value channel* is a channel that is *bound* to an asynchronous value.

A value channel can be created with the {ref}`channel.value <channel-value>` factory, certain operators (e.g., {ref}`operator-collect` and {ref}`operator-reduce`), and processes (under {ref}`certain conditions <process-out-singleton>`).

The value in a value channel cannot be accessed directly -- it can only be accessed by passing the channel as input to an operator or process. For example:

```nextflow
channel.value(1).view { v -> "value channel is ${v}" }
```

```console
value channel is 1
```

## Channel factories

Channel factories are functions that create channels from regular values.

The `channel.fromPath()` factory creates a channel from a file name or glob pattern, similar to the `files()` function:

```nextflow
channel.fromPath('input/*.txt').view()
```

See {ref}`channel-factory` for the full list of channel factories.

## Operators

Channel operators, or *operators* for short, are functions that consume and produce channels. Because channels are asynchronous, operators are necessary to manipulate the values in a channel. Operators are particularly useful for implementing glue logic between processes.

Commonly used operators include:

- {ref}`operator-combine`: emit the combinations of two channels

- {ref}`operator-collect`: collect the values from a channel into a list

- {ref}`operator-filter`: select the values in a channel that satisfy a condition

- {ref}`operator-flatMap`: transform each value from a channel into a list and emit each list 
element separately

- {ref}`operator-grouptuple`: group the values from a channel based on a grouping key

- {ref}`operator-join`: join the values from two channels based on a matching key

- {ref}`operator-map`: transform each value from a channel with a mapping function

- {ref}`operator-mix`: emit the values from multiple channels

- {ref}`operator-view`: print each value in a channel to standard output

See {ref}`operator-page` for the full list of operators.

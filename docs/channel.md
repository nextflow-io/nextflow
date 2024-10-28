(channel-page)=

# Channels

Nextflow is based on the dataflow programming model in which processes communicate through channels.

A channel has two major properties:

1. Sending a message is an *asynchronous* (i.e. non-blocking) operation, which means the sender doesn't have to wait for the receiving process.
2. Receiving a message is a *synchronous* (i.e. blocking) operation, which means the receiving process must wait until a message has arrived.

(channel-types)=

## Channel types

In Nextflow there are two kinds of channels: *queue channels* and *value channels*.

(channel-type-queue)=

### Queue channel

A *queue channel* is a non-blocking unidirectional FIFO queue connecting a *producer* process (i.e. outputting a value)
to a consumer process, or an operators.

A queue channel can be created by factory methods ({ref}`channel-of`, {ref}`channel-path`, etc), operators ({ref}`operator-map`, {ref}`operator-flatmap`, etc), and processes (see {ref}`Process outputs <process-output>`).

(channel-type-value)=

### Value channel

A *value channel* can be bound (i.e. assigned) with one and only one value, and can be consumed any number of times by
a process or an operator.

A value channel can be created with the {ref}`channel-value` factory method or by any operator that produces a single value
({ref}`operator-first`, {ref}`operator-collect`, {ref}`operator-reduce`, etc). Additionally, a process will emit value
channels if it is invoked with all value channels, including simple values which are implicitly wrapped in a value channel.

For example:

```nextflow
process foo {
  input:
  val x

  output:
  path 'x.txt'

  """
  echo $x > x.txt
  """
}

workflow {
  result = foo(1)
  result.view { file -> "Result: ${file}" }
}
```

In the above example, since the `foo` process is invoked with a simple value instead of a channel, the input is implicitly
wrapped in a value channel, and the output is also emitted as a value channel.

See also: {ref}`process-multiple-input-channels`.

## Channel factories

Channel factories are functions that can create channels.

For example, the `Channel.of()` factory can be used to create a channel from an arbitrary list of arguments:

```nextflow
Channel.of(1, 2, 3).view()
```

:::{versionadded} 20.07.0
`channel` was introduced as an alias of `Channel`, allowing factory methods to be specified as `channel.of()` or `Channel.of()`, and so on.
:::

See {ref}`channel-factory` for the full list of channel factories.

## Operators

Operators are methods that consume and produce channels. Because channels are asynchronous, operators are necessary to manipulate the values in a channel, without using a process. As a result, operators are useful for implementing the "glue logic" between processes.

See {ref}`operator-page` for the full list of operators. If you are new to Nextflow, here are some commonly-used operators to learn first:

Filtering:

- {ref}`operator-filter`: select all values in a channel that satisfy a condition
- {ref}`operator-first`: select the first value in a channel
- {ref}`operator-take`: select the first *n* values in a channel
- {ref}`operator-unique`: select the unique values in a channel (i.e. remove duplicates)

Transforming:

- {ref}`operator-collect`: collect the values from a channel into a list
- {ref}`operator-grouptuple`: group the values from a channel based on a grouping key
- {ref}`operator-map`: transform each value from a channel with a mapping function
- {ref}`operator-reduce`: accumulate each value from a channel into a single value

Combining multiple channels:

- {ref}`operator-combine`: emit the combinations of two channels
- {ref}`operator-concat`: emit the values from multiple channels (in the order in which the channels were given)
- {ref}`operator-join`: join the values from two channels based on a matching key
- {ref}`operator-mix`: emit the values from multiple channels (in the order in which items arrive)

Miscellaneous:

- {ref}`operator-ifempty`: emit a channel, or a default value if the channel is empty
- {ref}`operator-set`: assign a channel to a variable
- {ref}`operator-view`: print each value in a channel to standard output

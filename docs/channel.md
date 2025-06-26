# Channels

Nextflow is based on the dataflow programming model in which processes communicate through channels.

A channel has two major properties:

1. Sending a message is an *asynchronous* (i.e. non-blocking) operation, which means the sender doesn't have to wait for the receiving process.
2. Receiving a message is a *synchronous* (i.e. blocking) operation, which means the receiving process must wait until a message has arrived.

## Channel types

In Nextflow there are two kinds of channels: *queue channels* and *value channels*.

### Queue channel

A *queue channel* is a non-blocking unidirectional FIFO queue connecting a *producer* process (i.e. outputting a value)
to a consumer process or an operator.

A queue channel can be created by factory methods ([of][channel-of], [fromPath][channel-path], etc), operators ([map][operator-map], [flatMap][operator-flatmap], etc), and processes (see [Outputs][process-output]).

### Value channel

A *value channel* can be bound (i.e. assigned) with one and only one value, and can be consumed any number of times by
a process or an operator.

A value channel can be created with the [value][channel-value] factory method or by any operator that produces a single value
([first][operator-first], [collect][operator-collect], [reduce][operator-reduce], etc). Additionally, a process will emit value
channels if it is invoked with all value channels, including simple values which are implicitly wrapped in a value channel.

For example:

```nextflow
process echo {
  input:
  val x

  output:
  path 'x.txt'

  script:
  """
  echo $x > x.txt
  """
}

workflow {
  result = echo(1)
  result.view { file -> "Result: ${file}" }
}
```

In the above example, since the `echo` process is invoked with a simple value instead of a channel, the input is implicitly
wrapped in a value channel, and the output is also emitted as a value channel.

See also: [Multiple input channels][process-multiple-input-channels].

## Channel factories

Channel factories are functions that can create channels.

For example, the `channel.of()` factory can be used to create a channel from an arbitrary list of arguments:

```nextflow
channel.of(1, 2, 3).view()
```

See [Channel factories][channel-factory] for the full list of channel factories.

## Operators

Channel operators, or _operators_ for short, are functions that consume and produce channels. Because channels are asynchronous, operators are necessary to manipulate the values in a channel, aside from using a process. As a result, operators are useful for implementing the _glue logic_ between processes.

Commonly used operators include:

- [combine][operator-combine]: emit the combinations of two channels
- [collect][operator-collect]: collect the values from a channel into a list
- [filter][operator-filter]: select the values in a channel that satisfy a condition
- [flatMap][operator-flatMap]: transform each value from a channel into a list and emit each list element separately
- [groupTuple][operator-grouptuple]: group the values from a channel based on a grouping key
- [join][operator-join]: join the values from two channels based on a matching key
- [map][operator-map]: transform each value from a channel with a mapping function
- [mix][operator-mix]: emit the values from multiple channels
- [view][operator-view]: print each value in a channel to standard output

See [Operators][operator-page] for the full list of operators.

[channel-of]: /nextflow_docs/nextflow_repo/docs/reference/channel#of
[channel-path]: /nextflow_docs/nextflow_repo/docs/reference/channel#frompath
[operator-map]: /nextflow_docs/nextflow_repo/docs/reference/operator#map
[operator-flatmap]: /nextflow_docs/nextflow_repo/docs/reference/operator#flatmap
[process-output]: /nextflow_docs/nextflow_repo/docs/process#outputs
[channel-value]: /nextflow_docs/nextflow_repo/docs/reference/channel#value
[operator-first]: /nextflow_docs/nextflow_repo/docs/reference/operator#first
[operator-collect]: /nextflow_docs/nextflow_repo/docs/reference/operator#collect
[operator-reduce]: /nextflow_docs/nextflow_repo/docs/reference/operator#reduce
[process-multiple-input-channels]: /nextflow_docs/nextflow_repo/docs/process#multiple-input-channels
[channel-factory]: /nextflow_docs/nextflow_repo/docs/channel
[operator-combine]: /nextflow_docs/nextflow_repo/docs/reference/operator#combine
[operator-filter]: /nextflow_docs/nextflow_repo/docs/reference/operator#filter
[operator-flatMap]: /nextflow_docs/nextflow_repo/docs/reference/operator#flatmap
[operator-grouptuple]: /nextflow_docs/nextflow_repo/docs/reference/operator#grouptuple
[operator-join]: /nextflow_docs/nextflow_repo/docs/reference/operator#join
[operator-mix]: /nextflow_docs/nextflow_repo/docs/reference/operator#mix
[operator-view]: /nextflow_docs/nextflow_repo/docs/reference/operator#view
[operator-page]: /nextflow_docs/nextflow_repo/docs/reference/operator#page

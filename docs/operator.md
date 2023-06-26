(operator-page)=

# Operators

Nextflow **operators** are methods that allow you to manipulate channels. Every operator, with the exception of [set](#set) and [subscribe](#subscribe), produces one or more new channels, allowing you to chain operators to fit your needs.

This page is a comprehensive reference for all Nextflow operators. However, if you are new to Nextflow, here are some suggested operators to learn for common use cases:

- Filtering: [filter](#filter), [randomSample](#randomsample), [take](#take), [unique](#unique)
- Reduction: [collect](#collect), [groupTuple](#grouptuple), [reduce](#reduce)
- Parsing text data: [splitCsv](#splitcsv), [splitJson](#splitjson), [splitText](#splittext)
- Combining channels: [combine](#combine), [concat](#concat), [join](#join), [mix](#mix)
- Forking channels: [branch](#branch), [multiMap](#multimap)
- Maths: [count](#count), [max](#max), [min](#min), [sum](#sum)
- Other: [ifEmpty](#ifempty), [map](#map), [set](#set), [view](#view)

(operator-branch)=

## branch

:::{versionadded} 19.08.0-edge
:::

*Returns: map of queue channels*

The `branch` operator allows you to forward the items emitted by a source channel to one or more output channels, choosing one out of them at a time.

The selection criteria is defined by specifying a {ref}`closure <script-closure>` that provides one or more boolean expression, each of which is identified by a unique label. On the first expression that evaluates to a *true* value, the current item is bound to a named channel as the label identifier. For example:

```groovy
Channel
    .of(1, 2, 3, 40, 50)
    .branch {
        small: it < 10
        large: it > 10
    }
    .set { result }

 result.small.view { "$it is small" }
 result.large.view { "$it is large" }
```

```
1 is small
2 is small
3 is small
40 is large
50 is large
```

:::{note}
The above *small* and *large* strings may be printed in any order due to the asynchronous execution of the `view` operator.
:::

A default fallback condition can be specified using `true` as the last branch condition:

```groovy
Channel
    .from(1, 2, 3, 40, 50)
    .branch {
        small: it < 10
        large: it < 50
        other: true
    }
```

The value returned by each branch condition can be customised by specifying an optional expression statement(s) just after the condition expression. For example:

```groovy
Channel
    .from(1, 2, 3, 40, 50)
    .branch {
        foo: it < 10
            return it+2

        bar: it < 50
            return it-2

        other: true
            return 0
    }
```

:::{tip}
When the `return` keyword is omitted, the value of the last expression statement is implicitly returned.
:::

To create a branch criteria as variable that can be passed as an argument to more than one `branch` operator use the `branchCriteria` built-in method as shown below:

```groovy
def criteria = branchCriteria {
    small: it < 10
    large: it > 10
}

Channel.of(1, 2, 30).branch(criteria).set { ch1 }
Channel.of(10, 20, 1).branch(criteria).set { ch2 }
```

## buffer

*Returns: queue channel*

The `buffer` operator gathers the items emitted by the source channel into subsets and emits these subsets separately.

There are a number of ways you can regulate how `buffer` gathers the items from the source channel into subsets:

- `buffer( closingCondition )`: starts to collect the items emitted by the channel into a subset until the `closingCondition` is verified. After that the subset is emitted to the resulting channel and new items are gathered into a new subset. The process is repeated until the last value in the source channel is sent. The `closingCondition` can be specified either as a {ref}`regular expression <script-regexp>`, a Java class, a literal value, or a boolean predicate that has to be satisfied. For example:

  ```groovy
  Channel
      .of( 1, 2, 3, 1, 2, 3 )
      .buffer { it == 2 }
      .view()

  // emitted values
  [1,2]
  [3,1,2]
  ```

- `buffer( openingCondition, closingCondition )`: starts to gather the items emitted by the channel as soon as one of the them verify the `openingCondition` and it continues until there is one item which verify the `closingCondition`. After that the subset is emitted and it continues applying the described logic until the last channel item is emitted. Both conditions can be defined either as a {ref}`regular expression <script-regexp>`, a literal value, a Java class, or a boolean predicate that need to be satisfied. For example:

  ```groovy
  Channel
      .of( 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2 )
      .buffer( 2, 4 )
      .view()

  // emits bundles starting with '2' and ending with'4'
  [2,3,4]
  [2,3,4]
  ```

- `buffer( size: n )`: transform the source channel in such a way that it emits tuples made up of `n` elements. An incomplete tuple is discarded. For example:

  ```groovy
  Channel
      .of( 1, 2, 3, 1, 2, 3, 1 )
      .buffer( size: 2 )
      .view()

  // emitted values
  [1, 2]
  [3, 1]
  [2, 3]
  ```

  If you want to emit the last items in a tuple containing less than `n` elements, simply add the parameter `remainder` specifying `true`, for example:

  ```groovy
  Channel
      .of( 1, 2, 3, 1, 2, 3, 1 )
      .buffer( size: 2, remainder: true )
      .view()

  // emitted values
  [1, 2]
  [3, 1]
  [2, 3]
  [1]
  ```

- `buffer( size: n, skip: m )`: as in the previous example, it emits tuples containing `n` elements, but skips `m` values before starting to collect the values for the next tuple (including the first emission). For example:

  ```groovy
  Channel
      .of( 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2 )
      .buffer( size:3, skip:2 )
      .view()

  // emitted values
  [3, 4, 5]
  [3, 4, 5]
  ```

  If you want to emit the remaining items in a tuple containing less than `n` elements, simply add the parameter `remainder` specifying `true`, as shown in the previous example.

See also: [collate](#collate) operator.

## collate

*Returns: queue channel*

The `collate` operator transforms a channel in such a way that the emitted values are grouped in tuples containing `n` items. For example:

```groovy
Channel
    .of(1,2,3,1,2,3,1)
    .collate( 3 )
    .view()
```

```
[1, 2, 3]
[1, 2, 3]
[1]
```

As shown in the above example the last tuple may be incomplete e.g. contain fewer elements than the specified size. If you want to avoid this, specify `false` as the second parameter. For example:

```groovy
Channel
    .of(1,2,3,1,2,3,1)
    .collate( 3, false )
    .view()
```

```
[1, 2, 3]
[1, 2, 3]
```

A second version of the `collate` operator allows you to specify, after the `size`, the `step` by which elements are collected in tuples. For example:

```groovy
Channel
    .of(1,2,3,4)
    .collate( 3, 1 )
    .view()
```

```
[1, 2, 3]
[2, 3, 4]
[3, 4]
[4]
```

As before, if you don't want to emit the last items which do not complete a tuple, specify `false` as the third parameter.

See also: [buffer](#buffer) operator.

(operator-collect)=

## collect

*Returns: value channel*

The `collect` operator collects all the items emitted by a channel to a `List` and return the resulting object as a sole emission. For example:

```groovy
Channel
    .of( 1, 2, 3, 4 )
    .collect()
    .view()
```

```
[1,2,3,4]
```

An optional {ref}`closure <script-closure>` can be specified to transform each item before adding it to the resulting list. For example:

```groovy
Channel
    .of( 'hello', 'ciao', 'bonjour' )
    .collect { it.length() }
    .view()
```

```
[5,4,7]
```

Available options:

`flat`
: When `true` nested list structures are normalised and their items are added to the resulting list object (default: `true`).

`sort`
: When `true` the items in the resulting list are sorted by their natural ordering. It is possible to provide a custom ordering criteria by using either a {ref}`closure <script-closure>` or a [Comparator](https://docs.oracle.com/javase/8/docs/api/java/util/Comparator.html) object (default: `false`).

See also: [toList](#tolist) and [toSortedList](#tosortedlist) operator.

## collectFile

*Returns: queue channel*

The `collectFile` operator allows you to gather the items emitted by a channel and save them to one or more files. The operator returns a new channel that emits the collected file(s).

In the simplest case, just specify the name of a file where the entries have to be stored. For example:

```groovy
Channel
    .of('alpha', 'beta', 'gamma')
    .collectFile(name: 'sample.txt', newLine: true)
    .subscribe {
        println "Entries are saved to file: $it"
        println "File content is: ${it.text}"
    }
```

A second version of the `collectFile` operator allows you to gather the items emitted by a channel and group them together into files whose name can be defined by a dynamic criteria. The grouping criteria is specified by a {ref}`closure <script-closure>` that must return a pair in which the first element defines the file name for the group and the second element the actual value to be appended to that file. For example:

```groovy
Channel
    .of('Hola', 'Ciao', 'Hello', 'Bonjour', 'Halo')
    .collectFile() { item ->
        [ "${item[0]}.txt", item + '\n' ]
    }
    .subscribe {
        println "File ${it.name} contains:"
        println it.text
    }
```

```
File 'B.txt' contains:
Bonjour

File 'C.txt' contains:
Ciao

File 'H.txt' contains:
Halo
Hola
Hello
```

:::{tip}
When the items emitted by the source channel are files, the grouping criteria can be omitted. In this case the items content will be grouped into file(s) having the same name as the source items.
:::

Available options:

`cache`
: Controls the caching ability of the `collectFile` operator when using the *resume* feature. It follows the same semantic of the {ref}`process-cache` directive (default: `true`).

`keepHeader`
: Prepend the resulting file with the header fetched in the first collected file. The header size (ie. lines) can be specified by using the `skip` parameter (default: `false`), to determine how many lines to remove from all collected files except for the first (where no lines will be removed).

`name`
: Name of the file where all received values are stored.

`newLine`
: Appends a `newline` character automatically after each entry (default: `false`).

`seed`
: A value or a map of values used to initialise the files content.

`skip`
: Skip the first `n` lines e.g. `skip: 1`.

`sort`
: Defines sorting criteria of content in resulting file(s). Can be one of the following values:

  - `false`: Disable content sorting. Entries are appended as they are produced.
  - `true`: Order the content by the entry's natural ordering i.e. numerical for number, lexicographic for string, etc. See the [Java documentation](http://docs.oracle.com/javase/tutorial/collections/interfaces/order.html) for more information.
  - `'index'`: Order the content by the incremental index number assigned to each entry while they are collected.
  - `'hash'`: (default) Order the content by the hash number associated to each entry
  - `'deep'`: Similar to the previous, but the hash number is created on actual entries content e.g. when the entry is a file the hash is created on the actual file content.
  - A custom sorting criteria can be specified by using either a {ref}`Closure <script-closure>` or a [Comparator](http://docs.oracle.com/javase/7/docs/api/java/util/Comparator.html) object.

  The file content is sorted in such a way that it does not depend on the order in which entries were added to it, which guarantees that it is consistent (i.e. does not change) across different executions with the same data.

`storeDir`
: Folder where the resulting file(s) are be stored.

`tempDir`
: Folder where temporary files, used by the collecting process, are stored.

The following snippet shows how sort the content of the result file alphabetically:

```groovy
Channel
    .of('Z'..'A')
    .collectFile(name:'result', sort: true, newLine: true)
    .view { it.text }
```

```
A
B
C
:
Z
```

The following example shows how use a `closure` to collect and sort all sequences in a FASTA file from shortest to longest:

```groovy
Channel
    .fromPath('/data/sequences.fa')
    .splitFasta( record: [id: true, sequence: true] )
    .collectFile( name:'result.fa', sort: { it.size() } ) {
        it.sequence
    }
    .view { it.text }
```

:::{warning}
The `collectFile` operator needs to store files in a temporary folder that is automatically deleted on workflow completion. For performance reasons this folder is located in the machine's local storage, and it will require as much free space as the data that is being collected. Optionally, a different temporary data folder can be specified by using the `tempDir` parameter.
:::

(operator-combine)=

## combine

*Returns: queue channel*

The `combine` operator combines (cartesian product) the items emitted by two channels or by a channel and a `Collection` object (as right operand). For example:

```groovy
numbers = Channel.of(1, 2, 3)
words = Channel.of('hello', 'ciao')
numbers
    .combine(words)
    .view()
```

```
[1, hello]
[2, hello]
[3, hello]
[1, ciao]
[2, ciao]
[3, ciao]
```

A second version of the `combine` operator allows you to combine items that share a common matching key. The index of the key element is specified by using the `by` parameter (zero-based index, multiple indices can be specified as a list of integers). For example:

```groovy
left = Channel.of(['A', 1], ['B', 2], ['A', 3])
right = Channel.of(['B', 'x'], ['B', 'y'], ['A', 'z'], ['A', 'w'])

left
    .combine(right, by: 0)
    .view()
```

```
[A, 1, z]
[A, 3, z]
[A, 1, w]
[A, 3, w]
[B, 2, x]
[B, 2, y]
```

See also [join](#join) and [cross](#cross).

(operator-concat)=

## concat

*Returns: queue channel*

The `concat` operator allows you to *concatenate* the items emitted by two or more channels to a new channel. The items emitted by the resulting channel are in the same order as specified in the operator arguments.

In other words, given *N* channels, the items from the *i+1 th* channel are emitted only after all of the items from the *i th* channel have been emitted.

For example:

```groovy
a = Channel.of('a', 'b', 'c')
b = Channel.of(1, 2, 3)
c = Channel.of('p', 'q')

c.concat( b, a ).view()
```

```
p
q
1
2
3
a
b
c
```

(operator-count)=

## count

*Returns: value channel*

The `count` operator creates a channel that emits a single item: a number that represents the total number of items emitted by the source channel. For example:

```groovy
Channel
    .of(9,1,7,5)
    .count()
    .view()
// -> 4
```

An optional parameter can be provided to select which items are to be counted. The selection criteria can be specified either as a {ref}`regular expression <script-regexp>`, a literal value, a Java class, or a boolean predicate that needs to be satisfied. For example:

```groovy
Channel
    .of(4,1,7,1,1)
    .count(1)
    .view()
// -> 3

Channel
    .of('a','c','c','q','b')
    .count ( ~/c/ )
    .view()
// -> 2

Channel
    .of('a','c','c','q','b')
    .count { it <= 'c' }
    .view()
// -> 4
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

The `cross` operator allows you to combine the items of two channels in such a way that the items of the source channel are emitted along with the items emitted by the target channel for which they have a matching key.

The key is defined, by default, as the first entry in an array, a list or map object, or the value itself for any other data type. For example:

```groovy
source = Channel.of( [1, 'alpha'], [2, 'beta'] )
target = Channel.of( [1, 'x'], [1, 'y'], [1, 'z'], [2, 'p'], [2, 'q'], [2, 't'] )

source.cross(target).view()
```

```
[ [1, alpha], [1, x] ]
[ [1, alpha], [1, y] ]
[ [1, alpha], [1, z] ]
[ [2, beta],  [2, p] ]
[ [2, beta],  [2, q] ]
[ [2, beta],  [2, t] ]
```

The above example shows how the items emitted by the source channels are associated to the ones emitted by the target channel (on the right) having the same key.

There are two important caveats when using the `cross` operator:

1. The operator is not `commutative`, i.e. the result of `a.cross(b)` is different from `b.cross(a)`
2. The source channel should emits items for which there's no key repetition i.e. the emitted items have an unique key identifier.

Optionally, a mapping function can be specified in order to provide a custom rule to associate an item to a key.

## distinct

*Returns: queue channel*

The `distinct` operator allows you to remove *consecutive* duplicated items from a channel, so that each emitted item is different from the preceding one. For example:

```groovy
Channel
    .of( 1,1,2,2,2,3,1,1,2,2,3 )
    .distinct()
    .subscribe onNext: { println it }, onComplete: { println 'Done' }
```

```
1
2
3
1
2
3
Done
```

You can also specify an optional {ref}`closure <script-closure>` that customizes the way it distinguishes between distinct items. For example:

```groovy
Channel
    .of( 1,1,2,2,2,3,1,1,2,4,6 )
    .distinct { it % 2 }
    .subscribe onNext: { println it }, onComplete: { println 'Done' }
```

```
1
2
3
2
Done
```

(operator-dump)=

## dump

*Returns: queue channel or value channel, depending on the input*

The `dump` operator prints the items emitted by the channel to which is applied only when the option `-dump-channels` is specified on the `run` command line, otherwise it is ignored.

This is useful to enable the debugging of one or more channel content on-demand by using a command line option instead of modifying your script code.

An optional `tag` parameter allows you to select which channel to dump. For example:

```groovy
Channel
    .of(1,2,3)
    .map { it+1 }
    .dump(tag: 'foo')

Channel
    .of(1,2,3)
    .map { it^2 }
    .dump(tag: 'bar')
```

Then you will be able to specify the tag `foo` or `bar` as an argument of the `-dump-channels` option to print either the content of the first or the second channel. Multiple tag names can be specified separating them with a `,` character.

:::{versionadded} 22.10.0
The output can be formatted by enabling the `pretty` option:

```groovy
Channel
    .fromSRA('SRP043510')
    .dump(tag: 'foo', pretty: true)
```
:::

## filter

*Returns: queue channel*

The `filter` operator allows you to get only the items emitted by a channel that satisfy a condition and discarding all the others. The filtering condition can be specified by using either a {ref}`regular expression <script-regexp>`, a literal value, a type qualifier (i.e. a Java class) or any boolean predicate.

The following example shows how to filter a channel by using a regular expression that returns only strings that begin with `a`:

```groovy
Channel
    .of( 'a', 'b', 'aa', 'bc', 3, 4.5 )
    .filter( ~/^a.*/ )
    .view()
```

```
a
aa
```

The following example shows how to filter a channel by specifying the type qualifier `Number` so that only numbers are returned:

```groovy
Channel
    .of( 'a', 'b', 'aa', 'bc', 3, 4.5 )
    .filter( Number )
    .view()
```

```
3
4.5
```

Finally, a filtering condition can be defined by using any a boolean predicate. A predicate is expressed by a {ref}`closure <script-closure>` returning a boolean value. For example the following fragment shows how filter a channel emitting numbers so that the odd values are returned:

```groovy
Channel
    .of( 1, 2, 3, 4, 5 )
    .filter { it % 2 == 1 }
    .view()
```

```
1
3
5
```

:::{tip}
In the above example the filter condition is wrapped in curly brackets, instead of parentheses, because it specifies a {ref}`closure <script-closure>` as the operator's argument. In reality it is just syntactic sugar for `filter({ it % 2 == 1 })`
:::

(operator-first)=

## first

*Returns: value channel*

The `first` operator creates a channel that returns the first item emitted by the source channel, or eventually the first item that matches an optional condition. The condition can be specified by using a {ref}`regular expression<script-regexp>`, a Java `class` type or any boolean predicate. For example:

```groovy
// no condition is specified, emits the very first item: 1
Channel
    .of( 1, 2, 3 )
    .first()
    .view()

// emits the first String value: 'a'
Channel
    .of( 1, 2, 'a', 'b', 3 )
    .first( String )
    .view()

// emits the first item matching the regular expression: 'aa'
Channel
    .of( 'a', 'aa', 'aaa' )
    .first( ~/aa.*/ )
    .view()

// emits the first item for which the predicate evaluates to true: 4
Channel
    .of( 1,2,3,4,5 )
    .first { it > 3 }
    .view()
```

(operator-flatmap)=

## flatMap

*Returns: queue channel*

The `flatMap` operator applies a function of your choosing to every item emitted by a channel, and returns the items so obtained as a new channel. Whereas the mapping function returns a list of items, this list is flattened so that each single item is emitted on its own.

For example:

```groovy
// create a channel of numbers
numbers = Channel.of( 1, 2, 3 )

// map each number to a tuple (array), which items are emitted separately
results = numbers.flatMap { n -> [ n*2, n*3 ] }

// print the final results
results.subscribe onNext: { println it }, onComplete: { println 'Done' }
```

```
2
3
4
6
6
9
Done
```

Associative arrays are handled in the same way, so that each array entry is emitted as a single key-value pair. For example:

```groovy
Channel
    .of ( 1, 2, 3 )
    .flatMap { it -> [ number: it, square: it*it ] }
    .view { it.key + ': ' + it.value }
```

```
number: 1
square: 1
number: 2
square: 4
number: 3
square: 9
```

(operator-flatten)=

## flatten

*Returns: queue channel*

The `flatten` operator transforms a channel in such a way that every item of type `Collection` or `Array` is flattened so that each single entry is emitted separately by the resulting channel. For example:

```groovy
Channel
    .of( [1,[2,3]], 4, [5,[6]] )
    .flatten()
    .view()
```

```
1
2
3
4
5
6
```

See also: [flatMap](#flatmap) operator.

(operator-grouptuple)=

## groupTuple

*Returns: queue channel*

The `groupTuple` operator collects tuples (or lists) of values emitted by the source channel grouping together the elements that share the same key. Finally it emits a new tuple object for each distinct key collected.

In other words, the operator transforms a sequence of tuple like *(K, V, W, ..)* into a new channel emitting a sequence of *(K, list(V), list(W), ..)*

For example:

```groovy
Channel
    .of( [1, 'A'], [1, 'B'], [2, 'C'], [3, 'B'], [1, 'C'], [2, 'A'], [3, 'D'] )
    .groupTuple()
    .view()
```

```
[1, [A, B, C]]
[2, [C, A]]
[3, [B, D]]
```

By default the first entry in the tuple is used as grouping key. A different key can be chosen by using the `by` parameter and specifying the index of the entry to be used as key (the index is zero-based). For example, grouping by the second value in each tuple:

```groovy
Channel
    .of( [1, 'A'], [1, 'B'], [2, 'C'], [3, 'B'], [1, 'C'], [2, 'A'], [3, 'D'] )
    .groupTuple(by: 1)
    .view()
```

```
[[1, 2], A]
[[1, 3], B]
[[2, 1], C]
[[3], D]
```

By default, if you don't specify a size, the `groupTuple` operator will not emit any groups until *all* inputs have been received. If possible, you should always try to specify the number of expected elements in each group using the `size` option, so that each group can be emitted as soon as it's ready. In cases where the size of each group varies based on the grouping key, you can use the built-in `groupKey` function, which allows you to create a special grouping key with an associated size:

```groovy
chr_frequency = [ "chr1": 2, "chr2": 3 ]

Channel.of(
        [ 'region1', 'chr1', '/path/to/region1_chr1.vcf' ],
        [ 'region2', 'chr1', '/path/to/region2_chr1.vcf' ],
        [ 'region1', 'chr2', '/path/to/region1_chr2.vcf' ],
        [ 'region2', 'chr2', '/path/to/region2_chr2.vcf' ],
        [ 'region3', 'chr2', '/path/to/region3_chr2.vcf' ]
    )
    .map { region, chr, vcf -> tuple( groupKey(chr, chr_frequency[chr]), vcf ) }
    .groupTuple()
    .view()
```

```
[chr1, [/path/to/region1_chr1.vcf, /path/to/region2_chr1.vcf]]
[chr2, [/path/to/region1_chr2.vcf, /path/to/region2_chr2.vcf, /path/to/region3_chr2.vcf]]
```

Available options:

`by`
: The index (zero based) of the element to be used as grouping key. A key composed by multiple elements can be defined specifying a list of indices e.g. `by: [0,2]`

`remainder`
: When `false` incomplete tuples (i.e. with less than `size` grouped items) are discarded (default). When `true` incomplete tuples are emitted as the ending emission. Only valid when a `size` parameter is specified.

`size`
: The number of items the grouped list(s) has to contain. When the specified size is reached, the tuple is emitted.

`sort`
: Defines the sorting criteria for the grouped items. Can be one of the following values:

  - `false`: No sorting is applied (default).
  - `true`: Order the grouped items by the item's natural ordering i.e. numerical for number, lexicographic for string, etc. See the [Java documentation](http://docs.oracle.com/javase/tutorial/collections/interfaces/order.html) for more information.
  - `hash`: Order the grouped items by the hash number associated to each entry.
  - `deep`: Similar to the previous, but the hash number is created on actual entries content e.g. when the item is a file, the hash is created on the actual file content.
  - A custom sorting criteria used to order the tuples element holding list of values. It can be specified by using either a {ref}`Closure <script-closure>` or a [Comparator](http://docs.oracle.com/javase/7/docs/api/java/util/Comparator.html) object.

(operator-ifempty)=

## ifEmpty

*Returns: value channel*

The `ifEmpty` operator creates a channel which emits a default value, specified as the operator parameter, when the channel to which is applied is *empty* i.e. doesn't emit any value. Otherwise it will emit the same sequence of entries as the original channel.

Thus, the following example prints:

```groovy
Channel .of(1,2,3) .ifEmpty('Hello') .view()
```

```
1
2
3
```

Instead, this one prints:

```groovy
Channel .empty() .ifEmpty('Hello') .view()
```

```
Hello
```

The `ifEmpty` value parameter can be defined with a {ref}`closure <script-closure>`. In this case the result value of the closure evaluation will be emitted when the empty condition is satisfied.

See also: {ref}`channel-empty` method.

(operator-join)=

## join

*Returns: queue channel*

The `join` operator creates a channel that joins together the items emitted by two channels for which exists a matching key. The key is defined, by default, as the first element in each item emitted.

For example:

```groovy
left  = Channel.of(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
right = Channel.of(['Z', 6], ['Y', 5], ['X', 4])
left.join(right).view()
```

```
[Z, 3, 6]
[Y, 2, 5]
[X, 1, 4]
```

The `index` of a different matching element can be specified by using the `by` parameter.

The `join` operator can emit all the pairs that are incomplete, i.e. the items for which a matching element is missing, by specifying the optional parameter `remainder` as shown below:

```groovy
left  = Channel.of(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
right = Channel.of(['Z', 6], ['Y', 5], ['X', 4])
left.join(right, remainder: true).view()
```

```
[Y, 2, 5]
[Z, 3, 6]
[X, 1, 4]
[P, 7, null]
```

Available options:

`by`
: The index (zero based) of the element to be used as grouping key. A key composed by multiple elements can be defined specifying a list of indices e.g. `by: [0,2]`.

`failOnDuplicate`
: An error is reported when the same key is found more than once.

`failOnMismatch`
: An error is reported when a channel emits a value for which there isn't a corresponding element in the joining channel. This option cannot be used with `remainder`.

`remainder`
: When `false` incomplete tuples (i.e. with less than `size` grouped items) are discarded (default). When `true` incomplete tuples are emitted as the ending emission.

(operator-last)=

## last

*Returns: value channel*

The `last` operator creates a channel that only returns the last item emitted by the source channel. For example:

```groovy
Channel
    .of( 1,2,3,4,5,6 )
    .last()
    .view()
```

```
6
```

(operator-map)=

## map

*Returns: queue channel*

The `map` operator applies a function of your choosing to every item emitted by a channel, and returns the items so obtained as a new channel. The function applied is called the mapping function and is expressed with a {ref}`closure <script-closure>` as shown in the example below:

```groovy
Channel
    .of( 1, 2, 3, 4, 5 )
    .map { it * it }
    .subscribe onNext: { println it }, onComplete: { println 'Done' }
```

```
1
4
9
16
25
Done
```

(operator-max)=

## max

*Returns: value channel*

The `max` operator waits until the source channel completes, and then emits the item that has the greatest value. For example:

```groovy
Channel
    .of( 8, 6, 2, 5 )
    .max()
    .view { "Max value is $it" }
```

```
Max value is 8
```

An optional {ref}`closure <script-closure>` parameter can be specified in order to provide a function that returns the value to be compared. The example below shows how to find the string item that has the maximum length:

```groovy
Channel
    .of("hello","hi","hey")
    .max { it.size() }
    .view()
```

```
"hello"
```

Alternatively it is possible to specify a comparator function i.e. a {ref}`closure <script-closure>` taking two parameters that represent two emitted items to be compared. For example:

```groovy
Channel
    .of("hello","hi","hey")
    .max { a,b -> a.size() <=> b.size() }
    .view()
```

(operator-merge)=

## merge

*Returns: queue channel*

The `merge` operator lets you join items emitted by two (or more) channels into a new channel.

For example, the following code merges two channels together: one which emits a series of odd integers and the other which emits a series of even integers:

```groovy
odds  = Channel.of(1, 3, 5, 7, 9)
evens = Channel.of(2, 4, 6)

odds
    .merge( evens )
    .view()
```

```
[1, 2]
[3, 4]
[5, 6]
```

An optional closure can be provided to customise the items emitted by the resulting merged channel. For example:

```groovy
odds  = Channel.of(1, 3, 5, 7, 9)
evens = Channel.of(2, 4, 6)

odds
    .merge( evens ) { a, b -> tuple(b*b, a) }
    .view()
```

:::{danger}
In general, the use of the `merge` operator is discouraged. Processes and channel operators are not guaranteed to emit items in the order that they were received, due to their parallel and asynchronous nature. Therefore, if you try to merge output channels from different processes, the resulting channel may be different on each run, which will cause resumed runs to not work properly.

You should always use a matching key (e.g. sample ID) to merge multiple channels, so that they are combined in a deterministic way. For this purpose, you can use the [join](#join) operator.
:::

(operator-min)=

## min

*Returns: value channel*

The `min` operator waits until the source channel completes, and then emits the item that has the lowest value. For example:

```groovy
Channel
    .of( 8, 6, 2, 5 )
    .min()
    .view { "Min value is $it" }
```

```
Min value is 2
```

An optional {ref}`closure <script-closure>` parameter can be specified in order to provide a function that returns the value to be compared. The example below shows how to find the string item that has the minimum length:

```groovy
Channel
    .of("hello","hi","hey")
    .min { it.size() }
    .view()
```

```
"hi"
```

Alternatively it is possible to specify a comparator function i.e. a {ref}`closure <script-closure>` taking two parameters that represent two emitted items to be compared. For example:

```groovy
Channel
    .of("hello","hi","hey")
    .min { a,b -> a.size() <=> b.size() }
    .view()
```

(operator-mix)=

## mix

*Returns: queue channel*

The `mix` operator combines the items emitted by two (or more) channels into a single channel.

For example:

```groovy
c1 = Channel.of( 1, 2, 3 )
c2 = Channel.of( 'a', 'b' )
c3 = Channel.of( 'z' )

c1.mix(c2,c3)
    .subscribe onNext: { println it }, onComplete: { println 'Done' }
```

```
1
2
3
'a'
'b'
'z'
```

:::{note}
The items emitted by the resulting mixed channel may appear in any order, regardless of which source channel they came from. Thus, the following example could also be a possible result of the above example:

```
'z'
1
'a'
2
'b'
3
```
:::

(operator-multimap)=

## multiMap

:::{versionadded} 19.11.0-edge
:::

*Returns: map of queue channels*

The `multiMap` operator allows you to forward the items emitted by a source channel to two or more output channels, mapping each input value as a separate element.

The mapping criteria is defined with a {ref}`closure <script-closure>` that specifies the target channels (labelled with a unique identifier) followed by an expression that maps each item from the input channel to the target channel.

For example:

```groovy
Channel.of(1, 2, 3, 4)
    .multiMap { it ->
        foo: it + 1
        bar: it * it
    }
    .set { result }

result.foo.view { "foo $it" }
result.bar.view { "bar $it" }
```

```
foo 2
foo 3
foo 4
foo 5
bar 1
bar 4
bar 9
bar 16
```

The mapping expression can be omitted when the value to be emitted is the same as the following one. If you just need to forward the same value to multiple channels, you can use the following shorthand:

```groovy
Channel
    .of(1,2,3)
    .multiMap { it -> foo: bar: it }
    .set { result }
```

As before, this creates two channels, but now both of them receive the same source items.

You can use the `multiMapCriteria` method to create a multi-map criteria as a variable that can be passed as an argument to one or more `multiMap` operations, as shown below:

```groovy
def criteria = multiMapCriteria {
    small: it < 10
    large: it > 10
}

Channel.of(1, 2, 30).multiMap(criteria).set { ch1 }
Channel.of(10, 20, 1).multiMap(criteria).set { ch2 }
```

:::{note}
If you use `multiMap` to split a tuple or map into multiple channels, it is recommended that you retain a matching key (e.g. sample ID) with *each* new channel, so that you can re-combine these channels later on if needed. In general, you should not expect to be able to merge channels correctly without a matching key, due to the parallel and asynchronous nature of Nextflow pipelines.
:::

(operator-randomsample)=

## randomSample

*Returns: queue channel*

The `randomSample` operator allows you to create a channel emitting the specified number of items randomly taken from the channel to which is applied. For example:

```groovy
Channel
    .of( 1..100 )
    .randomSample( 10 )
    .view()
```

The above snippet will print 10 numbers in the range from 1 to 100.

The operator supports a second parameter that allows you to set the initial `seed` for the random number generator. By setting it, the `randomSample` operator will always return the same pseudo-random sequence. For example:

```groovy
Channel
    .of( 1..100 )
    .randomSample( 10, 234 )
    .view()
```

The above example will print 10 random numbers in the range between 1 and 100. At each run of the script, the same sequence will be returned.

(operator-reduce)=

## reduce

*Returns: value channel*

The `reduce` operator applies a function of your choosing to every item emitted by a channel. Each time this function is invoked it takes two parameters: the accumulated value and the *i-th* emitted item. The result is passed as the accumulated value to the next function call, along with the *i+1 th* item, until all the items are processed.

Finally, the `reduce` operator emits the result of the last invocation of your function as the sole output.

For example:

```groovy
Channel
    .of( 1, 2, 3, 4, 5 )
    .reduce { a, b -> println "a: $a b: $b"; return a+b }
    .view { "result = $it" }
```

```
a: 1 b: 2
a: 3 b: 3
a: 6 b: 4
a: 10 b: 5
result = 15
```

:::{tip}
A common use case for this operator is to use the first parameter as an accumulator and the second parameter as the `i-th` item to be processed.
:::

Optionally you can specify an initial value for the accumulator as shown below:

```groovy
myChannel.reduce( initialValue ) { a, b -> ... }
```

(operator-set)=

## set

*Returns: nothing*

The `set` operator assigns the channel to a variable whose name is specified as a closure parameter. For example:

```groovy
Channel.of(10, 20, 30).set { my_channel }
```

This is semantically equivalent to the following assignment:

```groovy
my_channel = Channel.of(10, 20, 30)
```

However the `set` operator is more idiomatic in Nextflow scripting, since it can be used at the end of a chain of operator transformations, thus resulting in a more fluent and readable operation.

(operator-splitcsv)=

## splitCsv

*Returns: queue channel*

The `splitCsv` operator allows you to parse text items emitted by a channel, that are formatted using the [CSV format](http://en.wikipedia.org/wiki/Comma-separated_values), and split them into records or group them into list of records with a specified length.

In the simplest case just apply the `splitCsv` operator to a channel emitting a CSV formatted text files or text entries. For example:

```groovy
Channel
    .of( 'alpha,beta,gamma\n10,20,30\n70,80,90' )
    .splitCsv()
    .view { row -> "${row[0]} - ${row[1]} - ${row[2]}" }
```

The above example shows hows CSV text is parsed and is split into single rows. Values can be accessed by its column index in the row object.

When the CSV begins with a header line defining the column names, you can specify the parameter `header: true` which allows you to reference each value by its name, as shown in the following example:

```groovy
Channel
    .of( 'alpha,beta,gamma\n10,20,30\n70,80,90' )
    .splitCsv(header: true)
    .view { row -> "${row.alpha} - ${row.beta} - ${row.gamma}" }
```

It will print

```
10 - 20 - 30
70 - 80 - 90
```

Alternatively you can provide custom header names by specifying a the list of strings in the `header` parameter as shown below:

```groovy
Channel
    .of( 'alpha,beta,gamma\n10,20,30\n70,80,90' )
    .splitCsv(header: ['col1', 'col2', 'col3'], skip: 1 )
    .view { row -> "${row.col1} - ${row.col2} - ${row.col3}" }
```

Available options:

`by`
: The number of rows in each `chunk`

`charset`
: Parse the content by using the specified charset e.g. `UTF-8`

`decompress`
: When `true` decompress the content using the GZIP format before processing it (note: files whose name ends with `.gz` extension are decompressed automatically)

`elem`
: The index of the element to split when the operator is applied to a channel emitting list/tuple objects (default: first file object or first element)

`header`
: When `true` the first line is used as columns names. Alternatively it can be used to provide the list of columns names.

`limit`
: Limits the number of retrieved records for each file to the specified value.

`quote`
: Values may be quoted by single or double quote characters.

`sep`
: The character used to separate the values (default: `,`)

`skip`
: Number of lines since the file beginning to ignore when parsing the CSV content.

`strip`
: Removes leading and trailing blanks from values (default: `false`)

(operator-splitfasta)=

## splitFasta

*Returns: queue channel*

The `splitFasta` operator allows you to split the entries emitted by a channel, that are formatted using the [FASTA format](http://en.wikipedia.org/wiki/FASTA_format). It returns a channel which emits text item for each sequence in the received FASTA content.

The number of sequences in each text chunk produced by the `splitFasta` operator can be set by using the `by` parameter. The following example shows how to read a FASTA file and split it into chunks containing 10 sequences each:

```groovy
Channel
     .fromPath('misc/sample.fa')
     .splitFasta( by: 10 )
     .view()
```

:::{warning}
Chunks are stored in memory by default. When splitting large files, specify the parameter `file: true` to save the chunks into files in order to avoid an `OutOfMemoryException`. See the parameter table below for details.
:::

A second version of the `splitFasta` operator allows you to split a FASTA content into record objects, instead of text chunks. A record object contains a set of fields that let you access and manipulate the FASTA sequence information with ease.

In order to split a FASTA content into record objects, simply use the `record` parameter specifying the map of required the fields, as shown in the example below:

```groovy
Channel
     .fromPath('misc/sample.fa')
     .splitFasta( record: [id: true, seqString: true ])
     .filter { record -> record.id =~ /^ENST0.*/ }
     .view { record -> record.seqString }
```

In this example, the file `misc/sample.fa` is split into records containing the `id` and the `seqString` fields (i.e. the sequence id and the sequence data). The following `filter` operator only keeps the sequences whose ID starts with the `ENST0` prefix, finally the sequence content is printed by using the `subscribe` operator.

Available options:

`by`
: Defines the number of sequences in each `chunk` (default: `1`)

`charset`
: Parse the content by using the specified charset e.g. `UTF-8`.

`compress`
: When `true` resulting file chunks are GZIP compressed. The `.gz` suffix is automatically added to chunk file names.

`decompress`
: When `true`, decompress the content using the GZIP format before processing it (note: files whose name ends with `.gz` extension are decompressed automatically).

`elem`
: The index of the element to split when the operator is applied to a channel emitting list/tuple objects (default: first file object or first element).

`file`
: When `true` saves each split to a file. Use a string instead of `true` value to create split files with a specific name (split index number is automatically added). Finally, set this attribute to an existing directory, in order to save the split files into the specified folder.

`limit`
: Limits the number of retrieved sequences for each file to the specified value.

`record`
: Parse each entry in the FASTA file as record objects. The following fields are available:

  - `id`: The FASTA sequence identifier i.e. the word following the `>` symbol up to the first `blank` or `newline` character
  - `header`: The first line in a FASTA sequence without the `>` character
  - `desc`: The text in the FASTA header following the ID value
  - `text`: The complete FASTA sequence including the header
  - `seqString`: The sequence data as a single line string i.e. containing no `newline` characters
  - `sequence`: The sequence data as a multi-line string (always ending with a `newline` character)
  - `width`: Define the length of a single line when the `sequence` field is used, after that the sequence data continues on a new line.

`size`
: Defines the size in memory units of the expected chunks e.g. `1.MB`.

See also: [countFasta](#countfasta)

(operator-splitfastq)=

## splitFastq

*Returns: queue channel*

The `splitFastq` operator allows you to split the entries emitted by a channel, that are formatted using the [FASTQ format](http://en.wikipedia.org/wiki/FASTQ_format). It returns a channel which emits a text chunk for each sequence in the received item.

The number of sequences in each text chunk produced by the `splitFastq` operator is defined by the parameter `by`. The following example shows you how to read a FASTQ file and split it into chunks containing 10 sequences each:

```groovy
Channel
    .fromPath('misc/sample.fastq')
    .splitFastq( by: 10 )
    .view()
```

:::{warning}
Chunks are stored in memory by default. When splitting large files, specify the parameter `file: true` to save the chunks into files in order to avoid an `OutOfMemoryException`. See the parameter table below for details.
:::

A second version of the `splitFastq` operator allows you to split a FASTQ formatted content into record objects, instead of text chunks. A record object contains a set of fields that let you access and manipulate the FASTQ sequence data with ease.

In order to split FASTQ sequences into record objects simply use the `record` parameter specifying the map of the required fields, or just specify `record: true` as in the example shown below:

```groovy
Channel
    .fromPath('misc/sample.fastq')
    .splitFastq( record: true )
    .view { record -> record.readHeader }
```

Finally the `splitFastq` operator is able to split paired-end read pair FASTQ files. It must be applied to a channel which emits tuples containing at least two elements that are the files to be split. For example:

```groovy
Channel
    .fromFilePairs('/my/data/SRR*_{1,2}.fastq', flat: true)
    .splitFastq(by: 100_000, pe: true, file: true)
    .view()
```

:::{note}
The `fromFilePairs` requires the `flat: true` option in order to emit the file pairs as separate elements in the produced tuples.
:::

:::{note}
This operator assumes that the order of the paired-end reads correspond with each other and both files contain the same number of reads.
:::

Available options:

`by`
: Defines the number of *reads* in each `chunk` (default: `1`)

`charset`
: Parse the content by using the specified charset e.g. `UTF-8`

`compress`
: When `true` resulting file chunks are GZIP compressed. The `.gz` suffix is automatically added to chunk file names.

`decompress`
: When `true` decompress the content using the GZIP format before processing it (note: files whose name ends with `.gz` extension are decompressed automatically)

`elem`
: The index of the element to split when the operator is applied to a channel emitting list/tuple objects (default: first file object or first element)

`file`
: When `true` saves each split to a file. Use a string instead of `true` value to create split files with a specific name (split index number is automatically added). Finally, set this attribute to an existing directory, in order to save the split files into the specified folder.

`limit`
: Limits the number of retrieved *reads* for each file to the specified value.

`pe`
: When `true` splits paired-end read files, therefore items emitted by the source channel must be tuples in which at least two elements are the read-pair files to be split.

`record`
: Parse each entry in the FASTQ file as record objects. The following fields are available:

  - `readHeader`: Sequence header (without the `@` prefix)
  - `readString`: The raw sequence data
  - `qualityHeader`: Base quality header (it may be empty)
  - `qualityString`: Quality values for the sequence

See also: [countFastq](#countfastq)

(operator-splitjson)=

## splitJson

*Returns: queue channel*

The `splitJson` operator allows you to split a JSON document from a source channel into individual records. If the document is a JSON array, each element of the array will be emitted. If the document is a JSON object, each key-value pair will be emitted as a map with the properties `key`  and `value`.

An example with a JSON array:

```groovy
Channel.of('[1,null,["A",{}],true]')
    .splitJson()
    .view{"Item: ${it}"}
```

Produces the following output:

```
Item: 1
Item: null
Item: [A, [:]]
Item: true
```

An example with a JSON object:

```groovy
Channel.of('{"A":1,"B":[1,2,3],"C":{"D":null}}')
    .splitJson()
    .view{"Item: ${it}"}
```

Produces the following output:

```
Item: [value:1, key:A]
Item: [value:[1, 2, 3], key:B]
Item: [value:[D:null], key:C]
```

You can optionally query a section of the JSON document to parse and split, using the `path` option:

```groovy
Channel.of('{"A":1,"B":[2,3,{"C":{"D":null,"E":4,"F":5}}]}')
    .splitJson(path: "B[2].C")
    .view{"Item: ${it}"}
```

Produces the following output:

```
Item: [value:null, key:D]
Item: [value:4, key:E]
Item: [value:5, key:F]
```

Available options:

`limit`
: Limits the number of retrieved lines for each file to the specified value.

`path`
: Define the section of the JSON document that you want to extract. The expression is a set of paths separated by a dot, similar to [JSONPath](https://goessner.net/articles/JsonPath/). The empty string is the document root (default). An integer in brackets is the 0-based index in a JSON array. A string preceded by a dot `.` is the key in a JSON object.

See also: [countJson](#countjson)

(operator-splittext)=

## splitText

*Returns: queue channel*

The `splitText` operator allows you to split multi-line strings or text file items, emitted by a source channel into chunks containing `n` lines, which will be emitted by the resulting channel.

For example:

```groovy
Channel
    .fromPath('/some/path/*.txt')
    .splitText()
    .view()
```

It splits the content of the files with suffix `.txt`, and prints it line by line.

By default the `splitText` operator splits each item into chunks of one line. You can define the number of lines in each chunk by using the parameter `by`, as shown in the following example:

```groovy
Channel
    .fromPath('/some/path/*.txt')
    .splitText( by: 10 )
    .subscribe {
        print it;
        print "--- end of the chunk ---\n"
    }
```

An optional {ref}`closure <script-closure>` can be specified in order to transform the text chunks produced by the operator. The following example shows how to split text files into chunks of 10 lines and transform them to capital letters:

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
: Parse the content by using the specified charset e.g. `UTF-8`.

`compress`
: When `true` resulting file chunks are GZIP compressed. The `.gz` suffix is automatically added to chunk file names.

`decompress`
: When `true`, decompress the content using the GZIP format before processing it (note: files whose name ends with `.gz` extension are decompressed automatically).

`elem`
: The index of the element to split when the operator is applied to a channel emitting list/tuple objects (default: first file object or first element).

`file`
: When `true` saves each split to a file. Use a string instead of `true` value to create split files with a specific name (split index number is automatically added). Finally, set this attribute to an existing directory, in order to save the split files into the specified folder.

`keepHeader`
: Parses the first line as header and prepends it to each emitted chunk.

`limit`
: Limits the number of retrieved lines for each file to the specified value.

See also: [countLines](#countlines)

(operator-subscribe)=

## subscribe

*Returns: nothing*

The `subscribe` operator allows you to execute a user defined function each time a new value is emitted by the source channel.

The emitted value is passed implicitly to the specified function. For example:

```groovy
// define a channel emitting three values
source = Channel.of( 'alpha', 'beta', 'delta' )

// subscribe a function to the channel printing the emitted values
source.subscribe { println "Got: $it" }
```

```
Got: alpha
Got: beta
Got: delta
```

:::{note}
In Groovy, the language on which Nextflow is based, the user defined function is called a **closure**. Read the {ref}`script-closure` section to learn more about closures.
:::

If needed the closure parameter can be defined explicitly, using a name other than `it` and, optionally, specifying the expected value type, as shown in the following example:

```groovy
Channel
    .of( 'alpha', 'beta', 'lambda' )
    .subscribe { String str ->
        println "Got: ${str}; len: ${str.size()}"
    }
```

```
Got: alpha; len: 5
Got: beta; len: 4
Got: lambda; len: 6
```

The `subscribe` operator may accept one or more of the following event handlers:

- `onNext`: function that is invoked whenever the channel emits a value. Equivalent to using the `subscribe` with a plain closure as described in the examples above.
- `onComplete`: function that is invoked after the last value is emitted by the channel.
- `onError`: function that it is invoked when an exception is raised while handling the `onNext` event. It will not make further calls to `onNext` or `onComplete`. The `onError` method takes as its parameter the `Throwable` that caused the error.

For example:

```groovy
Channel
    .of( 1, 2, 3 )
    .subscribe onNext: { println it }, onComplete: { println 'Done' }
```

```
1
2
3
Done
```

(operator-sum)=

## sum

*Returns: value channel*

The `sum` operator creates a channel that emits the sum of all the items emitted by the channel itself. For example:

```groovy
Channel
    .of( 8, 6, 2, 5 )
    .sum()
    .view { "The sum is $it" }
```

```
The sum is 21
```

An optional {ref}`closure <script-closure>` parameter can be specified in order to provide a function that, given an item, returns the value to be summed. For example:

```groovy
Channel
    .of( 4, 1, 7, 5 )
    .sum { it * it }
    .view { "Square: $it" }
```

```
Square: 91
```

## take

*Returns: queue channel*

The `take` operator allows you to filter only the first `n` items emitted by a channel. For example:

```groovy
Channel
    .of( 1, 2, 3, 4, 5, 6 )
    .take( 3 )
    .subscribe onNext: { println it }, onComplete: { println 'Done' }
```

```
1
2
3
Done
```

:::{tip}
Specifying a size of `-1` causes the operator to take all values.
:::

See also [until](#until).

## tap

*Returns: queue channel*

The `tap` operator is like the [set](#set) operator in that it assigns a source channel to a new target channel.
but it also emits the source channel for downstream use. This operator is a useful way to extract intermediate
output channels from a chain of operators. For example:

```groovy
Channel
    .of ( 'a', 'b', 'c' )
    .tap { log1 }
    .map { it * 2 }
    .tap { log2 }
    .map { it.toUpperCase() }
    .view { "Result: $it" }

log1.view { "Log 1: $it" }
log2.view { "Log 2: $it" }
```

```
Result: AA
Result: BB
Result: CC

Log 1: a
Log 1: b
Log 1: c

Log 2: aa
Log 2: bb
Log 2: cc
```

## toInteger

*Returns: queue channel*

The `toInteger` operator allows you to convert the string values emitted by a channel to `Integer` values. For example:

```groovy
Channel
    .of( '1', '7', '12' )
    .toInteger()
    .sum()
    .view()
```

:::{tip}
You can also use `toLong`, `toFloat`, and `toDouble` to convert to other numerical types.
:::

## toList

*Returns: value channel*

The `toList` operator collects all the items emitted by a channel to a `List` object and emits the resulting collection as a single item. For example:

```groovy
Channel
    .of( 1, 2, 3, 4 )
    .toList()
    .subscribe onNext: { println it }, onComplete: { println 'Done' }
```

```
[1,2,3,4]
Done
```

:::{note}
There are two differences between `toList` and `collect`:

- When there is no input, `toList` emits an empty list whereas `collect` emits nothing.
- By default, `collect` flattens list items by one level.

In other words, `toList` is equivalent to:

```groovy
collect(flat: false).ifEmpty([])
```
:::

See also: [collect](#collect) operator.

## toSortedList

*Returns: value channel*

The `toSortedList` operator collects all the items emitted by a channel to a `List` object where they are sorted and emits the resulting collection as a single item. For example:

```groovy
Channel
    .of( 3, 2, 1, 4 )
    .toSortedList()
    .subscribe onNext: { println it }, onComplete: { println 'Done' }
```

```
[1,2,3,4]
Done
```

You may also pass a comparator closure as an argument to the `toSortedList` operator to customize the sorting criteria. For example, to sort by the second element of a tuple in descending order:

```groovy
Channel
    .of( ["homer", 5], ["bart", 2], ["lisa", 10], ["marge", 3], ["maggie", 7] )
    .toSortedList( { a, b -> b[1] <=> a[1] } )
    .view()
```

```
[[lisa, 10], [maggie, 7], [homer, 5], [marge, 3], [bart, 2]]
```

See also: [collect](#collect) operator.

## transpose

*Returns: queue channel*

The `transpose` operator transforms a channel in such a way that the emitted items are the result of a transposition of all tuple elements in each item. For example:

```groovy
Channel.of(
        [1, ['A', 'B', 'C']],
        [2, ['C', 'A']],
        [3, ['B', 'D']]
    )
    .transpose()
    .view()
```

The above snippet prints:

```
[1, A]
[1, B]
[1, C]
[2, C]
[2, A]
[3, B]
[3, D]
```

Available options:

` by`
: The index (zero based) of the element to be transposed. Multiple elements can be defined specifying as list of indices e.g. `by: [0,2]`

` remainder`
: When `false` incomplete tuples are discarded (default). When `true` incomplete tuples are emitted containing a `null` in place of a missing element.

## unique

*Returns: queue channel*

The `unique` operator allows you to remove duplicate items from a channel and only emit single items with no repetition.

For example:

```groovy
Channel
    .of( 1, 1, 1, 5, 7, 7, 7, 3, 3 )
    .unique()
    .view()
```

```
1
5
7
3
```

You can also specify an optional {ref}`closure <script-closure>` that customizes the way it distinguishes between unique items. For example:

```groovy
Channel
    .of(1, 3, 4, 5)
    .unique { it % 2 }
    .view()
```

```
1
4
```

## until

*Returns: queue channel*

The `until` operator creates a channel that returns the items emitted by the source channel and stop when the condition specified is verified. For example:

```groovy
Channel
    .of( 3, 2, 1, 5, 1, 5 )
    .until { it == 5 }
    .view()
```

```
3
2
1
```

See also [take](#take).

(operator-view)=

## view

*Returns: queue channel*

The `view` operator prints the items emitted by a channel to the console standard output. For example:

```groovy
Channel.of(1, 2, 3).view()

1
2
3
```

Each item is printed on a separate line unless otherwise specified by using the `newLine: false` optional parameter.

How the channel items are printed can be controlled by using an optional closure parameter. The closure must return the actual value of the item to be printed:

```groovy
Channel.of(1, 2, 3)
    .map { it -> [it, it*it] }
    .view { num, sqr -> "Square of: $num is $sqr" }
```

```
Square of: 1 is 1
Square of: 2 is 4
Square of: 3 is 9
```

The `view` operator also emits every item that it receives, allowing it to be chained with other operators.

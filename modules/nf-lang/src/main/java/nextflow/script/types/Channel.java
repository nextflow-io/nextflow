/*
 * Copyright 2024-2025, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nextflow.script.types;

import java.nio.file.Path;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import groovy.lang.Closure;
import nextflow.script.dsl.Description;
import nextflow.script.dsl.Operator;

@Description("""
    A `Channel` is an asynchronous collection that is produced by a process, operator, or channel factory.

    [Read more](https://nextflow.io/docs/latest/reference/channel.html)
""")
public abstract class Channel<E> {

    protected static ChannelFactory instance;

    @Description("""
        Create a channel that emits nothing.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#empty)
    """)
    public static Channel empty() {
        return instance.empty();
    }

    @Deprecated
    @Description("""
        Create a channel that emits each argument.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#from)
    """)
    public static <E> Channel<E> from(E... values) {
        return instance.from(values);
    }

    @Deprecated
    @Description("""
        Create a channel that emits each element in a collection.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#from)
    """)
    public static <E> Channel<E> from(Collection<E> values) {
        return instance.from(values);
    }

    @Description("""
        Create a channel that emits all file pairs matching a glob pattern.

        An optional closure can be used to customize the grouping strategy.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
    """)
    public static Channel fromFilePairs(Map<String,?> opts, String pattern, Closure grouping) {
        return instance.fromFilePairs(opts, pattern, grouping);
    }

    @Description("""
        Create a channel that emits each element in a collection.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#fromlist)
    """)
    public static <E> Channel<E> fromList(Collection<E> values) {
        return instance.fromList(values);
    }

    @Description("""
        Create a channel that emits all paths matching a name or glob pattern.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#frompath)
    """)
    public static Channel<Path> fromPath(Map<String,?> opts, String pattern) {
        return instance.fromPath(opts, pattern);
    }

    @Description("""
        Create a channel that queries the [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) database and emits all FASTQ files matching the given project or accession ids.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#fromsra)
    """)
    public static Channel fromSRA(Map<String,?> opts, String query) {
        return instance.fromSRA(opts, query);
    }

    @Description("""
        Create a channel that emits each argument.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#of)
    """)
    public static <E> Channel<E> of(E... values) {
        return instance.of(values);
    }

    @Description("""
        Create a channel that emits all values in the given topic.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#topic)
    """)
    public static Channel topic(String name) {
        return instance.topic(name);
    }

    @Description("""
        Create a value channel.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#value)
    """)
    public static <E> Channel<E> value(E value) {
        return instance.value(value);
    }

    @Description("""
        Create a channel that watches for filesystem events for all files matching the given pattern.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#watchpath)
    """)
    public static Channel<Path> watchPath(String filePattern, String events) {
        return instance.watchPath(filePattern, events);
    }

    @Operator
    @Description("""
        The `branch` operator forwards each value from a source channel to one of multiple output channels, based on a selection criteria.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#branch)
    """)
    public abstract Object branch(Closure action);

    @Operator
    @Description("""
        The `buffer` operator collects values from a source channel into subsets and emits each subset separately.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#buffer)
    """)
    public abstract Channel buffer(Closure openingCondition, Closure closingCondition);

    @Operator
    @Description("""
        The `collate` operator collects values from a source channel into groups of *N* values.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#collate)
    """)
    public abstract Channel collate(int size, int step, boolean remainder);

    @Operator
    @Description("""
        The `collect` operator collects all values from a source channel into a list and emits it as a single value.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#collect)
    """)
    public abstract Channel collect(Closure action);

    @Operator
    @Description("""
        The `collectFile` operator collects the values from a source channel and saves them to one or more files, emitting the collected file(s).

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#collectfile)
    """)
    public abstract Channel collectFile(Map<String,?> opts, Closure closure);

    @Operator
    @Description("""
        The `combine` operator produces the combinations (i.e. cross product, “Cartesian” product) of two source channels, or a channel and a list (as the right operand), emitting each combination separately.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#combine)
    """)
    public abstract Channel combine(Map<String,?> opts, Object right);

    @Operator
    @Description("""
        The `concat` operator emits the values from two or more source channels into a single output channel. Each source channel is emitted in the order in which it was specified.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#concat)
    """)
    public abstract Channel concat(Channel... others);

    @Operator
    @Description("""
        The `count` operator computes the total number of values from a source channel and emits it.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#count)
    """)
    public abstract Channel count();

    @Operator
    @Description("""
        The `cross` operator emits every pairwise combination of two channels for which the pair has a matching key.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#cross)
    """)
    public abstract Channel cross(Channel right, Closure mapper);

    @Operator
    @Description("""
        The `distinct` operator forwards a source channel with consecutively repeated values removed, such that each emitted value is different from the preceding one.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#distinct)
    """)
    public abstract Channel distinct();

    @Operator
    @Description("""
        When the pipeline is executed with the `-dump-channels` command-line option, the `dump` operator prints each value in a source channel, otherwise it does nothing.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#dump)
    """)
    public abstract Channel dump(Map<String,?> opts);

    @Operator
    @Description("""
        The `filter` operator emits the values from a source channel that satisfy a condition, discarding all other values. The filter condition can be a literal value, a regular expression, a type qualifier, or a boolean predicate.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#filter)
    """)
    public abstract Channel filter(Closure<Boolean> closure);

    @Operator
    @Description("""
        The `first` operator emits the first value from a source channel, or the first value that satisfies a condition. The condition can be a regular expression, a type qualifier (i.e. Java class), or a boolean predicate.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#first)
    """)
    public abstract Channel first(Object criteria);

    @Operator
    @Description("""
        The `flatMap` operator applies a mapping function to each value from a source channel.
        
        When the mapping function returns a list, each element in the list is emitted separately. When the mapping function returns a map, each key-value pair in the map is emitted separately.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#flatmap)
    """)
    public abstract Channel flatMap(Closure closure);

    @Operator
    @Description("""
        The `flatten` operator flattens each value from a source channel that is a list or other collection, such that each element in each collection is emitted separately. Deeply nested collections are also flattened.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#flatten)
    """)
    public abstract Channel flatten();

    @Operator
    @Description("""
        The `groupTuple` operator collects tuples from a source channel into groups based on a grouping key. A new tuple is emitted for each distinct key.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#grouptuple)
    """)
    public abstract Channel groupTuple(Map<String,?> opts);

    @Operator
    @Description("""
        The `ifEmpty` operator emits a source channel, or a default value if the source channel is empty.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#ifempty)
    """)
    public abstract Channel ifEmpty(Object value);

    @Operator
    @Description("""
        The `join` operator emits the inner product of two source channels using a matching key.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#join)
    """)
    public abstract Channel join(Channel right);

    @Operator
    @Description("""
        The `last` operator emits the last value from a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#last)
    """)
    public abstract Channel last();

    @Operator
    @Description("""
        The `map` operator applies a mapping function to each value from a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#map)
    """)
    public abstract Channel map(Closure closure);

    @Operator
    @Description("""
        The `max` operator emits the item with the greatest value from a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#max)
    """)
    public abstract Channel max(Closure comparator);

    @Deprecated
    @Operator
    @Description("""
        The `merge` operator joins the values from two or more channels into a new channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#merge)
    """)
    public abstract Channel merge(Channel... others);

    @Operator
    @Description("""
        The `min` operator emits the item with the lowest value from a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#min)
    """)
    public abstract Channel min(Closure comparator);

    @Operator
    @Description("""
        The `mix` operator emits the values from two or more source channels into a single output channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#mix)
    """)
    public abstract Channel mix(Channel... others);

    @Operator
    @Description("""
        The `multiMap` operator applies a set of mapping functions to a source channel, producing a separate output channel for each mapping function.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#multimap)
    """)
    public abstract Object multiMap(Closure action);

    @Operator
    @Description("""
        The `randomSample` operator emits a randomly-selected subset of values from a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#randomsample)
    """)
    public abstract Channel randomSample(int n, Long seed);

    @Operator
    @Description("""
        The `reduce` operator applies an accumulator function sequentially to each value from a source channel, and emits the accumulated value. The accumulator function takes two parameters -- the accumulated value and the *i*-th emitted value -- and it should return the accumulated result, which is passed to the next invocation with the *i+1*-th value. This process is repeated for each value in the source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#reduce)
    """)
    public abstract Channel reduce(Object seed, Closure closure);

    @Operator
    @Description("""
        The `set` operator assigns a source channel to a variable, whose name is specified in a closure.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#set)
    """)
    public abstract void set(Closure holder);

    @Operator
    @Description("""
        The `splitCsv` operator parses and splits [CSV-formatted](http://en.wikipedia.org/wiki/Comma-separated_values) text from a source channel into records, or groups of records with a given size.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#splitcsv)
    """)
    public abstract Channel splitCsv(Map<String,?> opts);

    @Operator
    @Description("""
        The `splitFasta` operator splits [FASTA formatted](http://en.wikipedia.org/wiki/FASTA_format) text from a source channel into individual sequences.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#splitfasta)
    """)
    public abstract Channel splitFasta(Map<String,?> opts);

    @Operator
    @Description("""
        The `splitFastq` operator splits [FASTQ formatted](http://en.wikipedia.org/wiki/FASTQ_format) text from a source channel into individual sequences.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#splitfastq)
    """)
    public abstract Channel splitFastq(Map<String,?> opts);

    @Operator
    @Description("""
        The `splitText` operator splits multi-line text content from a source channel into chunks of *N* lines.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#splittext)
    """)
    public abstract Channel splitText(Map<String,?> opts, Closure action);

    @Operator
    @Description("""
        The `subscribe` operator invokes a custom function for each value in a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#subscribe)
    """)
    public abstract void subscribe(Closure closure);

    @Operator
    @Description("""
        The `sum` operator emits the sum of all values in a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#sum)
    """)
    public abstract Channel sum(Closure closure);

    @Operator
    @Description("""
        The `take` operator takes the first *N* values from a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#take)
    """)
    public abstract Channel take(int n);

    @Operator
    @Description("""
        The `tap` operator assigns a source channel to a variable, whose name is specified in a closure.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#tap)
    """)
    public abstract Channel tap(Closure holder);

    @Operator
    @Description("""
        The `toList` operator collects all the values from a source channel into a list and emits the list as a single value.

   public abstract      [Read more](https://nextflow.io/docs/latest/reference/operator.html#to;ist)
    """)
    public abstract Channel toList();

    @Operator
    @Description("""
        The `toSortedList` operator collects all the values from a source channel into a sorted list and emits the list as a single value.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#tosortedlist)
    """)
    public abstract Channel toSortedList();

    @Operator
    @Description("""
        The `transpose` operator transposes each tuple from a source channel by flattening any nested list in each tuple, emitting each nested value separately.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#transpose)
    """)
    public abstract Channel transpose(Map<String,?> opts);

    @Operator
    @Description("""
        The `unique` operator emits the unique values from a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#unique)
    """)
    public abstract Channel unique(Closure comparator);

    @Operator
    @Description("""
        The `until` operator emits each value from a source channel until a stopping condition is satisfied.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#until)
    """)
    public abstract Channel until(Closure<Boolean> closure);

    @Operator
    @Description("""
        The `view` operator prints each value from a source channel to standard output.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#view)
    """)
    public abstract Channel view(Closure closure);

}

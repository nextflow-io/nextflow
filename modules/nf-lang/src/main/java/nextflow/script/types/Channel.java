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
import java.util.function.BiFunction;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;

import groovy.lang.Closure;
import nextflow.script.dsl.Description;
import nextflow.script.dsl.Operator;

@Description("""
    A channel is an asynchronous sequence of values. It is used to facilitate dataflow logic in a workflow.

    [Read more](https://nextflow.io/docs/latest/reference/stdlib-types.html#channel-e)
""")
public interface Channel<E> {

    @Operator
    @Description("""
        The `branch` operator forwards each value from a source channel to one of multiple output channels, based on a selection criteria.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#branch)
    """)
    Record branch(Closure closure);

    @Operator
    @Description("""
        The `buffer` operator collects values from a source channel into subsets and emits each subset separately.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#buffer)
    """)
    Channel buffer(Closure openingCondition, Closure closingCondition);

    @Operator
    @Description("""
        The `collate` operator collects values from a source channel into groups of *N* values.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#collate)
    """)
    Channel collate(int size, int step, boolean remainder);

    @Operator
    @Description("""
        The `collect` operator collects all values from a source channel into a list and emits it as a single value.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#collect)
    """)
    Value<Bag<E>> collect();

    @Operator
    @Description("""
        The `collectFile` operator collects the values from a source channel and saves them to one or more files, emitting the collected file(s).

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#collectfile)
    """)
    Channel collectFile(Map<String,?> opts, Closure closure);

    @Operator
    @Description("""
        The `combine` operator produces the combinations (i.e. cross product, “Cartesian” product) of two source channels, or a channel and a list (as the right operand), emitting each combination separately.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#combine)
    """)
    Channel combine(Map<String,?> opts, Channel right);
    Channel combine(Channel right);
    Channel combine(Collection right);

    @Operator
    @Description("""
        The `concat` operator emits the values from two or more source channels into a single output channel. Each source channel is emitted in the order in which it was specified.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#concat)
    """)
    Channel concat(Channel... others);

    @Operator
    @Description("""
        The `count` operator computes the total number of values from a source channel and emits it.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#count)
    """)
    Value<Integer> count();

    @Operator
    @Description("""
        The `cross` operator emits every pairwise combination of two channels for which the pair has a matching key.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#cross)
    """)
    Channel cross(Channel right);

    @Operator
    @Description("""
        The `distinct` operator forwards a source channel with consecutively repeated values removed, such that each emitted value is different from the preceding one.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#distinct)
    """)
    Channel distinct();

    @Operator
    @Description("""
        When the pipeline is executed with the `-dump-channels` command-line option, the `dump` operator prints each value in a source channel, otherwise it does nothing.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#dump)
    """)
    Channel<E> dump(Map<String,?> opts);

    @Operator
    @Description("""
        The `filter` operator emits the values from a source channel that satisfy a condition, discarding all other values. The filter condition can be a literal value, a regular expression, a type qualifier, or a boolean predicate.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#filter)
    """)
    Channel<E> filter(Predicate<E> condition);

    @Operator
    @Description("""
        The `first` operator emits the first value from a source channel, or the first value that satisfies a condition. The condition can be a regular expression, a type qualifier (i.e. Java class), or a boolean predicate.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#first)
    """)
    Value<E> first(Predicate<E> criteria);
    Value<E> first();

    @Operator
    @Description("""
        The `flatMap` operator applies a mapping function to each value from a source channel.
        
        When the mapping function returns a list, each element in the list is emitted separately. When the mapping function returns a map, each key-value pair in the map is emitted separately.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#flatmap)
    """)
    <R> Channel<R> flatMap(Function<E,Iterable<R>> transform);

    @Operator
    @Description("""
        The `flatten` operator flattens each value from a source channel that is a list or other collection, such that each element in each collection is emitted separately. Deeply nested collections are also flattened.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#flatten)
    """)
    Channel flatten();

    @Operator
    @Description("""
        The `groupTuple` operator collects tuples from a source channel into groups based on a grouping key. A new tuple is emitted for each distinct key.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#grouptuple)
    """)
    Channel groupTuple(Map<String,?> opts);

    @Operator
    @Description("""
        The `ifEmpty` operator emits a source channel, or a default value if the source channel is empty.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#ifempty)
    """)
    Channel ifEmpty(Object value);

    @Operator
    @Description("""
        The `join` operator emits the inner product of two source channels using a matching key.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#join)
    """)
    Channel join(Map<String,?> opts, Channel right);
    Channel join(Channel right);

    @Operator
    @Description("""
        The `last` operator emits the last value from a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#last)
    """)
    Value<E> last();

    @Operator
    @Description("""
        The `map` operator applies a mapping function to each value from a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#map)
    """)
    <R> Channel<R> map(Function<E,R> transform);

    @Operator
    @Description("""
        The `max` operator emits the item with the greatest value from a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#max)
    """)
    Value<E> max(Closure comparator);

    @Operator
    @Description("""
        The `merge` operator joins the values from two or more channels into a new channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#merge)
    """)
    Channel merge(Channel... others);

    @Operator
    @Description("""
        The `min` operator emits the item with the lowest value from a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#min)
    """)
    Value<E> min(Closure comparator);

    @Operator
    @Description("""
        The `mix` operator emits the values from two or more source channels into a single output channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#mix)
    """)
    Channel<E> mix(Channel<E>... others);

    @Operator
    @Description("""
        The `multiMap` operator applies a set of mapping functions to a source channel, producing a separate output channel for each mapping function.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#multimap)
    """)
    Record multiMap(Closure closure);

    @Operator
    @Description("""
        The `randomSample` operator emits a randomly-selected subset of values from a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#randomsample)
    """)
    Channel randomSample(int n, Long seed);

    @Operator
    @Description("""
        The `reduce` operator applies an accumulator function sequentially to each value from a source channel, and emits the accumulated value. The accumulator function takes two parameters -- the accumulated value and the *i*-th emitted value -- and it should return the accumulated result, which is passed to the next invocation with the *i+1*-th value. This process is repeated for each value in the source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#reduce)
    """)
    <R> Value<R> reduce(R seed, BiFunction<R,E,R> accumulator);
    <R> Value<R> reduce(BiFunction<R,E,R> accumulator);

    @Operator
    @Description("""
        The `set` operator assigns a source channel to a variable, whose name is specified in a closure.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#set)
    """)
    void set(Closure holder);

    @Operator
    @Description("""
        The `splitCsv` operator parses and splits [CSV-formatted](http://en.wikipedia.org/wiki/Comma-separated_values) text from a source channel into records, or groups of records with a given size.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#splitcsv)
    """)
    Channel splitCsv(Map<String,?> opts);
    Channel splitCsv();

    @Operator
    @Description("""
        The `splitFasta` operator splits [FASTA formatted](http://en.wikipedia.org/wiki/FASTA_format) text from a source channel into individual sequences.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#splitfasta)
    """)
    Channel splitFasta(Map<String,?> opts);
    Channel splitFasta();

    @Operator
    @Description("""
        The `splitFastq` operator splits [FASTQ formatted](http://en.wikipedia.org/wiki/FASTQ_format) text from a source channel into individual sequences.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#splitfastq)
    """)
    Channel splitFastq(Map<String,?> opts);
    Channel splitFastq();

    @Operator
    @Description("""
        The splitJson operator splits JSON files or text from a source channel into individual records.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#splitjson)
    """)
    Channel splitJson(Map<String,?> opts);
    Channel splitJson();

    @Operator
    @Description("""
        The `splitText` operator splits multi-line text content from a source channel into chunks of *N* lines.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#splittext)
    """)
    Channel splitText(Map<String,?> opts);
    Channel splitText();

    @Operator
    @Description("""
        The `subscribe` operator invokes a custom function for each value in a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#subscribe)
    """)
    void subscribe(Consumer<E> action);

    @Operator
    @Description("""
        The `sum` operator emits the sum of all values in a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#sum)
    """)
    Value sum(Closure transform);

    @Operator
    @Description("""
        The `take` operator takes the first *N* values from a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#take)
    """)
    Channel take(int n);

    @Operator
    @Description("""
        The `tap` operator assigns a source channel to a variable, whose name is specified in a closure.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#tap)
    """)
    Channel tap(Closure holder);

    @Operator
    @Description("""
        The `toList` operator collects all the values from a source channel into a list and emits the list as a single value.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#to;ist)
    """)
    Value toList();

    @Operator
    @Description("""
        The `toSortedList` operator collects all the values from a source channel into a sorted list and emits the list as a single value.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#tosortedlist)
    """)
    Value toSortedList();

    @Operator
    @Description("""
        The `transpose` operator transposes each tuple from a source channel by flattening any nested list in each tuple, emitting each nested value separately.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#transpose)
    """)
    Channel transpose(Map<String,?> opts);

    @Operator
    @Description("""
        The `unique` operator emits the unique values from a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#unique)
    """)
    Channel unique(Closure comparator);
    Channel unique();

    @Operator
    @Description("""
        The `until` operator emits each value from a source channel until a stopping condition is satisfied.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#until)
    """)
    Channel<E> until(Predicate<E> condition);

    @Operator
    @Description("""
        The `view` operator prints each value from a source channel to standard output.

        [Read more](https://nextflow.io/docs/latest/reference/operator.html#view)
    """)
    Channel<E> view(Function<E,String> transform);
    Channel<E> view();

}

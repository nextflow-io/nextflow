/*
 * Copyright 2013-2026, Seqera Labs
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

import java.util.Map;
import java.util.function.BiFunction;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;

import groovy.lang.Closure;
import groovy.transform.NamedParam;
import groovy.transform.NamedParams;
import nextflow.script.dsl.Description;
import nextflow.script.dsl.Operator;
import nextflow.script.types.Record;

@Description("""
    A channel is an asynchronous sequence of values. It is used to facilitate dataflow logic in a workflow.

    [Read more](https://nextflow.io/docs/latest/reference/stdlib-types.html#channel-e)
""")
public interface Channel<E> {

    @Operator
    @Description("""
        The `collect` operator collects all values from a source channel into a list and emits it as a single value.

        [Read more](https://nextflow.io/docs/latest/reference/operator-typed.html#collect)
    """)
    Value<Bag<E>> collect();

    @Operator
    @Description("""
        The `combine` operator emits every pairwise combination of two channels.

        [Read more](https://nextflow.io/docs/latest/reference/operator-typed.html#combine)
    """)
    Channel<Tuple> combine(Channel right);
    Channel<Tuple> combine(Value right);

    @Operator
    @Description("""
        When the `combine` operator is called on a channel of records with named arguments, the named arguments are appended to each record in the source channel. Each named argument can be a value or dataflow value.

        [Read more](https://nextflow.io/docs/latest/reference/operator-typed.html#combine)
    """)
    Channel<Record> combine(Map<String,?> fields);

    @Operator
    @Description("""
        The `filter` operator emits the values from a source channel that satisfy a condition, discarding all other values.

        [Read more](https://nextflow.io/docs/latest/reference/operator-typed.html#filter)
    """)
    Channel<E> filter(Predicate<E> condition);

    @Operator
    @Description("""
        The `flatMap` operator applies a mapping function to each value from a source channel. The mapping function should return a collection, and each element in the collection is emitted separately.

        [Read more](https://nextflow.io/docs/latest/reference/operator-typed.html#flatmap)
    """)
    <R> Channel<R> flatMap(Function<E,Iterable<R>> transform);
    <R> Channel<R> flatMap();

    @Operator
    @Description("""
        The `groupBy` operator collects values from a source channel into groups based on a grouping key. A tuple is emitted for each group containing the grouping key and collection of values.

        [Read more](https://nextflow.io/docs/latest/reference/operator-typed.html#groupby)
    """)
    Channel<Tuple> groupBy();

    @Operator
    @Description("""
        The `join` operator emits the relational join of two source channels using a matching key.

        [Read more](https://nextflow.io/docs/latest/reference/operator-typed.html#join)
    """)
    Channel<Record> join(
        @NamedParams({
            @NamedParam(value = "by", type = String.class),
            @NamedParam(value = "remainder", type = Boolean.class)
        })
        Map<String,?> opts,
        Channel<Record> right
    );

    @Operator
    @Description("""
        The `map` operator applies a mapping function to each value from a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator-typed.html#map)
    """)
    <R> Channel<R> map(Function<E,R> transform);

    @Operator
    @Description("""
        The `mix` operator emits the values from two source channels into a single output channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator-typed.html#mix)
    """)
    Channel<E> mix(Channel<E> other);
    Channel<E> mix(Value<E> other);

    @Operator
    @Description("""
        The `reduce` operator applies an accumulator function sequentially to each value from a source channel, and emits the accumulated value.

        The accumulator function takes two parameters -- the accumulated value and the *i*-th emitted value -- and it should return the accumulated result, which is passed to the next invocation with the *i+1*-th value. This process is repeated for each value in the source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator-typed.html#reduce)
    """)
    <R> Value<R> reduce(R seed, BiFunction<R,E,R> accumulator);
    Value<E> reduce(BiFunction<E,E,E> accumulator);

    @Operator
    @Description("""
        The `subscribe` operator performs an action for each value in a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator-typed.html#subscribe)
    """)
    void subscribe(Consumer<E> action);
    void subscribe(Map<String,Consumer<E>> events);

    @Operator
    @Description("""
        The `unique` operator emits the unique values from a source channel.

        [Read more](https://nextflow.io/docs/latest/reference/operator-typed.html#unique)
    """)
    Channel<E> unique(Function<E,?> transform);
    Channel<E> unique();

    @Operator
    @Description("""
        The `until` operator emits each value from a source channel until a stopping condition is satisfied.

        [Read more](https://nextflow.io/docs/latest/reference/operator-typed.html#until)
    """)
    Channel<E> until(Predicate<E> condition);

    @Operator
    @Description("""
        The `view` operator prints each value from a source channel to standard output.

        [Read more](https://nextflow.io/docs/latest/reference/operator-typed.html#view)
    """)
    Channel<E> view(Map<String,?> opts, Function<E,String> transform);
    Channel<E> view(Function<E,String> transform);
    Channel<E> view(Map<String,?> opts);
    Channel<E> view();

    // legacy operators

    @Deprecated @Operator Record branch(Closure closure);
    @Deprecated @Operator Channel collectFile(Map<String,?> opts, Closure closure);
    @Deprecated @Operator Channel collectFile(Closure closure);
    @Deprecated @Operator Value count();
    @Deprecated @Operator Channel cross(Channel right);
    @Deprecated @Operator Channel dump(Map<String,?> opts);
    @Deprecated @Operator Channel dump();
    @Deprecated @Operator Channel flatten();
    @Deprecated @Operator Channel groupTuple(Map<String,?> opts);
    @Deprecated @Operator Channel ifEmpty(Object value);
    @Deprecated @Operator Value<E> max(Closure comparator);
    @Deprecated @Operator Value<E> min(Closure comparator);
    @Deprecated @Operator Record multiMap(Closure closure);
    @Deprecated @Operator Channel splitCsv(Map<String,?> opts);
    @Deprecated @Operator Channel splitCsv();
    @Deprecated @Operator Channel splitText(Map<String,?> opts);
    @Deprecated @Operator Channel splitText();
    @Deprecated @Operator Value sum(Closure transform);
    @Deprecated @Operator Channel take(int n);
    @Deprecated @Operator Value toList();
    @Deprecated @Operator Value toSortedList();
    @Deprecated @Operator Channel transpose(Map<String,?> opts);
    @Deprecated @Operator Channel transpose();

}

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
import java.util.function.Consumer;
import java.util.function.Function;

import nextflow.script.dsl.Description;
import nextflow.script.dsl.Operator;

@Description("""
    A dataflow value is an asynchronous value. It is used to facilitate dataflow logic in a workflow.

    [Read more](https://nextflow.io/docs/latest/reference/stdlib-types.html#value-v)
""")
public interface Value<V> {

    @Operator
    @Description("""
        The `combine` operator combines two dataflow values.
    """)
    Value<Tuple> combine(Value right);

    @Operator
    @Description("""
        When the `combine` operator is called with named arguments and the dataflow value is a record, the named arguments are appended to the record. Each named argument can be a value or dataflow value.
    """)
    Value<Record> combine(Map<String,?> fields);

    @Operator
    @Description("""
        Transforms the dataflow value into a collection with the given closure and emits the resulting values in a dataflow channel.
    """)
    <R> Channel<R> flatMap(Function<V,Iterable<R>> transform);
    <R> Channel<R> flatMap();

    @Operator
    @Description("""
        Transforms the dataflow value into another dataflow value with the given closure.
    """)
    <R> Value<R> map(Function<V,R> transform);

    @Operator
    @Description("""
        The `mix` operator emits the values from two source channels into a single output channel.
    """)
    Channel<V> mix(Channel<V> other);

    @Operator
    @Description("""
        Invokes the given closure on the dataflow value.
    """)
    void subscribe(Consumer<V> action);

    @Operator
    @Description("""
        Transforms the dataflow value using the given closure and print the result to standard output.
    """)
    Value<V> view(Map<String,?> opts, Function<V,String> transform);
    Value<V> view(Function<V,String> transform);
    Value<V> view(Map<String,?> opts);
    Value<V> view();

}

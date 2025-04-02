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
package nextflow.script.types.shim;

import java.util.function.BiFunction;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;

import nextflow.script.dsl.Description;

@Description("""
    An iterable is a trait shared by collections that support iteration.

    [Read more](https://nextflow.io/docs/latest/reference/stdlib.html#iterable-e)
""")
public interface Iterable<E> {

    @Description("""
        Returns `true` if any element in the iterable satisfies the given condition.
    """)
    boolean any(Predicate<E> predicate);

    @Description("""
        Returns a new iterable with each element transformed by the given closure.
    """)
    <R> Iterable<R> collect(Function<E,R> transform);

    @Description("""
        Transforms each element in the iterable into a collection with the given closure and concatenates the resulting collections into a list.
    """)
    <R> Iterable<R> collectMany(Function<E,Iterable<R>> transform);

    @Description("""
        Returns `true` if the iterable contains the given value.
    """)
    boolean contains(E value);

    @Description("""
        Invoke the given closure for each element in the iterable.
    """)
    void each(Consumer<E> action);

    @Description("""
        Returns `true` if every element in the iterable satisfies the given condition.
    """)
    boolean every(Predicate<E> predicate);

    @Description("""
        Returns the elements in the iterable that satisfy the given condition.
    """)
    Iterable<E> findAll(Predicate<E> predicate);

    @Description("""
        Collect the elements of an iterable into groups based on a matching key. The closure should return the key for a given element.
    """)
    <K> Map<K,Iterable<E>> groupBy(Function<E,K> transform);

    @Description("""
        Apply the given accumulator to each element in the iterable and return the final accumulated value. The closure should accept two parameters, corresponding to the current accumulated value and the current iterable element, and return the next accumulated value. The first element from the iterable is used as the initial accumulated value.
    """)
    E inject(BiFunction<E,E,E> accumulator);

    @Description("""
        Apply the given accumulator to each element in the iterable and return the final accumulated value. The closure should accept two parameters, corresponding to the current accumulated value and the current iterable element, and return the next accumulated value.
    """)
    <R> R inject(R initialValue, BiFunction<R,E,R> accumulator);

    @Description("""
        Returns `true` if the iterable is empty.
    """)
    boolean isEmpty();

    @Description("""
        Concatenates the string representation of each element in the iterable.
    """)
    String join();

    @Description("""
        Concatenates the string representation of each element in the iterable, with the given string as the separator between each element.
    """)
    String join(String separator);

    @Description("""
        Returns the maximum element in the iterable.
    """)
    E max();

    @Description("""
        Returns the maximum element in the iterable according to the given closure. The closure should follow the same semantics as the closure parameter of `toSorted()`.
    """)
    <R> E max(Function<E,R> comparator);

    @Description("""
        Returns the maximum element in the iterable.
    """)
    E min();

    @Description("""
        Returns the maximum element in the iterable according to the given closure. The closure should follow the same semantics as the closure parameter of `toSorted()`.
    """)
    <R> E min(Function<E,R>  comparator);

    @Description("""
        Returns the number of elements in the iterable.
    """)
    int size();

    @Description("""
        Returns the sum of the elements in the iterable. The elements should support addition (`+`).
    """)
    E sum();

    @Description("""
        Transforms each element in the iterable with the given closure and returns the sum. The values returned by the closure should support addition (`+`).
    """)
    <R> R sum(Function<E,R> transform);

    @Description("""
        Converts the iterable to a set. Duplicate elements are excluded.
    """)
    Set<E> toSet();

    @Description("""
        Returns a sorted list of the iterable's elements.
    """)
    List<E> toSorted();

    @Description("""
        Returns the iterable as a list sorted according to the given closure. The closure should accept one parameter and transform each element into the value that will be used for comparisons.
    """)
    <R> List<E> toSorted(Function<E,R> comparator);

    @Description("""
        Returns a shallow copy of the iterable with duplicate elements excluded.
    """)
    Iterable<E> toUnique();

    @Description("""
        Returns a shallow copy of the iterable with duplicate elements excluded. The closure should follow the same semantics as the closure parameter of `toSorted()`.
    """)
    <R> Iterable<E> toUnique(Function<E,R> comparator);

}

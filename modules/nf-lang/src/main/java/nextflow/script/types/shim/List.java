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

import java.util.function.Predicate;

import groovy.lang.Tuple2;
import nextflow.script.dsl.Description;

@Description("""
    A list is an unordered collection of elements.

    [Read more](https://nextflow.io/docs/latest/reference/stdlib.html#list-e)
""")
@ShimType(java.util.List.class)
public interface List<E> extends Iterable<E> {

    @Description("""
        Collates the list into a list of sub-lists of length `size`. If `keepRemainder` is `true`, any remaining elements are included as a partial sub-list, otherwise they are excluded.
    """)
    List<List<E>> collate(int size, boolean keepRemainder);

    @Description("""
        Collates the list into a list of sub-lists of length `size`, stepping through the list `step` elements for each sub-list. If `keepRemainder` is `true`, any remaining elements are included as a partial sub-list, otherwise they are excluded.
    """)
    List<List<E>> collate(int size, int step, boolean keepRemainder);

    @Description("""
        Returns the first value in the list that satisfies the given condition.
    """)
    E find(Predicate<E> condition);

    @Description("""
        Returns the first element in the list. Raises an error if the list is empty.
    """)
    E first();

    @Description("""
        Returns the list of integers from 0 to *n - 1*, where *n* is the number of elements in the list.
    """)
    List<Integer> getIndices();

    @Description("""
        Equivalent to `first()`.
    """)
    E head();

    @Description("""
        Returns the index of the first occurrence of the given value in the list, or -1 if the list does not contain the value.
    """)
    int indexOf(E value);

    @Description("""
        Returns a shallow copy of the list with the last element excluded.
    """)
    List<E> init();

    @Description("""
        Returns the last element in the list. Raises an error if the list is empty.
    """)
    E last();

    @Description("""
        Returns a shallow copy of the list with the elements reversed.
    """)
    List<E> reverse();

    @Description("""
        Returns the portion of the list between the given `fromIndex` (inclusive) and `toIndex` (exclusive).
    """)
    List<E> subList(int fromIndex, int toIndex);

    @Description("""
        Returns a shallow copy of the list with the first element excluded.
    """)
    List<E> tail();

    @Description("""
        Returns the first *n* elements of the list.
    """)
    List<E> take(int n);

    @Description("""
        Returns the longest prefix of the list where each element satisfies the given condition.
    """)
    List<E> takeWhile(Predicate<E> condition);

    @Description("""
        Returns a list of 2-tuples corresponding to the value and index of each element in the list.
    """)
    List<Tuple2<E,Integer>> withIndex();

}

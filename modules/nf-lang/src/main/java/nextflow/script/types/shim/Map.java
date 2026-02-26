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

import java.util.function.BiConsumer;
import java.util.function.BiPredicate;

import nextflow.script.dsl.Description;

@Description("""
    A map associates or "maps" keys to values. Each key can map to at most one value -- a map cannot contain duplicate keys.

    [Read more](https://nextflow.io/docs/latest/reference/stdlib.html#map-k-v)
""")
@ShimType(java.util.Map.class)
public interface Map<K,V> {

    @Description("""
        Returns `true` if any key-value pair in the map satisfies the given condition. The closure should accept two parameters corresponding to the key and value of an entry.
    """)
    boolean any(BiPredicate<K,V> condition);

    @Description("""
        Returns `true` if the map contains a mapping for the given key.
    """)
    boolean containsKey(K key);

    @Description("""
        Returns `true` if the map maps one or more keys to the given value.
    """)
    boolean containsValue(V value);

    @Description("""
        Invoke the given closure for each key-value pair in the map. The closure should accept two parameters corresponding to the key and value of an entry.
    """)
    void each(BiConsumer<K,V> action);

    @Description("""
        Returns a set of the key-value pairs in the map.
    """)
    Set<Entry<K,V>> entrySet();

    @Description("""
        Returns `true` if every key-value pair in the map satisfies the given condition. The closure should accept two parameters corresponding to the key and value of an entry.
    """)
    boolean every(BiPredicate<K,V> condition);

    @Description("""
        Returns `true` if the map is empty.
    """)
    boolean isEmpty();

    @Description("""
        Returns a set of the keys in the map.
    """)
    Set<K> keySet();

    @Description("""
        Returns the number of key-value pairs in the map.
    """)
    int size();

    @Description("""
        Returns a sub-map containing the given keys.
    """)
    Map<K,V> subMap(Iterable<K> keys);

    @Description("""
        Returns a collection of the values in the map.
    """)
    Bag<V> values();

    @Description("""
        A map entry is a key-value pair.
    """)
    @ShimType(java.util.Map.Entry.class)
    interface Entry<K,V> {
        K getKey();
        V getValue();
    }

}

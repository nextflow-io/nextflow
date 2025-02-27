/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.util

import groovy.transform.CompileStatic

/**
 * A {@link Map} that handles keys in a case insensitive manner
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class InsensitiveMap<K,V> implements Map<K,V> {

    @Delegate
    private Map<K,V> target

    private InsensitiveMap(Map<K,V> map) {
        this.target = map
    }

    @Override
    boolean containsKey(Object key) {
        target.any( it -> key?.toString()?.toLowerCase() == it.key?.toString()?.toLowerCase())
    }

    @Override
    V get(Object key) {
        target.find(it -> key?.toString()?.toLowerCase() == it.key?.toString()?.toLowerCase())?.value
    }

    static <K,V> Map<K,V> of(Map<K,V> target) {
        new InsensitiveMap(target)
    }
}

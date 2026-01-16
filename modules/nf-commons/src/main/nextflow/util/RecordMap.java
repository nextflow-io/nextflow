/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.util;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import nextflow.script.types.Record;
import org.codehaus.groovy.runtime.DefaultGroovyMethods;

/**
 * Implements Record as an immutable map.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class RecordMap extends HashMap<String,Object> implements Record {

    public RecordMap() {}

    public RecordMap(Map<String,Object> props) {
        super(props);
    }

    public RecordMap(int initialCapacity) {
        super(initialCapacity);
    }

    @Override
    public void clear() {
        throw new UnsupportedOperationException();
    }

    // @Override
    // public Object put(String key, Object value) {
    //     throw new UnsupportedOperationException();
    // }

    // @Override
    // public void putAll(Map<? extends String, ? extends Object> m) {
    //     throw new UnsupportedOperationException();
    // }

    @Override
    public Object remove(Object key) {
        throw new UnsupportedOperationException();
    }

    public Record subMap(Collection<String> keys) {
        return new RecordMap(DefaultGroovyMethods.subMap(this, keys));
    }

}

/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.util
import groovy.transform.CompileStatic
/**
 * A class for which the initial values cannot be changed
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ReadOnlyMap implements Map {

    /** The list of field names which value cannot change */
    private List<String> readOnlyNames

    /** The target map holding the values */
    @Delegate
    private Map target = new LinkedHashMap()

    /**
     * Create a map with values cannot be modified
     *
     * @param map Initial map values
     */
    ReadOnlyMap( Map map ) {
        assert map != null
        this.target.putAll(map)
        this.readOnlyNames = new ArrayList(map.keySet())
    }

    /**
     * @param map Initial map values
     * @param readOnly The list of values that cannot be changed after obj construction,
     */
    ReadOnlyMap( Map map, List<String> readOnly  ) {
        assert map != null
        assert readOnly != null
        this.target.putAll(map)
        this.readOnlyNames = new ArrayList(readOnly)
    }

    /**
     * Set a map entry only if it's not declared as read-only
     *
     * @param name The entry key
     * @param value The entry new value
     * @return the previous value associated with <tt>key</tt>, or
     *         <tt>null</tt> if there was no mapping for <tt>key</tt>.
     */
    def put(Object name, Object value) {
        final read_only = readOnlyNames.contains(name)
        if ( !read_only ) {
            target.put(name,value)
        }
    }

    /**
     * Overrides a map value even if it's declared as read-only
     *
     * @param name The entry key
     * @param value The entry new value
     * @return the previous value associated with <tt>key</tt>, or
     *         <tt>null</tt> if there was no mapping for <tt>key</tt>.
     */
    def force( Object name, Object value ) {
        target.put(name,value)
    }

}

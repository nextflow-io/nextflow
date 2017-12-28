/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.util
import groovy.transform.CompileStatic
/**
 * A class for which the initial values cannot be changed
 *
 * @author Paolo Di Tommaso
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

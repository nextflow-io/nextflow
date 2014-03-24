/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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


/**
 * A class for which the initial values cannot be changed
 *
 * @author Paolo Di Tommaso
 */
class ReadOnlyMap extends LinkedHashMap {

    List<String> readOnlyNames

    /**
     * When no values are specified it behave as a normal map
     */
    ReadOnlyMap() {

    }

    /**
     * @param map Initial map values
     * @param readOnlyNames1 The list of values that cannot be changed after obj construction,
     *          when {@code null} all the value in the initial map cannot be changed
     *
     */
    ReadOnlyMap( Map map, List<String> readOnlyNames = null  ) {
        super(map)
        this.readOnlyNames = new ArrayList( readOnlyNames != null ? readOnlyNames : map.keySet() )
    }

    def put(Object name, Object value) {
        final read_only = readOnlyNames.contains(name)
        if ( !read_only ) {
            super.put(name,value)
        }
    }

    def force( Object name, Object value ) {
        super.put(name,value)
    }

}

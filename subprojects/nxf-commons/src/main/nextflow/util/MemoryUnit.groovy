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

import java.text.DecimalFormat
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode

/**
 * Represent a memory unit
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@EqualsAndHashCode(includes = 'size')
class MemoryUnit implements Comparable<MemoryUnit>, Serializable, Cloneable {

    final static public MemoryUnit ZERO = new MemoryUnit(0)

    final static private Pattern FORMAT = ~/([0-9\.]+)\s*(\S)?B?/

    final static public List UNITS = [ "B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB" ]

    final long size

    /**
     * Default constructor is required by Kryo serializer
     * Do not remove of use directly
     */
    private MemoryUnit() { this.size=0 }

    /**
     * Create a memory unit instance
     *
     * @param value The number of bytes it represent
     */
    MemoryUnit( long value ) {
        assert value>=0, "Memory unit cannot be a negative number"
        this.size = value
    }

    /**
     * Create a memory unit instance with the given semantic string
     *
     * @param str A string using the following of of the following units: B, KB, MB, GB, TB, PB, EB, ZB
     */
    MemoryUnit( String str ) {

        def matcher = FORMAT.matcher(str)
        if( !matcher.matches() ) {
            throw new IllegalArgumentException("Not a valid FileSize value: '$str'")
        }

        final value = matcher.group(1)
        final unit = matcher.group(2)?.toUpperCase()

        if ( !unit || unit == "B" ) {
            size = Long.parseLong(value)
        }
        else {
            int p = UNITS.indexOf(unit)
            if ( p == -1 ) {
                // try adding a 'B' specified
                p = UNITS.indexOf(unit+'B')
                if( p == -1 ) {
                    throw new IllegalArgumentException("Not a valid file size unit: ${str}")
                }
            }

            size = Math.round( Double.parseDouble(value) * Math.pow(1024, p) )
        }

    }

    long toBytes() {
        size
    }

    long toKilo() {
        size >> 10
    }

    long toMega() {
        size >> 20
    }

    long toGiga() {
        size >> 30
    }

    MemoryUnit plus( MemoryUnit value )  {
        return value != null ? new MemoryUnit( size + value.size ) : this
    }

    MemoryUnit minus( MemoryUnit value )  {
        return value != null ? new MemoryUnit( size - value.size ) : this
    }

    MemoryUnit multiply( Number value ) {
        return new MemoryUnit( (long)(size * value) )
    }

    MemoryUnit div( Number value ) {
        return new MemoryUnit( Math.round((double)(size / value)) )
    }

    String toString() {
        if(size <= 0) {
            return "0"
        }

        // see http://stackoverflow.com/questions/2510434/format-bytes-to-kilobytes-megabytes-gigabytes
        int digitGroups = (int) (Math.log10(size) / Math.log10(1024))
        new DecimalFormat("#,##0.#").format(size / Math.pow(1024, digitGroups)) + " " + UNITS[digitGroups]
    }

    @Override
    int compareTo(MemoryUnit that) {
        return this.size <=> that.size
    }

    static MemoryUnit of( String value ) {
        new MemoryUnit(value)
    }

    static MemoryUnit of( long value ) {
        new MemoryUnit(value)
    }

    boolean asBoolean() {
        return size != 0
    }
}

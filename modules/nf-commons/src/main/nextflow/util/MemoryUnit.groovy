/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import java.text.DecimalFormat
import java.text.DecimalFormatSymbols
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
/**
 * Represent a memory unit
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@EqualsAndHashCode(includes = 'size', includeFields = true)
class MemoryUnit implements Comparable<MemoryUnit>, Serializable, Cloneable {

    final static public MemoryUnit ZERO = new MemoryUnit(0)

    final static private Pattern FORMAT = ~/([0-9\.]+)\s*(\S)?B?/

    final static public List UNITS = [ "B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB" ]

    private long size

    static final private DecimalFormatSymbols formatSymbols

    static {
        formatSymbols = new DecimalFormatSymbols()
        formatSymbols.setDecimalSeparator('.' as char)
    }

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

    long getBytes() { size }

    long toKilo() { size >> 10 }

    long getKilo() { size >> 10 }

    long toMega() { size >> 20 }

    long getMega() { size >> 20 }

    long toGiga() { size >> 30 }

    long getGiga() { size >> 30 }

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
        final formatter = new DecimalFormat("0.#", formatSymbols)
        formatter.setGroupingUsed(false)
        formatter.format(size / Math.pow(1024, digitGroups)) + " " + UNITS[digitGroups]
    }

    @Override
    int compareTo(MemoryUnit that) {
        return this.size <=> that.size
    }

    static int compareTo(MemoryUnit left, Object right) {
        assert left

        if( right==null )
            throw new IllegalArgumentException("Not a valid memory value: null")

        if( right instanceof MemoryUnit )
            return left <=> (MemoryUnit)right

        if( right instanceof Number )
            return left.size <=> right.toLong()

        if( right instanceof CharSequence )
            return left <=> MemoryUnit.of(right.toString())

        throw new IllegalArgumentException("Not a valid memory value: $right")
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

    /**
     * Function to parse/convert given memory unit
     *
     * @param unit String expressing memory unit in bytes, e.g. KB, MB, GB
     */
    long toUnit(String unit){
        int p = UNITS.indexOf(unit)
        if (p==-1)
            throw new IllegalArgumentException("Not a valid memory unit: $unit")
        return size.intdiv(Math.round(Math.pow(1024,p)))
    }
}

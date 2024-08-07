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
 */

package nextflow.util;

import java.io.Serializable;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import org.codehaus.groovy.runtime.DefaultGroovyMethods;
/**
 * Represent a memory unit
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class MemoryUnit implements Comparable<MemoryUnit>, Serializable, Cloneable {

    public static final List<String> UNITS = List.of("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB");

    public static final MemoryUnit ZERO = new MemoryUnit(0);

    private static final Pattern FORMAT = Pattern.compile("([0-9\\.]+)\\s*(\\S)?B?");

    private static final DecimalFormatSymbols formatSymbols = new DecimalFormatSymbols();

    static {
        formatSymbols.setDecimalSeparator('.');
    }

    private long size;

    /**
     * Default constructor is required by Kryo serializer
     * Do not remove of use directly
     */
    private MemoryUnit() {
        this.size = 0;
    }

    /**
     * Create a memory unit instance
     *
     * @param value The number of bytes it represent
     */
    public MemoryUnit(long value) {
        if( value < 0 )
            throw new IllegalArgumentException("Memory unit cannot be a negative number");
        this.size = value;
    }

    /**
     * Create a memory unit instance with the given semantic string
     *
     * @param str A string using the following of of the following units: B, KB, MB, GB, TB, PB, EB, ZB
     */
    public MemoryUnit(String str) {
        var matcher = FORMAT.matcher(str);
        if( !matcher.matches() ) {
            throw new IllegalArgumentException("Not a valid FileSize value: '" + str + "'");
        }

        var value = matcher.group(1);
        var unit = matcher.group(2);
        if( unit != null )
            unit = unit.toUpperCase();

        if ( !DefaultGroovyMethods.asBoolean(unit) || "B".equals(unit) ) {
            this.size = Long.parseLong(value);
        }
        else {
            int p = UNITS.indexOf(unit);
            if( p == -1 ) {
                // try adding a 'B' specified
                p = UNITS.indexOf(unit + "B");
                if( p == -1 ) {
                    throw new IllegalArgumentException("Not a valid file size unit: " + str);
                }
            }

            this.size = Math.round( Double.parseDouble(value) * Math.pow(1024, p) );
        }
    }

    public long toBytes() {
        return size;
    }

    public long getBytes() {
        return size;
    }

    public long toKilo() {
        return size >> 10;
    }

    public long getKilo() {
        return size >> 10;
    }

    public long toMega() {
        return size >> 20;
    }

    public long getMega() {
        return size >> 20;
    }

    public long toGiga() {
        return size >> 30;
    }

    public long getGiga() {
        return size >> 30;
    }

    public MemoryUnit plus(MemoryUnit value) {
        return value != null ? new MemoryUnit( size + value.size ) : this;
    }

    public MemoryUnit minus(MemoryUnit value) {
        return value != null ? new MemoryUnit( size - value.size ) : this;
    }

    public MemoryUnit multiply(Number value) {
        return new MemoryUnit( (long)(size * value.doubleValue()) );
    }

    public MemoryUnit div(Number value) {
        return new MemoryUnit( Math.round((double)size / value.doubleValue()) );
    }

    @Override
    public String toString() {
        if( size <= 0 ) {
            return "0";
        }

        // see http://stackoverflow.com/questions/2510434/format-bytes-to-kilobytes-megabytes-gigabytes
        var digitGroups = (int)(Math.log10(size) / Math.log10(1024));
        var formatter = new DecimalFormat("0.#", formatSymbols);
        formatter.setGroupingUsed(false);
        return formatter.format(size / Math.pow(1024, digitGroups)) + " " + UNITS.get(digitGroups);
    }

    @Override
    public int compareTo(MemoryUnit that) {
        return Long.compare(this.size, that.size);
    }

    public static int compareTo(MemoryUnit left, Object right) {
        assert left != null;

        if( right == null )
            throw new IllegalArgumentException("Not a valid memory value: null");

        if( right instanceof MemoryUnit )
            return left.compareTo((MemoryUnit) right);

        if( right instanceof Number )
            return Long.compare(left.size, ((Number) right).longValue());

        if( right instanceof CharSequence )
            return left.compareTo(MemoryUnit.of(right.toString()));

        throw new IllegalArgumentException("Not a valid memory value: " + right.toString());
    }

    public static MemoryUnit of(String value) {
        return new MemoryUnit(value);
    }

    public static MemoryUnit of(long value) {
        return new MemoryUnit(value);
    }

    public boolean asBoolean() {
        return size != 0;
    }

    /**
     * Function to parse/convert given memory unit
     *
     * @param unit String expressing memory unit in bytes, e.g. KB, MB, GB
     */
    public long toUnit(String unit) {
        int p = UNITS.indexOf(unit);
        if( p == -1 )
            throw new IllegalArgumentException("Not a valid memory unit: " + unit);
        return size / Math.round(Math.pow(1024, p));
    }

    @Override
    public boolean equals(Object that) {
        return compareTo(this, that) == 0;
    }

    @Override
    public int hashCode() {
        return Long.hashCode(size);
    }
}

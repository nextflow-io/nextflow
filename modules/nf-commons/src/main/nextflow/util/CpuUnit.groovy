/*
 * Copyright 2020-2021, Seqera Labs
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

import java.math.RoundingMode
import java.text.DecimalFormat
import java.text.DecimalFormatSymbols

import com.google.common.math.IntMath
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode

/**
 * Model a CPU requirement measured in cpu units. 1 CPU unit is equivalent to 1 physical CPU core, or 1 virtual core, depending
 * on whether the node is a physical host or a virtual machine running inside a physical machine.
 *
 * The following formats are allow:
 * - integer number: represents a full CPU unit eg. `1` => 1 CPUs
 * - decimal number: represents a percentage of the CPU 0.1 => 10% CPU or 100m, etc
 * - number with `m` suffix: represents milli cpus notation ie 100m == 0.1 CPU, 1000m == 1 CPU, etc
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@EqualsAndHashCode(includes = 'millis', includeFields = true)
class CpuUnit implements Comparable<CpuUnit>, Serializable, Cloneable {

    static public CpuUnit ONE_CORE = new CpuUnit(1_000)

    static private final DecimalFormat DECIMAL_FMT
    static {
        final formatSymbols = new DecimalFormatSymbols()
        formatSymbols.setDecimalSeparator('.' as char)
        DECIMAL_FMT = new DecimalFormat("##0.0", formatSymbols)
        DECIMAL_FMT.setRoundingMode(RoundingMode.CEILING)
    }

    final private int millis

    @Override
    int compareTo(CpuUnit that) {
        return this.millis <=> that.millis
    }

    static int compareTo(CpuUnit left, Object right) {
        assert left

        if( right==null )
            throw new IllegalArgumentException("Not a valid cpu value: null")

        if( right instanceof CpuUnit )
            return left <=> (CpuUnit)right

        if( right instanceof Number )
            return left <=> CpuUnit.of(right)

        if( right instanceof CharSequence )
            return left <=> CpuUnit.of(right)

        throw new IllegalArgumentException("Not a valid cpu value: $right")
    }

    static CpuUnit of(CharSequence value) {
        return new CpuUnit(parseCpuMillis(value) )
    }

    static CpuUnit of(Number value) {
        return new CpuUnit(parseCpuMillis(value) )
    }

    static CpuUnit of(Object value) {
        if( value instanceof Number )
            return new CpuUnit(parseCpuMillis((Number)value) )
        if( value instanceof CharSequence )
            return new CpuUnit(parseCpuMillis((CharSequence)value) )
        if( value==null )
            throw new IllegalArgumentException("Not a valid cpus value: null")
        else
            throw new IllegalArgumentException("Not a valid cpus value: $value [${value.getClass().getName()}]")
    }

    private CpuUnit(int value) {
        this.millis = value
    }

    /*
     * Default constructor is required by Kryo serializer
     * Do not remove or use directly
     */
    private CpuUnit() { this.millis=0 }

    static protected int parseCpuMillis(Number value) {
        if( value instanceof Integer )
            return value.toInteger() * 1_000
        if( value instanceof Number )
            return Math.ceil(value.toDouble() * 1_000)
        throw new IllegalArgumentException("Not a valid cpus value: '$value'")
    }

    static protected int parseCpuMillis(CharSequence value) {
        // turn millis into an integer unit
        final str = value.toString()
        if( str.endsWith('m') ) {
            try {
                return str.substring(0,str.length()-1).toInteger()
            }
            catch (NumberFormatException e) {
                throw new IllegalArgumentException("Not a valid cpus value: '$value'", e)
            }
        }
        // turn float to unit
        if( str.contains('.') ) {
            try {
                return Math.ceil(str.toDouble() * 1_000)
            }
            catch (NumberFormatException e) {
                throw new IllegalArgumentException("Not a valid cpus value: '$value'", e)
            }
        }
        // it should be an integer value
        try {
            return str.toInteger() * 1_000
        }
        catch (NumberFormatException e) {
            throw new IllegalArgumentException("Not a valid cpus value: '$value'", e)
        }
    }

    int toCores() {
        return IntMath.divide(millis, 1_000, RoundingMode.CEILING);
    }

    BigDecimal toDecimal() {
        return millis / 1_000
    }

    int toMillis() {
        return millis;
    }

    String toDecimalString() {
        return DECIMAL_FMT.format(millis / 1_000)
    }

    String toString() {
        return millis + 'm'
    }
}

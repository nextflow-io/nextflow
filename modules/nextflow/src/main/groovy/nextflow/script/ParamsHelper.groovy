/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.script

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.exception.ScriptRuntimeException
import nextflow.script.dsl.Types
import nextflow.script.types.Record
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import nextflow.util.RecordMap
import nextflow.util.TypeHelper
import nextflow.util.VersionNumber
import org.codehaus.groovy.runtime.typehandling.GroovyCastException
/**
 * Resolves a declared param against a given value.
 *
 * Used by pipeline params, typed workflows, and typed processes, so
 * that they all map a given value to a declared type the same way.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ParamsHelper {

    /**
     * Resolve a value given on the command line. Command-line values are
     * always strings, so they are parsed according to the declared type.
     *
     * @param decl
     * @param value
     */
    static Object resolveFromCli(Param decl, Object value) {
        if( value == null )
            return null

        if( value instanceof Collection || value instanceof Map )
            return asType(value, decl)

        final number = asNumberType(decl, value)
        if( number != null )
            return number

        if( value !instanceof CharSequence )
            return value

        final str = value.toString()

        if( decl.type == Boolean ) {
            if( str.toLowerCase() == 'true' ) return Boolean.TRUE
            if( str.toLowerCase() == 'false' ) return Boolean.FALSE
        }

        return resolveFromString(decl, str, value)
    }

    /**
     * Resolve a value given in a params file or the config. Such values are
     * already structured, so they only need to be converted where the
     * declared type is more specific than the source syntax.
     *
     * @param decl
     * @param value
     */
    static Object resolveFromCode(Param decl, Object value) {
        if( value == null )
            return null

        if( value instanceof Collection || value instanceof Map )
            return asType(value, decl)

        final number = asNumberType(decl, value)
        if( number != null )
            return number

        if( value !instanceof CharSequence )
            return value

        return resolveFromString(decl, value.toString(), value)
    }

    /**
     * Convert a value to a declared numeric type. Integer and Float are the
     * only numeric types in the Nextflow type system, but a value can be
     * given as any number (e.g. a Double or BigDecimal in the config), and a
     * command-line value is always a string, so the value is normalized to
     * the declared type.
     *
     * The value is converted to the narrowest representation that can hold
     * it -- Integer, Long, or BigInteger for an integral value, Float,
     * Double, or BigDecimal for a fractional one -- so that no precision is
     * lost for a value that is too large for the declared type.
     *
     * Returns null if the declared type is not numeric, or if the value
     * cannot be represented by it -- an Integer accepts an integral value in
     * any notation (e.g. `3`, `3.0`, `3e2`) but rejects one with a fractional
     * part, rather than silently truncating it. The value is then reported by
     * the caller as not assignable to the declared type.
     *
     * @param decl
     * @param value
     */
    private static Number asNumberType(Param decl, Object value) {
        if( decl.type != Integer && decl.type != Float )
            return null
        try {
            final number = new BigDecimal(value.toString().trim())
            return decl.type == Float
                ? asFloatType(number)
                : asIntegerType(number)
        }
        catch( NumberFormatException | ArithmeticException e ) {
            return null
        }
    }

    /**
     * Convert a number to a Float, widening it to a Double or BigDecimal
     * if it is too large to be represented by a Float.
     *
     * @param number
     */
    private static Number asFloatType(BigDecimal number) {
        final floatValue = number.floatValue()
        if( !floatValue.isInfinite() )
            return floatValue
        final doubleValue = number.doubleValue()
        return !doubleValue.isInfinite()
            ? doubleValue
            : number
    }

    /**
     * Convert a number to an Integer, widening it to a Long or BigInteger if
     * it is too large to be represented by an Integer. A fractional value is
     * rejected rather than being truncated.
     *
     * @param number
     */
    private static Number asIntegerType(BigDecimal number) {
        final integer = number.toBigIntegerExact()
        if( integer.bitLength() < Integer.SIZE )
            return integer.intValue()
        return integer.bitLength() < Long.SIZE
            ? integer.longValue()
            : integer
    }

    /**
     * Convert a string to a declared type that is always expressed as a
     * string, regardless of where the value came from. Returns the given
     * fallback if the declared type is not one of these.
     *
     * @param decl
     * @param str
     * @param fallback
     */
    private static Object resolveFromString(Param decl, String str, Object fallback) {
        if( decl.type == Path )
            return TypeHelper.asPathType(str)

        if( decl.type == Duration )
            return Duration.of(str)

        if( decl.type == MemoryUnit )
            return MemoryUnit.of(str)

        if( decl.type == VersionNumber )
            return new VersionNumber(str)

        return fallback
    }

    /**
     * Convert a composite value (a collection, map, or record) to the
     * declared type, reporting a conversion failure in terms of the param as
     * well as the underlying element -- a NumberFormatException from an
     * element of a `List<Integer>`, for example, names neither the param nor
     * the declared type on its own.
     *
     * @param value
     * @param decl
     */
    private static Object asType(Object value, Param decl) {
        try {
            return TypeHelper.asType(value, decl.type)
        }
        catch( GroovyCastException | UnsupportedOperationException | IllegalArgumentException e ) {
            final actualType = value.getClass()
            final detail = e.message ? " -- ${e.message}" : ''
            throw new ScriptRuntimeException("Parameter `${decl.name}` with type ${Types.getName(decl.type)} cannot be assigned to ${value} [${Types.getName(actualType)}]${detail}")
        }
    }

    static boolean isAssignableFrom(Class target, Class source) {
        // any numeric value can be assigned to Float
        if( target == Float.class )
            return Number.class.isAssignableFrom(source)

        // any integer value can be assigned to Integer
        if( target == Integer.class )
            return source == BigInteger.class || source == Long.class || source == Integer.class

        // any record can be assigned to a record type (validation is handled by asType())
        if( Record.class.isAssignableFrom(target) )
            return source == RecordMap.class

        return target.isAssignableFrom(source)
    }

}

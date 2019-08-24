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

import java.time.temporal.Temporal
import java.util.concurrent.TimeUnit

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.util.logging.Slf4j
import org.apache.commons.lang.time.DurationFormatUtils
/**
 * A simple time duration representation
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@EqualsAndHashCode(includes = 'durationInMillis')
class Duration implements Comparable<Duration>, Serializable, Cloneable {

    static private final FORMAT = ~/^(\d+\.?\d*)\s*([a-zA-Z]+)/

    static private final LEGACY = ~/^(\d{1,2}):(\d{1,2}):(\d{1,2})$/

    static private final List<String> MILLIS = ['ms','milli','millis']

    static private final List<String> SECONDS = ['s','sec','second','seconds']

    static private final List<String> MINUTES = ['m','min','minute','minutes']

    static private final List<String> HOURS = ['h','hour','hours']

    static private final List<String> DAYS = ['d','day','days']

    static public final List<String> UNITS

    static {
        UNITS = []
        UNITS.addAll(MILLIS)
        UNITS.addAll(SECONDS)
        UNITS.addAll(MINUTES)
        UNITS.addAll(HOURS)
        UNITS.addAll(DAYS)
    }

    /**
     * Duration in millis
     */
    final long durationInMillis

    /**
     * Create e a duration object having the specified number of millis
     *
     * @param duration The duration as milliseconds
     */
    Duration(long duration) {
        assert duration>=0, "Duration unit cannot be a negative number"
        this.durationInMillis = duration
    }


    /**
     * Default constructor is required by Kryo serializer
     * Do not removed or use it directly
     */
    private Duration() { durationInMillis=0 }

    /**
     * Create the object using a string 'duration' format.
     * Accepted prefix are:
     * <li>{@code ms}, {@code milli}, {@code millis}: for milliseconds
     * <li>{@code s}, {@code second}, {@code seconds}: for seconds
     * <li>{@code m}, {@code minute}, {@code minutes}: for minutes
     * <li>{@code h}, {@code hour}, {@code hours}: for hours
     * <li>{@code d}, {@code day}, {@code days}: for days
     *
     *
     * @param str
     */
    Duration(String str) {

        try {
            try {
                durationInMillis = parseSimple(str)
            }
            catch( IllegalArgumentException e ) {
                durationInMillis = parseLegacy(str)
            }
        }
        catch( IllegalArgumentException e ) {
            throw e
        }
        catch( Exception e ) {
            throw new IllegalArgumentException("Not a valid duration value: ${str}", e)
        }
    }

    /**
     * Parse a duration string in legacy format i.e. hh:mm:ss
     *
     * @param str The string to be parsed e.g. {@code 05:10:30} (5 hours, 10 min, 30 seconds)
     * @return The duration in millisecond
     */
    private long parseLegacy( String str ) {
        def matcher = (str =~ LEGACY)
        if( !matcher.matches() )
            new IllegalArgumentException("Not a valid duration value: ${str}")

        def groups = (List<String>)matcher[0]
        def hh = groups[1].toInteger()
        def mm = groups[2].toInteger()
        def ss = groups[3].toInteger()

        return TimeUnit.HOURS.toMillis(hh) + TimeUnit.MINUTES.toMillis(mm) + TimeUnit.SECONDS.toMillis(ss)
    }

    /**
     * Parse a duration string
     *
     * @param str A duration string containing one or more component e.g. {@code 1d 3h 10mins}
     * @return  The duration in millisecond
     */
    private long parseSimple( String str ) {

        long result=0
        for( int i=0; true; i++ ) {
            def matcher = (str =~ FORMAT)
            if( matcher.find() ) {
                def groups = (List<String>)matcher[0]
                def all = groups[0]
                def digit = groups[1]
                def unit = groups[2]

                result += convert( digit.toFloat(), unit )
                str = str.substring(all.length()).trim()
                continue
            }

            if( i == 0 || str )
                throw new IllegalArgumentException("Not a valid duration value: ${str}")
            break
        }

        return result
    }

    /**
     * Parse a single duration component
     *
     * @param digit
     * @param unit A valid duration unit e.g. {@code d}, {@code d}, {@code h}, {@code hour}, etc
     * @return The duration in millisecond
     */
    private long convert( float digit, String unit ) {

        if( unit in MILLIS ) {
            return Math.round(digit)
        }
        if ( unit in SECONDS ) {
            return Math.round(digit * 1_000)
        }
        if ( unit in MINUTES ) {
            return Math.round(digit * 60 * 1_000)
        }
        if ( unit in HOURS ) {
            return Math.round(digit * 60 * 60 * 1_000)
        }
        if ( unit in DAYS ) {
            return Math.round(digit * 24 * 60 * 60 * 1_000)
        }

        throw new IllegalStateException()
    }

    Duration(long value, TimeUnit unit) {
        assert value>=0, "Duration unit cannot be a negative number"
        assert unit, "Time unit cannot be null"
        this.durationInMillis = unit.toMillis(value)
    }

    static Duration of( long value ) {
        new Duration(value)
    }

    static Duration of( String str ) {
        new Duration(str)
    }

    static Duration of( String str, Duration fallback ) {
        try {
            return new Duration(str)
        }
        catch( IllegalArgumentException e ) {
            log.debug "Not a valid duration value: $str -- Fallback on default value: $fallback"
            return fallback
        }
    }

    static Duration between( Temporal start, Temporal end ) {
        new Duration(java.time.Duration.between(start, end).toMillis())
    }

    long toMillis() {
        durationInMillis
    }

    long getMillis() {
        durationInMillis
    }

    long toSeconds() {
        TimeUnit.MILLISECONDS.toSeconds(durationInMillis)
    }

    long getSeconds() {
        toSeconds()
    }

    long toMinutes() {
        TimeUnit.MILLISECONDS.toMinutes(durationInMillis)
    }

    long getMinutes() {
        toMinutes()
    }

    long toHours() {
        TimeUnit.MILLISECONDS.toHours(durationInMillis)
    }

    long getHours() {
        toHours()
    }

    long toDays() {
        TimeUnit.MILLISECONDS.toDays(durationInMillis)
    }

    long getDays() {
        toDays()
    }

    /**
     * Duration formatting utilities and constants. The following table describes the tokens used in the pattern language for formatting.
     * <p>
     * <pre>
     *   character	duration element
     *   y	        years
     *   d	        days
     *   H	        hours
     *   m	        minutes
     *   s	        seconds
     * </pre>
     *
     * @param fmt
     * @return
     */
    String format( String fmt ) {
        DurationFormatUtils.formatDuration(durationInMillis, fmt)
    }

    String toString() {

        // just prints the milliseconds
        if( durationInMillis < 1_000 ) {
            return durationInMillis + 'ms'
        }

        // when less than 60 seconds round up to 100th of millis
        if( durationInMillis < 60_000 ) {
            return String.valueOf( Math.round(durationInMillis / 1_000 * 10 as float) / 10 ) + 's'
        }

        def secs
        def mins
        def hours
        def days
        def result = []

        // round up to seconds
        secs = Math.round( (double)(durationInMillis / 1_000) )

        mins = secs.intdiv(60)
        secs = secs % 60
        if( secs )
            result.add( secs+'s' )

        hours = mins.intdiv(60)
        mins = mins % 60
        if( mins )
            result.add(0, mins+'m' )

        days = hours.intdiv(24)
        hours = hours % 24
        if( hours )
            result.add(0, hours+'h' )

        if( days )
            result.add(0, days+'d')

        return result.join(' ')
    }

    def plus( Duration value )  {
        return new Duration( durationInMillis + value.durationInMillis )
    }

    def minus( Duration value )  {
        return new Duration( durationInMillis - value.durationInMillis )
    }

    def multiply( Number value ) {
        return new Duration( (long)(durationInMillis * value) )
    }

    def div( Number value ) {
        return new Duration( Math.round((double)(durationInMillis / value)) )
    }

    boolean asBoolean() {
        return durationInMillis != 0
    }

    @Override
    int compareTo(Duration that) {
        return this.durationInMillis <=> that.durationInMillis
    }

    static int compareTo(Duration left, Object right) {
        assert left

        if( right==null )
            throw new IllegalArgumentException("Not a valid duration value: null")

        if( right instanceof Duration )
            return left <=> (Duration)right

        if( right instanceof Number )
            return left.durationInMillis <=> right.toLong()

        if( right instanceof CharSequence )
            return left <=> Duration.of(right.toString())

        throw new IllegalArgumentException("Not a valid duration value: $right")
    }

}

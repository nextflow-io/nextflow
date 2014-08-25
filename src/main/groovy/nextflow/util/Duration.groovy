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

import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.TimeUnit

import groovy.transform.Canonical
import groovy.transform.EqualsAndHashCode
import groovy.util.logging.Slf4j
import org.apache.commons.lang.time.DurationFormatUtils
/**
 * A simple time duration representation
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@EqualsAndHashCode(includes = 'durationInMillis')
class Duration implements Comparable<Duration> {

    static private final FORMAT = /(\d+)\s*(\S+)/

    static private final MILLIS = ['ms','milli','millis']

    static private final SECONDS = ['s','sec','second','seconds']

    static private final MINUTES = ['m','min','minute','minutes']

    static private final HOURS = ['h','hour','hours']

    static private final DAYS = ['d','day','days']

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
        assert duration>=0
        this.durationInMillis = duration
    }

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
        def matcher = (str =~~ FORMAT)
        if( !matcher.matches() ) {
            throw new IllegalArgumentException("Not a valid duration value: '$str'")
        }

        final val = Long.parseLong(matcher[0][1]?.toString())
        final unit = matcher[0][2]

        if( unit in MILLIS ) {
            this.durationInMillis = val
        }
        else if ( unit in SECONDS ) {
            this.durationInMillis = TimeUnit.SECONDS.toMillis(val)
        }
        else if ( unit in MINUTES ) {
            this.durationInMillis = TimeUnit.MINUTES.toMillis(val)
        }
        else if ( unit in HOURS ) {
            this.durationInMillis = TimeUnit.HOURS.toMillis(val)
        }
        else if ( unit in DAYS ) {
            this.durationInMillis = TimeUnit.DAYS.toMillis(val)
        }
        else {
            throw new IllegalArgumentException("Not a valid duration value: ${str}")
        }

    }

    Duration(long value0, TimeUnit unit) {
        assert unit
        this.durationInMillis = unit.toMillis(value0)
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

    long toMillis() {
        durationInMillis
    }

    long toSeconds() {
        TimeUnit.MILLISECONDS.toSeconds(durationInMillis)
    }

    long toMinutes() {
        TimeUnit.MILLISECONDS.toMinutes(durationInMillis)
    }

    long toHours() {
        TimeUnit.MILLISECONDS.toHours(durationInMillis)
    }

    long toDays() {
        TimeUnit.MILLISECONDS.toDays(durationInMillis)
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

        if( durationInMillis < 1000 ) {
            return durationInMillis + MILLIS[0]
        }

        def value = format("d:H:m:s").split(':').collect { Integer.parseInt(it) }
        def result = []

        // -- day / days
        if( value[0] >= 1 ) {
            result << value[0] + DAYS[0]
        }

        // hour / hours
        if( value[1] >= 1 ) {
            result << value[1] + HOURS[0]
        }

        // -- minute / minutes
        if( value[2] > 0 ) {
            result << value[2] + MINUTES[0]
        }

        // -- second / seconds
        if( value[3] > 0 ) {
            result << value[3] + SECONDS[0]
        }

        result.join(' ')
    }


    @Override
    int compareTo(Duration that) {
        return this.durationInMillis <=> that.durationInMillis
    }

    @Canonical
    static class ThrottleObj {
        Object result
        long timestamp
    }

    def throttle( Closure closure ) {
        throttle0( durationInMillis, null, closure)
    }

    def throttle( seed, Closure closure ) {
        def initialValue = new ThrottleObj( seed, System.currentTimeMillis() )
        throttle0( durationInMillis, initialValue, closure)
    }

    static final Map<Integer,ThrottleObj> throttleMap = new ConcurrentHashMap<>()

    private static throttle0( long delayMillis, ThrottleObj initialValue, Closure closure ) {
        assert closure != null

        def key = 17
        key  = 31 * key + closure.class.hashCode()
        key  = 31 * key + closure.owner.hashCode()
        key  = 31 * key + closure.delegate?.hashCode() ?: 0

        ThrottleObj obj = throttleMap.get(key)
        if( obj == null ) {
            obj = initialValue ?: new ThrottleObj()
            throttleMap.put(key,obj)
        }

        if( System.currentTimeMillis() - obj.timestamp > delayMillis ) {
            obj.timestamp = System.currentTimeMillis()
            obj.result = closure.call()
        }

        obj.result
    }

}

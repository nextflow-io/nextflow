/*
 * Copyright (c) 2012, the authors.
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

import java.util.concurrent.TimeUnit
import java.util.concurrent.TimeoutException

import org.apache.commons.lang.time.DurationFormatUtils

/**
 * A simple time duration representation
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class Duration {

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

    static Duration create( String str ) {
        new Duration(str)
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
        def value = format("d:H:m:s").split(':').collect { Integer.parseInt(it) }
        def result = []

        // -- day / days
        if( value[0] == 1 ) {
            result << value[0] + 'day'
        }
        else if( value[0] > 1 ) {
            result << value[0] + 'days'
        }

        // hour / hours
        if( value[1] == 1 ) {
            result << value[1] + 'hour'
        }
        else if( value[1] > 1 ) {
            result << value[1] + 'hour'
        }

        // -- minute / minutes
        if( value[2] > 0 ) {
            result << value[2] + 'min'
        }

        // -- second / seconds
        if( value[3] > 0 ) {
            result << value[3] + 'sec'
        }

        result.join(' ')
    }

    /**
     * Waits indefinitely the specified condition is verified
     *
     * @param condition A closure returning a condition to be satisfied, if waits until the condition is {@code false}
     *          and exit when {@code true}
     */
    static final protected  waitFor( Closure<Boolean> condition ) {
        waitFor((Duration)null, condition)
    }


    /**
     * Waits until the specified condition is verified or throws an exception of type {@code TimeoutException}
     * when the specified duration is exceeded
     *
     * @param maxDuration
     * @param condition
     */
    static final waitFor( Duration maxDuration, Closure<Boolean> condition ) {
        long time = System.currentTimeMillis()

        while( true ) {
            if( condition.call() ) {
                break
            }

            if( maxDuration && (System.currentTimeMillis()-time > maxDuration.durationInMillis ) ) {
                throw new TimeoutException("Wait condition exceed time-out: $maxDuration")
            }

            sleep(200)
        }
    }

    static final waitFor( String duration, Closure<Boolean> condition ) {
        waitFor(new Duration(duration), condition)
    }

    static final waitFor( long millis, Closure<Boolean> condition ) {
        waitFor(new Duration(millis), condition)
    }

}

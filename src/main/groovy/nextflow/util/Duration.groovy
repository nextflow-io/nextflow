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

/**
 * A simple time duration representation
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class Duration {

    static private FORMAT = /(\d+)\s*(\S+)/

    static private MILLIS = ['ms','milli','millis']

    static private SECONDS = ['s','sec','second','seconds']

    static private MINUTES = ['m','min','minute','minutes']

    static private HOURS = ['h','hour','hours']

    static private DAYS = ['d','day','days']


    /**
     * Duration in millis
     */
    long value

    /**
     * Create e a duration object having the specified number of millis
     *
     * @param duration The duration as milliseconds
     */
    Duration(long duration) {
        assert duration>=0
        this.value = duration
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
     *
     *
     * For example<p>:
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
            this.value = val
        }
        else if ( unit in SECONDS ) {
            this.value = TimeUnit.SECONDS.toMillis(val)
        }
        else if ( unit in MINUTES ) {
            this.value = TimeUnit.MINUTES.toMillis(val)
        }
        else if ( unit in HOURS ) {
            this.value = TimeUnit.HOURS.toMillis(val)
        }
        else if ( unit in DAYS ) {
            this.value = TimeUnit.DAYS.toMillis(val)
        }
        else {
            throw new IllegalArgumentException("Not a valid duration value: ${str}")
        }

    }

    Duration(long value0, TimeUnit unit) {
        assert unit
        this.value = unit.toMillis(value0)
    }

    static Duration create( String str ) {
        new Duration(str)
    }

    long toMillis() {
        value
    }

    long toSeconds() {
        TimeUnit.MILLISECONDS.toSeconds(value)
    }

    long toMinutes() {
        TimeUnit.MILLISECONDS.toMinutes(value)
    }

    long toHours() {
        TimeUnit.MILLISECONDS.toHours(value)
    }

    long toDays() {
        TimeUnit.MILLISECONDS.toDays(value)
    }

    String toString() {

    }

}

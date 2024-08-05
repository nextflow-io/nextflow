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
import java.time.temporal.Temporal;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;

import groovy.transform.EqualsAndHashCode;
import org.codehaus.groovy.runtime.DefaultGroovyMethods;
import org.codehaus.groovy.runtime.StringGroovyMethods;
import org.apache.commons.lang.time.DurationFormatUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
/**
 * A simple time duration representation
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@EqualsAndHashCode(includes = "durationInMillis")
public class Duration implements Comparable<Duration>, Serializable, Cloneable {

    private static final Logger log = LoggerFactory.getLogger(Duration.class);

    private static final Pattern FORMAT = Pattern.compile("^(\\d+\\.?\\d*)\\s*([a-zA-Z]+)");

    private static final Pattern LEGACY = Pattern.compile("^(\\d{1,2}):(\\d{1,2}):(\\d{1,2})$");

    private static final List<String> MILLIS = List.of("ms","milli","millis");

    private static final List<String> SECONDS = List.of("s","sec","second","seconds");

    private static final List<String> MINUTES = List.of("m","min","minute","minutes");

    private static final List<String> HOURS = List.of("h","hour","hours");

    private static final List<String> DAYS = List.of("d","day","days");

    public static final List<String> UNITS = new ArrayList<>();

    static {
        UNITS.addAll(MILLIS);
        UNITS.addAll(SECONDS);
        UNITS.addAll(MINUTES);
        UNITS.addAll(HOURS);
        UNITS.addAll(DAYS);
    }

    /**
     * Duration in millis
     */
    private final long durationInMillis;

    /**
     * Create e a duration object having the specified number of millis
     *
     * @param duration The duration as milliseconds
     */
    public Duration(long duration) {
        if( duration < 0 )
            throw new IllegalArgumentException("Duration unit cannot be a negative number");
        this.durationInMillis = duration;
    }

    /**
     * Default constructor is required by Kryo serializer
     * Do not removed or use it directly
     */
    private Duration() {
        this.durationInMillis = 0;
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
     * @param str
     */
    public Duration(String str) {
        try {
            long millis;
                millis = parseSimple(str);
            // try {
            // }
            // catch( IllegalArgumentException e ) {
            //     millis = parseLegacy(str);
            // }
            this.durationInMillis = millis;
        }
        catch( IllegalArgumentException e ) {
            throw e;
        }
        catch( Exception e ) {
            throw new IllegalArgumentException("Not a valid duration value: " + str, e);
        }
    }

    /**
     * Parse a duration string in legacy format i.e. hh:mm:ss
     *
     * @param str The string to be parsed e.g. {@code 05:10:30} (5 hours, 10 min, 30 seconds)
     * @return The duration in milliseconds
     */
    private long parseLegacy(String str) {
        var matcher = LEGACY.matcher(str);
        if( !matcher.matches() )
            new IllegalArgumentException("Not a valid duration value: " + str);

        var hh = StringGroovyMethods.toInteger(matcher.group(1));
        var mm = StringGroovyMethods.toInteger(matcher.group(2));
        var ss = StringGroovyMethods.toInteger(matcher.group(3));

        return TimeUnit.HOURS.toMillis(hh) + TimeUnit.MINUTES.toMillis(mm) + TimeUnit.SECONDS.toMillis(ss);
    }

    /**
     * Parse a duration string
     *
     * @param str A duration string containing one or more component e.g. {@code 1d 3h 10mins}
     * @return The duration in milliseconds
     */
    private long parseSimple(String str) {
        long result = 0;
        for( int i = 0; true; i++ ) {
            var matcher = FORMAT.matcher(str);
            if( matcher.find() ) {
                var all = matcher.group();
                var digit = matcher.group(1);
                var unit = matcher.group(2);

                result += convert(StringGroovyMethods.toDouble(digit), unit);
                str = str.substring(all.length()).trim();
                continue;
            }

            if( i == 0 || !str.isEmpty() )
                throw new IllegalArgumentException("Not a valid duration value: " + str);
            break;
        }

        return result;
    }

    /**
     * Parse a single duration component
     *
     * @param digit
     * @param unit A valid duration unit e.g. {@code d}, {@code d}, {@code h}, {@code hour}, etc
     * @return The duration in milliseconds
     */
    private long convert(double digit, String unit) {

        if( MILLIS.contains(unit) ) {
            return Math.round(digit);
        }
        if( SECONDS.contains(unit) ) {
            return Math.round(digit * 1_000);
        }
        if( MINUTES.contains(unit) ) {
            return Math.round(digit * 60 * 1_000);
        }
        if( HOURS.contains(unit) ) {
            return Math.round(digit * 60 * 60 * 1_000);
        }
        if( DAYS.contains(unit) ) {
            return Math.round(digit * 24 * 60 * 60 * 1_000);
        }

        throw new IllegalStateException();
    }

    public Duration(long value, TimeUnit unit) {
        if( value < 0 )
            throw new IllegalArgumentException("Duration unit cannot be a negative number");
        if( unit == null )
            throw new IllegalArgumentException("Time unit cannot be null");
        this.durationInMillis = unit.toMillis(value);
    }

    public static Duration of(long value) {
        return new Duration(value);
    }

    public static Duration of(String str) {
        return new Duration(str);
    }

    public static Duration of(String str, Duration fallback) {
        try {
            return new Duration(str);
        }
        catch( IllegalArgumentException e ) {
            log.debug("Not a valid duration value: " + str + " -- Fallback on default value: " + fallback.toString());
            return fallback;
        }
    }

    public static Duration between(Temporal start, Temporal end) {
        return new Duration(java.time.Duration.between(start, end).toMillis());
    }

    public long toMillis() {
        return durationInMillis;
    }

    public long getMillis() {
        return durationInMillis;
    }

    public long toSeconds() {
        return TimeUnit.MILLISECONDS.toSeconds(durationInMillis);
    }

    public long getSeconds() {
        return toSeconds();
    }

    public long toMinutes() {
        return TimeUnit.MILLISECONDS.toMinutes(durationInMillis);
    }

    public long getMinutes() {
        return toMinutes();
    }

    public long toHours() {
        return TimeUnit.MILLISECONDS.toHours(durationInMillis);
    }

    public long getHours() {
        return toHours();
    }

    public long toDays() {
        return TimeUnit.MILLISECONDS.toDays(durationInMillis);
    }

    public long getDays() {
        return toDays();
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
     */
    public String format(String fmt) {
        return DurationFormatUtils.formatDuration(durationInMillis, fmt);
    }

    public String toString() {

        // when less than 1 second just print the milliseconds
        if( durationInMillis < 1_000 ) {
            return String.valueOf( durationInMillis ) + "ms";
        }

        // when less than 60 seconds round to 100th of seconds
        if( durationInMillis < 60_000 ) {
            return String.valueOf( Math.round((double)durationInMillis / 1_000.0 * 10.0) / 10 ) + "s";
        }

        long secs;
        long mins;
        long hours;
        long days;
        var result = new ArrayList<String>();

        // round to seconds
        secs = Math.round( (double)durationInMillis / 1_000 );

        mins = secs / 60;
        secs = secs % 60;
        if( secs > 0 )
            result.add( String.valueOf(secs) + "s" );

        hours = mins / 60;
        mins = mins % 60;
        if( mins > 0 )
            result.add(0, String.valueOf(mins) + "m");

        days = hours / 24;
        hours = hours % 24;
        if( hours > 0 )
            result.add(0, String.valueOf(hours) + "h");

        if( days > 0 )
            result.add(0, String.valueOf(days) + "d");

        return DefaultGroovyMethods.join(result, " ");
    }

    public Duration plus(Duration value)  {
        return new Duration(durationInMillis + value.durationInMillis);
    }

    public Duration minus(Duration value)  {
        return new Duration(durationInMillis - value.durationInMillis);
    }

    public Duration multiply(Number value) {
        return new Duration(durationInMillis * value.longValue());
    }

    public Duration div(Number value) {
        return new Duration(Math.round((double)durationInMillis / value.longValue()));
    }

    public boolean asBoolean() {
        return durationInMillis != 0;
    }

    @Override
    public int compareTo(Duration that) {
        return Long.compare(this.durationInMillis, that.durationInMillis);
    }

    public static int compareTo(Duration left, Object right) {
        assert left != null;

        if( right == null )
            throw new IllegalArgumentException("Not a valid duration value: null");

        if( right instanceof Duration )
            return left.compareTo((Duration) right);

        if( right instanceof Number )
            return Long.compare(left.durationInMillis, ((Number) right).longValue());

        if( right instanceof CharSequence )
            return left.compareTo(Duration.of(right.toString()));

        throw new IllegalArgumentException("Not a valid duration value: " + right.toString());
    }

}

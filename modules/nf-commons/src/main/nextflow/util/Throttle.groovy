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

import java.util.concurrent.ConcurrentHashMap

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode

/**
 * Limit the execution of a code block on a specified time basis
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class Throttle {

    static final Map<Integer,ThrottleObj> throttleMap = new ConcurrentHashMap<>()

    @EqualsAndHashCode
    static class ThrottleObj {
        Object result
        long timestamp

        ThrottleObj() {}

        ThrottleObj( value, long timestamp ) {
            this.result = value
            this.timestamp = timestamp
        }
    }

    /**
     * Executes the specified code block and returns the resulting value. All following
     * invocations within the specified time {@code period} are skipped and the cached
     * result value is returned instead.
     *
     * @param period
     *          The time period in milliseconds, in which the throttler will allow at most to invoke
     *          the specified closure
     * @param closure
     *          The code block to execute
     * @return
     *          The value returned by the closure execution
     */
    static Object every( long period, Closure closure ) {
        throttle0( period, null, closure)
    }

    /**
     * Executes the specified code block and returns the resulting value. All following
     * invocations within the specified time {@code period} are skipped and the cached
     * result value is returned instead.
     *
     * @param period
     *          The time duration, in which the throttler will allow at most to invoke
     *          the specified closure
     * @param closure
     *          The code block to execute
     * @return
     *          The value returned by the closure execution
     */
    static Object every( Duration period, Closure closure ) {
        throttle0( period.millis, null, closure)
    }

    /**
     * Executes the specified code block and returns the resulting value. All following
     * invocations within the specified time {@code period} are skipped and the cached
     * result value is returned instead.
     *
     * @param period
     *          The duration string, in which the throttler will allow at most to invoke
     *          the specified closure
     * @param closure
     *          The code block to execute
     * @return
     *          The value returned by the closure execution
     */
    static Object every( String period, Closure closure ) {
        throttle0( Duration.of(period).millis, null, closure)
    }

    /**
     * Allows at most one execution after the specified time {@code delay} returning
     * {@code null} before the first execution, and the cached result value the
     * following times
     *
     * @param delay
     * @param closure
     * @return
     */
    static Object after( long delay, Closure closure ) {
        final obj = new ThrottleObj( null, System.currentTimeMillis() )
        throttle0( delay, obj, closure)
    }

    /**
     * Allows at most one execution after the specified time {@code delay} returning
     * {@code null} before the first execution, and the cached result value the
     * following times
     *
     * @param delay
     * @param closure
     * @return
     */
    static Object after( Duration delay, Closure closure ) {
        def obj = new ThrottleObj( null, System.currentTimeMillis() )
        throttle0( delay.millis, obj, closure)
    }

    /**
     * Allows at most one execution after the specified time {@code delay} returning
     * {@code null} before the first execution, and the cached result value the
     * following times
     *
     * @param delay
     * @param closure
     * @return
     */
    static Object after( String delay, Closure closure ) {
        def obj = new ThrottleObj( null, System.currentTimeMillis() )
        throttle0( Duration.of(delay).millis, obj, closure)
    }

    /**
     * Allows at most one execution after the specified time {@code delay} returning
     * {@code initialValue} before the first execution, and the cached result value the
     * following times
     *
     * @param delay
     * @param initialValue
     * @param closure
     * @return
     */
    static Object after( long delay, initialValue, Closure closure ) {
        final obj = new ThrottleObj( initialValue, System.currentTimeMillis() )
        throttle0( delay, obj, closure)
    }

    /**
     * Allows at most one execution after the specified time {@code delay} returning
     * {@code initialValue} before the first execution, and the cached result value the
     * following times
     *
     * @param delay
     * @param initialValue
     * @param closure
     * @return
     */
    static Object after( Duration delay, initialValue, Closure closure ) {
        def obj = new ThrottleObj( initialValue, System.currentTimeMillis() )
        throttle0( delay.millis, obj, closure)
    }

    /**
     * Allows at most one execution after the specified time {@code delay} returning
     * {@code initialValue} before the first execution, and the cached result value the
     * following times
     *
     * @param delay
     * @param initialValue
     * @param closure
     * @return
     */
    static Object after( String delay, initialValue, Closure closure ) {
        def obj = new ThrottleObj( initialValue, System.currentTimeMillis() )
        throttle0( Duration.of(delay).millis, obj, closure)
    }

    private static throttle0( long timeout, ThrottleObj initialValue, Closure closure ) {
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

        if( System.currentTimeMillis() - obj.timestamp > timeout ) {
            obj.timestamp = System.currentTimeMillis()
            obj.result = closure.call()
        }

        obj.result
    }

    static <V> V cache( Object key, Duration eviction, Closure<V> action ) {
        cache(key, eviction.toMillis(), action)
    }

    static <V> V cache( Object key, long eviction, Closure<V> closure ) {
        final int hash = key?.hashCode() ?: 0

        ThrottleObj obj = throttleMap.get(hash)
        if( obj == null ) {
            obj = new ThrottleObj()
            obj.result = closure.call()
            obj.timestamp = System.currentTimeMillis()
            throttleMap.put(hash,obj)
        }

        else if( System.currentTimeMillis() - obj.timestamp > eviction ) {
            obj.timestamp = System.currentTimeMillis()
            obj.result = closure.call()
        }

        (V)obj.result
    }

}

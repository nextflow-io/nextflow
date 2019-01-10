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

package nextflow.processor

import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic
/**
 * TaskRun unique identifier
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TaskId extends Number implements Comparable, Serializable, Cloneable {

    /**
     * Global count of all task instances
     */
    static final private AtomicInteger allCount = new AtomicInteger()

    static TaskId next() {
        new TaskId(allCount.incrementAndGet())
    }

    private final int value

    static TaskId of( value ) {
        if( value instanceof Integer )
            return new TaskId(value)
        else if( value != null )
            return new TaskId( Integer.valueOf(value.toString()))
        else
            throw new IllegalArgumentException("TaskId cannot be null")
    }

    TaskId( int id ) {
        value = id
    }

    TaskId( TaskId id ) {
        value = id.value
    }

    /**
     * Default constructor is only required to enable class serialisation
     * by Kryo library. Do not remove.
     */
    protected TaskId() { }

    @Override
    boolean equals( Object obj ) {
        if( obj instanceof TaskId ) {
            return value == obj.value
        }
        if( obj instanceof Integer ) {
            return value == obj
        }
        return false
    }

    @Override
    int hashCode() {
        value
    }

    @Override
    int compareTo( Object other ) {
        if( other instanceof TaskId ) {
            return Integer.compare(value,other.value)
        }
        if( other instanceof Integer  ) {
            return Integer.compare(value,other)
        }
        if( other == null ) {
            return 1
        }

        throw new IllegalArgumentException("Not a valid TaskId: $other [${other.getClass().getName()}]")
    }


    @Override
    String toString() { String.valueOf(value) }

    int intValue() { value }

    long longValue() { (long)value }

    float floatValue() { (float)value }

    double doubleValue() { (double)value }

 //    private static long toUnsignedLong(int x) {
//        return ((long) x) & 0xffffffffL
//    }
//
//    private static int toUnsignedInt(int x) {
//        return (x) & 0xffffffff
//    }
//
//    private static String toUnsignedString(int i) {
//        return Long.toString(toUnsignedLong(i));
//    }



}

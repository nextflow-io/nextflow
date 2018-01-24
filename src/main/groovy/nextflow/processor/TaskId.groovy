/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

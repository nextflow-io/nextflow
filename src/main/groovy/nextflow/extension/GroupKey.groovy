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

package nextflow.extension

import org.codehaus.groovy.runtime.InvokerHelper

/**
 * Helper class used to wrap a generic key object and to attach it
 * a size attribute to implement a dynamic `groupTuple` size rule
 *
 * See {@link nextflow.Nextflow#groupKey(java.lang.Object)}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GroupKey {

    private final Object target

    private final int size

    /** This constructor is needed by the Kryo deserializer */
    private GroupKey() { }

    /**
     * Create a grouping key object
     *
     * @param key
     * @param size
     */
    GroupKey(key, int size) {
        this.target = key
        this.size = size
    }

    int getGroupSize() { size }

    /**
     * Delegate any method invocation to the target key object
     *
     * @param name The method name
     * @param args The method arguments (if any)
     * @return The invoked method result (if any)
     */
    def methodMissing( String name, def args ) {
        InvokerHelper.invokeMethod(target, name, args)
    }

    /**
     * Delegate any property invocation to the target object
     *
     * @param name The property name
     * @return The resulting property value
     */
    def propertyMissing(String name) {
        InvokerHelper.getProperty(target,name)
    }

    @Override
    boolean equals(Object obj) {
        obj instanceof GroupKey ? target.equals(obj.target) : target.equals(obj)
    }

    @Override
    int hashCode() {
        target.hashCode()
    }

    @Override
    String toString() {
        target.toString()
    }
}
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

package nextflow.extension

import com.google.common.hash.Hasher
import nextflow.util.CacheFunnel
import nextflow.util.CacheHelper
import nextflow.util.CacheHelper.HashMode
import org.codehaus.groovy.runtime.InvokerHelper

/**
 * Helper class used to wrap a generic key object and to attach it
 * a size attribute to implement a dynamic `groupTuple` size rule
 *
 * See {@link nextflow.Nextflow#groupKey(java.lang.Object)}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GroupKey implements CacheFunnel {

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

    @Override
    Hasher funnel(Hasher hasher, HashMode mode) {
        return CacheHelper.hasher(hasher, target, mode)
    }
}

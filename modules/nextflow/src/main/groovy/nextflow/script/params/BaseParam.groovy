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

package nextflow.script.params

import groovy.util.logging.Slf4j

/**
 * Base class for input/output parameters
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class BaseParam implements Cloneable {

    final protected Binding binding

    final protected List<BaseParam> holder

    final short index

    final short mapIndex

    private boolean initialized

    BaseParam ( Binding binding, List holder, int ownerIndex = -1 ) {
        this.binding = binding
        this.holder = holder

        /*
         * by default the index is got from 'holder' current size
         * and the mapIndex is =1 (not defined)
         */
        if( ownerIndex == -1 ) {
            index = holder.size()
            mapIndex = -1
        }

        /*
         * when the owner index is provided (not -1) it is used as
         * the main index and the map index is got from the 'holder' size
         */
        else {
            index = ownerIndex
            mapIndex = holder.size()
        }

        // add the the param to the holder list
        holder.add(this)
    }

    String toString() {
        def p = mapIndex == -1 ? index : "$index:$mapIndex"
        return "${this.class.simpleName.toLowerCase()}<$p>"
    }

    /**
     * Lazy initializer
     */
    protected abstract void lazyInit()

    /**
     * Initialize the parameter fields if needed
     */
    final protected void init() {
        if( initialized ) return
        lazyInit()

        // flag as initialized
        initialized = true
    }


    /**
     * Get the value of variable {@code name} in the script context
     *
     * @param name The variable name
     * @param strict If {@code true} raises a {@code MissingPropertyException} when the specified variable does not exist
     * @return The variable object
     */
    protected getScriptVar(String name, boolean strict ) {
        if( binding.hasVariable(name) ) {
            return binding.getVariable(name)
        }

        if( strict )
            throw new MissingPropertyException(name,this.class)

        return null
    }

    protected getScriptVar( String name ) {
        getScriptVar(name,true)
    }

}

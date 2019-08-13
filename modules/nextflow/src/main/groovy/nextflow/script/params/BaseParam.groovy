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

    /**
     * The binding context to resolve param variables
     */
    final protected Binding binding

    protected List<BaseParam> holder

    /**
     * The param declaration index in the input/output block
     * Note the index do not change for nested parameters ie. declared in the same tuple param
     */
    final short index

    /**
     * The nested index of tuple composed parameters or -1 when it's a top level param ie. not a tuple element
     */
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

    @Override
    Object clone() {
        final copy = (BaseParam)super.clone()
        copy.holder = this.holder!=null ? new ArrayList<BaseParam>(holder) : new ArrayList<BaseParam>()
        return copy
    }

    String toString() {
        def p = mapIndex == -1 ? index : "$index:$mapIndex"
        return "${getTypeSimpleName()}<$p>"
    }

    String getTypeSimpleName() {
        this.class.simpleName.toLowerCase()
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

    protected BaseParam setOptions(Map<String,?> opts) {
        if( !opts )
            return this

        for( Map.Entry<String,?> entry : opts ) {
            setProperty(entry.key, entry.value)
        }
        return this
    }

    boolean isNestedParam() {
        return mapIndex >= 0
    }

}

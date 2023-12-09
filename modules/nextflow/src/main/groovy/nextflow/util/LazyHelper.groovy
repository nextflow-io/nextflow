/*
 * Copyright 2013-2023, Seqera Labs
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

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
/**
 * Helper methods for lazy binding and resolution.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class LazyHelper {

    /**
     * Resolve a lazy value against a given binding.
     *
     * @param binding
     * @param value
     */
    static Object resolve(Object binding, Object value) {
        if( value instanceof LazyAware )
            return value.resolve(binding)

        if( value instanceof Closure )
            return value.cloneWith(binding).call()

        if( value instanceof GString )
            return value.cloneAsLazy(binding).toString()

        return value
    }

}

/**
 * Interface for types that can be lazily resolved
 */
interface LazyAware {
    Object resolve(Object binding)
}

/**
 * A list that can be lazily resolved
 */
@CompileStatic
class LazyList implements LazyAware, List {

    @Delegate
    private List target

    LazyList() {
        target = []
    }

    LazyList(int size) {
        target = new ArrayList(size)
    }

    LazyList(Collection items) {
        target = new ArrayList(items)
    }

    @Override
    Object resolve(Object binding) {
        final result = new ArrayList(target.size())
        for( def item : target )
            result.add(LazyHelper.resolve(binding, item))
        return result
    }

}

/**
 * A map whose values can be lazily resolved
 */
@CompileStatic
class LazyMap implements Map<String,Object> {

    /** The target map holding the values */
    @Delegate
    private Map<String,Object> target

    /** The context map against which dynamic properties are resolved */
    private Map binding

    private boolean dynamic

    boolean isDynamic() { dynamic }

    protected void setDynamic(boolean val) { dynamic = val }

    protected Map getBinding() { binding }

    void setBinding(Map map) { this.binding = map }

    protected Map<String,Object> getTarget() { target }

    protected void setTarget(Map<String,Object> obj) { this.target = obj }

    LazyMap() {
        target = new HashMap<>()
    }

    LazyMap( Map<String,Object> entries ) {
        assert entries != null
        target = new HashMap<>()
        putAll(entries)
    }

    /**
     * Resolve a directive *dynamic* value i.e. defined with a closure or lazy string
     *
     * @param name The directive name
     * @param value The value to be resolved
     * @return The resolved value
     */
    protected resolve( String name, value ) {

        /*
         * directive with one value and optional named parameter are converted
         * to a list object in which the first element is a map holding the named parameters
         * and the second is the directive value
         */
        if( value instanceof LazyList ) {
            def copy = new ArrayList(value.size())
            for( Object item : value ) {
                if( item instanceof Map )
                    copy.add( resolveParams(name, item as Map) )
                else
                    copy.add( resolveImpl(name, item) )
            }
            return copy
        }

        /*
         * resolve the values in a map object, preserving
         * lazy maps as they are
         */
        else if( value instanceof Map && value !instanceof LazyMap ) {
            return resolveParams(name, value)
        }

        /*
         * simple value
         */
        else {
            return resolveImpl(name, value)
        }

    }

    /**
     * Resolve directive *dynamic* named params
     *
     * @param name The directive name
     * @param value The map holding the named params
     * @return A map in which dynamic params are resolved to the actual value
     */
    private Map resolveParams( String name, Map value ) {

        final copy = new LinkedHashMap()
        final attr = (value as Map)
        for( Entry entry : attr.entrySet() ) {
            copy[entry.key] = resolveImpl(name, entry.value, true)
        }
        return copy
    }

    /**
     * Resolve a directive dynamic value
     *
     * @param name The directive name
     * @param value The value to be resolved
     * @param param When {@code true} points that it is a named parameter value, thus closure are only cloned
     * @return The resolved directive value
     */
    private resolveImpl( String name, value, boolean param=false ) {

        if( value instanceof LazyVar ) {
            return binding.get(value.name)
        }

        else if( value instanceof Closure ) {
            def copy = value.cloneWith(getBinding())
            if( param ) {
                return copy
            }

            try {
                return copy.call()
            }
            catch( MissingPropertyException e ) {
                if( getBinding() == null ) throw new IllegalStateException("Directive `$name` doesn't support dynamic value (or context not yet initialized)")
                else throw e
            }
        }

        else if( value instanceof GString ) {
            return value.cloneAsLazy(getBinding()).toString()
        }

        return value
    }

    /**
     * Override the get method in such a way that {@link Closure} values are resolved against
     * the {@link #binding} map
     *
     * @param key The map entry key
     * @return The associated value
     */
    Object get( key ) {
        return getValue(key)
    }

    Object getValue(Object key) {
        final value = target.get(key)
        return resolve(key as String, value)
    }

    Object put( String key, Object value ) {
        if( value instanceof Closure ) {
            dynamic |= true
        }
        else if( value instanceof GString ) {
            for( int i=0; i<value.valueCount; i++ )
                if (value.values[i] instanceof Closure)
                    dynamic |= true
        }
        return target.put(key, value)
    }

    @Override
    void putAll( Map entries ) {
        entries.each { k, v -> put(k as String, v) }
    }

    @Override
    String toString() {
        final allKeys = keySet()
        final result = new ArrayList<String>(allKeys.size())
        for( String key : allKeys ) { result << "$key: ${getProperty(key)}".toString() }
        result.join('; ')
    }

}

/**
 * A variable that can be lazily resolved
 */
@CompileStatic
@EqualsAndHashCode
@ToString
class LazyVar implements LazyAware {
    String name

    LazyVar(String name) {
        this.name = name
    }

    @Override
    Object resolve(Object binding) {
        if( binding !instanceof Map )
            throw new IllegalArgumentException("Can't resolve lazy var `$name` because the given binding is not a map")

        return ((Map)binding).get(name)
    }
}

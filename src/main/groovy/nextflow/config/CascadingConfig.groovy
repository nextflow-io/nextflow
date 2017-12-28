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

package nextflow.config
import java.lang.reflect.Method

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.Memoized
/**
 * Configuration object which fallback to a parent {@link CascadingConfig} object
 * when an attribute value is missing.
 *
 * @see ConfigField
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
abstract class CascadingConfig<K,V> {

    protected Map<K,V> config = [:]

    protected CascadingConfig<K,V> parent

    CascadingConfig() {}

    CascadingConfig(Map<K,V> config, CascadingConfig<K,V> parent) {
        this.config = config
        this.parent = parent
        config.keySet().each { K it -> checkField(it) }
    }

    @CompileDynamic
    private Set<String> discoverFields(Closure<Boolean> accept) {
        def result = new HashSet<String>()
        def clazz = this.getClass()
        while( clazz != null ) {
            def methods = clazz.getMethods()
            def names = (methods as List<Method>)?.findResults { field ->
                def annotation = field.getAnnotation(ConfigField)
                if( annotation && accept(annotation) )
                    return annotation.value() ?: field.getName().replaceFirst('^get','').uncapitalize()
                return null
            }
            result.addAll(names)
            clazz = clazz.getSuperclass()
        }

        return result ?: Collections.emptySet()
    }

    @Memoized
    protected Set<String> validFields() {
        return discoverFields({ true })
    }

    protected Set<String> privateFields() {
        return discoverFields({ ConfigField field -> field._private() })
    }

    final protected void checkField(K name) throws IllegalArgumentException {
        def fields = validFields()

        if( !fields.contains(name.toString()) )
            throw new IllegalArgumentException("Not a valid config attribute: `$name`")
    }

    protected Map<String,Object> copyPublicAttributes() {
        def skip = privateFields()
        def copy = [:]
        this.config.each { k,v -> if(!skip.contains(k)) copy[k]=v }
        return copy
    }

    boolean isEmpty() { config.isEmpty() }

    boolean containsAttribute( K key ) {
        config.containsKey(key) ?: ( parent ? parent.containsAttribute(key) : false )
    }

    V getAttribute(K key) {
        getAttribute(key, null)
    }

    V getAttribute(K key, V defValue) {
        checkField(key)
        return config.containsKey(key) ? config.get(key) : ( parent && parent.containsAttribute(key) ? parent.getAttribute(key) : defValue )
    }

    V getOrCreateAttribute(K key, Closure<V> missing) {
        def result = getAttribute(key)
        if( result == null ) {
            result = missing.call()
            config.put(key, result)
        }
        return result
    }

    void setAttribute( K key, V value ) {
        checkField(key)
        config.put(key,value)
    }

    Set<K> getAttributeNames() {
        new HashSet(config.keySet())
    }


    /**
     * Convert this object to an equivalent {@link ConfigObject}
     *
     * @return A {@link ConfigObject} holding teh same data
     */
    ConfigObject toConfigObject() {
        toConfigObject0(this.config)
    }

    @CompileDynamic
    protected ConfigObject toConfigObject0( Map map ) {

        def result = new ConfigObject()
        map.each { key, value ->
            if( value instanceof Map ) {
                result.put( key, toConfigObject0((Map)value) )
            }
            else if( value instanceof CascadingConfig ) {
                result.put( key, value.toConfigObject() )
            }
            else {
                result.put( key, value )
            }
        }

        return result
    }

}

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

        return result ?: Collections.<String>emptySet()
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
        def copy = new LinkedHashMap()
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

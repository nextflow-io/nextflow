/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.cloud.aws
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.cloud.LaunchConfig
import nextflow.config.CascadingConfig

/**
 * Holds the cloud configuration setting needed to fill-up the
 * cloud-init template.
 *
 * See {@link AmazonCloudDriver#cloudInitScript(nextflow.cloud.LaunchConfig)}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@PackageScope
@CompileStatic
class CloudBootTemplateBinding implements Map<String,Object> {

    private LaunchConfig config

    private Map<String,Object> local = [:]

    CloudBootTemplateBinding(LaunchConfig config) {
        this.config = config
    }

    private valueFor( String name ) {

        if( local.containsKey(name) )
            return local.get(name)

        def meta = config.metaClass.getMetaProperty(name)
        if( meta )
            return meta.getProperty(config)

        if( config.getAttributeNames().contains(name) )
            return config.getAttribute(name)

        throw new MissingPropertyException(name, CloudBootTemplateBinding)
    }

    private normalise(String name) {

        final value = valueFor(name)

        if( value == null ) {
            return ''
        }

        if( value instanceof Collection ) {
            return name == 'dockerPull' ? value.join(' ') : value.join(',')
        }

        if( value instanceof Map ) {
            def copy = [:]
            ((Map<String,Object>)value).each {key, val -> copy.put(key, normalise(key)) }
            return copy
        }

        if( value instanceof CascadingConfig ) {
            return value
        }

        return value.toString()
    }

    @Override
    int size() {
        throw new UnsupportedOperationException()
    }

    @Override
    boolean isEmpty() {
        throw new UnsupportedOperationException()
    }

    @Override
    boolean containsKey(Object key) {
        return local.containsKey(key) || config.getAttributeNames().contains(key) || config.metaClass.getMetaProperty(key.toString())
    }

    @Override
    boolean containsValue(Object value) {
        throw new UnsupportedOperationException()
    }

    @Override
    Object get(Object key) {
        normalise(key.toString())
    }

    @Override
    Object put(String key, Object value) {
        local.put(key,value)
    }

    @Override
    Object remove(Object key) {
        throw new UnsupportedOperationException()
    }

    @Override
    void putAll(Map<? extends String, ?> m) {
        throw new UnsupportedOperationException()
    }

    @Override
    void clear() {
        throw new UnsupportedOperationException()
    }

    @Override
    Set<String> keySet() {
        throw new UnsupportedOperationException()
    }

    @Override
    Collection<Object> values() {
        throw new UnsupportedOperationException()
    }

    @Override
    Set<Map.Entry<String, Object>> entrySet() {
        throw new UnsupportedOperationException()
    }
}

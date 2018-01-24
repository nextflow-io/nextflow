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

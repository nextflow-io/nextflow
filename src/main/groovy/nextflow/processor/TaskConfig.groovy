/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.exception.AbortOperationException
import nextflow.executor.BashWrapperBuilder
import nextflow.util.CmdLineHelper
import nextflow.util.Duration
import nextflow.util.MemoryUnit
/**
 * Task local configuration properties
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TaskConfig implements Map<String,Object> {

    /** The target map holding the values */
    @Delegate
    final private Map<String,Object> target

    /** The context map against which dynamic properties are resolved */
    private Map context

    private boolean dynamic

    TaskConfig() {
        target = new HashMap()
    }

    TaskConfig( Map<String,Object> entries ) {
        assert entries != null
        target = new HashMap<>()
        putAll(entries)
    }

    TaskConfig setContext( Map value ) {
        assert value != null
        context = value
        return this
    }

    private normalize( String key, value ) {
        if( value instanceof Closure ) {
            if( context )
                return context.with((Closure)value)

            try {
                return (value as Closure).call()
            }
            catch( MissingPropertyException e ) {
                throw new IllegalStateException("Directive `$key` doesn't support dynamic value")
            }
        }

        return value
    }

    /**
     * Override the get method in such a way that {@link Closure} values are resolved against
     * the {@link #context} map
     *
     * @param key The map entry key
     * @return The associated value
     */
    Object get( key ) {
        def value = target.get(key)
        normalize(key as String, value)
    }

    Object put( String key, Object value ) {
        if( value instanceof Closure ) {
            dynamic |= true
        }
        else if( value instanceof List && key == 'module' ) {
            // 'module' directive can be defined as a list of dynamic values
            value.each { if (it instanceof Closure) dynamic |= true }
        }
        target.put(key, value)
    }

    @Override
    void putAll( Map entries ) {
        entries.each { k, v -> put(k as String, v) }
    }

    @PackageScope
    boolean isDynamic() { dynamic }

    @PackageScope
    boolean hasContext() { context != null }

    def getProperty(String name) {

        def meta = metaClass.getMetaProperty(name)
        if( meta )
            return meta.getProperty(this)

        return this.get(name)
    }

    boolean getEcho() {
        def value = get('echo')
        ProcessConfig.toBool(value)
    }

    List<Integer> getValidExitStatus() {
        def result = get('validExitStatus')
        if( result instanceof List<Integer> )
            return result as List<Integer>

        if( result != null )
            return [result as Integer]

        return [0]
    }

    ErrorStrategy getErrorStrategy() {
        final strategy = get('errorStrategy')
        if( strategy instanceof CharSequence )
            return strategy.toString().toUpperCase() as ErrorStrategy

        if( strategy instanceof ErrorStrategy )
            return (ErrorStrategy)strategy

        if( strategy == null )
            return null

        throw new IllegalArgumentException("Not a valid `ErrorStrategy` value: ${strategy}")
    }


    MemoryUnit getMemory() {
        def value = get('memory')

        if( !value )
            return null

        if( value instanceof MemoryUnit )
            return (MemoryUnit)value

        try {
            new MemoryUnit(value.toString().trim())
        }
        catch( Exception e ) {
            throw new AbortOperationException("Not a valid 'memory' value in process definition: $value")
        }
    }

    MemoryUnit getDisk() {
        def value = get('disk')

        if( !value )
            return null

        if( value instanceof MemoryUnit )
            return (MemoryUnit)value

        try {
            new MemoryUnit(value.toString().trim())
        }
        catch( Exception e ) {
            throw new AbortOperationException("Not a valid 'disk' value in process definition: $value")
        }
    }

    Duration getTime() {
        def value = get('time')

        if( !value )
            return null

        if( value instanceof Duration )
            return (Duration)value

        if( value instanceof Number )
            return new Duration(value as long)

        try {
            new Duration(value.toString().trim())
        }
        catch( Exception e ) {
            throw new AbortOperationException("Not a valid `time` value in process definition: $value")
        }
    }

    int getCpus() {
        final value = get('cpus')
        value ? value as int : 1  // note: always return at least 1 cpus
    }

    int getMaxRetries() {
        def result = get('maxRetries')
        result ? result as int : 0
    }

    int getMaxErrors() {
        def result = get('maxErrors')
        result ? result as int : 0
    }

    List<String> getModule() {
        def value = get('module')
        if( value instanceof String ) {
            return ProcessConfig.parseModule(value)
        }

        if( value instanceof List ) {
            def result = []
            value.each {
                result = ProcessConfig.parseModule(normalize('module',it), result)
            }
            return result
        }

        if( value == null )
            return Collections.emptyList()

        throw new IllegalStateException("Not a valid `module` value: $value")
    }

    List<String> getShell() {
        final value = get('shell')
        if( !value )
            return BashWrapperBuilder.BASH

        if( value instanceof List )
            return (List)value

        if( value instanceof CharSequence )
            return [ value.toString() ]

        throw new IllegalArgumentException("Not a valid `shell` configuration value: ${value}")
    }

    Path getStoreDir() {
        def path = get('storeDir')
        if( !path )
            return null

        return (path as Path).complete()
    }


    /**
     * @return Parse the {@code clusterOptions} configuration option and return the entries as a list of values
     */
    List<String> getClusterOptionsAsList() {

        def opts = get('clusterOptions')
        if ( !opts ) {
            return Collections.emptyList()
        }

        if( opts instanceof Collection ) {
            return new ArrayList<String>(opts)
        }
        else {
            return CmdLineHelper.splitter( opts.toString() )
        }
    }

    @Override
    String toString() {
        def result = []
        keySet().each { key -> result << "$key: ${getProperty(key)}" }
        result.join('; ')
    }
}

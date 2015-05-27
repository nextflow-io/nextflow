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
import static nextflow.processor.TaskProcessor.TASK_CONFIG

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.exception.AbortOperationException
import nextflow.exception.FailedGuardException
import nextflow.executor.BashWrapperBuilder
import nextflow.script.TaskClosure
import nextflow.util.CmdLineHelper
import nextflow.util.Duration
import nextflow.util.MemoryUnit
/**
 * Task local configuration properties
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TaskConfig extends LazyMap {

    TaskConfig() { }

    TaskConfig( Map<String,Object> entries ) {
        super(entries)
    }

    /**
     * Assign the context map for dynamic evaluation of task config properties
     * @param context The context binding object
     * @return The {@code TaskConfig} instance itself
     */
    TaskConfig setContext( Map context ) {
        assert context != null

        // set the binding context for this map
        this.binding = context

        // set the binding context for 'ext' map
        if( target.ext instanceof LazyMap )
            (target.ext as LazyMap).binding = context

        // set the this object in the task context  in order to allow task properties to be resolved in process script
        context.put(TASK_CONFIG, this)

        return this
    }

    def getProperty(String name) {

        def meta = metaClass.getMetaProperty(name)
        if( meta )
            return meta.getProperty(this)

        return get(name)
    }

    def get( String key ) {
        if( key != 'ext' ) {
            return super.get(key)
        }

        if( !target.containsKey('ext') ) {
            target.ext = new LazyMap()
        }

        return target.ext
    }

    Object put( String key, Object value ) {
        if( key == 'module' && value instanceof List ) {
            // 'module' directive can be defined as a list of dynamic values
            value.each { if (it instanceof Closure) dynamic |= true }
            target.put(key, value)
        }
        else if( key == 'ext' && value instanceof Map ) {
            super.put( key, new LazyMap(value) )
        }
        else {
            super.put(key,value)
        }
    }


    @PackageScope
    boolean isDynamic() {
        if( dynamic )
            return true

        if( target.ext instanceof LazyMap )
            return (target.ext as LazyMap).dynamic

        return false
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
                result = ProcessConfig.parseModule(resolve('module',it), result)
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

    /**
     * Get a closure guard condition and evaluate to a boolean result
     *
     * @param name The name of the guard to test e.g. {@code when}
     * @return {@code true} when the condition is verified
     */
    protected boolean getGuard( String name, boolean defValue=true ) throws FailedGuardException {

        final code = target.get(name)
        if( code == null )
            return defValue

        String source = null
        try {
            if( code instanceof Closure ) {
                if( code instanceof TaskClosure ) source = code.getSource()
                return code.cloneWith(binding).call()
            }
            // try to convert to a boolean value
            return code as Boolean
        }
        catch( Throwable e ) {
            throw new FailedGuardException("Cannot evaluate `$name` expression", source, e)
        }

    }

}

/**
 * A map that resolve closure and gstring in a lazy manner
 *
 */
@CompileStatic
class LazyMap implements Map<String,Object> {

    /** The target map holding the values */
    @Delegate
    final protected Map<String,Object> target

    /** The context map against which dynamic properties are resolved */
    protected Map binding

    boolean dynamic

    LazyMap() {
        target = new HashMap<>()
    }

    LazyMap( Map<String,Object> entries ) {
        assert entries != null
        target = new HashMap<>()
        putAll(entries)
    }

    protected resolve( String key, value ) {
        if( value instanceof Closure ) {
            def copy = value.cloneWith(binding)
            try {
                return copy.call()
            }
            catch( MissingPropertyException e ) {
                throw new IllegalStateException("Directive `$key` doesn't support dynamic value")
            }
        }
        else if( value instanceof GString ) {
            return value.cloneWith(binding).toString()
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
        def value = target.get(key)
        resolve(key as String, value)
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
        target.put(key, value)
    }

    @Override
    void putAll( Map entries ) {
        entries.each { k, v -> put(k as String, v) }
    }

    @Override
    String toString() {
        def result = []
        keySet().each { key -> result << "$key: ${getProperty(key)}" }
        result.join('; ')
    }
}
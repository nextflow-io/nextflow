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

package nextflow.processor

import static nextflow.processor.TaskProcessor.*

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.Const
import nextflow.exception.AbortOperationException
import nextflow.exception.FailedGuardException
import nextflow.executor.BashWrapperBuilder
import nextflow.k8s.model.PodOptions
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
class TaskConfig extends LazyMap implements Cloneable {

    private transient Map cache = [:]

    TaskConfig() {  }

    TaskConfig( Map<String,Object> entries ) {
        super(entries)
    }

    TaskConfig clone() {
        def copy = (TaskConfig)super.clone()
        copy.setTarget(new HashMap<>(this.getTarget()))
        copy.newCache()
        return copy
    }

    private void newCache() {
        cache = [:]
    }

    /**
     * Assign the context map for dynamic evaluation of task config properties
     * @param context A {@link TaskContext} object that holds the task evaluation context
     * @return The {@code TaskConfig} instance itself
     */
    TaskConfig setContext( Map context ) {
        assert context != null

        // set the binding context for this map
        this.binding = context

        // clear cache to force re-compute dynamic entries
        this.cache.clear()

        // set the binding context for 'ext' map
        if( target.ext instanceof LazyMap )
            (target.ext as LazyMap).binding = context

        // set the this object in the task context in order to allow task properties to be resolved in process script
        context.put(TASK_CONTEXT_PROPERTY_NAME, this)

        return this
    }

    def getProperty(String name) {

        def meta = metaClass.getMetaProperty(name)
        if( meta )
            return meta.getProperty(this)

        return get(name)
    }

    def get( String key ) {
        if( cache.containsKey(key) )
            return cache.get(key)

        def result
        if( key == 'ext' ) {
            if( target.containsKey(key) )
                result = target.get(key)
            else {
                result = new LazyMap()
                target.put(key, result)
            }
        }
        else
            result = super.get(key)

        cache.put(key,result)
        return result
    }

    Object put( String key, Object value ) {
        if( cache != null )
            cache.remove(key)
        if( key == 'module' && value instanceof List ) {
            // 'module' directive can be defined as a list of dynamic values
            for( Object it : value ) {
                if (it instanceof Closure) {
                    final flag = super.isDynamic() | true
                    super.setDynamic(flag)
                }
            }
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
        if( super.isDynamic() )
            return true

        if( target.ext instanceof LazyMap )
            return (target.ext as LazyMap).isDynamic()

        return false
    }


    boolean getEcho() {
        def value = get('echo')
        toBool(value)
    }

    private static boolean toBool( value )  {
        if( value instanceof Boolean ) {
            return value.booleanValue()
        }

        return value != null && value.toString().toLowerCase() in Const.BOOL_YES
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
            return ErrorStrategy.TERMINATE

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
        def defResult = getErrorStrategy() == ErrorStrategy.RETRY ? 1 : 0
        result ? result as int : defResult
    }

    int getMaxErrors() {
        def result = get('maxErrors')
        result ? result as int : 0
    }

    List<String> getModule() {
        def value = get('module')

        if( value instanceof List ) {
            def result = []
            for( String name : value ) {
                result.addAll( name.tokenize(':') )
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

    List<PublishDir> getPublishDir() {
        def dirs = get('publishDir')
        if( !dirs ) {
            return Collections.emptyList()
        }

        if( dirs instanceof List ) {
            final List<PublishDir> result = new ArrayList<>(dirs.size())
            for( Object params : dirs ) {
                if( !params ) continue
                if( params instanceof Map )
                    result.add( PublishDir.create(params) )
                else
                    throw new IllegalArgumentException("Not a valid PublishDir entry [${params.getClass().getName()}] $params")
            }
            return result
        }

        throw new IllegalArgumentException("Not a valid PublishDir collection [${dirs.getClass().getName()}] $dirs")
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

    Integer getAttempt() {
        get('attempt') as Integer ?: 1
    }

    Integer getErrorCount() {
        get('errorCount') as Integer ?: 0
    }

    Integer getRetryCount() {
        get('retryCount') as Integer ?: 0
    }

    PodOptions getPodOptions() {
        new PodOptions((List)get('pod'))
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
                return code.cloneWith(getBinding()).call()
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
 */
@CompileStatic
class LazyMap implements Map<String,Object> {

    /** The target map holding the values */
    @Delegate
    private Map<String,Object> target

    /** The context map against which dynamic properties are resolved */
    private Map binding

    private boolean dynamic

    protected boolean isDynamic() { dynamic }

    protected void setDynamic(boolean val) { dynamic = val }

    protected Map getBinding() { binding }

    protected void setBinding(Map map) { this.binding = map }

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
        if( value instanceof ConfigList ) {
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
         * resolve the values in a map object
         */
        else if( value instanceof Map ) {
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
    private resolveParams( String name, Map value ) {

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

        if( value instanceof Closure ) {
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
            return value.cloneWith(getBinding()).toString()
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

@CompileStatic
class ConfigList implements List {

    @Delegate
    private List target

    ConfigList() {
        target = []
    }

    ConfigList(int size) {
        target = new ArrayList(size)
    }

    ConfigList(Collection items) {
        target = new ArrayList(items)
    }

}
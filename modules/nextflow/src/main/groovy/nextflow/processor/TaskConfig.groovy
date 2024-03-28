/*
 * Copyright 2013-2024, Seqera Labs
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
import nextflow.ast.NextflowDSLImpl
import nextflow.exception.AbortOperationException
import nextflow.exception.FailedGuardException
import nextflow.executor.BashWrapperBuilder
import nextflow.executor.res.AcceleratorResource
import nextflow.executor.res.DiskResource
import nextflow.k8s.model.PodOptions
import nextflow.script.TaskClosure
import nextflow.util.CmdLineHelper
import nextflow.util.CmdLineOptionMap
import nextflow.util.Duration
import nextflow.util.LazyMap
import nextflow.util.MemoryUnit
/**
 * Task local configuration properties
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TaskConfig extends LazyMap implements Cloneable {

    static public final int EXIT_ZERO = 0

    private transient Map cache = new LinkedHashMap(20)

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

    /**
     * Evaluate a task config attribute. The main difference of this method
     * is that does not cache the result value.
     *
     * @param path
     *      The task config to be evaluated e.g. `cpus`. Note it allows
     *      traversing nested object separating keys with a dot e.g. `ext.args`
     * @return
     *      The value associate with the config key. Dynamic value i.e. closure are
     *      automatically resolved to target value.
     */
    Object eval(String path) {
        return eval0(this, path.tokenize('.'), path)
    }

    private Object eval0(Object object, List<String> path, String key ) {
        assert path, "Missing task attribute name"
        def result = null
        if( object instanceof LazyMap ) {
            result = ((LazyMap)object).getValue(path.first())
        }
        else if( Object instanceof Map ) {
            result = ((Map)object).get(path.first())
        }
        else if( path.size()>1 ) {
            throw new IllegalArgumentException()
        }

        if( path.size()==1 || result==null )
            return result

        return eval0( result, path.subList(1,path.size()), key )
    }

    def getProperty(String name) {

        def meta = metaClass.getMetaProperty(name)
        if( meta )
            return meta.getProperty(this)

        return get(name)
    }

    final getRawValue(String key) {
        return target.get(key)
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

    boolean isDynamic() {
        if( super.isDynamic() )
            return true

        if( target.ext instanceof LazyMap )
            return (target.ext as LazyMap).isDynamic()

        return false
    }

    String getBeforeScript() {
        return get('beforeScript')
    }

    String getAfterScript() {
        return get('afterScript')
    }

    def getCleanup() {
        return get('cleanup')
    }

    String getStageInMode() {
        return get('stageInMode')
    }

    String getStageOutMode() {
        return get('stageOutMode')
    }

    boolean getDebug() {
        // check both `debug` and `echo` for backward
        // compatibility until `echo` is not removed
        def value = get('debug') || get('echo')
        return toBool(value)
    }

    private static boolean toBool( value )  {
        if( value instanceof Boolean ) {
            return value.booleanValue()
        }

        return value != null && value.toString().toLowerCase() in Const.BOOL_YES
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

    DiskResource getDiskResource() {
        def value = get('disk')

        if( value instanceof Map )
            return new DiskResource((Map)value)

        if( value != null )
            return new DiskResource(value)

        return null
    }

    MemoryUnit getDisk() {
        getDiskResource()?.getRequest()
    }

    private Duration getDuration0(String key) {
        def value = get(key)

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
            throw new AbortOperationException("Not a valid `$key` value in process definition: $value")
        }

    }

    Duration getTime() {
        return getDuration0('time')
    }

    Duration getMaxSubmitAwait() {
        return getDuration0('maxSubmitAwait')
    }

    boolean hasCpus() {
        get('cpus') != null
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
            List<String> result = []
            for( String name : value ) {
                result.addAll( name.tokenize(':') )
            }
            return result
        }

        if( value == null )
            return Collections.emptyList()

        throw new IllegalStateException("Not a valid `module` value: $value")
    }

    List<String> getSecret() {
        return (List<String>) get('secret')
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

    String getClusterOptions() {
        return get('clusterOptions')
    }
    
    def getContainer() {
        return get('container')
    }

    Architecture getArchitecture() {
        final value = get('arch')
        if( value instanceof CharSequence )
            return new Architecture(value.toString())
        if( value instanceof Map )
            return new Architecture(value)
        if( value != null )
            throw new IllegalArgumentException("Invalid `arch` directive value: $value [${value.getClass().getName()}]")
        return null
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

    Integer getSubmitAttempt() {
        get('submitAttempt') as Integer ?: 1
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

    AcceleratorResource getAccelerator() {
        final value = get('accelerator')
        if( value instanceof Number )
            return new AcceleratorResource(value)
        if( value instanceof Map )
            return new AcceleratorResource(value)
        if( value != null )
            throw new IllegalArgumentException("Invalid `accelerator` directive value: $value [${value.getClass().getName()}]")
        return null
    }

    String getMachineType() {
        return get('machineType')
    }

    String getContainerOptions() {
        def opts = get('containerOptions')
        return opts instanceof CharSequence ? opts.toString() : null
    }

    CmdLineOptionMap getContainerOptionsMap() {
        def opts = get('containerOptions')
        if( opts instanceof Map )
            return CmdLineOptionMap.fromMap(opts)
        if( opts instanceof CharSequence )
            return CmdLineHelper.parseGnuArgs(opts.toString())
        if( opts!=null )
            throw new IllegalArgumentException("Invalid `containerOptions` directive value: $opts [${opts.getClass().getName()}]")
        return CmdLineOptionMap.emptyOption()
    }

    Map<String, String> getResourceLabels() {
        return get('resourceLabels') as Map<String, String> ?: Collections.<String,String>emptyMap()
    }

    String getResourceLabelsAsString() {
        final res = getResourceLabels()
        final result = new StringBuilder()
        int c=0
        for( Map.Entry<String,String> it : res ) {
            if(c++>0) result.append(',')
            result.append(it.key).append('=').append(it.value)
        }
        return result
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


    protected TaskClosure getStubBlock() {
        final code = target.get(NextflowDSLImpl.PROCESS_STUB)
        if( !code )
            return null
        if( code instanceof TaskClosure )
            return code
        throw new IllegalStateException()
    }

}

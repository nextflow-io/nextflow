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

package nextflow.processor
import java.nio.file.Path
import java.nio.file.Paths
import java.util.concurrent.atomic.AtomicBoolean

import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.io.Output
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Global
import nextflow.exception.ProcessException
import nextflow.util.KryoHelper
/**
 * Map used to delegate variable resolution to script scope
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class TaskContext implements Map<String,Object>, Cloneable {

    @Delegate
    private Map<String,Object> holder

    /**
     * Used to show the override warning message only the very first time
     */
    private transient final overrideWarnShown = new AtomicBoolean()

    /**
     * The main script owning the process
     */
    private Script script

    /**
     * The process name
     */
    private String name

    /**
     * The name of the variables not hold in the task context, but referenced in the global script binding object
     */
    private transient Set<String> variableNames

    TaskContext( TaskProcessor processor, Map holder = [:]) {
        assert holder != null
        this.holder = holder
        this.script = processor.ownerScript
        this.name = processor.name

        // fetch all the variables names referenced by the script body and retain
        // only the ones not declared as input or output, because these are supposed to
        // to be the ones provided by the *external* script context
        variableNames = processor.getTaskBody().getValNames() ?: []
        if( variableNames ) {
            Set<String> declaredNames = []
            declaredNames.addAll( processor.config.getInputs().getNames() )
            declaredNames.addAll( processor.config.getOutputs().getNames()  )
            if( declaredNames )
                variableNames = variableNames - declaredNames
        }

        log.trace "Binding names for '$name' > $variableNames"
    }


    protected TaskContext(Script script, Map holder, String name) {
        this.script = script
        this.holder = holder
        this.name = name
        def names = script.getBinding()?.getVariables()?.keySet()
        this.variableNames = names ? new HashSet<>(names) : new HashSet<>()
        log.trace "Binding names for '$name' > $variableNames"
    }

    /** ONLY FOR TEST PURPOSE -- do not use */
    protected TaskContext() { }

    TaskContext clone() {
        def copy = (TaskContext)super.clone()
        copy.setHolder( (Map)holder.clone() )
        return copy
    }

    /**
     * @return The inner map holding the process variables
     */
    public Map<String,Object> getHolder() { holder }

    private void setHolder( Map<String,Object> holder ) {
        this.holder = holder
    }

    /**
     * @return The script instance to which this map reference i.e. the main script object
     */
    public Script getScript() { script }

    /**
     * @return
     *      The set of variable and properties referenced in the user script.
     *      NOTE: it includes properties in the form {@code object.propertyName}
     */
    public Set<String> getVariableNames() { variableNames }

    @Override
    String toString() {
        "DelegateMap[process: $name; script: ${script?.class?.name}; holder: ${holder}]"
    }

    @Override
    public Object get(Object property) {
        assert property

        if( holder.containsKey(property) ) {
            return holder.get(property)
        }

        if( script?.getBinding()?.hasVariable(property.toString()) ) {
            return script.getBinding().getVariable(property.toString())
        }

        throw new MissingPropertyException("Unknown variable '$property' -- Make sure you didn't misspell it or define somewhere in the script before use it", property as String, null)
    }

    Object invokeMethod(String name, Object args) {
        if( name == 'template' )
            template(args)
        else
            script.invokeMethod(name, args)
    }

    public getProperty( String name ) {
        get((String)name)
    }

    public void setProperty( String name, def value ) {
        put(name, value)
    }

    @Override
    public put(String property, Object newValue) {

        if( property == 'task' && !(newValue instanceof TaskConfig ) && !overrideWarnShown.getAndSet(true) ) {
            log.warn "Process $name overrides reserved variable `task`"
        }

        holder.put(property, newValue)
    }


    def byte[] serialize() {
        try {
            def map = holder
            if( map.get(TaskProcessor.TASK_CONTEXT_PROPERTY_NAME) instanceof TaskConfig ) {
                map = new LinkedHashMap<String, Object>(holder)
                map.remove(TaskProcessor.TASK_CONTEXT_PROPERTY_NAME)
            }

            return KryoHelper.serialize(map)
        }
        catch( Exception e ) {
            log.warn "Cannot serialize context map. Cause: ${e.cause} -- Resume will not work on this process"
            log.debug "Failed to serialize delegate map items: ${dumpMap(holder)}", e
            return null
        }
    }

    static TaskContext deserialize(TaskProcessor processor, byte[] buffer) {
        def map = (Map)KryoHelper.deserialize(buffer)
        new TaskContext(processor, map)
    }


    @PackageScope
    static String dumpMap( Map map ) {
        def result = []
        result << "[ "
        map.each { key, value -> result << "  '$key':[${value?.class?.name}] = ${value}" }
        result << "]"
        return result.join('\n')
    }

    /**
     * Serialize the {@code DelegateMap} instance to a byte array
     */
    def byte[] dehydrate() {
        def kryo = KryoHelper.kryo()
        def buffer = new ByteArrayOutputStream(5*1024)
        def out = new Output(buffer)
        out.writeString(name)
        kryo.writeClassAndObject(out,holder)

        // -- the script class
        kryo.writeObject(out, script.class)

        // -- only the binding values for which there's an entry in the holder map
        final copy = new Binding()
        variableNames.each { String it ->
            // name can be a property, in this case use the root object
            def p = it.indexOf('.')
            def var = ( p == -1 ? it : it.substring(0,p) )
            checkAndSet(var, copy)
        }
        log.trace "Delegate for $name > binding copy: ${copy.getVariables()}"
        kryo.writeObject(out, copy)

        out.flush()
        return buffer.toByteArray()
    }

    private void checkAndSet( String name, Binding target ) {

        final binding = this.script.getBinding()
        if( !binding.hasVariable(name) )
            return

        def val = binding.getVariable(name)
        if( val instanceof DataflowReadChannel || val instanceof DataflowWriteChannel )
            return

        if( val instanceof Path || val instanceof Serializable ) {
            target.setVariable(name, val)
        }

    }

    /**
     * Deserialize and create a new instance of the {@code DelegateMap} using the provided byte array serialized binary
     *
     * @param binary
     *          The binary output of a previous {@code #dehydrate} invocation
     * @param loader
     *          An optional class loader to be used to resolve script class when this object
     *          need to be reacted in a remote JVM
     * @return
     *      A {@code DelegateMap} object instantiated using the provided binary byte[]
     */
    static TaskContext rehydrate(byte[] binary, ClassLoader loader = null) {
        assert binary
        final kryo = KryoHelper.kryo()

        def ClassLoader prev = null
        if( loader ) {
            prev = kryo.getClassLoader()
            kryo.setClassLoader(loader)
        }

        try {
            def input = new Input(new ByteArrayInputStream(binary))
            def name = input.readString()
            Map holder = (Map)kryo.readClassAndObject(input)
            Class<Script> clazz = kryo.readObject(input,Class)
            Binding binding = kryo.readObject(input,Binding)

            Script script = clazz.newInstance()
            script.setBinding(binding)
            return new TaskContext(script, holder, name)
        }
        finally {
            // set back the original class loader
            if( prev ) kryo.setClassLoader(prev)
        }

    }

    /**
     * Helper method to be used inside a process to load a script template. For example:
     *
     *     <pre>
     *         process foo {
     *             input:
     *             file file_name
     *
     *             script:
     *             template('script/file.sh')
     *         }
     *     </pre>
     *
     * @param path The path to the template file. When relative the file is resolved against script {@code $baseDir/templates} folder
     * @return The resolved template path
     */
    protected Path template( path ) {
        if( !path )
            throw new ProcessException("Process `$name` missing template name")

        if( !(path instanceof Path) )
            path = Paths.get(path.toString())

        // if the path is already absolute just return it
        if( path.isAbsolute() )
            return path

        // otherwise make
        def base = Global.session.baseDir
        if( base )
            return base.resolve('templates').resolve(path)

        // if the base dir is not available just use as it is
        return path
    }

}

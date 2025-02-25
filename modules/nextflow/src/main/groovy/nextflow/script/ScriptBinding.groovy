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

package nextflow.script

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.NF
import nextflow.Session
import nextflow.exception.AbortOperationException
/**
 * Defines the script execution context. By default provided the following variables
 * <li>{@code __$session}: the current execution session
 * <li>{@code params}: the parameters map specified by the users on the program CLI using the '--' prefix
 * <li>{@code args}: the list of programs arguments specified on the program CLI
 *
 * <p>
 *     These values cannot be overridden by definition
 *
 * Read more about 'binding variables'
 * http://groovy.codehaus.org/Scoping+and+the+Semantics+of+%22def%22
 *
 */
@Slf4j
@CompileStatic
class ScriptBinding extends WorkflowBinding {

    private Session session

    private Path scriptPath

    private List<String> args

    private ParamsMap params

    private Map configEnv = Collections.emptyMap()

    private String entryName

    /**
     * Creates a new nextflow script binding object
     *
     * @param config Nextflow configuration object
     * @return A new {@link ScriptBinding} instance
     */
    ScriptBinding() {
        this(new LinkedHashMap(20))
    }

    ScriptBinding(Map vars) {
        super(vars)

        // create and populate args
        args = new ArrayList<>()
        if( vars.args ) {
            if( !(vars.args instanceof List<String>) ) throw new IllegalArgumentException("ScriptBinding 'args' must be a List value")
            args.addAll((List<String>)vars.args)
        }
        vars.put('args', args)
        
        // create and populate args
        params = new ParamsMap()
        if( vars.params ) {
            if( !(vars.params instanceof Map) ) throw new IllegalArgumentException("ScriptBinding 'params' must be a Map value")
            params.putAll((Map)vars.params)
        }
        vars.params = params
    }

    ScriptBinding(Session session) {
        this()
        this.session = session
    }

    ScriptBinding setSession(Session session ) {
        this.session = session
        setConfigEnv(session.configEnv)
        return this
    }

    ScriptBinding setScriptPath(Path path) {
        scriptPath = path
        return this
    }

    ScriptBinding setEntryName(String entry) {
        this.entryName = entry
        return this
    }

    ScriptBinding setModule(boolean value ) {
        module = value
        return this
    }

    private ScriptBinding setConfigEnv(Map map ) {
        this.configEnv = map != null ? map : Collections.emptyMap()
        return this
    }

    ParamsMap getParams() { params }

    Session getSession() { session }

    Path getScriptPath() { scriptPath }

    String getEntryName() { entryName }

    @Memoized
    protected Map<String,String> getSysEnv() {
        new HashMap(System.getenv())
    }

    /**
     * The map of the CLI named parameters
     *
     * @param values
     */
    ScriptBinding setParams(Map<String,Object> values ) {
        if( values )
            params.putAll(values)
        return this
    }

    /**
     * The list of CLI arguments (unnamed)
     * @param values
     */
    ScriptBinding setArgs(List<String> values ) {
        args.clear()
        if( values )
         args.addAll(values)
        return this
    }

    /**
     * Try to get a value in the current bindings, if does not exist try
     * to fallback on the session environment.
     *
     * @param name
     * @return
     */
    Object getVariable( String name ) {

        if( super.hasVariable(name) )
            return super.getVariable(name)

        if( configEnv.containsKey(name) )
            return configEnv.get(name)

        if( sysEnv.containsKey(name) )
            return sysEnv.get(name)

        super.getVariable(name)
    }

    /**
     * Override variable existence
     *
     * @param name
     * @return
     */
    boolean hasVariable( String name ) {
        super.hasVariable(name) || configEnv.containsKey(name) || sysEnv.containsKey(name)
    }

    @Override
    void setVariable( String name, Object value ) {
        if( name == 'channel' ) {
            throw new IllegalAccessException("The use of the identifier `$name` as variable name is not allowed")
        }
        if( name != 'args' && name != 'params' )
            super.setVariable(name, value)
    }

    /**
     * Lookup the name of a variable giving the value reference
     *
     * @param value The value for which the variable name is needed
     * @return The associated variable name
     */
    String getVariableName(value) {
        super.getVariables().find { entry -> entry.value?.is(value) }?.getKey()
    }

    /**
     * Implements immutable params map
     */
    @CompileStatic
    static class ParamsMap implements Map<String,Object> {

        private List<String> scriptAssignment = []

        @Delegate
        private Map<String,Object> target

        ParamsMap() {
            target = new LinkedHashMap<>()
        }

        ParamsMap(Map<String,Object> values) {
            this()
            putAll(values)
        }

        private ParamsMap(ParamsMap values, Map<String,?> overrides) {
            this(values)

            for( def entry : overrides )
                put0(entry.key, entry.value, true)
        }

        ParamsMap copyWith(Map<String,?> overrides) {
            return new ParamsMap(this, overrides)
        }

        @Override
        Object get(Object key) {
            if( !target.containsKey(key) ) {
                final msg = "Access to undefined parameter `$key` -- Initialise it to a default value eg. `params.$key = some_value`"
                if( NF.isStrictMode() )
                    throw new AbortOperationException(msg)
                log.warn1(msg, firstOnly: true)
                return null
            }
            return target.get(key)
        }

        /**
         * A name-value pair to the map object. Note this API it supposed only to be invoked indirectly
         * when a parameter is defined/assigned. It should not be invoked directly.
         *
         * @param name The pair name
         * @param value The pair value
         * @return the previous value associated with <tt>name</tt>
         */
        @Override
        String put(String name, Object value) {
            assert name
            if( name in scriptAssignment )
                log.warn("`params.$name` is defined multiple times -- Assignments following the first are ignored")
            else
                scriptAssignment << name
            return put0(name,value)
        }

        private String put0(String name, Object value, boolean allowOverride=false) {
            return allowOverride || !target.containsKey(name)
                ? target.put(name, value)
                : null
        }

        @Override
        void putAll(Map<? extends String,? extends Object> other) {
            for( def entry : other.entrySet() ) {
                put0(entry.getKey(), entry.getValue())
            }
        }

    }

}

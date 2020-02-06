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

package nextflow.script

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Session
import org.apache.commons.lang.StringUtils
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
            if( !(vars.args instanceof List) ) throw new IllegalArgumentException("ScriptBinding 'args' must be a List value")
            args.addAll((List)vars.args)
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
        if( name == 'channel' )
            log.warn 'The use of the identifier `channel` as variable name is discouraged and will be deprecated in a future version'
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
     * Holds parameter immutable values
     */
    @CompileStatic
    static class ParamsMap implements Map<String,Object> {

        private Set<String> readOnlyNames

        private List<String> realNames

        private List<String> scriptAssignment = []

        @Delegate
        private Map<String,Object> target

        ParamsMap() {
            readOnlyNames = []
            realNames = []
            target = new LinkedHashMap<>()
        }

        ParamsMap(Map<String,Object> values) {
            this()
            putAll(values)
        }

        private ParamsMap(ParamsMap template, Map overrides) {
            this(template)
            allowNames(overrides.keySet())
            putAll(overrides)
        }

        ParamsMap copyWith(Map overrides) {
            return new ParamsMap(this, overrides)
        }

        private ParamsMap allowNames(Set names) {
            readOnlyNames.removeAll(names)
            return this
        }

        @Override
        Object get(Object key) {
            if( !target.containsKey(key) ) {
                log.warn1("Access to undefined parameter `$key` -- Initialise it to a default value eg. `params.$key = some_value`", firstOnly: true)
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
            (name in scriptAssignment
                    ? log.warn("`params.$name` is defined multiple times -- Assignments following the first are ignored")
                    : scriptAssignment << name )
            put0(name,value)
        }

        private String put0(String name, Object value) {
            // keep track of the real name
            realNames << name

            // normalize the name
            def name2 = name.contains('-') ? hyphenToCamelCase(name) : camelCaseToHyphen(name)

            final readOnly = name in readOnlyNames || name2 in readOnlyNames

            def result = null
            if( !readOnly ) {
                readOnlyNames << name
                readOnlyNames << name2
                result = target.put(name, value)
                target.put(name2, value)
            }

            return result
        }

        @Override
        void putAll(Map other) {
            for( Entry<String,Object> entry : other.entrySet() ) {
                put0(entry.getKey(), entry.getValue())
            }
        }

        /**
         * Renders a string containing all name-value pairs as a command line string
         *
         * @param opts String formatting options:
         *   {@code sep}: the key-value pair separator
         *   {@code quote}: the character to be used to quote the parameter value
         *   {@code prefix}: the character(s) prepended to the parameter key
         *
         * @return A command line formatted string e.g. {@code --foo x --bar y}
         */
        String all(Map opts = null) {

            String sep = opts?.sep ?: ' '
            String quote = opts?.quote ?: ''
            String prefix = opts?.prefix ?: '--'

            def result = []
            realNames.each {
                result << ("$prefix$it$sep${wrap(this.get(it).toString(), quote)}".toString())
            }
            return result.join(' ')
        }

        /**
         * Wrap a string value with quote characters
         *
         * @param str The string with value to be wrapped by quote characters
         * @param quote The quote character
         * @return The quoted string
         */
        @PackageScope
        static String wrap( String str, String quote ) {
            if( !quote && (str.contains(' ') || str.contains('"') || str.contains("'") ))
                quote = '"'

            if( !quote )
                return str

            def result = new StringBuilder()
            result.append(quote)
            str.each {
                result.append( it == quote ? '\\'+quote : it )
            }
            result.append(quote)
            return result.toString()
        }

        /**
         *  Converts a string containing words separated by a hyphen character to a camelCase string
         *
         * @param str The string to be converted
         * @return
         */
        @PackageScope
        static String hyphenToCamelCase( String str ) {

            if( !str ) { return str }

            def result = new StringBuilder()
            str.split('-').eachWithIndex{ String entry, int i ->
                result << (i>0 ? StringUtils.capitalize(entry) : entry )
            }

            return result.toString()
        }

        /**
         * Converts a camel-case string to a string where words are separated by hyphen character
         *
         * @param str The string to be converted
         * @return A string where camel-case words are converted to words separated by hyphen character
         */
        @PackageScope
        static String camelCaseToHyphen( String str ) {

            def lower = 'a'..'z'
            def upper = 'A'..'Z'

            def result = new StringBuilder()
            if( !str ) {
                return str
            }

            result << str[0]
            for( int i=1; i<str.size(); i++ ) {
                if( str[i] in upper && str[i-1] in lower  ) {
                    result << '-'
                    if( i+1<str.size() && str[i+1] in lower ) {
                        result << str[i].toLowerCase()
                    }
                    else {
                        result << str[i]
                    }
                }
                else {
                    result << str[i]
                }
            }

            return result.toString()
        }

    }

}

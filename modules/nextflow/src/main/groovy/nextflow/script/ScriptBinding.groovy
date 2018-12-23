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

package nextflow.script
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.util.ReadOnlyMap
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
class ScriptBinding extends Binding {

    private Map config

    /**
     * Creates a new nextflow script binding object
     *
     * @param config Nextflow configuration object
     * @return A new {@link ScriptBinding} instance
     */
    def ScriptBinding(Map config) {
        super( new ReadOnlyMap( [args:[], params: new ParamsMap() ]) )
        this.config = config != null ? config : [:]
    }

    Map getConfig() {
        return config
    }

    @Memoized
    protected Map<String,String> getSysEnv() {
        new HashMap(System.getenv())
    }

    /**
     * The fallback session environment
     */
    @Memoized
    protected Map getConfigEnv() {
        if( config.env instanceof Map )
            return (Map)config.env
        else
            [:]
    }


    /**
     * The map of the CLI named parameters
     *
     * @param values
     */
    def void setParams( Map<String,Object> values ) {
        if( !values )
            return
        def params = (ParamsMap)super.getVariable('params')
        params.putAll(values)
    }

    /**
     * The list of CLI arguments (unnamed)
     * @param values
     */
    def void setArgs( List<String> values ) {
        (getVariables() as ReadOnlyMap).force( 'args', values  )
    }

    /**
     * Try to get a value in the current bindings, if does not exist try
     * to fallback on the session environment.
     *
     * @param name
     * @return
     */
    def getVariable( String name ) {

        if( super.hasVariable(name) ) {
            return super.getVariable(name)
        }
        else if( configEnv.containsKey(name) ) {
            configEnv.get(name)
        }
        else if( sysEnv.containsKey(name) ) {
            return sysEnv.get(name)
        }
        else {
            throw new MissingPropertyException(name, getClass())
        }
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

        private Set<String> readOnlyNames = []

        private List<String> realNames = []

        private List<String> scriptAssignment = []

        @Delegate
        private Map<String,Object> target = new LinkedHashMap<>()

        ParamsMap() {}

        ParamsMap(Map<String,Object> copy) {
            putAll(copy)
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
            for( Map.Entry<String,Object> entry : other.entrySet() ) {
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

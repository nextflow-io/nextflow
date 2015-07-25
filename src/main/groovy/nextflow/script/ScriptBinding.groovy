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

        def params = super.getVariable('params') as ParamsMap
        values ?. each { String name, Object value ->
            params.put(name, value)
        }
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

    /**
     * Holds parameter immutable values
     */
    @CompileStatic
    static class ParamsMap implements Map<String,Object> {

        private Set<String> readOnlyNames = []

        private List<String> realNames = []

        @Delegate
        private Map<String,Object> target = new LinkedHashMap<>()

        ParamsMap() {}

        ParamsMap(Map<String,Object> copy) {
            copy.each { entry -> this.put(entry.key, entry.value) }
        }

        /**
         * A name-value pair to the map object
         *
         * @param name The pair name
         * @param value The pair value
         * @return the previous value associated with <tt>name</tt>
         */
        @Override
        String put(String name, Object value) {
            assert name

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

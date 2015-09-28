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

package nextflow.scm
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.Const
import nextflow.config.ComposedConfigSlurper
import nextflow.exception.ConfigParseException
/**
 * Models a repository provider configuration attributes
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ProviderConfig {

    private String name

    private Map attributes

    ProviderConfig( String name, Map values ) {
        this.name = name
        this.attributes = new HashMap(values)

        // init defaults
        switch( name ) {
            case 'github':
                attributes.platform = name
                if( !attributes.host ) attributes.host = 'https://github.com'
                if( !attributes.endpoint ) attributes.endpoint = 'https://api.github.com'
                break

            case 'gitlab':
                attributes.platform = name
                if( !attributes.host ) attributes.host = 'https://gitlab.com'
                break

            case 'bitbucket':
                attributes.platform = name
                if( !attributes.host ) attributes.host = 'https://bitbucket.org'
        }

        if( attributes.path )
            attributes.platform = 'file'
    }

    ProviderConfig( String name ) {
        this(name, [:])
    }

    /**
     * @return The provider name as defined in the configuration file
     */
    String getName() { name }

    /**
     * @return
     *      Remote host including protocol and eventually port number e.g.
     *      {@code https://github.com}
     */
    String getHost() { attributes.host }

    /**
     * @return
     *      Remote host name e.g. {@code github.com}
     */
    String getDomain() {
        def result = host
        def p = result.indexOf('://')
        p != -1 ? result.substring(p+3) : result
    }

    /**
     * @return The provider identifier e.g. {@code github}, {@code gitlab} or {@code bitbucket}
     */
    String getPlatform() {
        attributes.platform
    }

    /**
     * @return The authentication token
     */
    String getAuth() { attributes.auth }

    String getAuthObfuscated() {
        if( !auth ) return null
        def p=auth.indexOf(':')
        if( p != -1 ) {
            def result = []
            result[0] = p>0 ? auth.substring(0,p) : '-'
            result[1] = p+1<auth.size() ? '*' * (auth.size()-p-1) : '-'
            return result.join(':')
        }
        else {
            return '*' * auth.size()
        }
    }

    /**
     * Sets the authentication token
     * @param value The auth token string
     * @return The {@link ProviderConfig} object itself
     */
    ProviderConfig setAuth( String value ) {
        attributes.auth = value
        return this
    }

    /**
     * Sets the authentication token putting together the user name and password strings
     *
     * @param userName The user name string
     * @param password The password string
     * @return The {@link ProviderConfig} object itself
     */
    ProviderConfig setAuth( String userName, String password ) {
        def pair = [ userName?:'' , password?:'' ]
        attributes.auth = pair.join(':')
        return this
    }

    /**
     * @return
     *      The API root endpoint url e.g. {@code http://api.github.com}. Fallback on
     *      the {@link #getHost()} if not defined
     */
    String getEndpoint() {
        attributes.endpoint ?: attributes.host
    }

    String getUser() {
        if(!auth)
            return null

        def p=auth.indexOf(':')
        return p>0 ? auth.substring(0,p) : null
    }

    String getPassword() {
        if(!auth)
            return null

        def p = auth.indexOf(':')
        if( p != -1 ) {
            return p+1<auth.size() ? auth.substring(p+1) : null
        }
        else {
            return auth
        }
    }

    String getPath() { attributes.path?.toString() }

    String toString() {
        "ProviderConfig[name: $name, platform: $platform, host: $host]"
    }

    @PackageScope
    static Map parse(String text) {
        def slurper = new ComposedConfigSlurper()
        def env = new HashMap( System.getenv() )
        slurper.setBinding(env)
        return slurper.parse(text)
    }

    static List<ProviderConfig> createDefault() {
        createFromMap(null)
    }

    @PackageScope
    static List<ProviderConfig> createFromMap(Map<String,?> config) {

        def providers = (Map<String,?>)config?.providers

        List<ProviderConfig> result = []
        if( providers ) {
            providers.keySet().each { String name ->
                def attrs = (Map)providers.get(name)
                result << new ProviderConfig(name, attrs)
            }
        }

        addDefaults(result)
        return result
    }

    @PackageScope
    static List<ProviderConfig> createFromText(String text) {
        def config = parse(text)
        def result = createFromMap(config)
        return result
    }

    @PackageScope
    static Map getFromFile(File file) {
        try {
            parse(file.text)
        }
        catch( Exception e ) {
            def message = "Failed to parse config file: $file -- cause: ${e.message?:e.toString()}"
            throw new ConfigParseException(message,e)
        }
    }

    @PackageScope
    static Map getDefault() {
        def file = Const.APP_HOME_DIR.resolve('scm').toFile()
        return file.exists() ? getFromFile(file) : [:]
    }

    static private void addDefaults(List<ProviderConfig> result) {
        if( !result.find{ it.name == 'github' })
            result << new ProviderConfig('github')

        if( !result.find{ it.name == 'gitlab' })
            result << new ProviderConfig('gitlab')

        if( !result.find{ it.name == 'bitbucket' })
            result << new ProviderConfig('bitbucket')
    }

}

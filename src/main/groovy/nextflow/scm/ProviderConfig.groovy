/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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
import nextflow.config.ConfigParser
import nextflow.exception.AbortOperationException
import nextflow.exception.ConfigParseException
/**
 * Models a repository provider configuration attributes
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ProviderConfig {

    @PackageScope
    static public File SCM_FILE = Const.APP_HOME_DIR.resolve('scm').toFile()

    private String name

    private Map attr

    ProviderConfig( String name, Map values ) {
        this.name = name
        this.attr = new HashMap(values)

        // init defaults
        switch( name ) {
            case 'github':
                attr.platform = name
                if( !attr.server ) attr.server = 'https://github.com'
                if( !attr.endpoint ) attr.endpoint = 'https://api.github.com'
                break

            case 'gitlab':
                attr.platform = name
                if( !attr.server ) attr.server = 'https://gitlab.com'
                break

            case 'bitbucket':
                attr.platform = name
                if( !attr.server ) attr.server = 'https://bitbucket.org'
        }

        if( attr.path )
            attr.platform = 'file'

        if( !attr.platform ) {
            throw new AbortOperationException("Missing `platform` attribute for `$name` scm provider configuration -- Check file: ${SCM_FILE}")
        }

        if( attr.auth ) {
            setAuth(attr.auth.toString())
        }
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
     *      Remote server including protocol and eventually port number e.g.
     *      {@code https://github.com}
     */
    String getServer() {
        def result = (String)attr.server
        if( !result ) return null

        if( !result.contains('://') )
            result = 'https://' + result

        // remove ending slash
        while( result.endsWith('/')  )
            result = result.substring(0,result.size()-1)
        return result
    }

    /**
     * @return
     *      Remote server name e.g. {@code github.com}
     */
    String getDomain() {
        def result = server ?: path
        def p = result.indexOf('://')
        p != -1 ? result.substring(p+3) : result
    }

    /**
     * @return The provider identifier e.g. {@code github}, {@code gitlab} or {@code bitbucket}
     */
    String getPlatform() {
        attr.platform
    }

    /**
     * @return The authentication token
     */
    @PackageScope
    String getAuth() {
        def result = new StringBuilder()
        if( user ) {
            result << user
            result << ':'
        }
        if( password )
            result << password

        return result ? result.toString() : null
    }

    @PackageScope
    String getAuthObfuscated() {
        "${user ?: '-'}:${password? '*' * password.size() : '-'}"
    }

    @PackageScope
    ProviderConfig setAuth( String str ) {
        if( str ) {
            def items = str.tokenize(':')
            if( items.size()==1 ) {
                setToken(items[0])
            }
            else if( items.size()==2 ) {
                setUser(items[0])
                setPassword(items[1])
            }
            else if( items.size()>2 ) {
                setUser(items[0])
                setPassword(items[1])
                setToken(items[2])
            }
        }
        return this
    }

    /**
     * @return
     *      The API root endpoint url e.g. {@code http://api.github.com}. Fallback on
     *      the {@link #getServer()} if not defined
     */
    String getEndpoint() {
        attr.endpoint ?: attr.server
    }

    ProviderConfig setUser(String user) {
        attr.user = user
        return this
    }

    ProviderConfig setPassword( String password ) {
        attr.password = password
        return this
    }

    ProviderConfig setToken( String tkn ) {
        attr.token = tkn
        return this
    }

    String getUser() { attr.user }

    String getPassword() { attr.password }

    String getPath() { attr.path?.toString() }

    String getToken() { attr.token }

    String toString() {
        "ProviderConfig[name: $name, platform: $platform, server: $server]"
    }

    @PackageScope
    static Map parse(String text) {
        def slurper = new ConfigParser()
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

    static Map getDefault() {
        def file = SCM_FILE
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

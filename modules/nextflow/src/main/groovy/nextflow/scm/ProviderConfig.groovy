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

            case 'gitea':
                attr.platform = name
                if( !attr.server ) attr.server = 'https://try.gitea.io'
                if( !attr.endpoint ) attr.endpoint = attr.server.toString().stripEnd('/') + '/api/v1'
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
        return result.stripEnd('/')
    }

    /**
     * @return
     *      Remote server name e.g. {@code github.com}
     */
    String getDomain() {
        def result = server ?: path
        def p = result.indexOf('://')
        if( p != -1 )
            result = result.substring(p+3)
        p = result.indexOf('/')
        if( p != -1 )
            result = result.substring(0,p)
        return result
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
        attr.endpoint ? attr.endpoint.toString().stripEnd('/') : server
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

        if( !result.find{ it.name == 'gitea' })
            result << new ProviderConfig('gitea')

        if( !result.find{ it.name == 'bitbucket' })
            result << new ProviderConfig('bitbucket')
    }

}

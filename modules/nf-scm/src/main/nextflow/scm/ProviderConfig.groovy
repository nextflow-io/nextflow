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

package nextflow.scm

import org.apache.commons.lang3.StringUtils
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j

/**
 * Models a repository provider configuration attributes
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ProviderConfig {

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
                if( !attr.server ) attr.server = 'https://gitea.com' // default to free tier
                if( !attr.endpoint ) attr.endpoint = StringUtils.stripEnd(attr.server.toString(),'/') + '/api/v1'
                break

            case 'bitbucket':
                attr.platform = name
                if( !attr.server ) attr.server = 'https://bitbucket.org'
                if( !attr.endpoint ) attr.endpoint = 'https://api.bitbucket.org'
                break

            case 'azurerepos':
                attr.platform = name
                if( !attr.server ) attr.server = 'https://dev.azure.com'
                if( !attr.endpoint ) attr.endpoint = 'https://dev.azure.com'
                break
        }

        if( attr.path )
            attr.platform = 'file'

        if( !attr.platform ) {
            throw new IllegalArgumentException("Missing `platform` attribute for `$name` SCM provider configuration -- Check scm file")
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
        return StringUtils.stripEnd(result, '/')
    }

    /**
     * @return
     *      Remote server name e.g. {@code github.com}
     */
    String getDomain() {
        def result = server ?: path
        if( !result )
            return null
        def p = result.indexOf('://')
        if( p != -1 )
            result = result.substring(p+3)
        // a local file url (e.g. file:///path/to/repo or /path/to/repo)
        // so we need to return the full path as the domain
        if ( result.startsWith('/') )
            return result
        // a server url so we look for the domain without subdirectories
        p = result.indexOf('/')
        if( p != -1 )
            result = result.substring(0, p)
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

    @Deprecated
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
        attr.endpoint ? StringUtils.stripEnd(attr.endpoint.toString(), '/') : server
    }

    ProviderConfig setServer(String serverUrl) {
        attr.server = serverUrl
        return this
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

    static List<ProviderConfig> createDefault() {
        final result = new LinkedList<ProviderConfig>()
        addDefaults(result)
        return result
    }

    @PackageScope
    static void addDefaults(List<ProviderConfig> result) {
        if( !result.find{ it.name == 'github' })
            result << new ProviderConfig('github')

        if( !result.find{ it.name == 'gitlab' })
            result << new ProviderConfig('gitlab')

        if( !result.find{ it.name == 'gitea' })
            result << new ProviderConfig('gitea')

        if( !result.find{ it.name == 'bitbucket' })
            result << new ProviderConfig('bitbucket')

        if( !result.find{ it.name == 'azurerepos' })
            result << new ProviderConfig('azurerepos')

    }
    protected String resolveProjectName(String path) {
        assert path
        assert !path.startsWith('/')

        String project = path
        // fetch prefix from the server url
        def prefix = StringUtils.stripStart(new URL(server).path , '/')
        if( prefix && path.startsWith(prefix) ) {
            project = path.substring(prefix.length())
        }

        if( server == 'https://dev.azure.com' ) {
            project = AzureRepositoryProvider.getUniformPath(project).join('/')
        }

        return StringUtils.stripStart(project, '/')
    }
}

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

import java.nio.file.NoSuchFileException
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.config.ConfigParser
import nextflow.exception.AbortOperationException
import nextflow.exception.ConfigParseException
import nextflow.file.FileHelper
import nextflow.util.StringUtils
/**
 * Models a repository provider configuration attributes
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ProviderConfig {

    @PackageScope
    static Path DEFAULT_SCM_FILE = Const.APP_HOME_DIR.resolve('scm')

    @PackageScope
    static Map<String,String> env = new HashMap<>(System.getenv())

    static Path getScmConfigPath() {
        def cfg = env.get('NXF_SCM_FILE')
        if( !cfg ) {
            log.debug "Using SCM config path: ${DEFAULT_SCM_FILE}"
            return DEFAULT_SCM_FILE
        }
        log.debug "Detected SCM custom path: $cfg"
        // check and return it if valid
        return FileHelper.asPath(cfg)
    }

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
            throw new AbortOperationException("Missing `platform` attribute for `$name` SCM provider configuration -- Check file: ${getScmConfigPath().toUriString()}")
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
        attr.endpoint ? attr.endpoint.toString().stripEnd('/') : server
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

    @PackageScope
    static Map parse(String text) {
        def slurper = new ConfigParser()
        slurper.setBinding(env)
        return slurper.parse(text)
    }

    static List<ProviderConfig> createDefault() {
        createFromMap(null)
    }

    @PackageScope
    static List<ProviderConfig> createFromMap(Map<String,?> config) {

        final providers = (Map<String,?>)config?.providers

        List<ProviderConfig> result = []
        if( providers ) {
            for( String name : providers.keySet() ) {
                def attrs = (Map)providers.get(name)
                result << RepositoryFactory.newProviderConfig(name, attrs)
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
    static Map getFromFile(Path file) {
        try {
            // note: since this can be a remote file via e.g. read via HTTP
            // it should not read more than once
            final content = file.text
            dumpScmContent(file, content)
            final result = parse(content)
            dumpConfig(result)
            return result
        }
        catch (NoSuchFileException | FileNotFoundException e) {
            if( file == DEFAULT_SCM_FILE ) {
                return new LinkedHashMap()
            }
            else
                throw new AbortOperationException("Missing SCM config file: ${file.toUriString()} - Check the env variable NXF_SCM_FILE")
        }
        catch ( UnknownHostException e ) {
            final message = "Unable to access config file '${file?.toUriString()}' -- Unknown host: ${e}"
            throw new ConfigParseException(message,e)
        }
        catch( IOException e ) {
            final message = "Unable to access config file '${file?.toUriString()}' -- Cause: ${e.message?:e.toString()}"
            throw new ConfigParseException(message,e)
        }
        catch( Exception e ) {
            final message = "Failed to parse config file '${file?.toUriString()}' -- Cause: ${e.message?:e.toString()}"
            throw new ConfigParseException(message,e)
        }
    }

    static private void dumpScmContent(Path file, String content) {
        try {
            log.trace "Parsing SCM config path: ${file.toUriString()}\n${StringUtils.stripSecrets(content)}\n"
        }catch(Exception e){
            log.debug "Error dumping configuration ${file.toUriString()}", e
        }
    }

    static private void dumpConfig(Map config) {
        try {
            log.debug "Detected SCM config: ${StringUtils.stripSecrets(config).toMapString()}"
        }
        catch (Throwable e) {
            log.debug "Failed to dump SCM config: ${e.message ?: e}", e
        }
    }

    static Map getDefault() {
        final file = getScmConfigPath()
        return getFromFile(file)
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

        if( !result.find{ it.name == 'azurerepos' })
            result << new ProviderConfig('azurerepos')

    }

    protected String resolveProjectName(String path) {
        assert path
        assert !path.startsWith('/')

        String project = path
        // fetch prefix from the server url
        def prefix = new URL(server).path?.stripStart('/')
        if( prefix && path.startsWith(prefix) ) {
            project = path.substring(prefix.length())
        }

        if( server == 'https://dev.azure.com' ) {
            project = AzureRepositoryProvider.getUniformPath(project).join('/')
        }

        return project.stripStart('/')
    }
}

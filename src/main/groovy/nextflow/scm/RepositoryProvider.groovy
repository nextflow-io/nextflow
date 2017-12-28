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
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
/**
 *
 * Base class for a generic source repository provider
 *
 * @author Maria Chatzou
 * @author Paolo Di Tommaso
 */
@Slf4j
@CompileStatic
abstract class RepositoryProvider {

    /**
     * The pipeline qualified name following the syntax {@code owner/repository}
     */
    protected String project

    /**
     * Holds the provider configuration
     */
    protected ProviderConfig config

    RepositoryProvider setCredentials(String userName, String password) {
        config.user = userName
        config.password = password
        return this
    }

    boolean hasCredentials() {
        config ? config.user && config.password : false
    }

    String getUser() { config?.user }

    String getPassword() { config?.password }

    /**
     * @return The name of the source hub service e.g. github or bitbucket
     */
    abstract String getName()

    /**
     * The hub service API url to access the remote repository e.g. https://api.github.com/repos/nextflow-io/hello
     *
     * @return The repository API url string
     */
    abstract String getEndpointUrl()

    /**
     * The hub service API rul to access a content in the remote repository
     * e.g. https://api.github.com/repos/nextflow-io/hello/contents/my/file
     *
     * @param path The path relative to the remote repository of the content to be retrieved
     * @return The API url string to access the remote content
     */
    abstract String getContentUrl(String path)

    /**
     * @return The repository URL used by Git to clone the remote repository
     *      e.g. https://github.com/nextflow-io/hello.git
     */
    abstract String getCloneUrl()

    /**
     * @return The project home page e.g. https://github.com/nextflow-io/hello
     */
    abstract String getRepositoryUrl()

    /**
     * Invoke the API request specified
     *
     * @param api A API request url e.g. https://api.github.com/repos/nextflow-io/hello
     * @return The remote service response as a text
     */
    protected String invoke( String api ) {
        assert api

        log.debug "Request [credentials ${config.getAuthObfuscated() ?: '-'}] -> $api"
        def connection = new URL(api).openConnection() as URLConnection
        connection.setConnectTimeout(5_000)

        auth(connection)

        if( connection instanceof HttpURLConnection ) {
            checkResponse(connection)
        }

        InputStream content = connection.getInputStream()
        try {
            return content.text
        }
        finally{
            content?.close()
        }
    }

    /**
     * Sets the authentication credential on the connection object
     *
     * @param connection The URL connection object to be authenticated
     */
    protected void auth( URLConnection connection ) {
        if( hasCredentials() ) {
            String authString = "${config.user}:${config.password}".bytes.encodeBase64().toString()
            connection.setRequestProperty("Authorization","Basic " + authString)
        }

    }

    /**
     * Check for response error status. Throws a {@link AbortOperationException} exception
     * when a 401 or 403 error status is returned.
     *
     * @param connection A {@link HttpURLConnection} connection instance
     */
    protected checkResponse( HttpURLConnection connection ) {
        def code = connection.getResponseCode()

        switch( code ) {
            case 401:
                log.debug "Response status: $code -- ${connection.getErrorStream().text}"
                throw new AbortOperationException("Not authorized -- Check that the ${name.capitalize()} user name and password provided are correct")

            case 403:
                log.debug "Response status: $code -- ${connection.getErrorStream().text}"
                def limit = connection.getHeaderField('X-RateLimit-Remaining')
                if( limit == '0' ) {
                    def message = config.auth ? "Check ${name.capitalize()}'s API rate limits for more details" : "Provide your ${name.capitalize()} user name and password to get a higher rate limit"
                    throw new AbortOperationException("API rate limit exceeded -- $message")
                }
                else {
                    def message = config.auth ? "Check that the ${name.capitalize()} user name and password provided are correct" : "Provide your ${name.capitalize()} user name and password to access this repository"
                    throw new AbortOperationException("Forbidden -- $message")
                }
        }

    }


    /**
     * Invoke the API request specified and parse the JSON response
     *
     * @param request A API request url e.g. https://request.github.com/repos/nextflow-io/hello
     * @return The remote service response parsed to a {@link Map} object
     */
    @Memoized
    protected Map invokeAndParseResponse( String request ) {

        def response = invoke(request)
        return new JsonSlurper().parseText(response) as Map

    }

    /**
     * Read the content of a file stored in the remote repository
     *
     * @param path The relative path of a file stored in the repository
     * @return The file content a
     */
    abstract protected byte[] readBytes( String path )

    String readText( String path ) {
        def bytes = readBytes(path)
        return bytes ? new String(bytes) : null
    }


    /**
     * Validate the repository for the specified file.
     * It tries to read the content of the specified file throwing an {@link AbortOperationException} if it does not exist
     *
     * @param scriptName The path of the script to check
     */
    void validateFor( String scriptName ) {

        try {
            readBytes(scriptName)
        }
        catch( IOException e1 ) {

            try {
                invokeAndParseResponse( getEndpointUrl() )
            }
            catch( IOException e2 ) {
                throw new AbortOperationException("Cannot find `$project` -- Make sure exists a ${name.capitalize()} repository at this address `${getRepositoryUrl()}`", e2)
            }

            throw new AbortOperationException("Not a valid Nextflow project -- The repository `${getRepositoryUrl()}` must contain a the script `${AssetManager.DEFAULT_MAIN_FILE_NAME}` or the file `${AssetManager.MANIFEST_FILE_NAME}`", e1)
        }

    }

    /**
     * Factory method
     *
     * @param provider
     * @return
     */
    static RepositoryProvider create( ProviderConfig config, String project ) {
        switch(config.platform) {
            case 'github':
                return new GithubRepositoryProvider(project, config)

            case 'bitbucket':
                return new BitbucketRepositoryProvider(project, config)

            case 'gitlab':
                return new GitlabRepositoryProvider(project, config)

            case 'file':
                // remove the 'local' prefix for the file provider
                def localName = project.tokenize('/').last()
                return new LocalRepositoryProvider(localName, config)
        }

        throw new AbortOperationException("Unkwnon project repository platform: ${config.platform}")
    }

}




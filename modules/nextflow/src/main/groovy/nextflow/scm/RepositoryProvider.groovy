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

import static nextflow.util.StringUtils.*

import groovy.json.JsonSlurper
import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.exception.AbortOperationException
import nextflow.exception.RateLimitExceededException
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.lib.Ref
import org.eclipse.jgit.transport.CredentialsProvider
import org.eclipse.jgit.transport.UsernamePasswordCredentialsProvider
/**
 *
 * Base class for a generic source repository provider
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
abstract class RepositoryProvider {

    @Canonical
    static class TagInfo {
        String name
        String commitId
    }

    @Canonical
    static class BranchInfo {
        String name
        String commitId
    }

    /**
     * The pipeline qualified name following the syntax {@code owner/repository}
     */
    protected String project

    /**
     * Holds the provider configuration
     */
    protected ProviderConfig config

    /**
     * The name of the commit/branch/tag
     */
    protected String revision

    RepositoryProvider setCredentials(String userName, String password) {
        config.user = userName
        config.password = password
        return this
    }

    RepositoryProvider setRevision(String revision) {
        this.revision = revision
        return this
    }

    boolean hasCredentials() {
        getUser() && getPassword()
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

    @Memoized
    protected Collection<Ref> fetchRefs() {
        /*
         * fetch repos tags & branches
         * see https://github.com/centic9/jgit-cookbook/
         */
        return Git.lsRemoteRepository()
                .setRemote(getEndpointUrl())
                .setCredentialsProvider(getGitCredentials())
                .call()
    }

    List<BranchInfo> getBranches() {
        final PREFIX = 'refs/heads/'
        final refs = fetchRefs()
        final result = new ArrayList<BranchInfo>()
        for( Ref it : refs ) {
            if( !it.name.startsWith(PREFIX) )
                continue
            result.add( new BranchInfo(it.name.substring(PREFIX.size()), it.objectId.name()) )
        }
        return result
    }

    List<TagInfo> getTags() {
        final PREFIX = 'refs/tags/'
        final refs = fetchRefs()
        final result = new ArrayList<TagInfo>()
        for( Ref it : refs ) {
            if( !it.name.startsWith(PREFIX) )
                continue
            result.add( new TagInfo(it.name.substring(PREFIX.size()), it.objectId.name()) )
        }
        return result
    }

    /**
     * @return a org.eclipse.jgit.transport.CredentialsProvider object for authenticating git operations
     * like clone, fetch, pull, and update
     **/
    CredentialsProvider getGitCredentials() {
        return new UsernamePasswordCredentialsProvider(getUser(), getPassword())
    }

    /**
     * Invoke the API request specified
     *
     * @param api A API request url e.g. https://api.github.com/repos/nextflow-io/hello
     * @return The remote service response as a text
     */
    protected String invoke( String api ) {
        assert api

        log.debug "Request [credentials ${getAuthObfuscated() ?: '-'}] -> $api"
        def connection = new URL(api).openConnection() as URLConnection
        connection.setConnectTimeout(5_000)

        auth(connection)

        if( connection instanceof HttpURLConnection ) {
            checkResponse(connection)
        }

        InputStream content = connection.getInputStream()
        try {
            final result = content.text
            log.trace "Git provider HTTP request: '$api' -- Response:\n${result}"
            return result
        }
        finally{
            content?.close()
        }
    }

    protected String getAuthObfuscated() {
        final usr = getUser()
        final pwd = getPassword()
        return "${usr ? redact(usr) : '-'}:${pwd ? redact(pwd) : '-'}"
    }

    /**
     * Sets the authentication credential on the connection object
     *
     * @param connection The URL connection object to be authenticated
     */
    protected void auth( URLConnection connection ) {
        if( hasCredentials() ) {
            String authString = "${getUser()}:${getPassword()}".bytes.encodeBase64().toString()
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
                log.debug "Response status: $code -- ${connection.getErrorStream()?.text}"
                throw new AbortOperationException("Not authorized -- Check that the ${name.capitalize()} user name and password provided are correct")

            case 403:
                log.debug "Response status: $code -- ${connection.getErrorStream()?.text}"
                def limit = connection.getHeaderField('X-RateLimit-Remaining')
                if( limit == '0' ) {
                    def message = config.auth ? "Check ${name.capitalize()}'s API rate limits for more details" : "Provide your ${name.capitalize()} user name and password to get a higher rate limit"
                    throw new RateLimitExceededException("API rate limit exceeded -- $message")
                }
                else {
                    def message = config.auth ? "Check that the ${name.capitalize()} user name and password provided are correct" : "Provide your ${name.capitalize()} user name and password to access this repository"
                    throw new AbortOperationException("Forbidden -- $message")
                }
            case 404:
                log.debug "Response status: $code -- ${connection.getErrorStream()?.text}"
                throw new AbortOperationException("Remote resource not found: ${connection.getURL()}")
        }
    }

    protected <T> List<T> invokeAndResponseWithPaging(String request, Closure<T> parse) {
        // this is needed because apparently bytebuddy used by testing framework is not able
        // to handle properly this method signature using both generics and `@Memoized` annotation.
        // therefore the `@Memoized` has been moved to the inner method invocation
        return invokeAndResponseWithPaging0(request, parse)
    }

    @Memoized
    protected List invokeAndResponseWithPaging0(String request, Closure parse) {
        int page = 0
        final result = new ArrayList()
        while( true ) {
            final url = request + (request.contains('?') ? "&page=${++page}": "?page=${++page}")
            final response = invoke(url)
            final list = (List) new JsonSlurper().parseText(response)
            if( !list )
                break

            for( def item : list ) {
                final entry = parse(item)
                if( result.contains(entry) ) {
                    log.debug("Duplicate entry detected on request '$request'")
                    return result
                }
                result.add(entry)
            }

            // prevent endless looping
            if( page==100 ) {
                log.warn("Too many requests '$request'")
                break
            }
        }
        return result
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
            validateRepo()
            throw new AbortOperationException("Not a valid Nextflow project -- The repository `${getRepositoryUrl()}` must contain a `${Const.DEFAULT_MAIN_FILE_NAME}` script or the file `${Const.MANIFEST_FILE_NAME}`", e1)
        }
    }

    void validateRepo() {
        try {
            invokeAndParseResponse( getEndpointUrl() )
        }
        catch( IOException e ) {
            throw new AbortOperationException("Cannot find `$project` -- Make sure exists a ${name.capitalize()} repository at this address `${getRepositoryUrl()}`", e)
        }
    }

}

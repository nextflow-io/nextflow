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

import java.nio.channels.UnresolvedAddressException
import java.time.temporal.ChronoUnit
import java.util.function.Predicate

import static nextflow.util.StringUtils.*

import dev.failsafe.Failsafe
import dev.failsafe.FailsafeException
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import dev.failsafe.function.CheckedSupplier
import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.time.Duration
import java.util.concurrent.Executors

import groovy.json.JsonSlurper
import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.exception.AbortOperationException
import nextflow.exception.RateLimitExceededException
import nextflow.util.RetryConfig
import nextflow.util.Threads
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

    static final public String[] EMPTY_ARRAY = new String[0]

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
     * The client used to carry out http requests
     */
    private HttpClient httpClient

    /**
     * The retry options to be used for http requests
     */
    private RetryConfig retryConfig

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

    String getRevision() {
        return this.revision
    }

    RepositoryProvider setRevision(String revision) {
        this.revision = revision
        return this
    }

    RepositoryProvider setRetryConfig(RetryConfig retryConfig) {
        this.retryConfig = retryConfig
        return this
    }

    String getProject() {
        return this.project
    }

    ProviderConfig getConfig() {
        return this.config
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

    protected HttpRequest newRequest(String api) {
        final builder = HttpRequest
            .newBuilder()
            .uri(new URI(api))
        final auth0 = getAuth()
        if( auth0 )
            builder.headers(auth0)
        return builder.GET().build()
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
        final request = newRequest(api)
        // submit the request
        final HttpResponse<String> resp = httpSend(request)
        // check the response code
        checkResponse(resp)
        // return the body as string
        return resp.body()
    }

    protected String getAuthObfuscated() {
        final usr = getUser()
        final pwd = getPassword()
        return "${usr ? redact(usr) : '-'}:${pwd ? redact(pwd) : '-'}"
    }

    /**
     * Define the credentials to be used to authenticate the http request
     *
     * @return
     *      A string array holding the authentication HTTP headers e.g.
     *      {@code [ "Authorization", "Bearer 1234567890"] } or an empty
     *      array when the credentials are not available or provided.
     *      Note: {@code null} is not a valid return value for this method.
     */
    protected String[] getAuth() {
        if( hasCredentials() ) {
            String authString = "${getUser()}:${getPassword()}".bytes.encodeBase64().toString()
            return new String[] { "Authorization", "Basic " + authString }
        }
        return EMPTY_ARRAY
    }

    /**
     * Check for response error status. Throws a {@link AbortOperationException} exception
     * when a 401 or 403 error status is returned.
     *
     * @param response A {@link HttpURLConnection} response instance
     */
    protected checkResponse( HttpResponse<String> response ) {
        final code = response.statusCode()
        if( code==401 ) {
            log.debug "Response status: $code -- ${response.body()}"
            throw new AbortOperationException("Not authorized -- Check that the ${name.capitalize()} user name and password provided are correct")
        }
        if( code==403 ) {
            log.debug "Response status: $code -- ${response.body()}"
            def limit = response.headers().firstValue('X-RateLimit-Remaining').orElse(null)
            if( limit == '0' ) {
                def message = config.auth ? "Check ${name.capitalize()}'s API rate limits for more details" : "Provide your ${name.capitalize()} user name and password to get a higher rate limit"
                throw new RateLimitExceededException("API rate limit exceeded -- $message")
            }
            else {
                def message = config.auth ? "Check that the ${name.capitalize()} user name and password provided are correct" : "Provide your ${name.capitalize()} user name and password to access this repository"
                throw new AbortOperationException("Forbidden -- $message")
            }
        }
        if( code==404 ) {
            log.debug "Response status: $code -- ${response.body()}"
            throw new AbortOperationException("Remote resource not found: ${response.uri()}")
        }
        if( code>=500 ) {
            String msg = "Unexpected HTTP request error '${response.uri()}' [${code}]"
            try {
                final body = response.body()
                if ( body )
                    msg += " - response: ${body}"
            }
            catch (Throwable t) {
                log.debug "Enable to read response body - ${t.message}"
            }
            throw new IOException(msg)
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
        final response = invoke(request)
        return new JsonSlurper().parseText(response) as Map
    }

    /**
     * Read the content of a file stored in the remote repository
     *
     * @param path The relative path of a file stored in the repository
     * @return The file content a
     */
    abstract byte[] readBytes( String path )

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

    /**
     * Creates a retry policy using the configuration specified by {@link RetryConfig}
     *
     * @param cond
     *      A predicate that determines when a retry should be triggered
     * @param handle
     *
     * @return
     *      The {@link dev.failsafe.RetryPolicy} instance
     */
    protected <T> RetryPolicy<T> retryPolicy(Predicate<? extends Throwable> cond, Predicate<T> handle) {
        final listener = new EventListener<ExecutionAttemptedEvent<?>>() {
            @Override
            void accept(ExecutionAttemptedEvent<?> event) throws Throwable {
                def msg = "Git provider connection failure - attempt: ${event.attemptCount}"
                if( event.lastResult !=null )
                    msg += "; response: ${event.lastResult}"
                if( event.lastFailure != null )
                    msg += "; exception: [${event.lastFailure.class.name}] ${event.lastFailure.message}"
                log.debug(msg)
            }
        }
        return RetryPolicy.<T>builder()
            .handleIf(cond)
            .handleResultIf(handle)
            .withBackoff(retryConfig.delay.toMillis(), retryConfig.maxDelay.toMillis(), ChronoUnit.MILLIS)
            .withMaxAttempts(retryConfig.maxAttempts)
            .withJitter(retryConfig.jitter)
            .onRetry(listener as EventListener)
            .build()
    }

    static private final List<Integer> HTTP_RETRYABLE_ERRORS = [429, 500, 502, 503, 504]

    /**
     * Carry out the invocation of the specified action using a retry policy.
     *
     * @param action A {@link dev.failsafe.function.CheckedSupplier} instance modeling the action to be performed in a safe manner
     * @return The result of the supplied action
     */
    protected <T> HttpResponse<T> safeApply(CheckedSupplier action) {
        final retryOnException = ((Throwable e) -> isRetryable(e) || isRetryable(e.cause)) as Predicate<? extends Throwable>
        final retryOnStatusCode = ((HttpResponse<T> resp) -> resp.statusCode() in HTTP_RETRYABLE_ERRORS) as Predicate<HttpResponse<T>>
        final policy = retryPolicy(retryOnException, retryOnStatusCode)
        return Failsafe.with(policy).get(action)
    }

    protected boolean isRetryable(Throwable t) {
        // only retry SocketException and ignore generic IOException
        return t instanceof SocketException && !isCausedByUnresolvedAddressException(t)
    }

    private boolean isCausedByUnresolvedAddressException(Throwable t) {
        if( t instanceof UnresolvedAddressException )
            return true
        if( t.cause==null )
            return false
        else
            return isCausedByUnresolvedAddressException(t.cause)
    }

    protected HttpResponse<String> httpSend(HttpRequest request) {
        if( httpClient==null )
            httpClient = newHttpClient()
        if( retryConfig==null )
            retryConfig = new RetryConfig()
        try {
            safeApply(()-> httpClient.send(request, HttpResponse.BodyHandlers.ofString()))
        }
        catch (FailsafeException e) {
            throw e.cause
        }
    }

    private HttpClient newHttpClient() {
        final builder = HttpClient.newBuilder()
            .version(HttpClient.Version.HTTP_1_1)
            .connectTimeout(Duration.ofSeconds(60))
        // use virtual threads executor if enabled
        if( Threads.useVirtual() )
            builder.executor(Executors.newVirtualThreadPerTaskExecutor())
        // build and return the new client
        return builder.build()
    }

}

/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
 * Author Maria Chatzou
 * Author Paolo Di Tommaso
 */
@Slf4j
@CompileStatic
abstract class RepositoryProvider {

    /**
     * The pipeline qualified name following the syntax {@code owner/repository}
     */
    protected String pipeline

    /**
     * The user name to access a private repository
     */
    protected String user

    /**
     * The password to access a private repository
     */
    protected String pwd

    /**
     * @return The name of the source hub service e.g. github or bitbucket
     */
    abstract String getName()

    /**
     * The hub service API url to access the remote repository e.g. https://api.github.com/repos/nextflow-io/hello
     *
     * @return The repository API url string
     */
    abstract String getRepoUrl()

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
    abstract String getHomePage()

    /**
     * Invoke the API request specified
     *
     * @param api A API request url e.g. https://api.github.com/repos/nextflow-io/hello
     * @return The remote service response as a text
     */
    protected String invoke( String api ) {
        assert api

        log.debug "Request [credentials ${user?:'-'}:${pwd ? pwd.replaceAll('.','*') : '-'}] -> $api"
        def connection = new URL(api).openConnection()

        if( user && pwd ){
            String authString=("$user:$pwd").bytes.encodeBase64().toString()
            connection.setRequestProperty("Authorization","Basic " + authString)
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
     * Invoke the API request specified and parse the JSON response
     *
     * @param api A API request url e.g. https://api.github.com/repos/nextflow-io/hello
     * @return The remote service response parsed to a {@link Map} object
     */
    protected Map invokeAndParseResponse( String api ) {

        def response = invoke(api)
        return new JsonSlurper().parseText(response) as Map

    }

    /**
     * Read the content of a file stored in the remote repository
     *
     * @param path The relative path of a file stored in the repository
     * @return The file content a
     */
    abstract protected byte[] readContent( String path )

    /**
     * Read the content of a manifest file stored in the remote repository
     *
     * @param manifestPath Path of the manifest file in the remote repository to be read
     * @return The manifest parsed as key-value map object
     */
    Map readManifest( String manifestPath ) {
        InputStream input = null
        try {
            def bytes = readContent(manifestPath)
            def result = new Properties()
            result.load( new ByteArrayInputStream(bytes) )
            return result
        }
        catch( IOException e ) {
            log.debug "Unable to read manifest file: $manifestPath"
            return null
        }
        finally {
            input?.close()
        }
    }

    /**
     * Validate the repository for the specified file.
     * It tries to read the content of the specified file throwing an {@link AbortOperationException} if it does not exist
     *
     * @param scriptName The path of the script to check
     */
    void validateFor( String scriptName ) {

        try {
            readContent(scriptName)
        }
        catch( IOException e1 ) {

            if( e1.message?.startsWith('Server returned HTTP response code: 401'))
                throw new AbortOperationException("Not authorized -- Check that user name and password are correct")

            try {
                invokeAndParseResponse( getRepoUrl() )
            }
            catch( IOException e2 ) {
                throw new AbortOperationException("Cannot find $pipeline pipeline -- Make sure exists a ${getName()} repository at ${getHomePage()}")
            }

            throw new AbortOperationException("Illegal pipeline repository ${getHomePage()} -- It must contain a script named '${AssetManager.DEFAULT_MAIN_FILE_NAME}' or a file '${AssetManager.MANIFEST_FILE_NAME}'")
        }

    }


}


/**
 * Implements a repository provider for GitHub service
 *
 * Author Maria Chatzou
 * Author Paolo Di Tommaso
 */
@CompileStatic
final class GithubRepositoryProvider extends RepositoryProvider{

    /** {@inheritDoc} */
    @Override
    String getName() { "GitHub" }

    /** {@inheritDoc} */
    @Override
    String getRepoUrl() {
        "https://api.github.com/repos/${pipeline}"
    }

    /** {@inheritDoc} */
    @Override
    String getContentUrl( String path ) {
        "https://api.github.com/repos/$pipeline/contents/$path"
    }

    /** {@inheritDoc} */
    @Override
    String getCloneUrl() {
        Map response = invokeAndParseResponse( getRepoUrl() )

        def result = response.get('clone_url')
        if( !result )
            throw new IllegalStateException("Missing clone URL for: $pipeline")

        return result
    }

    /** {@inheritDoc} */
    @Override
    String getHomePage() {
        "https://github.com/$pipeline"
    }

    /** {@inheritDoc} */
    @Override
    byte[] readContent(String path) {

        def url = getContentUrl(path)
        Map response  = invokeAndParseResponse(url)
        response.get('content')?.toString()?.decodeBase64()

    }
}


/**
 * Implements a repository provider for the BitBucket service
 *
 * Author Maria Chatzou
 * Author Paolo Di Tommaso
 */
@Slf4j
final class BitbucketRepositoryProvider extends RepositoryProvider {

    /** {@inheritDoc} */
    @Override
    String getName() { "BitBucket" }

    @Override
    String getRepoUrl() {
        return "https://bitbucket.org/api/2.0/repositories/${pipeline}"
    }

    @Override
    String getContentUrl( String path ) {
        return "https://bitbucket.org/api/1.0/repositories/$pipeline/src/${getMainBranch()}/$path"
    }

    private String getMainBranchUrl() {
        "https://bitbucket.org/api/1.0/repositories/$pipeline/main-branch/"
    }

    @Memoized
    String getMainBranch() {
        invokeAndParseResponse(getMainBranchUrl()) ?. name
    }

    @Override
    String getCloneUrl() {
        Map response = invokeAndParseResponse( getRepoUrl() )

        if( response?.scm != "git" ){
            throw new AbortOperationException("Bitbucket repository at ${getHomePage()} is not supporting Git")
        }

        def result = response?.links?.clone?.find{ it.name == "https" } as Map
        if( !result )
            throw new IllegalStateException("Missing clone URL for: $pipeline")

        return result.href
    }

    @Override
    String getHomePage() {
        return "https://bitbucket.org/$pipeline"
    }

    @Override
    byte[] readContent(String path) {

        def url = getContentUrl(path)
        Map response  = invokeAndParseResponse(url)
        response.get('data')?.toString()?.getBytes()

    }
}
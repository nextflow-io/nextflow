package nextflow.scm

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException

/**
 * Created by mchatzou on 8/14/14.
 */
@Slf4j
@CompileStatic
abstract class RepositoryProvider {

    protected String pipeline

    protected String user

    protected String password

    protected String branch = "master"

    abstract String getName()

    abstract String getRepoURL()

    abstract String getContentURL(String path)

    abstract String getCloneURL()

    abstract String getHomePage()


    protected String invoke( String api ) {
        assert api

        def connection = new URL(api).openConnection()

        if( user && password ){
            String authString=("$user:$password").bytes.encodeBase64().toString()
            connection.setRequestProperty("Authorization","Basic " + authString)
        }

        return connection.getInputStream().text
    }


    protected Map invokeAndParseResponse( String api ) {

        def response = invoke(api)
        return new JsonSlurper().parseText(response) as Map

    }

    abstract String readContent( String path )


    Map readManifest( String manifestFileName ) {
        InputStream input = null
        try {
            def bytes = readContent(manifestFileName).bytes
            def result = new Properties()
            result.load( new ByteArrayInputStream(bytes) )
            return result
        }
        catch( IOException e ) {
            log.debug "Unable to read manifest file: $manifestFileName"
            return null
        }
        finally {
            input?.close()
        }
    }

    void validateFor( String scriptName ) {

        try {
            readContent(scriptName)
        }
        catch( IOException e1 ) {

            if( e1.message?.startsWith('Server returned HTTP response code: 401'))
                throw new AbortOperationException("Not authorized -- Check that user name and password are correct")

            try {
                invokeAndParseResponse( getRepoURL() )
            }
            catch( IOException e2 ) {
                throw new AbortOperationException("Cannot find $pipeline pipeline -- Make sure exists a ${getName()} repository at http://github.com/$pipeline")
            }

            throw new AbortOperationException("Illegal pipeline repository ${getHomePage()} -- It must contain a script named '${AssetManager.DEFAULT_MAIN_FILE_NAME}' or a file '${AssetManager.MANIFEST_FILE_NAME}'")
        }

    }


}



@CompileStatic
final class GithubRepositoryProvider extends RepositoryProvider{

    @Override
    String getName() { "GitHub" }

    @Override
    String getRepoURL() {
        "https://api.github.com/repos/${pipeline}"
    }

    @Override
    String getContentURL( String path ) {
        "https://api.github.com/repos/$pipeline/contents/$path"
    }

    @Override
    String getCloneURL() {
        Map response = invokeAndParseResponse( getRepoURL() )

        def result = response.get('clone_url')
        if( !result )
            throw new IllegalStateException("Missing clone URL for: $pipeline")

        return result
    }

    @Override
    String getHomePage() {
        "https://github.com/$pipeline"
    }

    @Override
    String readContent(String path) {

        def url = getContentURL(path)
        Map response  = invokeAndParseResponse(url)
        new String(response.get('content')?.toString()?.decodeBase64())

    }
}



final class BitbucketRepositoryProvider extends RepositoryProvider {

    @Override
    String getName() { "BitBucket" }

    @Override
    String getRepoURL() {
        return "https://bitbucket.org/api/2.0/repositories/${pipeline}"
    }

    @Override
    String getContentURL( String path ) {
        return "https://bitbucket.org/api/1.0/repositories/$pipeline/src/$branch/$path"
    }

    @Override
    String getCloneURL() {
        Map response = invokeAndParseResponse( getRepoURL() )

        if( response?.scm != "git" ){
            throw new AbortOperationException("Bitbucket repository at ${getHomePage()} is not supporting Git")
        }

        def result = response?.links?.clone?.find{ it.name == "https" }
        if( !result )
            throw new IllegalStateException("Missing clone URL for: $pipeline")

        return result.href
    }

    @Override
    String getHomePage() {
        return "https://bitbucket.org/$pipeline"
    }

    @Override
    String readContent(String path) {

        def url = getContentURL(path)
        Map response  = invokeAndParseResponse(url)
        response.get('data')?.toString()

    }
}
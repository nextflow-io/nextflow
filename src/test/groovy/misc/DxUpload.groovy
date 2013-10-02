package misc

import com.dnanexus.DXAPI
import com.dnanexus.DXEnvironment
import com.dnanexus.DXJSON
import com.fasterxml.jackson.databind.node.ObjectNode
import org.apache.http.HttpVersion
import org.apache.http.client.HttpClient
import org.apache.http.client.methods.HttpPost
import org.apache.http.entity.FileEntity
import org.apache.http.impl.client.DefaultHttpClient
import org.apache.http.params.CoreProtocolPNames
import org.apache.http.util.EntityUtils

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DXUpload {

    public static void main(String[] args) {

        println "security: " +  DXEnvironment.create().securityContext.toString()

        def projectId = System.getenv('DX_PROJECT_CONTEXT_ID')
        println "current project >> ${projectId}"

        File fileToUpload = new File('my-upload')
        fileToUpload.text = 'Hello there'

        def newFileRequest = mapToJsonNode([name: 'my-upload', project: projectId ?: 'project-B7Z6B700j58P8Y43B1vQ03gf'] )
        println ">> request:  ${newFileRequest.toString()}"

        def result = DXAPI.fileNew(  newFileRequest )
        println ">> fileNew: ${result.toString()}"

        def fileId = result.get('id').textValue()
        def upload = DXAPI.fileUpload(fileId)
        println ">> fileUpload: ${upload.toString()}"

        def url = upload.get('url').textValue()
        def headers = upload.get('headers')
        println "headers>>>\n" + headers.toString()

        def auth = headers.get('Authorization').textValue()

        HttpClient client = new DefaultHttpClient();
        client.getParams().setParameter(CoreProtocolPNames.PROTOCOL_VERSION, HttpVersion.HTTP_1_1);

        HttpPost post = new HttpPost( url );

        post.setHeader('Authorization', auth)

        FileEntity entity = new FileEntity(fileToUpload)
        post.setEntity(entity)

        String response = EntityUtils.toString( client.execute( post ).getEntity(), "UTF-8" );
        println "post >> ${response}"

        def close = DXAPI.fileClose(fileId)
        println ">> fileClose ${close.toString()}"



        //client.getConnectionManager().shutdown();


    }


    static ObjectNode mapToJsonNode( Map params ) {
        assert params != null

        def node = DXJSON.getObjectBuilder()
        params.each { String name, value ->

            switch( value ) {
                case boolean:
                    node = node.put(name, value.asBoolean() );
                    break

                case Integer:
                    node = node.put(name, value as Integer)
                    break

                case Long:
                    node = node.put(name, value as Long)
                    break

                case double:
                    node = node.put(name, value as Double)
                    break

                case Map:
                    node = node.put(name, mapToJsonNode(value as Map))
                    break

                default:
                    node = node.put(name, value?.toString() )

            }

        }

        node.build()
    }



}
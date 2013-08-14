package nextflow.util

import com.dnanexus.DXAPI
import com.dnanexus.DXJSON
import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.core.JsonParser
import com.fasterxml.jackson.core.JsonProcessingException
import com.fasterxml.jackson.core.JsonToken
import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.SerializerProvider
import com.fasterxml.jackson.databind.node.ArrayNode
import com.fasterxml.jackson.databind.node.BigIntegerNode
import com.fasterxml.jackson.databind.node.BinaryNode
import com.fasterxml.jackson.databind.node.BooleanNode
import com.fasterxml.jackson.databind.node.DecimalNode
import com.fasterxml.jackson.databind.node.DoubleNode
import com.fasterxml.jackson.databind.node.IntNode
import com.fasterxml.jackson.databind.node.JsonNodeFactory
import com.fasterxml.jackson.databind.node.LongNode
import com.fasterxml.jackson.databind.node.NullNode
import com.fasterxml.jackson.databind.node.NumericNode
import com.fasterxml.jackson.databind.node.ObjectNode
import com.fasterxml.jackson.databind.node.POJONode
import com.fasterxml.jackson.databind.node.TextNode
import groovy.util.logging.Slf4j
import org.apache.http.HttpVersion
import org.apache.http.client.HttpClient
import org.apache.http.client.methods.HttpGet
import org.apache.http.client.methods.HttpPost
import org.apache.http.entity.FileEntity
import org.apache.http.impl.client.DefaultHttpClient
import org.apache.http.params.CoreProtocolPNames
import org.apache.http.util.EntityUtils

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
class DxHelper {

    static void downloadFile( String fileId, File targetFile ) {

        def download = DXAPI.fileDownload(fileId)
        def url = download.get('url').textValue()
        def headers = download.get('headers')

        log.debug "download headers>>>\n" + headers.toString()

        HttpClient client = new DefaultHttpClient();
        client.getParams().setParameter(CoreProtocolPNames.PROTOCOL_VERSION, HttpVersion.HTTP_1_1);

        def get = new HttpGet(url)
        for( Map.Entry<String,JsonNode> item : headers.fields() ) {
            log.debug "setting header > ${item.key}: ${item.value.textValue()}"
            get.setHeader( item.key, item.value.textValue())
        }

        def commandOutFile = targetFile
        BufferedOutputStream buffer = new BufferedOutputStream(new FileOutputStream(commandOutFile))
        try {
            client.execute(get).getEntity().writeTo(buffer)
        }
        finally {
            buffer.close()
        }

    }

    /**
     * Uploads the defined file to the Dnanexus' space
     *
     * @param fileToUpload
     * @param fileName
     * @return fileId
     */
    static String uploadFile( File fileToUpload, String fileName) {

        def projectId = System.getenv('DX_PROJECT_CONTEXT_ID')
        log.debug  "Current project >> ${projectId}"

        def newFileRequest = mapToJsonNode([name: fileName, project: projectId ?: 'project-B7fQ9vj0FqXB2z80y5FQ0JGG'] )
        println "Request >> ${newFileRequest.toString()}"

        def result = DXAPI.fileNew(  newFileRequest )
        log.debug "FileNew >> ${result.toString()}"

        def fileId = result.get('id').textValue()
        def upload = DXAPI.fileUpload(fileId)
        log.debug "FileUpload >> ${upload.toString()}"

        def url = upload.get('url').textValue()
        def headers = upload.get('headers')
        log.debug "Headers >>>\n" + headers.toString()

        def auth = headers.get('Authorization').textValue()

        HttpClient client = new DefaultHttpClient();
        client.getParams().setParameter(CoreProtocolPNames.PROTOCOL_VERSION, HttpVersion.HTTP_1_1);

        HttpPost post = new HttpPost( url );

        post.setHeader('Authorization', auth)

        FileEntity entity = new FileEntity(fileToUpload)
        post.setEntity(entity)

        String response = EntityUtils.toString( client.execute( post ).getEntity(), "UTF-8" );
        log.debug "Post >> ${response}"

        def close = DXAPI.fileClose(fileId)
        log.debug "FileClose >> ${close.toString()}"

        return fileId
    }


//    /**
//     * Converts a {@code Map} to a {@code ObjectNode} instance
//     *
//     * @param params
//     * @return
//     */
//    static ObjectNode mapToJsonNode( Map params ) {
//        assert params != null
//
//        def node = DXJSON.getObjectBuilder()
//        params.each { String name, value ->
//
//            switch( value ) {
//
//                case String:
//                    node = node.put(name, value as String)
//                    break
//
//                case BigDecimal:
//                    node = node.put(name, new DecimalNode(value as BigDecimal))
//                    break
//
//                case BigInteger:
//                    node = node.put(name, new BigIntegerNode(value as BigInteger))
//                    break
//
//                case byte[]:
//                    node = node.put( name, new BinaryNode(value as byte[]) )
//                    break
//
//                case Boolean:
//                    node = node.put(name, value as Boolean)
//                    break
//
//                case Float:
//                case Double:
//                    node = node.put(name, value as Double)
//                    break
//
//                case Integer:
//                    node = node.put(name, value as Integer)
//                    break
//
//                case Long:
//                    node = node.put(name, value as Long)
//                    break
//
//                case Map:
//                    node = node.put(name, mapToJsonNode(value as Map))
//                    break
//
//                case Collection:
//                case Object[]:
//
//                    def array = new ArrayNode(JsonNodeFactory.instance)
//                    (value as List) .each {
//                        array.add( it )
//                    }
//                    node = node.put(name, array)
//                    break
//
//
//                default:
//                    log.debug "Type not supported: ${value.class.name}"
//                    node = node.put(name, value?.toString() )
//
//            }
//
//        }
//
//        node.build()
//    }

    static JsonNode mapToJsonNode( def value ) {

        switch( value ) {

            case String:
                return TextNode.valueOf(value as String)

            case BigDecimal:
                return DecimalNode.valueOf(value as BigDecimal)


            case Integer:
                return IntNode.valueOf(value as int)

            case BigInteger:
                return BigIntegerNode.valueOf(value as BigInteger)

            case byte[]:
                return BinaryNode.valueOf(value as byte[])

            case Boolean:
                return BooleanNode.valueOf(value as Boolean)

            case Float:
            case Double:
                return DoubleNode.valueOf(value as Double)

            case Long:
                return LongNode.valueOf(value as Long)

            case Map:
                def result = new ObjectNode(JsonNodeFactory.instance)
                (value as Map).each { String name, Object item ->
                    result.put( name, mapToJsonNode(item))
                }
                return result

            case Collection:
            case Object[]:

                def result = new ArrayNode(JsonNodeFactory.instance)
                (value as List) .each {
                    result.add( mapToJsonNode(it) )
                }
                return result
        }
    }


    static Map jsonNodeToMap( JsonNode node ) {
        def result = [:]

        if( node.isObject() ) {
            for( Map.Entry<String,JsonNode> item : node.fields() ) {
                jsonNodeToMap(item.key, item.value, result)
            }
        }
        else if( node.isArray() ) {

        }
        else {
            throw new IllegalStateException("Not a valid json node: ${node.toString()}")
        }

        return result
    }


    static Map jsonNodeToMap( String field, JsonNode node, Map result ) {

        if( node.isArray()){
            for( def item : node.elements() ) {
                result[field] = item
            }
        }
        else if( node.isBigDecimal() ) {
            result[field] = node.decimalValue()
        }
        else if( node.isBigInteger() ) {
            result[field] = node.bigIntegerValue()
        }
        else if( node.isBinary() ) {
            result[field] = node.binaryValue()
        }
        else if( node.isBoolean() ) {
            result[field] = node.booleanValue()
        }
        else if( node.isDouble() ) {
            result[field] = node.doubleValue()
        }
        else if( node.isFloatingPointNumber()) {
            result[field] = node
        }
        else if( node.isInt() ) {
            result[field] = node.intValue()
        }
        else if( node.isLong()) {
            result[field] = node.longValue()
        }
        else if( node.isNumber()) {
            result[field] = node.numberValue()
        }
        else if( node.isObject()) {
            result[field] = jsonNodeToMap(node)
        }
        else if( node.isPojo() ) {
            result[field] = node
        }
        else if( node.isTextual()) {
            result[field] = node.textValue()
        }
        else if( node.isValueNode()) {
            if(node.isContainerNode() || node.isMissingNode()){
                result[field] = false
            }else{
                result[field] = true
            }
        }
        else if( node.isNull()) {
            result[field] = null
        }
        else {

        }

        return result
    }


}

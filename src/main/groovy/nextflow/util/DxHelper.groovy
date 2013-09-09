/*
 * Copyright (c) 2012-2013, the authors.
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
package nextflow.util

import java.nio.ByteBuffer
import java.util.concurrent.ConcurrentLinkedQueue
import java.util.concurrent.atomic.AtomicInteger

import com.dnanexus.DXAPI
import com.fasterxml.jackson.databind.JsonNode
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
import com.fasterxml.jackson.databind.node.ObjectNode
import com.fasterxml.jackson.databind.node.TextNode
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import org.apache.http.HttpVersion
import org.apache.http.client.HttpClient
import org.apache.http.client.methods.HttpGet
import org.apache.http.client.methods.HttpPost
import org.apache.http.entity.InputStreamEntity
import org.apache.http.impl.client.DefaultHttpClient
import org.apache.http.params.CoreProtocolPNames
import org.apache.http.util.EntityUtils
import sun.nio.ch.DirectBuffer

/**
 * Helper methods to interact with DNAnexus API
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Beatriz M. San Juan <bmsanjuan@gmail.com>
 */

@Slf4j
class DxHelper {

    /**
     * Given a DNAnexus 'fileId' download it to the local file system
     *
     * @param fileId The remote file to be downloaded
     * @param targetFile A local {@code File}
     */
    static void downloadFile( String fileId, File targetFile ) {
        assert fileId
        assert targetFile

        log.debug "Dx download: $fileId; target: $targetFile"

        def download = DXAPI.fileDownload(fileId)
        def url = download.get('url').textValue()
        def headers = download.get('headers')

        log.trace "download headers>>>\n" + headers.toString()

        HttpClient client = new DefaultHttpClient();
        client.getParams().setParameter(CoreProtocolPNames.PROTOCOL_VERSION, HttpVersion.HTTP_1_1);

        def get = new HttpGet(url)
        for( Map.Entry<String,JsonNode> item : headers.fields() ) {
            log.trace "setting header > ${item.key}: ${item.value.textValue()}"
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
     * @param fileToUpload The file to be uploaded. The size cannot exceed 5TB
     * @param targetName Name used to store the file in the remote storage. If not specified the source file name is used.
     * @param workspaceId The DNAnexus workspace-id where the file has to be stored. If not specified the environment variable {@code DX_WORKSPACE_ID} is used
     * @return fileId
     */
    static String uploadFile( File fileToUpload, String targetName = null, String workspaceId = null ) {
        assert fileToUpload

        if( !fileToUpload.exists() ) {
            throw new IllegalArgumentException("Missing file: $fileToUpload -- upload failed" )
        }

        if( fileToUpload.size() > new MemoryUnit('5TB').toBytes() ) {
            throw new IllegalArgumentException("File $fileToUpload exceed the maximum allowed upload size (5TB) -- upload failed" )
        }

        // fallback to the source file name
        if( !targetName ) {
            targetName = fileToUpload.name
        }
        def targetFile = new File(targetName)

        // try to access to the current 'DX_WORKSPACE_ID'
        if( !workspaceId ) {
            workspaceId = System.getenv('DX_WORKSPACE_ID')
        }

        def uploadInfo = [project: workspaceId, name: targetFile.name]
        if( targetFile.parent ) {
            // set the parent folder where it must be stored
            uploadInfo.folder = targetFile.parent
            uploadInfo.parents = true
        }
        log.debug "Dx upload info: $uploadInfo"

        /*
         * Create a new remote file
         */
        def newFileRequest = toJsonNode(uploadInfo)
        log.trace "Dx fileNew request > ${newFileRequest.toString()}"

        def result = DXAPI.fileNew(newFileRequest)
        log.trace "Dx fileNew reply > ${result.toString()}"

        /*
         * The file is splitted in multiple chunks which are uploaded in parallel
         */
        final fileId = result.get('id').textValue()

        // defines the max size of each chunk
        // By the API spec if must be a value between 5MB and 5GB
        // See https://wiki.dnanexus.com/API-Specification-v1.0.0/Files#API-method:-/file-xxxx/upload
        final chunkSize = 1024 * 1024 * 100 // 100 MB
        final maxParallel = 5
        final chunks = new DataflowQueue()
        def totChunks = 0
        while( (totChunks*chunkSize) <= fileToUpload.size() ) {
            chunks << totChunks++
        }

        // synchronization vars
        final counter = new AtomicInteger()
        final done = new DataflowVariable()
        final Queue<ByteBuffer> bufferPool = new ConcurrentLinkedQueue<ByteBuffer>()

        // error handler to redirect exception on the application log
        def errorHandler = new DataflowEventAdapter() {
            public boolean onException(final DataflowProcessor processor, final Throwable e) {
                DxHelper.log.error "Failure uploading: $fileToUpload", e
                done << false
                return true;
            }
        }

        /*
         * parallel uploader
         */
        Dataflow.operator( inputs:[chunks], outputs:[], maxForks: maxParallel, listeners: [errorHandler] ) { int current ->

            final fileName = fileToUpload.name
            log.trace "Uploading chunk [$current] for file: ${fileName}"

            // open the file in 'read' mode
            final file = new RandomAccessFile(fileToUpload, 'r')
            final offset = current * chunkSize
            final channel = file.getChannel()

            // get am available buffer or allocate it
            def buffer = bufferPool.poll()
            if( !buffer ) {
                log.trace "File: $fileName; chunk [$current] > allocating buffer"
                buffer = ByteBuffer.allocateDirect(chunkSize)
            }
            else {
                log.trace "File: $fileName; chunk [$current] > clearing buffer"
                buffer.clear()
            }

            // read the chunk to be uploaded using a direct buffer
            def len = channel.read(buffer, offset)
            log.trace "File: $fileName; chunk [$current] > read buffer len: $len "

            // request to upload a new chunk
            // note: dnanexus upload chunk index is 1-based
            def params = toJsonNode( index: current+1 )
            def upload = DXAPI.fileUpload(fileId, params)
            log.trace "File: $fileName; chunk [$current] > FileUpload: ${upload.toString()}"

            // the response provide the url when 'post' the chunk and the
            // 'authorization' code
            def url = upload.get('url').textValue()
            def auth = upload.get('headers')?.get('Authorization')?.textValue()

            // create a 'post' request to upload the stuff
            HttpPost post = new HttpPost(url);
            post.setHeader('Authorization', auth)

            buffer.flip()
            log.trace "File: $fileName; chunk [$current] > buffer remaining: ${buffer.remaining()} "

            def payload = new InputStreamEntity(new ByteBufferBackedInputStream(buffer), len)
            post.setEntity(payload)

            HttpClient client = new DefaultHttpClient();
            client.getParams().setParameter(CoreProtocolPNames.PROTOCOL_VERSION, HttpVersion.HTTP_1_1);
            log.trace "File: $fileName; chunk [$current] > Post starting: $post "

            def entity = client.execute( post ).getEntity()
            String response = EntityUtils.toString( entity, "UTF-8" );
            log.trace "File: $fileName; chunk [$current] > post response: ${response}"

            // close the client (maybe not really necessary)
            client.getConnectionManager().shutdown()
            // offer the buffer so that it can be reused in the next iteration
            bufferPool.offer(buffer)

            // when ALL the chunks have been uploaded, signal the termination
            if( counter.incrementAndGet() == totChunks )  {
                done << true
            }

        }

        /*
         * await the termination
         */
        try {
            if( done.get() )  {
                /*
                 * when finished, close the remote file
                 */
                def close = DXAPI.fileClose(fileId)
                log.debug "Dx fileClose > ${close.toString()}"
                return fileId
            }

            // return null in case of error (!done)
            return null
        }
        finally {
            log.trace "Uploader > release buffers"
            // release the direct buffer
            bufferPool.each { DirectBuffer it ->
                it?.cleaner() ?.clean()
            }
        }

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

    static JsonNode toJsonNode( def value ) {

        switch( value ) {
            case null:
                return NullNode.instance

            case String:
                return TextNode.valueOf(value as String)

            case BigDecimal:
                return DecimalNode.valueOf(value as BigDecimal)

            case Integer:
                return IntNode.valueOf(value as int)

            case BigInteger:
                return BigIntegerNode.valueOf(value as BigInteger)

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
                    result.put( name, toJsonNode(item))
                }
                return result

            case ObjectNode:
                return value as ObjectNode

            case byte[]:
                return BinaryNode.valueOf(value as byte[])

            case Collection:
            case ArrayNode:
            case Object[]:

                def result = new ArrayNode(JsonNodeFactory.instance)
                (value as List) .each {
                    result.add( toJsonNode(it) )
                }
                return result

            default:
                log.debug "Unknown json type: ${value.class.name} -- mapping as string"
                return TextNode.valueOf(value.toString())
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

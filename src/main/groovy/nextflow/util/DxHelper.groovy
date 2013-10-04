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
import com.fasterxml.jackson.databind.ObjectMapper
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

    static private ObjectMapper mapper = new ObjectMapper();

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
        def newFileRequest = objToJson(uploadInfo)
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
            def params = objToJson( index: current+1 )
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


    /**
     * Parse Json formatted string to a {@code JsonNode} object instance
     *
     * @param str The source JSON string
     * @return The resulting {@code JsonNode}
     */
    static JsonNode parseJsonString( String str ) {
        try {
            return mapper.readValue(str, JsonNode.class);
        }
        catch( IOException e ) {
            throw new IllegalArgumentException(String.format("Error parsing json string:\n%s", str), e);
        }

    }


    /**
     * Converts a generic Java object to a {@code JsoNode} object
     * @param obj The source Java object, usually it is a {@code Map}
     * @return The mapped {@code JsonNode} instance
     */
    static JsonNode objToJson( def value ) {

        return mapper.valueToTree(value)

    }


    /**
     * Converts a {@code JsonNode} to a generic Java {@code Map} object
     *
     * @param node The source Json node
     * @return The resulting Map object
     */
    static <T> T jsonToObj( JsonNode node ) {

        return (T)mapper.convertValue(node, Map.class)

    }

}

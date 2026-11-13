/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.cloud.aws.nio

import java.nio.file.Files
import java.nio.file.Path
import java.util.concurrent.Callable
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.Future
import java.util.concurrent.ThreadFactory
import java.util.concurrent.TimeUnit
import java.util.concurrent.TimeoutException
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream

import com.sun.net.httpserver.HttpExchange
import com.sun.net.httpserver.HttpServer
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session
import nextflow.cloud.aws.util.S3PathFactory
import nextflow.file.FileHelper
import software.amazon.awssdk.core.ResponseInputStream
import spock.lang.Shared
import spock.lang.Specification

/**
 * Regression test for {@link S3FileSystemProvider#newInputStream} close-on-partial-read behavior.
 *
 * Before the fix, {@code newInputStream()} returned the raw {@code ResponseInputStream} from the
 * AWS SDK. Closing it without reading to EOF would trigger Apache HTTP client's
 * {@code ContentLengthInputStream.close()}, which drains the remaining response body to release
 * the connection back to the pool. For a multi-GB object this blocked the caller for many
 * minutes. The fix wraps the stream so {@code close()} calls {@code ResponseInputStream.abort()}
 * instead.
 *
 * The object is served by a local HTTP server standing in for S3: it declares a ~1GB body and
 * streams it at a throttled rate, as a large object over a real network would arrive, so
 * draining it would take far longer than the wall-clock bound while aborting returns at once.
 *
 * A dedicated spec is used because the test requires non-trivial orchestration:
 *   - the wall-clock bound must be enforced on the caller side — Spock {@code @Timeout} relies
 *     on {@code Thread.interrupt()}, which does not unblock a thread parked in a socket read;
 *   - when the regression is present the worker thread cannot be stopped by interrupt; the
 *     spec captures the underlying {@link ResponseInputStream} so it can call {@code abort()}
 *     from the test thread on timeout to force-release the HTTP connection.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
class S3InputStreamAbortTest extends Specification {

    static final String BUCKET = 'nf-abort-test'

    static final String KEY = 'SRR7890919_WES_HCC1395BL-EA_normal_1.fastq.gz'

    static final long BODY_SIZE = 1L << 30

    static final long TIMEOUT_SECONDS = 10

    @Shared
    private ExecutorService executor

    @Shared
    private HttpServer server

    @Shared
    private volatile boolean stopped

    def setupSpec() {
        executor = Executors.newSingleThreadExecutor({ Runnable r ->
            def t = new Thread(r, 's3-abort-test-worker')
            t.daemon = true   // so a hung worker cannot keep the JVM alive
            return t
        } as ThreadFactory)

        final head = gzippedFastq()
        server = HttpServer.create(new InetSocketAddress(InetAddress.loopbackAddress, 0), 0)
        server.executor = Executors.newCachedThreadPool({ Runnable r ->
            def t = new Thread(r, 's3-abort-test-server')
            t.daemon = true
            return t
        } as ThreadFactory)
        server.createContext('/') { HttpExchange exchange -> serveObject(exchange, head) }
        server.start()
    }

    def cleanupSpec() {
        stopped = true
        server?.stop(0)
        executor?.shutdownNow()
    }

    def setup() {
        // Anonymous access to the local server, with path-style URLs so the bucket goes in the path
        def cfg = [aws: [client: [anonymous: true, endpoint: "http://localhost:${server.address.port}", s3PathStyleAccess: true], region: 'eu-west-1']]
        FileHelper.getOrCreateFileSystemFor(URI.create("s3:///$BUCKET/"), cfg.aws)
        Global.config = cfg
        Global.session = Mock(Session) { getConfig() >> cfg }
    }

    /**
     * Serve the object: a gzipped FASTQ head, then filler up to {@link #BODY_SIZE} at ~1MB/s,
     * until the client closes the connection.
     */
    private void serveObject(HttpExchange exchange, byte[] head) {
        try {
            if( exchange.requestMethod != 'GET' || exchange.requestURI.path != "/$BUCKET/$KEY" ) {
                exchange.sendResponseHeaders(404, -1)
                return
            }
            exchange.responseHeaders.add('Content-Type', 'application/octet-stream')
            exchange.responseHeaders.add('ETag', '"0123456789abcdef"')
            exchange.sendResponseHeaders(200, BODY_SIZE)
            final out = exchange.responseBody
            out.write(head)
            out.flush()
            final chunk = new byte[8 * 1024]
            long sent = head.length
            while( !stopped && sent < BODY_SIZE ) {
                final len = (int) Math.min(chunk.length, BODY_SIZE - sent)
                out.write(chunk, 0, len)
                out.flush()
                sent += len
                sleep(8)
            }
        }
        catch( IOException e ) {
            // the client aborted the connection
        }
        finally {
            exchange.close()
        }
    }

    private static byte[] gzippedFastq() {
        final buffer = new ByteArrayOutputStream()
        new GZIPOutputStream(buffer).withCloseable { gz ->
            for( int i=1; i<=1000; i++ )
                gz.write("@SRR7890919.${i} ${i}/1\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n".getBytes('ASCII'))
        }
        return buffer.toByteArray()
    }

    def 'close on a partially-consumed newInputStream should abort, not drain'() {
        given: 'an S3 path to a large (~1GB) gzipped object'
        final Path path = (Path) S3PathFactory.parse("s3://$BUCKET/$KEY")

        and: 'open the stream on the test thread'
        final InputStream raw = Files.newInputStream(path)

        when: 'read the first line and close on a background thread, bounded by a wall-clock timeout'
        final Future<String> future = executor.submit({
            String line = null
            raw.withCloseable { InputStream is ->    // close() here is the code path under test
                def gz = new GZIPInputStream(is)
                def reader = new BufferedReader(new InputStreamReader(gz, 'ASCII'))
                line = reader.readLine()
            }
            return line
        } as Callable<String>)

        String firstLine
        try {
            firstLine = future.get(TIMEOUT_SECONDS, TimeUnit.SECONDS)
        }
        catch (TimeoutException e) {
            // Thread.interrupt() cannot unblock the socket read — forcibly release the
            // HTTP connection by calling abort() on the underlying ResponseInputStream so the
            // worker thread can exit instead of lingering until the full body has drained.
            log.warn("Timed out after ${TIMEOUT_SECONDS}s waiting for close() — aborting underlying S3 stream")
            raw.abort()
            throw e
        }
        finally {
            future.cancel(true)
        }

        then: 'no timeout occurred and the first FASTQ record identifier was returned'
        noExceptionThrown()
        firstLine == '@SRR7890919.1 1/1'
    }
}

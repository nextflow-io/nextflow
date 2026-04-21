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

import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session
import nextflow.cloud.aws.util.S3PathFactory
import nextflow.file.FileHelper
import software.amazon.awssdk.core.ResponseInputStream
import spock.lang.IgnoreIf
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
 * A dedicated spec is used because the test requires non-trivial orchestration that would
 * clutter {@link AwsS3NioTest}:
 *   - the wall-clock bound must be enforced on the caller side — Spock {@code @Timeout} relies
 *     on {@code Thread.interrupt()}, which does not unblock a thread parked in
 *     {@code NioSocketImpl.timedRead()} on a native SSL read;
 *   - when the regression is present the worker thread cannot be stopped by interrupt; the
 *     spec captures the underlying {@link ResponseInputStream} so it can call {@code abort()}
 *     from the test thread on timeout to force-release the HTTP connection.
 *
 * The test reads the first line of a public ~1GB FASTQ in the {@code ngi-igenomes} bucket
 * (eu-west-1, anonymous). Without the fix the run blows the 30s wall-clock bound; with the fix
 * it completes in seconds.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@IgnoreIf({ System.getenv('NXF_SMOKE') })
class S3InputStreamAbortTest extends Specification {

    static final String PUBLIC_FASTQ =
            's3://ngi-igenomes/test-data/sarek/SRR7890919_WES_HCC1395BL-EA_normal_1.fastq.gz'

    static final long TIMEOUT_SECONDS = 30

    @Shared
    private ExecutorService executor

    def setupSpec() {
        executor = Executors.newSingleThreadExecutor({ Runnable r ->
            def t = new Thread(r, 's3-abort-test-worker')
            t.daemon = true   // so a hung worker cannot keep the JVM alive
            return t
        } as ThreadFactory)
    }

    def cleanupSpec() {
        executor?.shutdownNow()
    }

    def setup() {
        // Anonymous S3 access — ngi-igenomes is public, bucket lives in eu-west-1.
        def cfg = [aws: [client: [anonymous: true], region: 'eu-west-1']]
        FileHelper.getOrCreateFileSystemFor(URI.create('s3:///'), cfg.aws)
        Global.config = cfg
        Global.session = Mock(Session) { getConfig() >> cfg }
    }

    def 'close on a partially-consumed newInputStream should abort, not drain'() {
        given: 'an S3 path to a large (~1GB) gzipped object'
        final Path path = (Path) S3PathFactory.parse(PUBLIC_FASTQ)

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
            // Thread.interrupt() cannot unblock the native SSL read — forcibly release the
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
        firstLine?.startsWith('@')
    }
}

/*
 * Copyright 2020-2022, Seqera Labs
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
 *
 */

package io.seqera.tower.plugin

import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.time.temporal.ChronoUnit
import java.util.concurrent.ExecutorService
import java.util.function.Predicate

import dev.failsafe.Failsafe
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import dev.failsafe.function.CheckedSupplier
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.file.FileHelper
import nextflow.util.Duration
import nextflow.util.ThreadPoolBuilder
import nextflow.util.ThreadPoolHelper

/**
 * Helper class to resolve bucket archive paths
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TowerArchiver {

    private Map<String,String> env = System.getenv()

    private final Path baseDir
    private final Path targetDir

    private ExecutorService executor
    private Duration delay
    private Duration maxDelay
    private Integer maxAttempts
    private Double jitter
    private Duration maxAwait

    Path getBaseDir() { baseDir }

    Path getTargetDir() { targetDir }

    protected TowerArchiver(Path baseDir, Path targetDir, Session session, Map<String,String> env=null) {
        log.debug "Creating tower archiver for base-dir: '$baseDir'; target-dir: '$targetDir'"
        this.baseDir = baseDir
        this.targetDir = targetDir
        // retry settings
        this.delay = session.config.navigate('tower.archiver.delay', '100ms') as Duration
        this.maxDelay = session.config.navigate('tower.archiver.maxDelay', '1m') as Duration
        this.maxAttempts = session.config.navigate('tower.archiver.maxAttempts', '2') as Integer
        this.jitter = session.config.navigate('tower.archiver.jitter', '0.25') as Double
        this.maxAwait = session.config.navigate('tower.archiver.shutdown.maxAwait', '1h') as Duration
        executor = ThreadPoolBuilder.io(10,10,1000, 'tower-archiver')
        if( env!=null )
            this.env = env
    }

    static TowerArchiver create(Session session, Map<String,String> env) {
        def paths = parse(env.get('NXF_ARCHIVE_DIR'))
        if( !paths )
            return null
        new TowerArchiver(Path.of(paths[0]), FileHelper.asPath(paths[1]), session, env)
    }

    static protected List<String> parse(String archiveDef) {
        if( !archiveDef )
            return Collections.<String>emptyList()
        final paths = splitPaths(archiveDef)
        if( !paths )
            return Collections.<String>emptyList()

        if( paths.size()!=2 )
            throw new IllegalArgumentException("Invalid NXF_ARCHIVE_DIR format - expected exactly two paths separated by a command - offending value: ${System.getenv('NXF_ARCHIVE_DIR')}")
        if( !paths[0].startsWith('/') )
            throw new IllegalArgumentException("Invalid NXF_ARCHIVE_DIR base path - it must start with a slash character - offending value: '${paths[0]}'")
        final scheme = FileHelper.getUrlProtocol(paths[1])
        if ( !scheme && !paths[1].startsWith('/') )
            throw new IllegalArgumentException("Invalid NXF_ARCHIVE_DIR target path - it must start be a remote path - offending value: '${paths[1]}'")

        return paths
    }

    static List<String> splitPaths(String paths){
        // multiple paths should be separated by comma
        // allow to escape the separator using backslash
        paths.split(/(?<!\\),/).collect( it-> unescapeQuote(it))
    }

    static protected String unescapeQuote(String uri) {
        uri.replaceAll(/\\,/,',').trim()
    }

    Path archivePath(Path source) {
        if( baseDir==null )
            return null
        if( source==null )
            return null
        if( !source.startsWith(baseDir) )
            return null
        final delta = baseDir.relativize(source)
        // convert to string to prevent 'ProviderMismatchException'
        return targetDir.resolve(delta.toString())
    }

    void archiveLogs() {
        archiveFile(env.get('NXF_OUT_FILE'))
        archiveFile(env.get('NXF_LOG_FILE'))
        archiveFile(env.get('NXF_TML_FILE'))
        archiveFile(env.get('TOWER_CONFIG_FILE'))
        archiveFile(env.get('TOWER_REPORTS_FILE'))
    }

    void archiveTaskLogs(String workDir) {
        final base = Path.of(workDir)
        archiveFile(base.resolve('.command.out'))
        archiveFile(base.resolve('.command.err'))
        archiveFile(base.resolve('.command.log'))
        archiveFile(base.resolve('.command.run'))
        archiveFile(base.resolve('.command.sh'))
        archiveFile(base.resolve('.exitcode'))
    }

    protected void archiveFile(String name) {
        if( name )
            archiveFile(Path.of(name).toAbsolutePath())
    }

    protected void archiveFile(Path source) {
        try {
            final target = archivePath(source)
            if( target==null )
                return
            executor.submit( submitArchive(source,target))
        }
        catch (Throwable t) {
            log.warn("Unable to archive file: $source -- cause: ${t.message ?: t}", t)
        }
    }

    Runnable submitArchive(Path source, Path target) {
        new Runnable() {
            @Override
            void run() {
                try {
                    safeExecute(() -> FileHelper.copyPath(source, target, StandardCopyOption.REPLACE_EXISTING) )
                    log.trace("Archived file: '$source to: '$target'")
                }
                catch (Exception e) {
                    log.warn("Unable to archive file: $source -- cause: ${e.message ?: e}", e)
                }
            }
        }
    }


    protected <T> RetryPolicy<T> retryPolicy() {
        final listener = new EventListener<ExecutionAttemptedEvent<T>>() {
            @Override
            void accept(ExecutionAttemptedEvent<T> event) throws Throwable {
                log.debug("Archiver failed to transfer file - attempt: ${event.attemptCount}", event.lastFailure)
            }
        }
        return RetryPolicy.<T>builder()
                .handleIf(retryCondition())
                .withBackoff(delay.toMillis(), maxDelay.toMillis(), ChronoUnit.MILLIS)
                .withMaxAttempts(maxAttempts)
                .withJitter(jitter)
                .onRetry(listener)
                .build()
    }

    @Memoized
    protected Predicate<? extends Throwable> retryCondition() {
        return new Predicate<Throwable>() {
            @Override
            boolean test(Throwable failure) {
                return failure instanceof IOException
            }
        }
    }

    protected <T> T safeExecute(CheckedSupplier<T> action) {
        final policy = retryPolicy()
        return Failsafe.with(policy).get(action)
    }

    void shutdown() throws IOException {
        executor.shutdown()
        final waitMsg = "Waiting file archiver to complete (%d files)"
        final exitMsg = "Exiting before file archiver thread pool complete -- Some files maybe lost"
        ThreadPoolHelper.await(executor, maxAwait, waitMsg, exitMsg)
    }
}

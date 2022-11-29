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

import nextflow.exception.AbortOperationException

import java.nio.file.FileAlreadyExistsException
import java.nio.file.Path
import java.time.temporal.ChronoUnit
import java.util.concurrent.ExecutorService
import java.util.function.Predicate
import java.util.regex.Pattern

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
import nextflow.util.ThreadPoolManager
/**
 * This class stores all nextflow task '.command.*' files and pipeline reports
 *  into a storage path specified via the NXF_ARCHIVE_DIR env variable
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TowerArchiver {

    private static final String RETRY_REASON = 'slowdown|slow down|toomany|too many'
    private Map<String,String> env = System.getenv()

    private final Path baseDir
    private final Path targetDir

    private Duration delay
    private Duration maxDelay
    private Integer maxAttempts
    private Double jitter
    private Duration maxAwait
    private String retryReason
    private ExecutorService executor
    private ThreadPoolManager poolManager = new ThreadPoolManager('TowerArchiver')

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
        this.retryReason = session.config.navigate('tower.archiver.shutdown.retryReason', RETRY_REASON) as String
        this.executor = poolManager.withConfig(session.config).create()
        if( env!=null )
            this.env = env
    }

    static TowerArchiver create(Session session, Map<String,String> env) {
        final paths = parse(env.get('NXF_ARCHIVE_DIR'))
        if( !paths )
            return null
        final result = new TowerArchiver(Path.of(paths[0]), FileHelper.asPath(paths[1]), session, env)
        return result
    }

    static protected List<String> parse(String archiveDef) {
        if( !archiveDef )
            return Collections.<String>emptyList()
        final paths = splitPaths(archiveDef)
        if( !paths )
            return Collections.<String>emptyList()

        if( paths.size()!=2 )
            throw new IllegalArgumentException("Invalid NXF_ARCHIVE_DIR format - expected exactly two paths separated by a comma - offending value: ${System.getenv('NXF_ARCHIVE_DIR')}")
        if( !paths[0].startsWith('/') )
            throw new IllegalArgumentException("Invalid NXF_ARCHIVE_DIR base path - it must start with a slash character - offending value: '${paths[0]}'")
        final scheme = FileHelper.getUrlProtocol(paths[1])
        if ( !scheme && !paths[1].startsWith('/') )
            throw new IllegalArgumentException("Invalid NXF_ARCHIVE_DIR target path - it must be a remote path - offending value: '${paths[1]}'")

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
        // Nextflow log file are expected to be found at the working directory
        // when running from Tower launcher the cache command, if needed, has
        // copied them before starting the archive process

        final work = env.get('NXF_WORK') ?: env.get('NXF_TEST_WORK')
        if( !work ) {
            log.error("Missing target work dir - Tower archive cannot be performed")
            return
        }
        if( !work.startsWith('/') ) {
            log.error("Invalid NXF_WORK work directory - it must start with a slash character - offending value: '${work}'")
            return
        }

        final base = Path.of(work)
        archiveFile(base, env.get('NXF_OUT_FILE'))
        archiveFile(base, env.get('NXF_LOG_FILE'))
        archiveFile(base, env.get('NXF_TML_FILE'))
        archiveFile(base, env.get('TOWER_CONFIG_FILE'))
        archiveFile(base, env.get('TOWER_REPORTS_FILE'))
    }

    void archiveTaskLogs(String workDir) {
        final base = Path.of(workDir)
        archiveFile(base, '.command.out')
        archiveFile(base, '.command.err')
        archiveFile(base, '.command.log')
        archiveFile(base, '.command.run')
        archiveFile(base, '.command.sh')
        archiveFile(base, '.exitcode')
    }

    protected void archiveFile(Path base, String name) {
        if( name )
            archiveFile(base.resolve(name))
    }

    protected void archiveFile(Path source) {
        try {
            if( !source.exists() ) {
                log.debug "File does not exist: $source -- skipping archiving"
                return
            }
            // go ahead with archiving
            final target = archivePath(source)
            if( target==null )
                return
            log.debug "Submit file archive request: ${source.toUriString()}"
            executor.submit(submitArchive(source,target))
        }
        catch (Throwable t) {
            log.warn("Unable to archive file: $source -- cause: ${t.message ?: t}", t)
        }
    }

    protected void copyPath(Path source, Path target) {
        FileHelper.copyPath(source, target)
    }

    Runnable submitArchive(Path source, Path target) {
        new Runnable() {
            @Override
            void run() {
                try {
                    safeExecute(() -> copyPath(source,target) )
                    log.trace("Archived file: '$source to: '$target'")
                }
                catch (FileAlreadyExistsException e) {
                    log.debug "Skipping archive of file: $source; target alredy exists: $target"
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

    protected Pattern retryPattern() {
        log.debug "File archiver retry-reason='$retryReason'"
        return Pattern.compile(retryReason, Pattern.CASE_INSENSITIVE | Pattern.MULTILINE)
    }

    @Memoized
    protected Predicate<? extends Throwable> retryCondition() {

        return new Predicate<Throwable>() {
            @Override
            boolean test(Throwable failure) {
                final reason = failure.message ?: failure.toString()
                log.trace "Testing file archiver failing reason for retry: '$reason'"
                return retryPattern()
                        .matcher(reason)
                        .find()
            }
        }
    }

    protected <T> T safeExecute(CheckedSupplier<T> action) {
        final policy = retryPolicy()
        return Failsafe.with(policy).get(action)
    }

    void shutdown(Session session) {
        poolManager.shutdown(session)
    }
}

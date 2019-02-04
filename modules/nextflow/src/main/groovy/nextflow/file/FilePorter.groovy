/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.file

import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.util.concurrent.Callable
import java.util.concurrent.ExecutionException
import java.util.concurrent.ExecutorService
import java.util.concurrent.Future
import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.ThreadFactory
import java.util.concurrent.ThreadPoolExecutor
import java.util.concurrent.TimeUnit
import java.util.concurrent.TimeoutException
import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.ProcessStageException
import nextflow.extension.FilesEx
import nextflow.util.CacheHelper
import nextflow.util.Duration
/**
 * Move foreign (ie. remote) files to the staging work area
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ToString(includeFields = true, includeNames = true, includePackage = false, includes = 'coreThreads,maxThreads,maxRetries,keepAlive,pollTimeout')
@CompileStatic
class FilePorter {

    static final private Random RND = Random.newInstance()

    static final private int MAX_RETRIES = 3

    static final private int CORE_THREADS = Runtime.runtime.availableProcessors()

    static final private int MAX_THREADS = 2 * CORE_THREADS

    static final private Duration KEEP_ALIVE = Duration.of('60sec')

    static final private Duration POLL_TIMEOUT = Duration.of('2sec')

    final Map<FileStageAction,Future<NamePathPair>> stagingFutures = new WeakHashMap<>()

    @Lazy private ExecutorService stagingExecutor = createExecutor()

    private Duration keepAlive

    private Duration pollTimeout

    final int coreThreads

    final int maxThreads

    final int maxRetries

    final Session session

    FilePorter( Session session ) {
        this.session = session
        keepAlive = session.config.navigate('filePorter.keepAlive') as Duration ?: KEEP_ALIVE
        pollTimeout = session.config.navigate('filePorter.pollTimeout') as Duration ?: POLL_TIMEOUT
        coreThreads = session.config.navigate('filePorter.coreThreads') as Integer ?: CORE_THREADS
        maxThreads = session.config.navigate('filePorter.maxThreads') as Integer ?: MAX_THREADS
        maxRetries = session.config.navigate('filePorter.maxRetries') as Integer ?: MAX_RETRIES
    }

    /**
     * Given a map of files, copies all the ones stored in a foreign file system
     * and store them in the current working directory
     *
     * @param filesMap A map of files
     * @return A new files map in which all foreign {@link Path} are replaced with local paths
     */
    Map<String,Path> stageForeignFiles(Map<String, Path> filesMap, Path stageDir) {
        final List<FileStageAction> actions = []
        final List<Path> paths = []
        final scheme = stageDir.scheme

        // check for foreign file to copy
        for( Map.Entry<String,Path> entry : filesMap ) {
            def name = entry.getKey()
            def path = entry.getValue()
            if( path.scheme == scheme )
                continue
            // copy the path with a thread pool
            actions.add( new FileStageAction(name, path, stageDir, maxRetries) )
            paths.add(path)
        }

        // no foreign file to copy, just return the original map
        if( !actions )
            return filesMap

        log.trace "Stage foreign files: $paths"
        final result = new HashMap(filesMap)

        final futures = submitStagingActions(actions)
        log.trace "Stage foreign files completed: $paths"

        for( Future<NamePathPair> fut : futures ) {
            final pair = fut.get()
            result.put( pair.name, pair.path )
        }
        return result
    }

    protected List<Future<NamePathPair>> submitStagingActions(List<FileStageAction> actions) {

        String msg = null
        final result = new ArrayList<Future<NamePathPair>>(actions.size())
        for ( def action : actions ) {
            // here's the magic: use a Map to check if a future for the staging action already exist
            // - if exists take it to wait to the submit action termination
            // - otherwise create a new future submitting the action operation
            final future = stagingFutures.getOrCreate(action) {
                msg = "Staging foreign file: ${action.path.toUriString()}"
                return stagingExecutor.submit(action)
            }
            result.add(future)
        }

        // wait for staging actions completion
        // note: make a copy of the list to avoid a ConcurrentModificationException
        final futures = new ArrayList<Future<NamePathPair>>(result)
        while( futures.size() ) {
            final itr = futures.iterator()
            while( itr.hasNext() ) {
                try {
                    final fut = itr.next()
                    fut.get(pollTimeout.millis, TimeUnit.MILLISECONDS)
                    itr.remove()
                }
                catch( TimeoutException e ) {
                    // the timeout event is used to report an info message that
                    // a download is taking place only the very first time
                    if( msg ) {
                        log.info((String)msg)
                        msg = null
                    }
                }
                catch( ExecutionException e ) {
                    throw e.cause ?: e
                }
            }
        }

        return result
    }

    @ToString
    @CompileStatic
    @EqualsAndHashCode
    @PackageScope
    @TupleConstructor
    static class FileStageAction implements Callable<NamePathPair> {

        final String name
        final Path path
        final Path stageDir
        final int maxRetries

        @Override
        NamePathPair call() throws Exception {
            final stagedPath = stageForeignFile(path,stageDir)
            return new NamePathPair(name, stagedPath)
        }

        /**
         * Download a foreign file (ie. remote) storing it the current pipeline execution directory
         *
         * @param filePath The {@link Path} of the remote file to copy
         * @return The path of a local copy of the remote file
         */
        protected Path stageForeignFile(Path filePath, Path stageDir) {
            // create cache directory
            final hash = CacheHelper.hasher([filePath, stageDir]).hash().toString()
            final target = getCacheDir(stageDir, hash).resolve(filePath.getName())

            int count = 0
            while( true ) {
                try {
                    return stageForeignFile0(filePath, target)
                }
                catch( IOException e ) {
                    if( count++ < maxRetries && !(e instanceof NoSuchFileException )) {
                        def message = "Unable to stage foreign file: ${filePath.toUriString()} (try ${count}) -- Cause: $e.message"
                        log.isDebugEnabled() ? log.warn(message, e) : log.warn(message)

                        sleep (10 + RND.nextInt(300))
                        continue
                    }

                    throw new ProcessStageException(fmtError(filePath,e), e)
                }
            }
        }


        protected String getStagingMessage(List<Path> filePaths) {
            def msg = 'Staging foreign file'
            if (filePaths.size() > 1)
                msg += 's:\n'
            else
                msg += ': '
            msg += filePaths.collect() { it.toUriString() }.join('\n')
            return msg
        }

        private String fmtError(Path filePath, Exception e) {
            def message = "Can't stage file ${filePath.toUri().toString()}"
            if( e instanceof NoSuchFileException )
                message += " -- file does not exist"
            else if( e.message )
                message += " -- reason: ${e.message}"
            return message
        }

        private Path stageForeignFile0(Path source, Path target) {
            if( target.exists() ) {
                log.debug "Local cache found for foreign file ${source.toUriString()} at ${target.toUriString()}"
                return target
            }
            log.debug "Copying foreign file ${source.toUriString()} to work dir: ${target.toUriString()}"
            return FileHelper.copyPath(source, target)
        }

        private Path getCacheDir(Path workDir, String hash) {
            def bucket = hash.substring(0,2)
            def result = workDir.resolve( "$bucket/${hash.substring(2)}")

            if( !FilesEx.mkdirs(result) ) {
                throw new IOException("Unable to create cache directory: $result -- Verify file system access permissions or if a file having the same name exists")
            }
            return result.toAbsolutePath()
        }

    }


    @ToString
    @EqualsAndHashCode
    @TupleConstructor
    @PackageScope
    static class NamePathPair {
        String name
        Path path

        String toString() {
            "name=$name,path=$path"
        }
    }

    private ExecutorService createExecutor() {
        log.debug "Creating executor service for $this"

        final result = new ThreadPoolExecutor(
                coreThreads,
                maxThreads,
                keepAlive.millis,
                TimeUnit.MILLISECONDS,
                new LinkedBlockingQueue<Runnable>(),
                new DownloaderThreadFactory() )

        result.allowCoreThreadTimeOut(true);

        // register the shutdown on termination
        session.onShutdown {
            result.shutdown()
            result.awaitTermination(1, TimeUnit.MINUTES)
        }

        return result
    }

    /**
     * Custom thread factory
     */
    static private class DownloaderThreadFactory implements ThreadFactory {
        private final ThreadGroup group;
        private final AtomicInteger threadNumber = new AtomicInteger(1);
        private final String namePrefix = "FilePorter-"

        DownloaderThreadFactory() {
            SecurityManager s = System.getSecurityManager();
            group = (s != null) ? s.getThreadGroup() : Thread.currentThread().getThreadGroup();
        }

        Thread newThread(Runnable r) {
            Thread t = new Thread(group, r, namePrefix + threadNumber.getAndIncrement(), 0);
            if (t.isDaemon())
                t.setDaemon(false);
            if (t.getPriority() != Thread.NORM_PRIORITY)
                t.setPriority(Thread.NORM_PRIORITY);
            return t;
        }
    }
}

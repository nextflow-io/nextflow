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
import java.nio.file.Paths
import java.util.concurrent.Callable
import java.util.concurrent.CompletionService
import java.util.concurrent.ExecutionException
import java.util.concurrent.ExecutorCompletionService
import java.util.concurrent.Future
import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.ThreadFactory
import java.util.concurrent.ThreadPoolExecutor
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.Memoized
import groovy.transform.ToString
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session
import nextflow.exception.ProcessStageException
import nextflow.extension.FilesEx
import nextflow.util.CacheHelper

/**
 * Move foreign (ie. remote) files to the staging work area
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class FilePorter {

    static Random rnd = Random.newInstance()

    private Path stagingDir

    // Timeout for displaying the download message
    private long pollTimeout = 5L
    private TimeUnit pollTimeoutUnit = TimeUnit.SECONDS

    FilePorter() {
        getStagingDir()
    }

    FilePorter( Path stagingDir) {
        assert stagingDir, "Argument stagingDir cannot be null"
        this.stagingDir = stagingDir
    }

    /**
     * Given a map of files, copies all the ones stored in a foreign file system
     * and store them in the current working directory
     *
     * @param filesMap A map of files
     * @return A new files map in which all foreign {@link Path} are replaced with local paths
     */
    Map<String,Path> stageForeignFiles(Map<String, Path> filesMap) {
        List<Callable<NamePathPair>> actions = []
        List<Path> paths = []

        // check for foreign file to copy
        for( Map.Entry<String,Path> entry : filesMap ) {
            def name = entry.getKey()
            def path = entry.getValue()
            if( path.scheme == stagingDir.scheme )
                continue
            // copy the path with a thread pool
            actions << { new NamePathPair(name, stageForeignFile(path)) } as Callable<NamePathPair>
            paths << path
        }

        // no foreign file to copy, just return the original map
        if( !actions )
            return filesMap

        log.trace "Stage foreign files: $paths"
        def result = new HashMap(filesMap)

        def stageExecutor = getExecutor()
        def futures = submitStagingActions(stageExecutor, actions, paths)
        log.trace "Stage foreign files completed: $paths"

        try {
            for( Future<NamePathPair> fut : futures ) {
                final pair = fut.get()
                result.put( pair.name, pair.path )
            }
        }
        catch( ExecutionException e ) {
            throw e.cause ?: e
        }
        return result
    }

    protected ArrayList<Future> submitStagingActions(CompletionService executor, List<Callable> actions, List<Path> paths) {
        // submit actions
        for ( def action : actions ) {
            executor.submit(action)
        }
        def futures = new ArrayList<Future>(actions.size())
        // display a message if staging takes more than the defined timeout
        while ( futures.size() < actions.size() ) {
            def future = executor.poll(pollTimeout, pollTimeoutUnit)
            if ( !future ) {
                log.info1(getDownloadMessage(paths))
                continue
            }
            futures << future
        }
        return futures
    }

    /**
     * Download a foreign file (ie. remote) storing it the current pipeline execution directory
     *
     * @param filePath The {@link Path} of the remote file to copy
     * @return The path of a local copy of the remote file
     */
    protected Path stageForeignFile(Path filePath) {
        int max = getMaxRetries()
        int count = 0

        // create cache directory
        def hash = CacheHelper.hasher(filePath).hash().toString()
        final target = getCacheDir(stagingDir, hash).resolve(filePath.getName())

        // Synchronize download for multiple Nextflow instances running at the same time
        final lockFile = FileHelper.getLocalTempLockPath().resolve("${hash}.lock").toFile()
        final wait = "Another Nextflow instance is downloading the foreign file $filePath -- please wait the download completes"
        final err =  "Unable to acquire exclusive lock on file: $lockFile"
        final mutex = new FileMutex(target: lockFile, waitMessage: wait, errorMessage: err)

        while( true ) {
            try {
                mutex .lock { stageForeignFile0(filePath, target) }
            }
            catch( IOException e ) {
                if( count++ < max && !(e instanceof NoSuchFileException )) {
                    def message = "Unable to stage foreign file: ${filePath.toUriString()} (try ${count}) -- Cause: $e.message"
                    log.isDebugEnabled() ? log.warn(message, e) : log.warn(message)

                    sleep (10 + rnd.nextInt(300))
                    continue
                }

                throw new ProcessStageException(fmtError(filePath,e), e)
            }
            finally {
                lockFile.delete()
            }

            return target
        }
    }

    protected String getDownloadMessage(List<Path> filePaths) {
        def msg = 'Staging foreign file'
        if (filePaths.size() > 1)
            msg += 's:\n'
        else
            msg += ' '
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

    @ToString
    @EqualsAndHashCode
    @TupleConstructor
    static private class NamePathPair {
        String name
        Path path

        String toString() {
            "name=$name,path=$path"
        }
    }

    /**
     * @return
     *      The number of core threads used for parallel files copy.
     *      This value can be defined by using the {@code filePorter.coreThreads} configuration setting
     *      (default: <number of avail CPU cores>)
     */
    static protected int getCoreThreads() {
        Integer result = Global.session.config.navigate('filePorter.coreThreads') as Integer
        if( !result )
            result = Runtime.runtime.availableProcessors()
        log.debug "filePorter.coreThreads=$result"
        return result
    }


    /**
     * @return
     *      The maximum number of threads used for parallel files copy.
     *      This value can be defined by using the {@code filePorter.maxThreads} configuration setting
     *      (default: <2 * number of avail CPU cores>)
     */
    static protected int getMaxThreads() {
        Integer result = Global.session.config.navigate('filePorter.maxThreads') as Integer
        if( !result )
            result = 2 * Runtime.runtime.availableProcessors()
        log.debug "filePorter.maxThreads=$result"
        return result
    }

    /**
     * @return
     *      The maximum number of time the copy process is retried in case of an unexpected error.
     *      This value can be defined by using the {@code filePorter.maxRetries} configuration setting (default: 3).
     */
    static int getMaxRetries() {
        Integer result = Global.session.config.navigate('filePorter.maxRetries') as Integer
        if( !result )
            result = 3
        log.debug "filePorter.maxRetries=$result"
        return result
    }

    /**
     * Lazy getter for the stagingDir property
     * @return
     *      The path the foreign files are copied to.
     *      This value can be defined by using the {@code filePorter.stagingDir} configuration setting (default: <session_workDir/download>).
     */
     Path getStagingDir() {
         if (!stagingDir) {
             String result = Global.session.config.navigate('filePorter.stagingDir') as String
             if (!result)
                 result = System.getenv().get('NXF_DOWNLOAD_PATH')
             if (!result)
                 result = Global.session.workDir.resolve("download").toString()
             log.debug "filePorter.stagingDir=$result"
             stagingDir = Paths.get(result).toAbsolutePath()
         }
         return stagingDir
    }

    /**
     * @return Creates lazily the executor service used to stage remote files
     */
    static private CompletionService getExecutor() {
        createExecutor(Global.session as Session)
    }

    // note: memoized guarantees to use the same executor during the same session
    @Memoized
    static synchronized private CompletionService createExecutor(Session session) {
        def coreThreads = getCoreThreads()
        def maxThreads = getMaxThreads()

        final result = new ThreadPoolExecutor(
                coreThreads,
                maxThreads,
                60L,
                TimeUnit.SECONDS,
                new LinkedBlockingQueue<Runnable>(),
                new DownloaderThreadFactory() )

        result.allowCoreThreadTimeOut(true);

        // register the shutdown on termination
        session?.onShutdown {
            result.shutdown()
            result.awaitTermination(1, TimeUnit.MINUTES)
        }

        return new ExecutorCompletionService<>(result)
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

        public Thread newThread(Runnable r) {
            Thread t = new Thread(group, r, namePrefix + threadNumber.getAndIncrement(), 0);
            if (t.isDaemon())
                t.setDaemon(false);
            if (t.getPriority() != Thread.NORM_PRIORITY)
                t.setPriority(Thread.NORM_PRIORITY);
            return t;
        }
    }
}

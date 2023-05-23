/*
 * Copyright 2013-2023, Seqera Labs
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

import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.SimpleFileVisitor
import java.nio.file.attribute.BasicFileAttributes
import java.util.concurrent.ExecutionException
import java.util.concurrent.ExecutorService
import java.util.concurrent.Future
import java.util.concurrent.Semaphore
import java.util.concurrent.TimeUnit
import java.util.concurrent.TimeoutException
import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.ProcessStageException
import nextflow.extension.FilesEx
import nextflow.util.CacheHelper
import nextflow.util.Duration
import nextflow.util.ThreadPoolManager
import nextflow.util.Threads

/**
 * Move foreign (ie. remote) files to the staging work area
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ToString(includeFields = true, includeNames = true, includePackage = false, includes = 'maxRetries,pollTimeout')
@CompileStatic
class FilePorter {

    static final private Random RND = Random.newInstance()

    static final private int MAX_RETRIES = 3

    static final private int MAX_TRANSFERS = 50

    static final private Duration POLL_TIMEOUT = Duration.of('2sec')

    final Map<FileCopy,FileTransfer> stagingTransfers = new HashMap<>()

    private ExecutorService threadPool

    private Duration pollTimeout

    private final Semaphore semaphore

    final private Lock sync

    final int maxRetries

    final int maxTransfers

    final private Session session

    FilePorter( Session session ) {
        this.session = session
        maxRetries = session.config.navigate('filePorter.maxRetries') as Integer ?: MAX_RETRIES
        pollTimeout = session.config.navigate('filePorter.pollTimeout') as Duration ?: POLL_TIMEOUT
        maxTransfers = session.config.navigate('filePorter.maxTransfers') as Integer ?: MAX_TRANSFERS
        log.debug "File porter settings maxRetries=$maxRetries; maxTransfers=$maxTransfers; pollTimeout=$threadPool"
        sync = new ReentrantLock()
        // use a semaphore to cap the number of max transfer when using virtual thread
        // when using platform threads the max transfers are limited by the thread pool itself
        semaphore = Threads.useVirtual() ? new Semaphore(maxTransfers) : null
        threadPool = new ThreadPoolManager('FileTransfer')
                .withConfig(session.config)
                .createAndRegisterShutdownCallback(session)
    }

    Batch newBatch(Path stageDir) { new Batch(this, stageDir) }

    void transfer(Batch batch) {
        if( batch.size() ) {
            log.trace "Stage foreign files: $batch"
            submitStagingActions(batch.foreignPaths)
            log.trace "Stage foreign files completed: $batch"
        }
    }

    protected FileTransfer createFileTransfer(Path source, Path target) {
        return new FileTransfer(source, target, maxRetries, semaphore)
    }

    protected FileTransfer getOrSubmit(FileCopy copy) {
        sync.lock()
        try {
            FileTransfer transfer = stagingTransfers.get(copy)
            if( transfer == null ) {
                transfer = createFileTransfer(copy.source, copy.target)
                transfer.result = threadPool.submit(transfer)
                stagingTransfers.put(copy, transfer)
            }
            // increment the ref count
            transfer.refCount.incrementAndGet()

            return transfer
        }
        finally {
            sync.unlock()
        }
    }

    protected void decOrRemove(FileTransfer action) {
        //assert stagingTransfers.containsKey(key)
        if( action.refCount.decrementAndGet() == 0 ) {
            stagingTransfers.remove(action.source)
        }
    }

    /**
     * Stages i.e. copies the file from the remote source to a local staging path
     * using a thread pool
     * @param copies
     *      A map where each key-value pair represent a file to be copied.
     *      The key element is the file source path. The value element represent the target path
     * @return
     *      A list of {@link FileTransfer} operations
     */
    protected List<FileTransfer> submitStagingActions(List<FileCopy> copies) {

        final result = new ArrayList<FileTransfer>(copies.size())
        for ( FileCopy it : copies ) {
            // here's the magic: use a Map to check if a future for the staging action already exist
            // - if exists take it to wait to the submit action termination
            // - otherwise create a new future submitting the action operation
            result << getOrSubmit(it)
        }

        // wait for staging actions completion
        // note: make a copy of the list to avoid a ConcurrentModificationException
        final futures = new ArrayList<FileTransfer>(result)
        while( futures.size() ) {
            final itr = futures.iterator()
            while( itr.hasNext() ) {
                final action = itr.next()
                try {
                    action.result.get(pollTimeout.millis, TimeUnit.MILLISECONDS)
                    itr.remove()
                    decOrRemove(action)
                }
                catch( TimeoutException e ) {
                    // the timeout event is used to report an info message that
                    // a download is taking place only the very first time
                    final msg = action.getMessageAndClear()
                    if( msg ) {
                        log.info((String)msg)
                    }
                }
                catch( ExecutionException e ) {
                    throw e.cause ?: e
                }
            }
        }

        return result
    }

    /**
     * Model a file stage requirement
     */
    @Canonical
    static class FileCopy {
        final Path source
        final Path target
    }

    /**
     * Models a batch (collection) of foreign files that need to be transferred to
     * the process staging are co-located with the work directory
     */
    static class Batch  {

        final private FilePorter owner

        /**
         * Holds the list of foreign files to be transferred
         */
        private List<FileCopy> foreignPaths = new ArrayList<>(100)

        /**
         * The *local* directory where against where files need to be staged.
         */
        private Path stageDir

        /**
         * The stage directory scheme. A file is considered *foreign* when its scheme is different from
         * the stage directory scheme.
         */
        private String stageScheme

        Batch(FilePorter owner, Path stageDir) {
            this.owner = owner
            this.stageDir = stageDir
            this.stageScheme = stageDir.scheme
        }

        /**
         * Add the specified path to list of foreign files when
         * the path scheme is different from the stage directory scheme
         *
         * @param path
         *      The path to include in the foreign files batch
         * @return
         *
         */
        Path addToForeign(Path path) {
            // copy the path with a thread pool
            final copy = owner.getCachePathFor(path, stageDir)
            foreignPaths << copy
            return copy.target
        }

        /**
         * @return The number of foreign files in the batch
         */
        int size() { foreignPaths.size() }

        /**
         * @return {@code true} when it contains one or more files, {@code false} otherwise
         */
        boolean asBoolean() { foreignPaths.size()>0 }

        @Override
        String toString() {
            return "FilePorter.Batch[stageDir=${stageDir.toUriString()}; foreignPaths=${foreignPaths}]"
        }
    }

    /**
     * Model a foreign file transfer
     */
    @ToString
    @CompileStatic
    @PackageScope
    static class FileTransfer implements Runnable {

        /**
         * The source file (foreign) to be transfer
         */
        final Path source

        /**
         * The target path where the file need to be copied
         */
        final Path target

        /**
         * Max number of retries in case of error
         */
        final int maxRetries

        final private Semaphore semaphore
        final AtomicInteger refCount
        volatile Future result
        private String message
        private int debugDelay

        FileTransfer(Path foreignPath, Path stagePath, int maxRetries, Semaphore semaphore) {
            this.semaphore = semaphore
            this.source = foreignPath
            this.target = stagePath
            this.maxRetries = maxRetries
            this.message = "Staging foreign file: ${source.toUriString()}"
            this.refCount = new AtomicInteger(0)
            this.debugDelay = System.getProperty('filePorter.debugDelay') as Integer ?: 0
        }

        @Override
        void run() throws Exception {
            if( semaphore )
                semaphore.acquire()
            try {
                stageForeignFile(source, target)
            }
            finally {
                if( semaphore )
                    semaphore.release()
            }
        }

        /**
         * Download a foreign file (ie. remote) storing it the current pipeline execution directory
         *
         * @param filePath The {@link Path} of the remote file to copy
         * @return The path of a local copy of the remote file
         */
        protected Path stageForeignFile(Path filePath, Path stagePath) {

            int count = 0
            while( true ) {
                try {
                    return stageForeignFile0(filePath, stagePath)
                }
                catch( IOException e ) {
                    if( count++ < maxRetries && e !instanceof NoSuchFileException && e !instanceof InterruptedIOException && !Thread.currentThread().isInterrupted() ) {
                        def message = "Unable to stage foreign file: ${filePath.toUriString()} (try ${count}) -- Cause: $e.message"
                        log.isDebugEnabled() ? log.warn(message, e) : log.warn(message)

                        sleep (10 + RND.nextInt(300))
                        continue
                    }

                    throw new ProcessStageException(fmtError(filePath,e), e)
                }
            }
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
            if( debugDelay )
                sleep ( new Random().nextInt(debugDelay) )
            return FileHelper.copyPath(source, target)
        }

        synchronized String getMessageAndClear() {
            def result = message
            message = null
            return result
        }
    }

    protected FileCopy getCachePathFor(Path sourcePath, Path stageDir) {
        sync.lock()
        try {
            final dirPath = stageDir.toUriString() // <-- use a string to avoid changes in the dir to alter the hashing
            int i=0
            while( true ) {
                final uniq = List.of(sourcePath, dirPath, i++)
                final hash = CacheHelper.hasher(uniq).hash().toString()
                final targetPath = getCacheDir0(stageDir, hash).resolve(sourcePath.getName())
                final result = new FileCopy(sourcePath, targetPath)
                if( stagingTransfers.containsKey(result) )
                    return result
                final exist = targetPath.exists()
                if( !exist || checkPathIntegrity(sourcePath, targetPath) )
                    return result
            }
        }
        finally {
            sync.unlock()
        }
    }

    static private Path getCacheDir0(Path workDir, String hash) {
        def bucket = hash.substring(0,2)
        def result = workDir.resolve( "$bucket/${hash.substring(2)}")

        if( !FilesEx.mkdirs(result) ) {
            throw new IOException("Unable to create cache directory: $result -- Make sure a file with the same name doesn't exist and you have write permissions")
        }
        return result.toAbsolutePath()
    }

    static private boolean checkPathIntegrity(Path source, Path target) {
        try {
            // the file must have the same size. this is needed
            // to prevent re-using broken files left by a previous interrupted download
            final attrs = Files.readAttributes(source, BasicFileAttributes)
            final same = attrs.isDirectory()
                    ? checkDirIntegrity0(source, target)
                    : attrs.size() == Files.size(target)
            return same
        }
        catch (NoSuchFileException e) {
            log.warn "Path integrity check failed because the following file has been deleted: $e.message -- make sure to not run more than one nextflow instance using the same work directory"
            return false
        }
        catch (Exception e) {
            log.debug "Unable to determine stage file integrity: source=$source; target=$target", e
            return false
        }
    }

    static private boolean checkDirIntegrity0(Path sourceDir, Path targetDir) {

        // traverse the sourceDir directory and for each file check exists
        // a corresponding file in the target directory having the same size
        boolean same = true
        Files.walkFileTree(sourceDir, new SimpleFileVisitor<Path>() {
            FileVisitResult visitFile(Path sourceFile, BasicFileAttributes attrs) {
                final rel = sourceDir.relativize(sourceFile).toString()

                if( attrs.size() != Files.size(targetDir.resolve(rel)) ) {
                    same = false
                    return FileVisitResult.TERMINATE
                }
                else {
                    same = true
                    return FileVisitResult.CONTINUE
                }
            }} )
        return same
    }
}

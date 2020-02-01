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
import java.util.concurrent.ExecutionException
import java.util.concurrent.ExecutorService
import java.util.concurrent.Future
import java.util.concurrent.TimeUnit
import java.util.concurrent.TimeoutException

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.transform.ToString
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
@ToString(includeFields = true, includeNames = true, includePackage = false, includes = 'maxRetries,pollTimeout')
@CompileStatic
class FilePorter {

    static final private Random RND = Random.newInstance()

    static final private int MAX_RETRIES = 3

    static final private Duration POLL_TIMEOUT = Duration.of('2sec')

    final Map<Path,FileTransfer> stagingTransfers = new HashMap<>()

    @Lazy private ExecutorService threadPool = session.getFileTransferThreadPool()

    private Duration pollTimeout

    final int maxRetries

    final Session session

    FilePorter( Session session ) {
        this.session = session
        maxRetries = session.config.navigate('filePorter.maxRetries') as Integer ?: MAX_RETRIES
        pollTimeout = session.config.navigate('filePorter.pollTimeout') as Duration ?: POLL_TIMEOUT
    }

    Batch newBatch(Path stageDir) { new Batch(stageDir) }

    void transfer(Batch batch) {
        if( batch.size() ) {
            log.trace "Stage foreign files: $batch"
            submitStagingActions(batch.foreignPaths, batch.stageDir)
            log.trace "Stage foreign files completed: $batch"
        }
    }

    protected FileTransfer createFileTransfer(Path source, Path stageDir) {
        final stagePath = getCachePathFor(source,stageDir)
        return new FileTransfer(source, stagePath, maxRetries)
    }

    protected Future submitForExecution(FileTransfer transfer) {
        threadPool.submit(transfer)
    }

    protected FileTransfer getOrSubmit(Path source, Path stageDir) {
        synchronized (stagingTransfers) {
            FileTransfer transfer = stagingTransfers.get(source)
            if( transfer == null ) {
                transfer = createFileTransfer(source, stageDir)
                transfer.refCount = 1
                transfer.result = submitForExecution(transfer)
                stagingTransfers.put(source, transfer)
            }
            else
                transfer.refCount ++

            return transfer
        }
    }

    protected void decOrRemove(FileTransfer action) {
        synchronized (stagingTransfers) {
            final key = action.source
            assert stagingTransfers.containsKey(key)
            if( --action.refCount == 0 ) {
                stagingTransfers.remove(key)
            }
        }
    }

    protected List<FileTransfer> submitStagingActions(List<Path> paths, Path stageDir) {

        final result = new ArrayList<FileTransfer>(paths.size())
        for ( def file : paths ) {
            // here's the magic: use a Map to check if a future for the staging action already exist
            // - if exists take it to wait to the submit action termination
            // - otherwise create a new future submitting the action operation
            result << getOrSubmit(file,stageDir)
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
     * Models a batch (collection) of foreign files that need to be transferred to
     * the process staging are co-located with the work directory
     */
    static class Batch  {

        /**
         * Holds the list of foreign files to be transferred
         */
        private List<Path> foreignPaths = new ArrayList<>()

        /**
         * The *local* directory where against where files need to be staged.
         */
        private Path stageDir

        /**
         * The stage directory scheme. A file is considered *foreign* when its scheme is different from
         * the stage directory scheme.
         */
        private String stageScheme

        Batch(Path stageDir) {
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
            foreignPaths << path
            return getCachePathFor(path, stageDir)
        }

        /**
         * @return The number of foreign files in the batch
         */
        int size() { foreignPaths.size() }

        /**
         * @return {@code true} when it contains one or more files, {@code false} otherwise
         */
        boolean asBoolean() { foreignPaths.size()>0 }

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

        volatile int refCount
        volatile Future result
        private String message

        FileTransfer(Path foreignPath, Path stagePath, int maxRetries=0) {
            this.source = foreignPath
            this.target = stagePath
            this.maxRetries = maxRetries
            this.message = "Staging foreign file: ${source.toUriString()}"
        }

        @Override
        void run() throws Exception {
            stageForeignFile(source, target)
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

        synchronized String getMessageAndClear() {
            def result = message
            message = null
            return result
        }
    }


    static protected Path getCachePathFor(Path filePath, Path stageDir) {
        final dirPath = stageDir.toUriString() // <-- use a string to avoid changes in the dir to alter the hashing
        final hash = CacheHelper.hasher([filePath, dirPath]).hash().toString()
        final result = getCacheDir0(stageDir, hash).resolve(filePath.getName())
        return result
    }

    static private Path getCacheDir0(Path workDir, String hash) {
        def bucket = hash.substring(0,2)
        def result = workDir.resolve( "$bucket/${hash.substring(2)}")

        if( !FilesEx.mkdirs(result) ) {
            throw new IOException("Unable to create cache directory: $result -- Verify file system access permissions or if a file having the same name exists")
        }
        return result.toAbsolutePath()
    }


}

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

package nextflow.util

import java.nio.ByteBuffer
import java.nio.file.FileAlreadyExistsException
import java.nio.file.Files
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.StandardOpenOption

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * A cloud-based lock manager that uses atomic file creation for distributed locking
 * across multiple processes accessing cloud storage (S3, Azure Blob, GCS, etc.)
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class GlobalCacheLockManager {

    private Path workDir

    GlobalCacheLockManager(Path workDir){
        assert workDir
        this.workDir = workDir
    }

    /**
     * Attempts to acquire a distributed lock by atomically creating a lock file
     * in the cloud storage. The lock file is created at: workDir/.hash.lock
     *
     * @param hash The task hash code as a string
     * @return A GlobalCacheLockHandle if the lock was successfully acquired, null otherwise
     */
    LockHandle acquire( String hash ) {
        assert hash
        final lockPath = this.workDir.resolve(".${hash}.lock")

        try {
            // Try to create the lock file atomically using CREATE_NEW
            // This will fail with FileAlreadyExistsException if another process holds the lock
            final channel = Files.newByteChannel(
                lockPath,
                StandardOpenOption.CREATE_NEW,
                StandardOpenOption.WRITE
            )

            // Write some metadata to the lock file (optional, for debugging)
            final timestamp = System.currentTimeMillis()
            final pid = ProcessHandle.current().pid()
            final content = "locked_at=${timestamp}\npid=${pid}\n"
            final buffer = ByteBuffer.wrap(content.getBytes())
            channel.write(buffer)
            channel.close()

            log.trace "Acquired global cache lock for hash=${hash} at path=${lockPath}"
            return new GlobalCacheLockHandle(lockPath)
        }
        catch (FileAlreadyExistsException e) {
            // Another process holds the lock
            log.debug "Another process holds the global lock for hash=${hash} - already exists at path=${lockPath}"
            return null
        }
        catch (Exception e) {
            log.warn "Error while acquiring cloud lock for hash=${hash} at path=${lockPath}", e
            throw e
        }
    }

    /**
     * Handle for a cloud file-based lock that can be released
     */
    @CompileStatic
    static class GlobalCacheLockHandle implements LockHandle{
        private final Path lockPath
        private volatile boolean released = false

        GlobalCacheLockHandle(Path lockPath) {
            this.lockPath = lockPath
        }

        /**
         * Releases the lock by deleting the lock file
         */
        void release() {
            if (released) {
                return
            }

            try {
                Files.deleteIfExists(lockPath)
                log.trace "Released cloud lock at path=${lockPath}"
                released = true
            }
            catch (NoSuchFileException e) {
                // Lock file was already deleted, ignore
                log.trace "Cloud lock file already deleted: ${lockPath}"
            }
            catch (Exception e) {
                log.warn "Error while releasing cloud lock at path=${lockPath}", e
            }
        }
    }
}
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

package nextflow.executor

import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.file.FileHelper
import nextflow.util.Duration
/**
 * Reads the exit status from a task's {@code .exitcode} file in a way that tolerates
 * the file being momentarily missing or empty right after the job is reported terminated
 * by a remote API (e.g. Kubernetes, AWS/Google/Azure Batch), which can happen when the
 * task work directory is backed by an eventually-consistent or network filesystem.
 *
 * This is meant to be polled repeatedly from {@code TaskHandler#checkIfCompleted()}: a
 * {@code null} result means "keep waiting and call again on the next poll cycle" rather
 * than blocking, since {@code checkIfCompleted()} runs on a single shared polling thread.
 */
@Slf4j
@CompileStatic
class ExitStatusAwaiter {

    private final long timeoutMillis

    private long missingSinceMillis

    private long emptySinceMillis

    ExitStatusAwaiter(Duration timeout) {
        this.timeoutMillis = timeout.toMillis()
    }

    /**
     * @param exitFile The path to the task's {@code .exitcode} file
     * @return
     *      the exit status read from the file once it's present and has content,
     *      {@code null} if the file is missing or empty but still within the timeout
     *      window (the caller should treat the task as not-yet-complete and retry on
     *      the next poll cycle),
     *      or {@link Integer#MAX_VALUE} once the timeout has elapsed without finding
     *      usable content.
     */
    Integer read(Path exitFile) {

        final attrs = exitFile ? FileHelper.readAttributes(exitFile) : null
        if( !attrs || !attrs.lastModifiedTime()?.toMillis() ) {
            if( !missingSinceMillis ) {
                log.trace "Exit file does not exist yet -- waiting: ${exitFile}"
                missingSinceMillis = System.currentTimeMillis()
                return null
            }
            final delta = System.currentTimeMillis() - missingSinceMillis
            if( delta < timeoutMillis )
                return null
            log.warn "Unable to find exit status file after ${delta}ms: ${exitFile?.toUriString()}"
            return Integer.MAX_VALUE
        }

        final status = exitFile.text?.trim()
        if( status ) {
            try {
                return status.toInteger()
            }
            catch( Exception e ) {
                log.warn "Unable to parse process exit file: ${exitFile.toUriString()} -- bad value: '$status'"
                return Integer.MAX_VALUE
            }
        }

        // file exists but it's empty -- this can happen because of write propagation
        // delays on shared/eventually-consistent filesystems; wait before giving up
        if( !emptySinceMillis ) {
            log.debug "Exit file is empty -- waiting: ${exitFile}"
            emptySinceMillis = System.currentTimeMillis()
            return null
        }
        final delta = System.currentTimeMillis() - emptySinceMillis
        if( delta < timeoutMillis )
            return null
        log.warn "Exit status file is still empty after ${delta}ms: ${exitFile.toUriString()}"
        return Integer.MAX_VALUE
    }
}

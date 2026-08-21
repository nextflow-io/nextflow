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

package nextflow.file

import java.io.IOException
import java.io.OutputStream

/**
 * A handle to an in-progress resumable upload, produced by {@link ResumableFileSystem#resumeUpload}.
 */
interface ResumableUpload {

    /**
     * @return The number of bytes already committed to the upload.
     */
    long committedBytes()

    /**
     * @return An output stream that writes the remaining bytes of the upload.
     */
    OutputStream outputStream()

    /**
     * Finish the upload, making the target object visible.
     */
    void complete() throws IOException

    /**
     * Abandon the upload, leaving it in-progress so a later attempt can resume it.
     */
    void abandon() throws IOException

    /**
     * Abort the upload, discarding the in-progress upload entirely.
     */
    void abort() throws IOException

}

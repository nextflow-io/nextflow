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
import java.nio.file.CopyOption
import java.nio.file.Path

/**
 * A file system provider whose targets can resume an upload after an interruption.
 *
 * Implemented by providers that support resumable writes (e.g. S3 multipart upload). Used by a
 * range-capable source provider to continue staging into a partially-uploaded target.
 */
interface ResumableFileSystem {

    /**
     * Start a fresh upload of {@code target}, returning a handle that can write and complete it.
     */
    ResumableUpload newUpload(Path target, CopyOption... options) throws IOException

    /**
     * Recover an in-progress upload of {@code target}, returning a handle that can continue it,
     * or {@code null} when there is no resumable in-progress upload.
     */
    ResumableUpload resumeUpload(Path target) throws IOException

}

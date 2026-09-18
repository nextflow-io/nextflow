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

package nextflow.cloud.aws.nio

import java.nio.file.Path
import java.nio.file.ProviderMismatchException

import groovy.transform.CompileStatic
import nextflow.file.AtomicLockProvider
import org.pf4j.Extension

/**
 * {@link AtomicLockProvider} for {@code s3://} paths. Creates the lock object via a
 * conditional PUT ({@code If-None-Match: *}); a 412 (object exists) maps to a lost claim.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Extension
@CompileStatic
class S3AtomicLockProvider extends AtomicLockProvider {

    @Override
    boolean canHandle(String scheme) {
        return scheme == 's3'
    }

    @Override
    boolean tryCreate(Path lockPath) {
        // a path from another provider can never be locked here: fail fast, so the caller can tell
        // "cannot lock" apart from "lost the race" (a null lock) and does not retry forever
        if( !(lockPath instanceof S3Path) )
            throw new ProviderMismatchException("Not a valid S3 path -- cannot create the atomic lock object: `${lockPath}` [${lockPath?.class?.name ?: '-'}]")
        final s3 = (S3Path) lockPath
        // an IOException from a non-412 error propagates — no silent fallback
        return s3.getFileSystem().getClient().putObjectIfAbsent(s3.getBucket(), s3.getKey())
    }


}

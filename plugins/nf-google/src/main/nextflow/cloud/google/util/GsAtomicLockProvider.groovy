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

package nextflow.cloud.google.util

import java.nio.file.Path
import java.nio.file.ProviderMismatchException

import com.google.cloud.storage.BlobId
import com.google.cloud.storage.BlobInfo
import com.google.cloud.storage.Storage
import com.google.cloud.storage.StorageException
import com.google.cloud.storage.contrib.nio.CloudStoragePath
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import nextflow.file.AtomicLockProvider
import org.pf4j.Extension

/**
 * {@link AtomicLockProvider} for {@code gs://} paths.
 *
 * Google's {@code CloudStorageFileSystemProvider} cannot be extended, so this handler uses
 * the GCS {@link Storage} client directly, creating the lock object with a precondition
 * ({@code doesNotExist}); a 412 (object exists) maps to a lost claim.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Extension
@CompileStatic
class GsAtomicLockProvider extends AtomicLockProvider {

    @Override
    boolean canHandle(String scheme) {
        return scheme == 'gs'
    }

    /** Carries the configured credentials, project id, timeouts and retry policy -- see {@link GsStorageOptions}. */
    @Memoized
    protected Storage storage() {
        return GsStorageOptions.sharedClientFor(GsStorageOptions.sessionOpts())   // one client per JVM, not one per SPI
    }

    /** The project to bill on a requester-pays bucket, or {@code null}. */
    @Memoized
    protected String userProject() {
        return GsStorageOptions.userProject(GsStorageOptions.sessionOpts())
    }

    @Override
    boolean tryCreate(Path lockPath) {
        // a path from another provider can never be locked here: fail fast, so the caller can tell
        // "cannot lock" apart from "lost the race" (a null lock) and does not retry forever
        if( !(lockPath instanceof CloudStoragePath) )
            throw new ProviderMismatchException("Not a valid Google Cloud Storage path -- cannot create the atomic lock object: `${lockPath}` [${lockPath?.class?.name ?: '-'}]")
        final gs = (CloudStoragePath) lockPath
        final bucket = gs.getFileSystem().bucket()
        final str = gs.toString()
        final name = str.startsWith('/') ? str.substring(1) : str
        final blobInfo = BlobInfo.newBuilder(BlobId.of(bucket, name)).build()
        final opts = new ArrayList<Storage.BlobTargetOption>()
        opts.add(Storage.BlobTargetOption.doesNotExist())
        final billTo = userProject()
        if( billTo )
            opts.add(Storage.BlobTargetOption.userProject(billTo))
        try {
            storage().create(blobInfo, new byte[0], opts as Storage.BlobTargetOption[])
            return true
        }
        catch( StorageException e ) {
            if( e.getCode() == 412 )
                return false
            throw e
        }
    }


}

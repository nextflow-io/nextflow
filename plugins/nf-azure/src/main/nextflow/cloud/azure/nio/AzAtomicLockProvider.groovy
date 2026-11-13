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

package nextflow.cloud.azure.nio

import java.nio.file.Path
import java.nio.file.ProviderMismatchException

import com.azure.core.util.BinaryData
import com.azure.core.util.Context
import com.azure.storage.blob.models.BlobErrorCode
import com.azure.storage.blob.models.BlobRequestConditions
import com.azure.storage.blob.models.BlobStorageException
import com.azure.storage.blob.options.BlockBlobSimpleUploadOptions
import groovy.transform.CompileStatic
import nextflow.file.AtomicLockProvider
import org.pf4j.Extension

/**
 * {@link AtomicLockProvider} for {@code az://} paths. Creates the lock blob with a
 * conditional upload ({@code If-None-Match: *}); the error CODE {@code BlobAlreadyExists} or
 * {@code ConditionNotMet} maps to a lost claim. Not the HTTP status: Azure answers 409 and 412 for
 * several unrelated conditions too, and reading those as "someone else claimed it" would hide a
 * real failure as a routine race (see {@link #uploadIfAbsent}).
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Extension
@CompileStatic
class AzAtomicLockProvider extends AtomicLockProvider {

    @Override
    boolean canHandle(String scheme) {
        return scheme == 'az'
    }

    @Override
    boolean tryCreate(Path lockPath) {
        // a path from another provider can never be locked here: fail fast, so the caller can tell
        // "cannot lock" apart from "lost the race" (a null lock) and does not retry forever
        if( !(lockPath instanceof AzPath) )
            throw new ProviderMismatchException("Not a valid Azure blob storage path -- cannot create the atomic lock object: `${lockPath}` [${lockPath?.class?.name ?: '-'}]")
        try {
            uploadIfAbsent((AzPath) lockPath)
            return true
        }
        catch( BlobStorageException e ) {
            // Match the ERROR CODE, not the HTTP status. Azure answers 409 for BlobAlreadyExists but
            // also for ContainerBeingDeleted / LeaseIdMismatch, and 412 for ConditionNotMet but also
            // for LeaseIdMissing / LeaseLost -- none of which mean "the blob already exists". The SPI
            // contract is that `false` means one thing only, because the caller bumps the hash and
            // retries on a lost race: mapping an unrelated failure to `false` sends it walking hashes
            // to the 100-conflict abort, whose message points at cache configuration.
            final code = e.getErrorCode()
            if( code == BlobErrorCode.BLOB_ALREADY_EXISTS || code == BlobErrorCode.CONDITION_NOT_MET )
                return false
            throw e
        }
    }

    /**
     * The conditional upload itself. Declared as a protected method (rather than inlined) so a test
     * can drive the error mapping above without a real Azure connection -- the SDK's client classes
     * are final and {@code BlobStorageException} cannot be subclassed.
     */
    protected void uploadIfAbsent(AzPath az) {
        final opts = new BlockBlobSimpleUploadOptions(BinaryData.fromBytes(new byte[0]))
                .setRequestConditions(new BlobRequestConditions().setIfNoneMatch('*'))
        az.blobClient().getBlockBlobClient().uploadWithResponse(opts, null, Context.NONE)
    }


}

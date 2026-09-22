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

import com.azure.core.http.HttpHeaderName
import com.azure.core.http.HttpHeaders
import com.azure.core.http.HttpResponse
import com.azure.storage.blob.models.BlobStorageException
import spock.lang.Specification
import spock.lang.Unroll

class AzAtomicLockProviderTest extends Specification {

    def 'handles only the az scheme'() {
        given:
        def provider = new AzAtomicLockProvider()
        expect:
        provider.canHandle('az')
        !provider.canHandle('s3')
        !provider.canHandle('gs')
    }

    def 'tryCreate throws for a path of another provider, instead of reporting a lost race'() {
        when:
        new AzAtomicLockProvider().tryCreate(Mock(Path))
        then: 'a false return would be read as "lost the race" and retried forever'
        thrown(ProviderMismatchException)
    }

    /** Drives the error mapping without a connection: the SDK clients are final and
     *  BlobStorageException cannot be subclassed, so the upload itself is the seam. */
    static class Failing extends AzAtomicLockProvider {
        BlobStorageException error
        @Override protected void uploadIfAbsent(AzPath az) { throw error }
    }

    /** `getErrorCode()` reads the `x-ms-error-code` response header, so the header IS the fixture. */
    private BlobStorageException azError(String code, int status) {
        final headers = new HttpHeaders().set(HttpHeaderName.fromString('x-ms-error-code'), code)
        final resp = Stub(HttpResponse) { getHeaders() >> headers; getStatusCode() >> status }
        return new BlobStorageException("simulated ${code}", resp, null)
    }

    @Unroll
    def 'a #errorCode response means the object already existed -- a lost race'() {
        expect:
        new Failing(error: azError(errorCode, status)).tryCreate(Mock(AzPath)) == false
        where:
        errorCode            | status
        'BlobAlreadyExists'  | 409
        'ConditionNotMet'    | 412
    }

    @Unroll
    def 'an unrelated #errorCode failure propagates instead of looking like a lost race'() {
        when: 'Azure reuses 409 and 412 for conditions that have nothing to do with existence'
        new Failing(error: azError(errorCode, status)).tryCreate(Mock(AzPath))

        then: 'it must throw -- a false here sends the caller bumping hashes to the conflict abort,'
        and: 'with an error message pointing at cache configuration rather than at the real cause'
        thrown(BlobStorageException)

        where:
        errorCode                | status
        'ContainerBeingDeleted'  | 409
        'LeaseIdMismatchWithBlobOperation' | 409
        'LeaseIdMissing'         | 412
        'LeaseLost'              | 412
    }

}

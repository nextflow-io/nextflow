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

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.util.Locale

import com.azure.storage.blob.BlobServiceClient
import com.azure.storage.blob.BlobServiceClientBuilder
import com.azure.storage.common.StorageSharedKeyCredential
import nextflow.Global
import nextflow.Session
import spock.lang.IgnoreIf
import spock.lang.Requires
import spock.lang.Shared
import spock.lang.Specification

/**
 * Credential-gated real-Azure test for {@link AzAtomicLockProvider}: skipped without Azure
 * credentials, runs in CI/CD. Verifies the conditional-upload atomic semantics.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@IgnoreIf({System.getenv('NXF_SMOKE')})
@Requires({System.getenv('AZURE_STORAGE_ACCOUNT_NAME') && System.getenv('AZURE_STORAGE_ACCOUNT_KEY')})
class AzAtomicLockProviderAzTest extends Specification implements AzBaseSpec {

    @Shared
    BlobServiceClient storageClient

    def setupSpec() {
        def accountName = System.getenv('AZURE_STORAGE_ACCOUNT_NAME')
        def accountKey = System.getenv('AZURE_STORAGE_ACCOUNT_KEY')
        def credential = new StorageSharedKeyCredential(accountName, accountKey)
        def endpoint = String.format(Locale.ROOT, "https://%s.blob.core.windows.net", accountName)
        storageClient = new BlobServiceClientBuilder().endpoint(endpoint).credential(credential).buildClient()
        and:
        def cfg = [azure: [storage: [accountName: accountName, accountKey: accountKey]]]
        Global.session = Mock(Session) { getConfig() >> cfg }
    }

    def cleanupSpec() {
        Global.session = null
    }

    def 'tryCreate is atomic: the first caller wins, the second conflicts'() {
        given:
        def bucket = createBucket()
        def lock = Paths.get(new URI("az://$bucket/.nf-lock-${UUID.randomUUID()}"))
        def provider = new AzAtomicLockProvider()

        expect: 'first create wins'
        provider.tryCreate(lock)

        and: 'a second create conflicts'
        !provider.tryCreate(lock)

        cleanup:
        // the SPI has no release: in production the marker is never removed, so a test that
        // needs the object gone deletes it directly (what the dropped impls did anyway)
        Files.deleteIfExists(lock)
        tryDeleteBucket(bucket)
    }

}

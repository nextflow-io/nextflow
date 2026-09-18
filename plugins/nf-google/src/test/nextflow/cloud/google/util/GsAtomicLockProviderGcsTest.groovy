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

import java.nio.file.Files
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import nextflow.Global
import nextflow.Session
import spock.lang.IgnoreIf
import spock.lang.Requires
import spock.lang.Specification

/**
 * Credential-gated real-GCS test for {@link GsAtomicLockProvider}: skipped without GCS
 * credentials and a writable test bucket, runs in CI/CD. Verifies the
 * {@code doesNotExist} precondition atomic semantics.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@IgnoreIf({System.getenv('NXF_SMOKE')})
@Requires({System.getenv('GOOGLE_APPLICATION_CREDENTIALS') && System.getenv('NXF_GCS_TEST_BUCKET')})
class GsAtomicLockProviderGcsTest extends Specification {

    def setup() {
        def cfg = [google: [project: System.getenv('GOOGLE_PROJECT_ID')]]
        Global.session = Mock(Session) { getConfig() >> cfg }
    }

    def cleanup() {
        Global.session = null
    }

    def 'tryCreate is atomic: the first caller wins, the second conflicts'() {
        given:
        def bucket = System.getenv('NXF_GCS_TEST_BUCKET')
        def lock = CloudStorageFileSystem.forBucket(bucket).getPath("/.nf-lock-${UUID.randomUUID()}")
        def provider = new GsAtomicLockProvider()

        expect: 'first create wins'
        provider.tryCreate(lock)

        and: 'a second create conflicts (412)'
        !provider.tryCreate(lock)

        cleanup:
        // the SPI has no release: in production the marker is never removed, so a test that
        // needs the object gone deletes it directly (what the dropped impls did anyway)
        Files.deleteIfExists(lock)
    }

}

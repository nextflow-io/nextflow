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

import java.nio.file.Files
import nextflow.Global
import nextflow.Session
import nextflow.file.FileHelper
import software.amazon.awssdk.services.s3.S3Client
import spock.lang.IgnoreIf
import spock.lang.Requires
import spock.lang.Specification

/**
 * Credential-gated real-S3 test for {@link S3AtomicLockProvider}: skipped without AWS
 * credentials, runs in CI/CD. Verifies the conditional-PUT atomic semantics.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@IgnoreIf({System.getenv('NXF_SMOKE')})
@Requires({System.getenv('AWS_S3FS_ACCESS_KEY') && System.getenv('AWS_S3FS_SECRET_KEY')})
class S3AtomicLockProviderAwsTest extends Specification implements AwsS3BaseSpec {

    private S3Client s3Client0

    S3Client getS3Client() { s3Client0 }

    static private Map config0() {
        def accessKey = System.getenv('AWS_S3FS_ACCESS_KEY')
        def secretKey = System.getenv('AWS_S3FS_SECRET_KEY')
        return [aws: [accessKey: accessKey, secretKey: secretKey]]
    }

    def setup() {
        def fs = (S3FileSystem) FileHelper.getOrCreateFileSystemFor(URI.create("s3:///"), config0().aws)
        s3Client0 = fs.client.getClient()
        and:
        def cfg = config0()
        Global.config = cfg
        Global.session = Mock(Session) { getConfig() >> cfg }
    }

    def 'tryCreate is atomic: the first caller wins, the second conflicts'() {
        given:
        def bucket = createBucket()
        def lock = s3path("s3://$bucket/.nf-lock-${UUID.randomUUID()}")
        def provider = new S3AtomicLockProvider()

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

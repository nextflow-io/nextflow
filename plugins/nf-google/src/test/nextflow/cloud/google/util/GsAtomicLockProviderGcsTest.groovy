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

import com.google.cloud.storage.BucketInfo
import com.google.cloud.storage.Storage
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session
import nextflow.cloud.google.GoogleOpts
import org.junit.Assume
import spock.lang.IgnoreIf
import spock.lang.Requires
import spock.lang.Shared
import spock.lang.Specification

/**
 * Credential-gated real-GCS test for {@link GsAtomicLockProvider}, verifying that the
 * {@code doesNotExist} precondition is genuinely atomic: only a live service can answer that,
 * so a mock cannot stand in for this.
 *
 * <p>It creates its own bucket and deletes it, the way {@code AwsS3BaseSpec} and
 * {@code AzBaseSpec} do for S3 and Azure, so the only thing it needs from the environment is
 * credentials. The project comes from the credentials too --
 * {@link GoogleOpts#getProjectIdFromCreds} reads {@code project_id} out of the service-account
 * JSON that {@code GOOGLE_APPLICATION_CREDENTIALS} points at -- so there is nothing else to
 * configure and the test runs wherever those credentials are present.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@IgnoreIf({System.getenv('NXF_SMOKE')})
@Requires({System.getenv('GOOGLE_APPLICATION_CREDENTIALS')})
class GsAtomicLockProviderGcsTest extends Specification {

    @Shared String bucket
    @Shared Storage storage

    def setupSpec() {
        // a user ADC (`gcloud auth application-default login`) carries no project_id -- only a
        // service-account key does, and getProjectIdFromCreds THROWS rather than returning null
        // for one that does not. Skip rather than fail, so a developer with a user ADC is not
        // stopped by a test CI can still run with its service account.
        String project = null
        try {
            project = GoogleOpts.getProjectIdFromCreds(System.getenv('GOOGLE_APPLICATION_CREDENTIALS'))
        }
        catch( Exception e ) {
            log.warn "Skipping: cannot read a project_id from the credentials -- ${e.message}"
        }
        Assume.assumeTrue('no project_id in the credentials file', project != null)

        Global.session = Mock(Session) { getConfig() >> [google: [project: project]] }
        storage = GsStorageOptions.sharedClientFor(GsStorageOptions.sessionOpts())
        bucket = "nf-gs-test-${UUID.randomUUID()}"
        log.debug "Creating gcs bucket '$bucket'"
        storage.create(BucketInfo.of(bucket))
    }

    def cleanupSpec() {
        try {
            if( bucket && storage )
                storage.delete(bucket)
        }
        catch( Exception e ) {
            log.warn "Unable to remove gcs bucket '$bucket' -- ${e.message}"
        }
        Global.session = null
    }

    def 'tryCreate is atomic: the first caller wins, the second conflicts'() {
        given:
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

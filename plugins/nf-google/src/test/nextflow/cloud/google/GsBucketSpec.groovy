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

package nextflow.cloud.google

import com.google.cloud.storage.BucketInfo
import com.google.cloud.storage.Storage
import com.google.cloud.storage.StorageOptions
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Creates and disposes a throw-away Google Storage bucket, so that tests needing real storage
 * only require credentials. Mirrors {@code AwsS3BaseSpec} in the nf-amazon plugin.
 *
 * @author Ben Sherman <ben.sherman@seqera.io>
 */
trait GsBucketSpec {

    static final Logger log = LoggerFactory.getLogger(GsBucketSpec)

    private Storage storage0

    Storage getStorage() {
        if( !storage0 )
            storage0 = StorageOptions.getDefaultInstance().getService()
        return storage0
    }

    String createBucket() {
        final name = "nf-gcsfs-test-${UUID.randomUUID()}".toString()
        log.debug "Creating gs bucket '$name'"
        storage.create(BucketInfo.newBuilder(name).build())
        return name
    }

    void deleteBucket(String bucketName) {
        if( !bucketName )
            return
        log.debug "Deleting gs bucket '$bucketName'"
        // a bucket can only be removed once empty, and the gcsfuse placeholder objects
        // deliberately survive the code under test
        storage.list(bucketName).iterateAll().each { it.delete() }
        storage.delete(bucketName)
    }
}

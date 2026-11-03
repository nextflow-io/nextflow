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

package nextflow.extension

import spock.lang.IgnoreIf
import spock.lang.Requires
import spock.lang.Specification
import spock.lang.Unroll

import java.nio.file.Files
import java.nio.file.Path

import com.google.cloud.storage.BlobId
import com.google.cloud.storage.BlobInfo
import com.google.cloud.storage.Storage
import com.google.cloud.storage.StorageOptions
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import com.google.cloud.storage.contrib.nio.CloudStoragePath
import nextflow.Global
import nextflow.Session

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FilesExTest2 extends Specification {

    /** Shared test bucket, the same one used by validation/google.config */
    static final String TEST_BUCKET = 'rnaseq-nf'

    @Unroll
    def 'should return uri string for #PATH' () {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [google:[project:'foo', region:'x']]
        }

        when:
        def path = PATH as Path
        then:
        path instanceof CloudStoragePath
        println FilesEx.toUriString(path)
        FilesEx.toUriString(path) == PATH

        where:
        PATH                    | _
        'gs://foo/bar'          | _
        'gs://foo'              | _
        'gs://foo/'             | _
        'gs://foo/bar/baz'      | _
        'gs://foo/bar/baz/'     | _
        'gs://foo/bar - baz/'   | _
    }

    /**
     * A directory created by gcsfuse leaves behind a zero-byte placeholder object that the NIO
     * provider refuses to delete, which makes the parent look non-empty and trips the unchecked
     * CloudStoragePseudoDirectoryException.
     *
     * See https://github.com/nextflow-io/nextflow/issues/5647
     */
    @IgnoreIf({System.getenv('NXF_SMOKE')})
    @Requires({System.getenv('GOOGLE_APPLICATION_CREDENTIALS')})
    def 'should delete a directory holding a gcsfuse placeholder object' () {
        given:
        def storage = StorageOptions.getDefaultInstance().getService()
        def prefix = "nf-test-${UUID.randomUUID()}"
        def base = CloudStorageFileSystem.forBucket(TEST_BUCKET).getPath("/$prefix")
        and:
        Files.write(base.resolve('output/reads.txt'), 'data'.bytes)
        storage.create(BlobInfo.newBuilder(BlobId.of(TEST_BUCKET, "$prefix/output/")).build(), new byte[0])

        when:
        def result = FilesEx.deleteDir(base)

        then:
        noExceptionThrown()
        result
        and:
        // every real object is gone, only the placeholder may survive
        storage.list(TEST_BUCKET, Storage.BlobListOption.prefix(prefix)).iterateAll().count { !it.name.endsWith('/') } == 0

        cleanup:
        storage?.list(TEST_BUCKET, Storage.BlobListOption.prefix(prefix))?.iterateAll()?.each { it.delete() }
    }

}

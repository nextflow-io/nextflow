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

package nextflow.processor

import java.nio.file.Path

import com.google.cloud.storage.BlobId
import com.google.cloud.storage.BlobInfo
import com.google.cloud.storage.Storage
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import com.google.cloud.storage.contrib.nio.testing.LocalStorageHelper
import nextflow.exception.MissingFileException
import spock.lang.Specification
/**
 * Test the collection of task outputs from a Google Cloud Storage work directory
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TaskFileCollectorGsTest extends Specification {

    Storage storage
    Path workDir

    def setup() {
        final options = LocalStorageHelper.getOptions()
        storage = options.getService()
        workDir = CloudStorageFileSystem.forBucket('bucket', CloudStorageConfiguration.DEFAULT, options).getPath('/work/ab/cdef')
    }

    private void createBlob(String name) {
        storage.create(BlobInfo.newBuilder(BlobId.of('bucket', name)).build(), 'hello'.bytes)
    }

    private TaskFileCollector collector(List<String> patterns, Map opts=[:]) {
        final task = Mock(TaskRun) {
            getTargetDir() >> workDir
            getStagedInputs() >> []
        }
        return new TaskFileCollector(patterns, opts, task)
    }

    def 'should strip the trailing slash of a directory output' () {
        given:
        createBlob('work/ab/cdef/outdir/hello.txt')

        when:
        def result = collector(['outdir/']).collect()
        then:
        result.size() == 1
        result[0].toString() == '/work/ab/cdef/outdir'
        workDir.relativize(result[0]).toString() == 'outdir'
    }

    def 'should strip multiple trailing slashes of a directory output' () {
        given:
        createBlob('work/ab/cdef/outdir/hello.txt')

        when:
        def result = collector(['outdir//']).collect()
        then:
        result*.toString() == ['/work/ab/cdef/outdir']
    }

    def 'should collect an empty directory backed by a placeholder object' () {
        // gcsfuse creates a placeholder object for each directory
        given:
        createBlob('work/ab/cdef/emptydir/')

        when:
        def result = collector(['emptydir/']).collect()
        then:
        result*.toString() == ['/work/ab/cdef/emptydir']
    }

    def 'should report a missing directory output' () {
        // gcs considers any path ending with a slash to exist,
        // see https://github.com/nextflow-io/nextflow/issues/7885
        when:
        collector(['outdir/']).collect()
        then:
        thrown(MissingFileException)

        when:
        def result = collector(['outdir/'], [optional: true]).collect()
        then:
        result == []
    }

}

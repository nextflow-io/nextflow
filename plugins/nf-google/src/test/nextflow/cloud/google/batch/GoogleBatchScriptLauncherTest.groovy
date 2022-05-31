/*
 * Copyright 2022, Google Inc.
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
 *
 */

package nextflow.cloud.google.batch

import java.nio.file.Paths

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GoogleBatchScriptLauncherTest extends Specification{

    @Unroll
    def 'should convert to container path' () {
        given:
        def launcher = new GoogleBatchScriptLauncher()

        expect:
        def path = CloudStorageFileSystem.forBucket(BUCKET).getPath(PATH)
        launcher.toContainerMount(path, PARENT) == EXPECTED
        and:
        launcher.getContainerMounts() == [MOUNTS] 
        where:
        BUCKET  | PATH          | PARENT    | EXPECTED                          | MOUNTS
        'foo'   | '/'           | false     | Paths.get('/mnt/foo')             | '/mnt/foo:/mnt/foo:rw'
        'foo'   | '/some/dir'   | false     | Paths.get('/mnt/foo/some/dir')    | '/mnt/foo/some/dir:/mnt/foo/some/dir:rw'
        'foo'   | '/some/dir'   | true      | Paths.get('/mnt/foo/some/dir')    | '/mnt/foo/some:/mnt/foo/some:rw'

    }


    def 'should compute vol mounts' () {
        given:
        def launcher = new GoogleBatchScriptLauncher()
        def PATH1 = CloudStorageFileSystem.forBucket('alpha').getPath('/data/sample1.bam')
        def PATH2 = CloudStorageFileSystem.forBucket('alpha').getPath('/data/sample2.bam')
        def PATH3 = CloudStorageFileSystem.forBucket('omega').getPath('/data/sample3.bam')

        expect:
        launcher.toContainerMount(PATH1) == Paths.get('/mnt/alpha/data/sample1.bam')
        launcher.toContainerMount(PATH2) == Paths.get('/mnt/alpha/data/sample2.bam')
        launcher.toContainerMount(PATH3) == Paths.get('/mnt/omega/data/sample3.bam')
        and:
        launcher.taskVolumes.size() == 2
        and:
        def vol0 = launcher.taskVolumes.get(0)
        vol0.gcs == [remotePath:'alpha']
        vol0.mountPath == '/mnt/alpha'
        vol0.mountOptions == ['-o rw,allow_other', '-implicit-dirs']
        and:
        def vol1 = launcher.taskVolumes.get(1)
        vol1.gcs == [remotePath:'omega']
        vol1.mountPath == '/mnt/omega'
        vol1.mountOptions == ['-o rw,allow_other', '-implicit-dirs']

        and:
        launcher.containerMounts.size() == 2
        and:
        launcher.containerMounts[0] == '/mnt/alpha/data:/mnt/alpha/data:rw'
        launcher.containerMounts[1] == '/mnt/omega/data/sample3.bam:/mnt/omega/data/sample3.bam:rw'
    }

}

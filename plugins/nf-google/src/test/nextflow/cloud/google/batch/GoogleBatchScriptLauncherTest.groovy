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
import nextflow.cloud.google.GoogleOpts
import nextflow.cloud.google.batch.client.BatchConfig
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
        BUCKET  | PATH          | PARENT    | EXPECTED                              | MOUNTS
        'foo'   | '/'           | false     | Paths.get('/mnt/disks/foo')           | '/mnt/disks/foo:/mnt/disks/foo:rw'
        'foo'   | '/some/dir'   | false     | Paths.get('/mnt/disks/foo/some/dir')  | '/mnt/disks/foo/some/dir:/mnt/disks/foo/some/dir:rw'
        'foo'   | '/some/dir'   | true      | Paths.get('/mnt/disks/foo/some/dir')  | '/mnt/disks/foo/some:/mnt/disks/foo/some:rw'
    }

    def 'should compute volume mounts' () {
        given:
        def launcher = new GoogleBatchScriptLauncher()
        launcher.config = Mock(BatchConfig) {
            getGoogleOpts() >> Mock(GoogleOpts) {
                getProjectId() >> 'my-project'
                getEnableRequesterPaysBuckets() >> true
            }
        }
        and:
        def PATH1 = CloudStorageFileSystem.forBucket('alpha').getPath('/data/sample1.bam')
        def PATH2 = CloudStorageFileSystem.forBucket('alpha').getPath('/data/sample2.bam')
        def PATH3 = CloudStorageFileSystem.forBucket('omega').getPath('/data/sample3.bam')

        expect:
        launcher.toContainerMount(PATH1) == Paths.get('/mnt/disks/alpha/data/sample1.bam')
        launcher.toContainerMount(PATH2) == Paths.get('/mnt/disks/alpha/data/sample2.bam')
        launcher.toContainerMount(PATH3) == Paths.get('/mnt/disks/omega/data/sample3.bam')

        and:
        def containerMounts = launcher.getContainerMounts()
        and:
        containerMounts.size() == 2
        containerMounts[0] == '/mnt/disks/alpha/data:/mnt/disks/alpha/data:rw'
        containerMounts[1] == '/mnt/disks/omega/data/sample3.bam:/mnt/disks/omega/data/sample3.bam:rw'

        and:
        def volumes = launcher.getVolumes()
        and:
        volumes.size() == 2
        volumes[0].getGcs().getRemotePath() == 'alpha'
        volumes[0].getMountPath() == '/mnt/disks/alpha'
        volumes[0].getMountOptionsList() == ['-o rw', '-implicit-dirs', '--billing-project my-project']
        volumes[1].getGcs().getRemotePath() == 'omega'
        volumes[1].getMountPath() == '/mnt/disks/omega'
        volumes[1].getMountOptionsList() == ['-o rw', '-implicit-dirs', '--billing-project my-project']
    }

}

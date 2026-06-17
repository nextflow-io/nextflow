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

package nextflow.cloud.google.batch

import java.nio.file.Paths

import nextflow.cloud.google.batch.client.BatchConfig
import nextflow.processor.TaskBean
import spock.lang.Specification

class GoogleBatchFileCopyStrategyTest extends Specification {

    def 'should map container mount path to gs uri' () {
        expect:
        GoogleBatchFileCopyStrategy.toGsUriFromContainerMount(Paths.get('/mnt/disks/mybucket/work/dir')) == 'gs://mybucket/work/dir'
        GoogleBatchFileCopyStrategy.toGsUriFromContainerMount(Paths.get('/mnt/disks/onlybucket')) == 'gs://onlybucket/'
    }

    def 'should build gs download staging line when gcloud copy transport' () {
        given:
        def batch = new BatchConfig([
            stageInCopyTransport: 'gcloud'
        ])
        def bean = Mock(TaskBean) {
            getWorkDir() >> Paths.get('/mnt/disks/w/x')
            getTargetDir() >> Paths.get('/mnt/disks/w/x')
            getStageInMode() >> 'copy'
        }
        def path = Paths.get('/mnt/disks/b/data/in.txt')
        def copy = new GoogleBatchFileCopyStrategy(bean, batch)

        expect:
        copy.stageInputFile(path, 'in.txt') == 'downloads+=("nxf_gs_download \'gs://b/data/in.txt\' in.txt")'
    }

    def 'should build gs download staging line with retry when maxTransferAttempts > 1' () {
        given:
        def batch = new BatchConfig([
            stageInCopyTransport: 'gcloud',
            maxTransferAttempts: 3
        ])
        def bean = Mock(TaskBean) {
            getWorkDir() >> Paths.get('/mnt/disks/w/x')
            getTargetDir() >> Paths.get('/mnt/disks/w/x')
            getStageInMode() >> 'copy'
        }
        def path = Paths.get('/mnt/disks/b/data/in.txt')
        def copy = new GoogleBatchFileCopyStrategy(bean, batch)

        expect:
        copy.stageInputFile(path, 'in.txt') == 'downloads+=("nxf_cp_retry nxf_gs_download \'gs://b/data/in.txt\' in.txt")'
    }

    def 'should build gs upload unstage script' () {
        given:
        def batch = new BatchConfig([stageOutCopyTransport: 'gcloud'])
        def wd = Paths.get('/mnt/disks/foo/wd')
        def bean = Mock(TaskBean) {
            getWorkDir() >> wd
            getTargetDir() >> wd
            getStageOutMode() >> null
        }
        def copy = new GoogleBatchFileCopyStrategy(bean, batch)

        when:
        def script = copy.getUnstageOutputFilesScript(['out.txt'], wd)
        then:
        script.trim() == '''
                    uploads=()
                    IFS=$'\\n'
                    for name in $(eval "ls -1d out.txt" | sort | uniq); do
                        uploads+=("nxf_gs_upload '$name' 'gs://foo/wd'")
                    done
                    unset IFS
                    nxf_parallel "${uploads[@]}"
                    '''
                    .stripIndent().trim()
    }

    def 'should use posix symlink when stageIn is not cli copy' () {
        given:
        def batch = new BatchConfig([stageOutCopyTransport: 'gcloud'])
        def wd = Paths.get('/mnt/disks/foo/wd')
        def bean = Mock(TaskBean) {
            getWorkDir() >> wd
            getTargetDir() >> wd
            getStageInMode() >> 'symlink'
        }
        def path = Paths.get('/mnt/disks/foo/bucket/data.txt')
        def copy = new GoogleBatchFileCopyStrategy(bean, batch)

        expect:
        copy.stageInputFile(path, 'data.txt') == 'ln -s /mnt/disks/foo/bucket/data.txt data.txt'
    }

    def 'should use posix move when stageOutMode is move' () {
        given:
        def batch = new BatchConfig([stageOutCopyTransport: 'gcloud'])
        def wd = Paths.get('/mnt/disks/b/work')
        def td = Paths.get('/mnt/disks/b/out')
        def bean = Mock(TaskBean) {
            getWorkDir() >> wd
            getTargetDir() >> td
            getStageOutMode() >> 'move'
        }
        def copy = new GoogleBatchFileCopyStrategy(bean, batch)

        when:
        def script = copy.getUnstageOutputFilesScript(['x.txt'], td)
        then:
        script.contains('nxf_fs_move')
        script.contains('/mnt/disks/b/out')
    }

    def 'should use posix copy when stageOutCopyTransport is posix' () {
        given:
        def batch = new BatchConfig([stageOutCopyTransport: 'posix'])
        def wd = Paths.get('/mnt/disks/b/work')
        def bean = Mock(TaskBean) {
            getWorkDir() >> wd
            getTargetDir() >> wd
            getStageOutMode() >> 'copy'
        }
        def copy = new GoogleBatchFileCopyStrategy(bean, batch)

        when:
        def script = copy.getUnstageOutputFilesScript(['out.txt'], wd)
        then:
        script.contains('nxf_fs_copy')
    }
}

/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.cloud.google.lifesciences


import nextflow.cloud.google.GoogleSpecification
import nextflow.processor.TaskBean
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GoogleLifeSciencesFileCopyStrategyTest extends GoogleSpecification {


    def 'create stage files' () {
        given:
        def bean = Mock(TaskBean) {
            getWorkDir() >> mockGsPath('gs://my-bucket/work/xx/yy')
        }
        def handler = Mock(GoogleLifeSciencesTaskHandler) {
            getExecutor() >> Mock(GoogleLifeSciencesExecutor) {
                getConfig() >> Mock(GoogleLifeSciencesConfig)
            }
        }
        and:
        def strategy = new GoogleLifeSciencesFileCopyStrategy(bean, handler)

        // file with the same name
        when:
        def inputs = ['foo.txt': mockGsPath('gs://my-bucket/bar/foo.txt')]
        def result = strategy.getStageInputFilesScript(inputs)
        then:
        result == '''\
                echo start | gsutil -q cp  -c - gs://my-bucket/work/xx/yy/.command.begin
                downloads=(true)
                downloads+=("nxf_gs_download gs://my-bucket/bar/foo.txt foo.txt ")
                nxf_parallel "${downloads[@]}"
                '''.stripIndent()

        // file with a different name
        when:
        inputs = ['hola.txt': mockGsPath('gs://my-bucket/bar/foo.txt')]
        result = strategy.getStageInputFilesScript(inputs)
        then:
        result == '''\
                echo start | gsutil -q cp  -c - gs://my-bucket/work/xx/yy/.command.begin
                downloads=(true)
                downloads+=("nxf_gs_download gs://my-bucket/bar/foo.txt hola.txt ")
                nxf_parallel "${downloads[@]}"
                '''.stripIndent()

        // file name with blanks
        when:
        inputs = ['f o o.txt': mockGsPath('gs://my-bucket/bar/f o o.txt')]
        result = strategy.getStageInputFilesScript(inputs)
        then:
        result == '''\
                echo start | gsutil -q cp  -c - gs://my-bucket/work/xx/yy/.command.begin
                downloads=(true)
                downloads+=("nxf_gs_download gs://my-bucket/bar/f\\ o\\ o.txt f\\ o\\ o.txt ")
                nxf_parallel "${downloads[@]}"
                '''.stripIndent()

        // file in a sub-directory
        when:
        inputs = ['subdir/foo.txt': mockGsPath('gs://my-bucket/bar/foo.txt')]
        result = strategy.getStageInputFilesScript(inputs)
        then:
        result == '''\
                echo start | gsutil -q cp  -c - gs://my-bucket/work/xx/yy/.command.begin
                downloads=(true)
                downloads+=("nxf_gs_download gs://my-bucket/bar/foo.txt subdir/foo.txt ")
                nxf_parallel "${downloads[@]}"
                '''.stripIndent()

    }

    def 'create stage files using Requester Pays' () {
        given:
        def bean = Mock(TaskBean) {
            getWorkDir() >> mockGsPath('gs://my-bucket/work/xx/yy')
        }
        def handler = Mock(GoogleLifeSciencesTaskHandler) {
            getExecutor() >> Mock(GoogleLifeSciencesExecutor) {
                getConfig() >> GoogleLifeSciencesConfig.fromSession0([google:[project:'foo', region:'x', enableRequesterPaysBuckets:true]])
            }
        }
        and:
        def strategy = new GoogleLifeSciencesFileCopyStrategy(bean, handler)

        // file with the same name
        when:
        def inputs = ['foo.txt': mockGsPath('gs://my-bucket/bar/foo.txt')]
        def result = strategy.getStageInputFilesScript(inputs)
        then:
        result == '''\
                echo start | gsutil -q cp  -c - gs://my-bucket/work/xx/yy/.command.begin
                downloads=(true)
                downloads+=("nxf_gs_download gs://my-bucket/bar/foo.txt foo.txt foo")
                nxf_parallel "${downloads[@]}"
                '''.stripIndent()

        // file a directory name
        when:
        inputs = ['dir1': mockGsPath('gs://my-bucket/foo/dir1', true)]
        result = strategy.getStageInputFilesScript(inputs)
        then:
        result == '''\
                echo start | gsutil -q cp  -c - gs://my-bucket/work/xx/yy/.command.begin
                downloads=(true)
                downloads+=("nxf_gs_download gs://my-bucket/foo/dir1 dir1 foo")
                nxf_parallel "${downloads[@]}"
                '''.stripIndent()

    }

    def 'should unstage files' () {
        given:
        def WORK = gsPath('gs://my-bucket/work/xx/yy')
        def bean = Mock(TaskBean) {
            getWorkDir() >> WORK
        }
        def handler = Mock(GoogleLifeSciencesTaskHandler) {
            getExecutor() >> Mock(GoogleLifeSciencesExecutor) {
                getConfig() >> Mock(GoogleLifeSciencesConfig)
            }
        }
        and:
        def strategy = new GoogleLifeSciencesFileCopyStrategy(bean, handler)

        when:
        def result = strategy.getUnstageOutputFilesScript(['foo.txt', 'dirx/'], gsPath('gs://other/dir'))
        then:
        result ==
                '''\
                uploads=()
                IFS=$'\\n'
                for name in $(eval "ls -1d foo.txt dirx" | sort | uniq); do
                    uploads+=("nxf_gs_upload '$name' gs://other/dir")
                done
                unset IFS
                nxf_parallel "${uploads[@]}"
                '''
                .stripIndent()

        when:
        result = strategy.getUnstageOutputFilesScript(['fo o.txt'], gsPath('gs://other/dir x/'))
        then:
        result ==
                '''\
                uploads=()
                IFS=$'\\n'
                for name in $(eval "ls -1d fo\\ o.txt" | sort | uniq); do
                    uploads+=("nxf_gs_upload '$name' gs://other/dir\\ x")
                done
                unset IFS
                nxf_parallel "${uploads[@]}"
                '''
                .stripIndent()

    }


    def 'should create task env' () {
        given:
        def bean = Mock(TaskBean) {
            getWorkDir() >> mockGsPath('gs://my-bucket/work/xx/yy')
        }
        def handler = Mock(GoogleLifeSciencesTaskHandler) {
            getExecutor() >> Mock(GoogleLifeSciencesExecutor) {
                getConfig() >> Mock(GoogleLifeSciencesConfig)
            }
        }
        and:
        def strategy = new GoogleLifeSciencesFileCopyStrategy(bean, handler)


        when:
        def result = strategy.getEnvScript([FOO:1, BAR: 2], false)
        then:
        result == '''\
                export FOO="1"
                export BAR="2"
                '''.stripIndent()
    }

    def 'should create task with remote bin' () {
        given:
        def bean = Mock(TaskBean) {
            getWorkDir() >> mockGsPath('gs://my-bucket/work/xx/yy')
        }
        def handler = Mock(GoogleLifeSciencesTaskHandler) {
            getExecutor() >> Mock(GoogleLifeSciencesExecutor) {
                getConfig() >> Mock(GoogleLifeSciencesConfig) {
                    getRemoteBinDir() >> mockGsPath('gs://my-bucket/bin/d i r')
                }
            }
        }
        and:
        def strategy = new GoogleLifeSciencesFileCopyStrategy(bean, handler)


        when:
        def envScript = strategy.getEnvScript([FOO:1, BAR: 2, PATH: 3], false)
        then:
        envScript == '''\
                chmod +x /work/xx/yy/nextflow-bin/* || true
                export PATH=/work/xx/yy/nextflow-bin:$PATH
                export FOO="1"
                export BAR="2"
                '''.stripIndent()

        when:
        def stageScript = strategy.getStageInputFilesScript([:])
        then:
        stageScript == '''\
            echo start | gsutil -q cp  -c - gs://my-bucket/work/xx/yy/.command.begin
            mkdir -p /work/xx/yy/nextflow-bin
            gsutil -m -q cp -P -r gs://my-bucket/bin/d\\ i\\ r/* /work/xx/yy/nextflow-bin
            '''.stripIndent()
    }

    
}

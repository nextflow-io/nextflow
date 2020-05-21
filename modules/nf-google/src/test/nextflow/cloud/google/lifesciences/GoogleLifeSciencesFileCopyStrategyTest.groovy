/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
                gsutil -m -q cp gs://my-bucket/bar/foo.txt /work/xx/yy/foo.txt
                '''.stripIndent()

        // file with a different name
        when:
        inputs = ['hola.txt': mockGsPath('gs://my-bucket/bar/foo.txt')]
        result = strategy.getStageInputFilesScript(inputs)
        then:
        result == '''\
                echo start | gsutil -q cp  -c - gs://my-bucket/work/xx/yy/.command.begin
                gsutil -m -q cp gs://my-bucket/bar/foo.txt /work/xx/yy/hola.txt
                '''.stripIndent()

        // file name with blanks
        when:
        inputs = ['f o o.txt': mockGsPath('gs://my-bucket/bar/f o o.txt')]
        result = strategy.getStageInputFilesScript(inputs)
        then:
        result == '''\
                echo start | gsutil -q cp  -c - gs://my-bucket/work/xx/yy/.command.begin
                gsutil -m -q cp gs://my-bucket/bar/f\\ o\\ o.txt /work/xx/yy/f\\ o\\ o.txt
                '''.stripIndent()

        // file in a sub-directory
        when:
        inputs = ['subdir/foo.txt': mockGsPath('gs://my-bucket/bar/foo.txt')]
        result = strategy.getStageInputFilesScript(inputs)
        then:
        result == '''\
                echo start | gsutil -q cp  -c - gs://my-bucket/work/xx/yy/.command.begin
                mkdir -p /work/xx/yy/subdir
                gsutil -m -q cp gs://my-bucket/bar/foo.txt /work/xx/yy/subdir/foo.txt
                '''.stripIndent()

        // file a directory name
        when:
        inputs = ['dir1': mockGsPath('gs://my-bucket/foo/dir1', true)]
        result = strategy.getStageInputFilesScript(inputs)
        then:
        result == '''\
                echo start | gsutil -q cp  -c - gs://my-bucket/work/xx/yy/.command.begin
                gsutil -m -q cp -R gs://my-bucket/foo/dir1/ /work/xx/yy
                '''.stripIndent()

        // stage file is a directory with a different name
        when:
        inputs = ['dir2': mockGsPath('gs://my-bucket/foo/dir1', true)]
        result = strategy.getStageInputFilesScript(inputs)
        then:
        result == '''\
                echo start | gsutil -q cp  -c - gs://my-bucket/work/xx/yy/.command.begin
                gsutil -m -q cp -R gs://my-bucket/foo/dir1/ /work/xx/yy
                mv /work/xx/yy/dir1 /work/xx/yy/dir2
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
                gsutil -m -q -u foo cp gs://my-bucket/bar/foo.txt /work/xx/yy/foo.txt
                '''.stripIndent()

        // file a directory name
        when:
        inputs = ['dir1': mockGsPath('gs://my-bucket/foo/dir1', true)]
        result = strategy.getStageInputFilesScript(inputs)
        then:
        result == '''\
                echo start | gsutil -q cp  -c - gs://my-bucket/work/xx/yy/.command.begin
                gsutil -m -q -u foo cp -R gs://my-bucket/foo/dir1/ /work/xx/yy
                '''.stripIndent()

    }

    def 'should unstage files' () {
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
        def result = strategy.getUnstageOutputFilesScript(['foo.txt'], mockGsPath('gs://other/dir'))
        then:
        result == '''\
                IFS=$'\\n'; for name in $(eval "ls -1d foo.txt" 2>/dev/null);do gsutil -m -q cp -R $name gs://other/dir/\$name; done; unset IFS
                '''
                .stripIndent()

        when:
        result = strategy.getUnstageOutputFilesScript(['fo o.txt'], mockGsPath('gs://other/dir x/'))
        then:
        result == '''\
                IFS=$'\\n'; for name in $(eval "ls -1d fo\\ o.txt" 2>/dev/null);do gsutil -m -q cp -R $name gs://other/dir\\ x/\$name; done; unset IFS
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
                chmod +x /work/xx/yy/nextflow-bin/*
                export PATH=/work/xx/yy/nextflow-bin:$PATH
                export BAR="2"
                export FOO="1"
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

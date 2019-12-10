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

import java.nio.file.Files

import nextflow.Session
import nextflow.cloud.google.GoogleSpecification
import nextflow.processor.TaskBean
import nextflow.util.MustacheTemplateEngine
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GoogleLifeSciencesScriptLauncherTest extends GoogleSpecification {


    def 'should create task env' () {

        given:
        def session = new Session()
        def bucket = mockGsPath('gs://bucket/work/xx/yy')
        def binDir = mockGsPath('gs://bucket/bin')
        def handler = Mock(GoogleLifeSciencesTaskHandler) {
            getExecutor() >> Mock(GoogleLifeSciencesExecutor) {
                getConfig() >> new GoogleLifeSciencesConfig(remoteBinDir: binDir)
            }
        }
        def bean = [
                name: 'Hello 1',
                workDir: bucket,
                targetDir: bucket,
                environment: [PATH:'/this:/that', FOO: 'xxx'],
                script: 'echo Hello world!' ] as TaskBean

        when:
        def binding = new GoogleLifeSciencesScriptLauncher(bean, handler).makeBinding()
        then:
        binding.touch_file == null
        binding.stage_cmd == null
        binding.unstage_cmd == null
        binding.task_env == '''\
                chmod +x /work/xx/yy/nextflow-bin/*
                export PATH=/work/xx/yy/nextflow-bin:$PATH
                export FOO="xxx"
                '''.stripIndent()
    }

    def 'should render launcher script' () {

        given:
        def WORK_DIR = mockGsPath('gs://bucket/work/dir')
        def folder = Files.createTempDirectory('test')
        def handler = Mock(GoogleLifeSciencesTaskHandler) {
            getExecutor() >> Mock(GoogleLifeSciencesExecutor) {
                getConfig() >> new GoogleLifeSciencesConfig()
            }
        }
        def bean = [name: 'Hello 1',
                    script: 'echo Hello world!',
                    workDir: WORK_DIR] as TaskBean
        /*
         * simple bash run
         */
        when:
        def wrapper = new GoogleLifeSciencesScriptLauncher(bean, handler) .buildNew0()

        then:
        wrapper == load('bash-wrapper-gcp.txt', [folder: WORK_DIR.toString()])

        cleanup:
        folder?.deleteDir()
    }

    private String load(String name, Map<String,String> binding=[:]) {
        def template = new File("src/test/nextflow/cloud/google/lifesciences/$name").text
        return binding ? new MustacheTemplateEngine().render(template, binding) : template
    }
}

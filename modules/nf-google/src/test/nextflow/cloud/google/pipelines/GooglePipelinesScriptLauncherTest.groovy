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

package nextflow.cloud.google.pipelines

import spock.lang.Specification

import java.nio.file.Path
import java.nio.file.Paths

import nextflow.processor.TaskBean
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GooglePipelinesScriptLauncherTest extends Specification {

    def 'should create task env' () {

        given:
        def bucket = Paths.get('/bucket/work')
        def config = new GooglePipelinesConfiguration(remoteBinDir: 'gs://bucket/bin' as Path)
        def handler = [:] as GooglePipelinesTaskHandler
        handler.pipelineConfiguration = config
        def bean = [
                name: 'Hello 1',
                workDir: bucket,
                targetDir: bucket,
                environment: [PATH:'/this:/that', FOO: 'xxx'],
                script: 'echo Hello world!' ] as TaskBean

        when:
        def binding = new GooglePipelinesScriptLauncher(bean, handler).makeBinding()
        then:
        binding.task_env == '''\
                chmod +x /bucket/work/nextflow-bin/*
                export PATH=/bucket/work/nextflow-bin:$PATH
                export FOO="xxx"
                '''.stripIndent()
    }
}

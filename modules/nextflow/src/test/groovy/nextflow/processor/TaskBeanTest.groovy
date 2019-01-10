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

package nextflow.processor

import java.nio.file.Paths

import nextflow.Session
import nextflow.container.ContainerConfig
import nextflow.util.MemoryUnit
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskBeanTest extends Specification {

    def 'should create a bean object' () {

        given:

        def session = Mock(Session)
        session.getStatsEnabled() >> true
        session.getWorkDir() >> Paths.get('/work/dir')
        session.getBinDir() >> Paths.get('/bin/dir')

        def process = Mock(TaskProcessor)
        process.getConfig() >> ([stageInMode: 'link', stageOutMode: 'rsync'] as ProcessConfig)
        process.getSession() >> session

        def config = new TaskConfig()
        config.module = ['blast/1.1']
        config.shell = ['bash', '-x']
        config.beforeScript = 'before do this'
        config.afterScript = 'after do that'
        config.memory = '1GB'

        def task = Mock(TaskRun)
        task.getId() >> '123'
        task.getName() >> 'Hello'
        task.getStdin() >> 'input from stdin'
        task.getScratch() >> '/tmp/x'
        task.getWorkDir() >> Paths.get('/work/dir')
        task.getTargetDir() >> Paths.get('/target/dir')
        task.getScript() >>  'echo Ciao mondo'

        task.getConfig() >> config
        task.getProcessor() >> process
        task.getInputFilesMap() >> [file_1: Paths.get('/file/one'), file_2: Paths.get('/file/two')]
        task.getOutputFilesNames() >> [ 'simple.txt', 'my/path/file.bam' ]
        task.getTargetDir() >> Paths.get('/target/work/dir')
        task.getEnvironment() >> [alpha: 'one', beta: 'xxx', gamma: 'yyy']
        task.getContainer() >> 'busybox:latest'
        task.getContainerConfig() >> [docker: true, registry: 'x']

        when:
        def bean = new TaskBean(task)

        then:
        bean.name == 'Hello'
        bean.input == 'input from stdin'
        bean.scratch == '/tmp/x'
        bean.workDir == Paths.get('/work/dir')
        bean.targetDir == Paths.get('/target/dir')

        bean.environment == [alpha: 'one', beta:'xxx', gamma: 'yyy']
        bean.moduleNames ==  ['blast/1.1']
        bean.shell ==  ['bash', '-x']
        bean.script == 'echo Ciao mondo'
        bean.beforeScript == 'before do this'
        bean.afterScript == 'after do that'

        bean.containerImage == 'busybox:latest'
        bean.containerConfig == [docker: true, registry: 'x'] as ContainerConfig
        bean.containerMemory == new MemoryUnit('1GB')
        bean.statsEnabled

        bean.inputFiles == [file_1: Paths.get('/file/one'), file_2: Paths.get('/file/two')]
        bean.outputFiles ==  [ 'simple.txt', 'my/path/file.bam' ]
        bean.workDir == Paths.get('/work/dir')
        bean.binDir == Paths.get('/bin/dir')
        bean.stageInMode == 'link'
        bean.stageOutMode == 'rsync'

    }

    def 'should clone task bean' () {

        given:
        def task = new TaskBean(
                name: 'Hello',
                environment: [A: 'Alpha', B: 'Beta'],
                moduleNames: ['x','y'],
                workDir: Paths.get('/a/b'),
                containerMemory: new MemoryUnit('2GB')
        )

        when:
        def copy = task.clone()

        then:
        copy.name == task.name
        copy.environment == task.environment
        copy.moduleNames == task.moduleNames
        copy.workDir == task.workDir
        copy.containerMemory == task.containerMemory
    }
}

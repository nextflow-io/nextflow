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

package nextflow.config.parser.v2

import nextflow.config.ConfigClosurePlaceholder
import nextflow.processor.TaskConfig
import nextflow.script.ProcessConfig
import nextflow.script.dsl.ProcessConfigBuilder
import nextflow.util.ConfigHelper
import spock.lang.Specification

/**
 * Verify the compile-time collection of the {@code task} directives referenced by a
 * dynamic directive value defined in the config file.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ConfigDirectiveRefsTest extends Specification {

    private TaskConfig taskConfig(String text) {
        final config = new ConfigParserV2().parse(text)
        return new TaskConfig(config.process as Map)
    }

    def 'should detect the memory reference in a dynamic ext value' () {
        given:
        def config = taskConfig('''
            process.ext.args = { "-Xmx${task.memory.toGiga()}g" }
            ''')

        expect:
        config.isDirectiveReferenced('memory')
        !config.isDirectiveReferenced('cpus')
    }

    def 'should detect the memory reference in a dynamic directive value' () {
        given:
        def config = taskConfig('''
            process.beforeScript = { "echo ${task.memory}" }
            ''')

        expect:
        config.isDirectiveReferenced('memory')
    }

    private static final String SELECTOR_CONFIG = '''
        process {
            withName: FOO {
                ext.args = { "-Xmx${task.memory.toGiga()}g" }
            }
        }
        '''

    private TaskConfig taskConfigFor(String text, String processName) {
        final config = new ConfigParserV2().parse(text)
        final processConfig = new ProcessConfig([:])
        new ProcessConfigBuilder(processConfig)
            .applyConfig(config.process as Map, processName, processName, processName)
        return processConfig.createTaskConfig()
    }

    def 'should detect the memory reference for a process matching the selector' () {
        expect:
        taskConfigFor(SELECTOR_CONFIG, 'FOO').isDirectiveReferenced('memory')
    }

    def 'should not detect the memory reference for a process not matching the selector' () {
        expect:
        // the un-applied `withName:` blocks are copied verbatim into the process config,
        // they must not be mistaken for a directive of this process
        !taskConfigFor(SELECTOR_CONFIG, 'BAR').isDirectiveReferenced('memory')
    }

    def 'should render a dynamic directive as a string placeholder' () {
        when:
        def config = new ConfigParserV2()
            .setRenderClosureAsString(true)
            .parse('''
                process.clusterOptions = { "--mem ${task.memory.toMega()}" }
                ''')

        then:
        // `nextflow config` and `kuberun` render the config back as text -- the directive
        // reference must not get in the way of the placeholder substitution
        config.process.clusterOptions instanceof ConfigClosurePlaceholder
        ConfigHelper.toCanonicalString(config).contains('task.memory.toMega()')
    }

    def 'should report a config directive reference through the task run' () {
        given:
        def task = new nextflow.processor.TaskRun(config: taskConfig('''
            process.ext.args = { "-Xmx${task.memory.toGiga()}g" }
            '''))

        expect:
        task.isDirectiveReferenced('memory')
        !task.isDirectiveReferenced('cpus')
    }

    def 'should not detect anything when no directive is referenced' () {
        given:
        def config = taskConfig('''
            process.ext.args = { "--threads ${task.cpus}" }
            process.cpus = 4
            ''')

        expect:
        !config.isDirectiveReferenced('memory')
        config.isDirectiveReferenced('cpus')
    }
}

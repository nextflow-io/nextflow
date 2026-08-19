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
import nextflow.processor.TaskRun
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

    def 'should detect a reference nested in a list or map value' () {
        given:
        def config = taskConfig('''
            process.publishDir = [[path: { "${task.memory}" }]]
            ''')

        expect:
        config.isDirectiveReferenced('memory')
    }

    def 'should ignore the task properties that are not directives' () {
        given:
        def config = taskConfig('''
            process.ext.args = { "--seed ${task.attempt} ${task.exitStatus}" }
            ''')

        expect:
        !config.isDirectiveReferenced('attempt')
        !config.isDirectiveReferenced('exitStatus')
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

    def 'should still resolve the wrapped directive value' () {
        given:
        // the wrapper must be transparent: LazyMap clones and calls the value with `task`
        // as delegate, exactly as it does for a closure the user wrote
        def config = taskConfig('''
            process.memory = '8 GB'
            process.beforeScript = { "echo -Xmx${task.memory.toGiga()}g" }
            ''')
        config.setContext([:])

        expect:
        config.getBeforeScript() == 'echo -Xmx8g'
    }

    // --- source directive scoping ---

    def 'should restrict the check to the given source directives' () {
        given:
        def config = taskConfig('''
            process.clusterOptions = { "-l h_vmem=${task.memory}" }
            ''')

        expect:
        // the union over every source directive still sees it
        config.isDirectiveReferenced('memory')
        and:
        // ..but a caller that only renders `ext`/`beforeScript` into the command must not,
        // `clusterOptions` is a grid-executor submit option that never reaches the command
        !config.isDirectiveReferenced('memory', ['ext','beforeScript','afterScript'])
    }

    def 'should detect a reference made by one of the given source directives' () {
        given:
        def config = taskConfig('''
            process.ext.args = { "-Xmx${task.memory.toGiga()}g" }
            ''')

        expect:
        config.isDirectiveReferenced('memory', ['ext','beforeScript','afterScript'])
        !config.isDirectiveReferenced('memory', ['beforeScript'])
    }

    // --- config selectors ---

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

    // --- config rendering ---

    def 'should render a dynamic directive as a string placeholder' () {
        when:
        def config = new ConfigParserV2()
            .setRenderClosureAsString(true)
            .parse('''
                process.clusterOptions = { "--mem ${task.memory.toMega()}" }
                ''')

        then:
        // `nextflow config` and `kuberun` render the config back as text -- the directive
        // reference collection must not get in the way of the placeholder substitution
        config.process.clusterOptions instanceof ConfigClosurePlaceholder
        ConfigHelper.toCanonicalString(config).contains('task.memory.toMega()')
    }

    // --- task run ---

    private TaskRun taskRun(String text) {
        // the script variable references are the concern of ProcessDirectiveRefsTest, here
        // they are stubbed empty so that only the config half of the union is exercised
        def task = Spy(TaskRun)
        task.config = taskConfig(text)
        task.getVariableNames() >> (Collections.<String>emptySet())
        return task
    }

    def 'should report a config directive reference through the task run' () {
        given:
        def task = taskRun('''
            process.ext.args = { "-Xmx${task.memory.toGiga()}g" }
            ''')

        expect:
        task.isDirectiveReferenced('memory')
        !task.isDirectiveReferenced('cpus')
    }

    def 'should not report a config reference made by a directive that does not reach the command' () {
        given:
        // `clusterOptions` is a grid-executor submit option: it never reaches the rendered
        // command, so a reference made there is not a command dependency
        def task = taskRun('''
            process.clusterOptions = { "-l h_vmem=${task.memory}" }
            ''')

        expect:
        !task.isDirectiveReferenced('memory')
    }
}

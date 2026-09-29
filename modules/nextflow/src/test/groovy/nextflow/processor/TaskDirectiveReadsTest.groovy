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

package nextflow.processor

import java.nio.file.Files

import ch.artecat.grengine.Grengine
import nextflow.Session
import nextflow.config.parser.v1.ConfigParserV1
import nextflow.config.parser.v2.ConfigParserV2
import nextflow.script.BodyDef
import nextflow.script.ProcessConfig
import nextflow.script.dsl.ProcessConfigBuilder
import spock.lang.Specification
import spock.lang.Unroll

/**
 * Verify the directives read while the task command is rendered are reported by
 * {@link TaskRun#isDirectiveReferenced}.
 *
 * @author Ben Sherman <ben.sherman@seqera.io>
 */
class TaskDirectiveReadsTest extends Specification {

    private TaskRun taskWith(Map directives, BodyDef body) {
        return taskWith(new TaskConfig(directives), body)
    }

    private TaskRun taskWith(TaskConfig config, BodyDef body) {
        final task = unresolvedTaskWith(config)
        task.resolve(body)
        return task
    }

    private TaskRun unresolvedTaskWith(TaskConfig config) {
        final processor = new TaskProcessor()
        processor.@session = new Session()
        processor.@grengine = new Grengine()

        final task = new TaskRun(processor: processor, config: config)
        task.context = new TaskContext(holder: [:])
        task.config.setContext(task.context)
        return task
    }

    /**
     * Build the task config as the engine does, by applying a parsed config file to the process
     * config -- so the directive values are the closures the parser produced.
     */
    private TaskConfig taskConfigFrom(parser, String text, String processName) {
        final config = parser.parse(text)
        final processConfig = new ProcessConfig([:])
        new ProcessConfigBuilder(processConfig)
            .applyConfig(config.process as Map, processName, processName, processName)
        return processConfig.createTaskConfig()
    }

    def 'should report a directive read by the task script' () {
        when:
        def task = taskWith(
            [memory: '8 GB', cpus: 4],
            new BodyDef({-> "java -Xmx${task.memory.toGiga()}g -jar app.jar"}, 'java ...', 'script') )

        then:
        task.script == 'java -Xmx8g -jar app.jar'
        and:
        task.isDirectiveReferenced('memory')
        !task.isDirectiveReferenced('cpus')
    }

    def 'should report a directive read behind a dynamic ext directive' () {
        when:
        // the script mentions `task.ext.args`, the memory reference is inside the closure --
        // this is the common nf-core config idiom `ext.args = { ... }`
        def task = taskWith(
            [memory: '8 GB', ext: [args: {"-Xmx${task.memory.toGiga()}g"}]],
            new BodyDef({-> "java ${task.ext.args} -jar app.jar"}, 'java ...', 'script') )

        then:
        task.script == 'java -Xmx8g -jar app.jar'
        and:
        task.isDirectiveReferenced('memory')
    }

    def 'should report a directive read by a shell block' () {
        when:
        def task = taskWith(
            [memory: '8 GB'],
            new BodyDef({-> 'java -Xmx!{task.memory.toGiga()}g -jar app.jar'}, 'java ...', 'shell') )

        then:
        task.script == 'java -Xmx8g -jar app.jar'
        and:
        task.isDirectiveReferenced('memory')
    }

    def 'should report a directive read by a template file' () {
        given:
        def file = Files.createTempDirectory('test').resolve('foo.sh')
        file.text = 'java -Xmx${task.memory.toGiga()}g -jar app.jar'

        when:
        def task = taskWith(
            [memory: '8 GB'],
            new BodyDef({-> template(file)}, 'template(file)', 'script') )

        then:
        task.script == 'java -Xmx8g -jar app.jar'
        and:
        task.isDirectiveReferenced('memory')
    }

    def 'should not report a directive accessed outside the tracked action' () {
        given:
        // the engine accesses the directives all the time e.g. the executor asking for the
        // memory to request -- only the accesses made while rendering the command count
        def task = taskWith(
            [memory: '8 GB'],
            new BodyDef({-> 'echo hello'}, 'echo hello', 'script') )

        when:
        task.config.getMemory()

        then:
        !task.isDirectiveReferenced('memory')
    }

    def 'should not report any directive when the command interpolates none' () {
        when:
        def task = taskWith(
            [memory: '8 GB', cpus: 4],
            new BodyDef({-> 'echo hello'}, 'echo hello', 'script') )

        then:
        // pins the absence of false positives: what the engine itself touches while the
        // command is rendered must not be reported as a dependency of it
        !task.isDirectiveReferenced('memory')
        !task.isDirectiveReferenced('cpus')
    }

    def 'should answer false before the command has been rendered' () {
        given:
        def task = unresolvedTaskWith(new TaskConfig(memory: '8 GB'))

        expect:
        !task.isDirectiveReferenced('memory')
    }

    def 'should carry the accessed directives over to a task copy' () {
        given:
        // a copy keeps the rendered command, so it depends on the same directives -- the
        // retryable/spot path in TaskProcessor copies the task without resolving it again
        def task = taskWith(
            [memory: '8 GB'],
            new BodyDef({-> "java -Xmx${task.memory.toGiga()}g -jar app.jar"}, 'java ...', 'script') )

        when:
        def copy = task.makeCopy()

        then:
        copy.script == 'java -Xmx8g -jar app.jar'
        and:
        copy.isDirectiveReferenced('memory')
    }

    @Unroll
    def 'should report a directive accessed behind a dynamic ext value from the config file [#parser.class.simpleName, #processName]' () {
        given:
        // the nf-core idiom: the script only mentions `task.ext.args`, the memory reference
        // sits inside the closure the config parser produced
        def config = taskConfigFrom(parser, '''
            process {
                memory = '8 GB'
                withName: FOO {
                    ext.args = { "-Xmx${task.memory.toGiga()}g" }
                }
            }
            ''', processName)

        when:
        def task = taskWith(config, new BodyDef({-> "java ${task.ext.args} -jar app.jar"}, 'java ...', 'script'))

        then:
        task.script == expected
        and:
        task.isDirectiveReferenced('memory') == referenced

        where:
        parser               | processName | expected                   | referenced
        new ConfigParserV1() | 'FOO'       | 'java -Xmx8g -jar app.jar' | true
        new ConfigParserV2() | 'FOO'       | 'java -Xmx8g -jar app.jar' | true
        new ConfigParserV1() | 'BAR'       | 'java null -jar app.jar'   | false
        new ConfigParserV2() | 'BAR'       | 'java null -jar app.jar'   | false
    }
}

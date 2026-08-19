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
import nextflow.script.BodyDef
import spock.lang.Specification

/**
 * Verify the directives read while the task command is rendered are reported by
 * {@link TaskRun#isDirectiveReferenced}.
 *
 * @author Ben Sherman <ben.sherman@seqera.io>
 */
class TaskDirectiveReadsTest extends Specification {

    private TaskRun taskWith(Map directives, BodyDef body) {
        final processor = new TaskProcessor()
        processor.@session = new Session()
        processor.@grengine = new Grengine()

        final task = new TaskRun(processor: processor, config: new TaskConfig(directives))
        task.context = new TaskContext(holder: [:])
        task.config.setContext(task.context)
        task.resolve(body)
        return task
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

    def 'should not report a directive read outside the rendering window' () {
        given:
        // the engine reads the directives all the time e.g. the executor asking for the
        // memory to request -- only the reads made while rendering the command count
        def task = taskWith(
            [memory: '8 GB'],
            new BodyDef({-> 'echo hello'}, 'echo hello', 'script') )

        when:
        task.config.getMemory()

        then:
        !task.isDirectiveReferenced('memory')
    }
}

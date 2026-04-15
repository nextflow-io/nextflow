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

import java.nio.file.Path

import nextflow.exception.IllegalArityException
import nextflow.script.ScriptType
import nextflow.script.params.FileOutParam
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskOutputResolverV1Test extends Specification {

    @Unroll
    def 'should collect output files' () {
        given:
        def context = new TaskContext(holder: new HashMap<String, Object>())
        def task = new TaskRun(
                name: 'foo',
                type: ScriptType.SCRIPTLET,
                context: context,
                config: new TaskConfig(),
                workDir: Path.of('/work')
        )
        and:
        def resolver = Spy(new TaskOutputResolverV1(task))

        when:
        def param = new FileOutParam(new Binding(), [])
                .setPathQualifier(true)
                .optional(OPTIONAL)
                .bind(FILE_NAME) as FileOutParam
        if( ARITY )
            param.setArity(ARITY)
        and:
        resolver.resolve(param)
        then:
        resolver.collectOutFiles0(_,_) >> RESULTS
        and:
        task.getOutputs().get(param) == EXPECTED

        where:
        FILE_NAME       | RESULTS                                   | OPTIONAL  | ARITY         | EXPECTED
        'file.txt'      | [Path.of('/work/file.txt')]               | false     | null          | Path.of('/work/file.txt')
        '*'             | [Path.of('/work/file.txt')]               | false     | null          | Path.of('/work/file.txt')
        '*'             | [Path.of('/work/A'), Path.of('/work/B')]  | false     | null          | [Path.of('/work/A'), Path.of('/work/B')]
        '*'             | []                                        | true      | null          | []
        and:
        'file.txt'      | [Path.of('/work/file.txt')]               | false     | '1'           | Path.of('/work/file.txt')
        '*'             | [Path.of('/work/file.txt')]               | false     | '1'           | Path.of('/work/file.txt')
        '*'             | [Path.of('/work/file.txt')]               | false     | '1..*'        | [Path.of('/work/file.txt')]
        '*'             | [Path.of('/work/A'), Path.of('/work/B')]  | false     | '2'           | [Path.of('/work/A'), Path.of('/work/B')]
        '*'             | [Path.of('/work/A'), Path.of('/work/B')]  | false     | '1..*'        | [Path.of('/work/A'), Path.of('/work/B')]
        '*'             | []                                        | false     | '0..*'        | []
    }

    @Unroll
    def 'should report output file arity error' () {
        given:
        def context = new TaskContext(holder: new HashMap<String, Object>())
        def task = new TaskRun(
                name: 'foo',
                type: ScriptType.SCRIPTLET,
                context: context,
                config: new TaskConfig(),
                workDir: Path.of('/work')
        )
        and:
        def resolver = Spy(new TaskOutputResolverV1(task))

        when:
        def param = new FileOutParam(new Binding(), [])
                .setPathQualifier(true)
                .optional(OPTIONAL)
                .bind(FILE_NAME) as FileOutParam
        if( ARITY )
            param.setArity(ARITY)
        and:
        resolver.resolve(param)
        then:
        resolver.collectOutFiles0(_,_) >> RESULTS
        and:
        def e = thrown(EXCEPTION)
        e.message == ERROR

        where:
        FILE_NAME       | RESULTS                                   | OPTIONAL  | ARITY         | EXCEPTION             | ERROR
        'file.txt'      | [Path.of('/work/file.txt')]               | false     | '2'           | IllegalArityException | "Incorrect number of output files for process `foo` -- expected 2, found 1"
        '*'             | [Path.of('/work/file.txt')]               | false     | '2'           | IllegalArityException | "Incorrect number of output files for process `foo` -- expected 2, found 1"
        '*'             | [Path.of('/work/file.txt')]               | false     | '2..*'        | IllegalArityException | "Incorrect number of output files for process `foo` -- expected 2..*, found 1"
        '*'             | []                                        | true      | '1..*'        | IllegalArityException | "Incorrect number of output files for process `foo` -- expected 1..*, found 0"
    }

}

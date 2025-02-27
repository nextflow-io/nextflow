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

package nextflow.extension

import nextflow.Channel
import nextflow.Session
import spock.lang.Specification
import spock.lang.Unroll
import test.OutputCapture

import java.nio.file.Path

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DumpOpTest extends Specification {

    /*
     * Read more http://mrhaki.blogspot.com.es/2015/02/spocklight-capture-and-assert-system.html
     */
    @org.junit.Rule
    OutputCapture capture = new OutputCapture()

    def 'should dump channel items'() {

        given:
        new Session(dumpChannels: ['*'])

        when:
        def result = Channel.of(1, 2, 3).dump()
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP
        capture.toString().contains('[DUMP] 1')
        capture.toString().contains('[DUMP] 2')
        capture.toString().contains('[DUMP] 3')

    }

    def 'should dump channel items with closure'() {

        given:
        new Session(dumpChannels: ['*'])

        when:
        def result = Channel.of(1, 2, 3).dump { it * it }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP
        capture.toString().contains('[DUMP] 1')
        capture.toString().contains('[DUMP] 4')
        capture.toString().contains('[DUMP] 9')

    }

    def 'should print channel items with a tag'() {

        given:
        new Session(dumpChannels: ['*'])

        when:
        def result = Channel.of(1, 2, 3).dump(tag: 'foo')
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP
        capture.toString().contains('[DUMP: foo] 1')
        capture.toString().contains('[DUMP: foo] 2')
        capture.toString().contains('[DUMP: foo] 3')

    }

    @Unroll
    def 'should validate isEnabled when tag=#tag and names=#names'() {

        given:
        def op = new DumpOp(tag: tag, dumpNames: names.tokenize(','))

        expect:
        op.isEnabled() == expected

        where:
        tag   | names   | expected
        'foo' | 'foo'   | true
        'foo' | '*'     | true
        'foo' | 'bar'   | false
        'foo' | 'f'     | false
        'foo' | 'f*'    | true
        'bar' | 'f*'    | false
        'bar' | 'f*,b*' | true
        null  | '*'     | true

    }

    @Unroll
    def 'should pretty print channel with items #items'() {

        given:
        new Session(dumpChannels: ['*'])

        when:
        def result = Channel.of(items).dump(tag: 'foo', pretty: true)

        then:
        result.val
        result.val == Channel.STOP
        capture.toString().contains(expected)

        where:
        items                    | expected
        'a string'               | "[DUMP: foo] 'a string'"
        123123                   | '[DUMP: foo] 123123'
        ['a', 'b']               | '"a",'
        ['a': 'b']               | '"a": "b"'
        ['a': Path.of('b')]      | '"a": "b"'
        ['a': [b: Path.of('c')]] | '"b": "c"'
    }

    def 'should pretty print a complex channel'() {

        given:
        new Session(dumpChannels: ['*'])

        when:
        def result = Channel.of([
                id     : 'test',
                samples: [
                        [id: 'S1', path: '/path/to/s1'],
                        [id: 'S2', path: '/path/to/s2'],
                        [id: 'S3', path: '/path/to/s3']
                ]
        ]).dump(tag: 'foo', pretty: true)

        then:
        result.val
        result.val == Channel.STOP
        capture.toString().contains('"path": "/path/to/s1"')
    }

}

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

package nextflow.extension

import java.nio.file.Path

import nextflow.Session
import spock.lang.Specification
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class PublishOpTest extends Specification {

    def 'should normalize the target directory' () {
        given:
        def session = Mock(Session) { getOutputDir() >> Path.of('/work/results') }

        when:
        def op = new PublishOp(session, 'foo', null, [path: PATH])
        then:
        op.getTargetDir(null) == Path.of(EXPECTED)

        where:
        PATH        | EXPECTED
        '.'         | '/work/results'
        './'        | '/work/results'
        'bam'       | '/work/results/bam'
        './bam'     | '/work/results/bam'
        'bam/../bam'| '/work/results/bam'
    }

    def 'should normalize the target directory returned by a closure' () {
        given:
        def session = Mock(Session) { getOutputDir() >> Path.of('/work/results') }
        def resolver = { v -> '.' }

        when:
        def op = new PublishOp(session, 'foo', null, [path: '.', pathResolver: resolver])
        then:
        op.getTargetDir(null) == Path.of('/work/results')
    }

    def 'should map source files to target paths returned by a closure' () {
        given:
        def session = Mock(Session) {
            getOutputDir() >> Path.of('/work/results')
            getWorkDir() >> Path.of('/work')
        }
        def foo = Path.of('/work/ab/cdef/foo.txt')
        def bar = Path.of('/work/ab/cdef/bar.txt')
        def resolver = { v -> [(foo): 'foo/', (bar): 'bar/renamed.txt'] }

        when:
        def op = new PublishOp(session, 'foo', null, [path: '.', pathResolver: resolver])
        def saveAs = op.getTargetDir(null)
        then:
        saveAs instanceof Closure
        saveAs.call('foo.txt') == Path.of('/work/results/foo/foo.txt')
        saveAs.call('bar.txt') == Path.of('/work/results/bar/renamed.txt')
        saveAs.call('baz.txt') == null
    }

}

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
        op.getTargetDir(PATH) == Path.of(EXPECTED)

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
        def session = Mock(Session) {
            getOutputDir() >> Path.of('/work/results')
            getWorkDir() >> Path.of('/work')
        }
        def file = Path.of('/work/ab/1234/file.txt')
        def resolver = { v -> './bam/../bam' }

        when:
        def op = new PublishOp(session, 'foo', null, [path: '.', pathResolver: resolver])
        then:
        op.resolveTargets(file) == [(file): Path.of('/work/results/bam/file.txt')]
    }

    def 'should resolve the target of each file declared by a publish statement' () {
        given:
        def session = Mock(Session) {
            getOutputDir() >> Path.of('/work/results')
            getWorkDir() >> Path.of('/work')
        }
        def file1 = Path.of('/work/ab/1234/report.txt')
        def file2 = Path.of('/work/cd/5678/report.txt')
        def resolver = { v ->
            publish(v.alpha, 'reports/alpha.txt')
            publish(v.beta, 'reports/./beta.txt')
        }

        when:
        def op = new PublishOp(session, 'foo', null, [path: '.', pathResolver: resolver])
        then:
        op.resolveTargets([alpha: file1, beta: file2]) == [
            (file1): Path.of('/work/results/reports/alpha.txt'),
            (file2): Path.of('/work/results/reports/beta.txt')
        ]
    }

    def 'should publish nothing when all publish statements are no-ops' () {
        given:
        def session = Mock(Session) { getOutputDir() >> Path.of('/work/results') }
        def resolver = { v -> publish(v, null) }

        when:
        def op = new PublishOp(session, 'foo', null, [path: '.', pathResolver: resolver])
        then:
        op.resolveTargets(Path.of('/work/ab/cdef/out.txt')) == [:]
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
        then:
        op.resolveTargets(null) == [
            (foo): Path.of('/work/results/foo/foo.txt'),
            (bar): Path.of('/work/results/bar/renamed.txt')
        ]
    }

    def 'should publish files outside the work directory' () {
        given:
        def session = Mock(Session) {
            getOutputDir() >> Path.of('/work/results')
            getWorkDir() >> Path.of('/work')
        }
        def input = Path.of('/data/input.txt')
        def output = Path.of('/work/ab/cdef/sub/out.txt')

        when:
        def op = new PublishOp(session, 'foo', null, [path: 'txt', includeInputs: true])
        then:
        op.resolveTargets([input, output]) == [
            (input): Path.of('/work/results/txt/input.txt'),
            (output): Path.of('/work/results/txt/sub/out.txt')
        ]
    }

    def 'should not publish files outside the work directory by default' () {
        given:
        def session = Mock(Session) {
            getOutputDir() >> Path.of('/work/results')
            getWorkDir() >> Path.of('/work')
        }
        def input = Path.of('/data/input.txt')
        def output = Path.of('/work/ab/cdef/out.txt')

        when:
        def op = new PublishOp(session, 'foo', null, [path: 'txt'])
        then:
        op.resolveTargets([input, output]) == [
            (output): Path.of('/work/results/txt/out.txt')
        ]

        when:
        def resolver = { v -> publish(input, 'txt/'); publish(output, 'txt/') }
        op = new PublishOp(session, 'foo', null, [path: '.', pathResolver: resolver])
        then:
        op.resolveTargets(null) == [
            (output): Path.of('/work/results/txt/out.txt')
        ]
    }

    def 'should publish files outside the work directory in publish statements' () {
        given:
        def session = Mock(Session) {
            getOutputDir() >> Path.of('/work/results')
            getWorkDir() >> Path.of('/work')
        }
        def input = Path.of('/data/input.txt')
        def output = Path.of('/work/ab/cdef/out.txt')
        def resolver = { v -> publish(input, 'txt/'); publish(output, 'txt/') }

        when:
        def op = new PublishOp(session, 'foo', null, [path: '.', pathResolver: resolver, includeInputs: true])
        then:
        op.resolveTargets(null) == [
            (input): Path.of('/work/results/txt/input.txt'),
            (output): Path.of('/work/results/txt/out.txt')
        ]
    }

}

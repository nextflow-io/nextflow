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
import java.nio.file.Path

import nextflow.exception.MissingFileException
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskFileCollectorTest extends Specification {

    def 'should filter staged inputs'() {

        given:
        def workDir = Path.of('/work/dir')
        def task = Spy(new TaskRun(
                config: new TaskConfig(),
                workDir: workDir ))
        def collector = new TaskFileCollector([], [:], task)

        def FILE1 = workDir.resolve('alpha.txt')
        def FILE2 = workDir.resolve('beta.txt')
        def FILE3 = workDir.resolve('out/beta.txt')
        def FILE4 = workDir.resolve('gamma.fasta')

        when:
        def result = collector.excludeStagedInputs([ FILE1, FILE2, FILE3, FILE4 ])
        then:
        1 * task.getStagedInputs() >> [ 'beta.txt' ]
        and:
        result == [ FILE1, FILE3, FILE4 ]

    }


    private List<String> fetchResultFiles(Map opts=[:], String namePattern, Path folder) {
        def collector = new TaskFileCollector([], opts, Mock(TaskRun))
        return collector
            .fetchResultFiles(namePattern, folder)
            .collect { Path it -> folder.relativize(it).toString() }
    }

    def 'should return the list of output files'() {

        given:
        def folder = Files.createTempDirectory('test')
        folder.resolve('file1.txt').text = 'file 1'
        folder.resolve('file2.fa').text = 'file 2'
        folder.resolve('.hidden.fa').text = 'hidden'
        folder.resolve('dir1').mkdir()
        folder.resolve('dir1').resolve('file3.txt').text = 'file 3'
        folder.resolve('dir1')
        folder.resolve('dir1').resolve('dir2').mkdirs()
        folder.resolve('dir1').resolve('dir2').resolve('file4.fa').text = 'file '
        Files.createSymbolicLink( folder.resolve('dir_link'), folder.resolve('dir1') )
        and:
        def result

        when:
        result = fetchResultFiles('*.fa', folder )
        then:
        result == ['file2.fa']

        when:
        result = fetchResultFiles('*.fa', folder, type: 'file')
        then:
        result == ['file2.fa']

        when:
        result = fetchResultFiles('*.fa', folder, type: 'dir')
        then:
        result == []

        when:
        result = fetchResultFiles('**.fa', folder)
        then:
        result == ['dir1/dir2/file4.fa', 'dir_link/dir2/file4.fa', 'file2.fa']

        when:
        result = fetchResultFiles('**.fa', folder, followLinks: false)
        then:
        result == ['dir1/dir2/file4.fa', 'file2.fa']

        when:
        result = fetchResultFiles('**.fa', folder, maxDepth: 1)
        then:
        result == ['file2.fa']

        when:
        result = fetchResultFiles('*', folder)
        then:
        result == ['dir1', 'dir_link', 'file1.txt', 'file2.fa']

        when:
        result = fetchResultFiles('*', folder, type: 'dir')
        then:
        result == ['dir1', 'dir_link']

        when:
        result = fetchResultFiles('*', folder, type: 'file')
        then:
        result == ['file1.txt', 'file2.fa']

        when:
        result = fetchResultFiles('*', folder, type: 'file', hidden: true)
        then:
        result == ['.hidden.fa', 'file1.txt', 'file2.fa']

        when:
        result = fetchResultFiles('.*', folder)
        then:
        result == ['.hidden.fa']

        when:
        result = fetchResultFiles('file{1,2}.{txt,fa}', folder)
        then:
        result == ['file1.txt', 'file2.fa']

        cleanup:
        folder?.deleteDir()

    }

    def defaultCollector(Map opts) {
        return new TaskFileCollector([], opts, Mock(TaskRun))
    }

    def 'should create the map of path visit options'() {

        given:
        def collector

        when:
        collector = defaultCollector([:])
        then:
        collector.visitOptions('file.txt') == [type:'any', followLinks: true, maxDepth: null, hidden: false, relative: false]
        collector.visitOptions('path/**') == [type:'file', followLinks: true, maxDepth: null, hidden: false, relative: false]
        collector.visitOptions('.hidden_file') == [type:'any', followLinks: true, maxDepth: null, hidden: true, relative: false]

        when:
        collector = defaultCollector([type: 'dir'])
        then:
        collector.visitOptions('dir-name') == [type:'dir', followLinks: true, maxDepth: null, hidden: false, relative: false]

        when:
        collector = defaultCollector([hidden: true])
        then:
        collector.visitOptions('dir-name') == [type:'any', followLinks: true, maxDepth: null, hidden: true, relative: false]

        when:
        collector = defaultCollector([followLinks: false])
        then:
        collector.visitOptions('dir-name') == [type:'any', followLinks: false, maxDepth: null, hidden: false, relative: false]

        when:
        collector = defaultCollector([maxDepth: 5])
        then:
        collector.visitOptions('dir-name') == [type:'any', followLinks: true, maxDepth: 5, hidden: false, relative: false]
    }

    def 'should collect output files' () {
        given:
        def task = new TaskRun(
                name: 'foo',
                config: new TaskConfig(),
                workDir: Path.of('/work') )
        and:
        def opts = [optional: OPTIONAL]
        def collector = Spy(new TaskFileCollector([FILE_NAME], opts, task))

        when:
        def result = collector.collect()
        then:
        collector.fetchResultFiles(_,_) >> RESULTS
        collector.checkFileExists(_) >> EXISTS
        and:
        result == EXPECTED

        where:
        FILE_NAME       | RESULTS                                   | EXISTS    | OPTIONAL  | EXPECTED
        'file.txt'      | []                                        | true      | false     | [Path.of('/work/file.txt')]
        '*'             | [Path.of('/work/file.txt')]               | true      | false     | [Path.of('/work/file.txt')]
        '*'             | [Path.of('/work/A'), Path.of('/work/B')]  | true      | false     | [Path.of('/work/A'), Path.of('/work/B')]
        '*'             | []                                        | true      | true      | []
    }

    @Unroll
    def 'should report missing output file error' () {
        given:
        def task = new TaskRun(
                name: 'foo',
                config: new TaskConfig(),
                workDir: Path.of('/work') )
        and:
        def opts = [optional: OPTIONAL]
        def collector = Spy(new TaskFileCollector([FILE_NAME], opts, task))

        when:
        collector.collect()
        then:
        collector.fetchResultFiles(_,_) >> RESULTS
        collector.checkFileExists(_) >> EXISTS
        and:
        def e = thrown(EXCEPTION)
        e.message == ERROR

        where:
        FILE_NAME       | RESULTS                                   | EXISTS    | OPTIONAL  | EXCEPTION             | ERROR
        'file.txt'      | null                                      | false     | false     | MissingFileException  | "Missing output file(s) `file.txt` expected by process `foo`"
        '*'             | []                                        | true      | false     | MissingFileException  | "Missing output file(s) `*` expected by process `foo`"

    }

}

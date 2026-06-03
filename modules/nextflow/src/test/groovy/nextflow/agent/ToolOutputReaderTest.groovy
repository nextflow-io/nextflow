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
package nextflow.agent

import java.nio.file.Files
import java.nio.file.Path

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ToolOutputReaderTest extends Specification {

    def 'should inline a small structured json file content'() {
        given:
        final dir = Files.createTempDirectory('test')
        final file = dir.resolve('stats.json')
        file.text = '{"n50":54321,"contigs":12}'

        when:
        final result = ToolOutputReader.readOrHandle(file, 1024 * 1024)

        then:
        result instanceof String
        result == '{"n50":54321,"contigs":12}'
    }

    def 'should return a path+note map when the file exceeds the cap'() {
        given:
        final dir = Files.createTempDirectory('test')
        final file = dir.resolve('big.json')
        file.text = '{"n50":54321,"contigs":12}'

        when:
        final result = ToolOutputReader.readOrHandle(file, 8)

        then:
        result instanceof Map
        (result as Map).path == file.toAbsolutePath().toString()
        ((result as Map).note as String).contains('not inlined')
    }

    def 'should return the path for a non text-like extension'() {
        given:
        final dir = Files.createTempDirectory('test')
        final file = dir.resolve('contigs.fa')
        file.text = '>contig1\nACGT'

        when:
        final result = ToolOutputReader.readOrHandle(file, 1024 * 1024)

        then:
        result instanceof String
        result == file.toAbsolutePath().toString()
    }

    def 'should return the path when the file looks binary'() {
        given:
        final dir = Files.createTempDirectory('test')
        final file = dir.resolve('data.txt')
        Files.write(file, [0x68, 0x00, 0x69] as byte[])

        when:
        final result = ToolOutputReader.readOrHandle(file, 1024 * 1024)

        then:
        result instanceof String
        result == file.toAbsolutePath().toString()
    }

    def 'should return the path for a file with no extension'() {
        given:
        final dir = Files.createTempDirectory('test')
        final file = dir.resolve('README')
        file.text = 'some readme content'

        when:
        final result = ToolOutputReader.readOrHandle(file, 1024 * 1024)

        then:
        result instanceof String
        result == file.toAbsolutePath().toString()
    }

    def 'should inline a small tsv file content'() {
        given:
        final dir = Files.createTempDirectory('test')
        final file = dir.resolve('report.tsv')
        file.text = 'name\tvalue\nn50\t54321'

        when:
        final result = ToolOutputReader.readOrHandle(file, 1024 * 1024)

        then:
        result instanceof String
        result == 'name\tvalue\nn50\t54321'
    }

    def 'should extract the lowercased extension'() {
        expect:
        ToolOutputReader.extensionOf(Path.of(NAME)) == EXT
        where:
        NAME        | EXT
        'x.JSON'    | 'json'
        'a.b.csv'   | 'csv'
        'README'    | ''
        'a.b/c'     | ''
    }
}

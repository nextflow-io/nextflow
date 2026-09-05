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

package nextflow.util

import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class CsvWriterTest extends Specification {

    def 'should write csv file'() {
        given:
        def file = TestHelper.createInMemTempFile()
        and:
        def records = [
            [id: 1, fastq_1: '1_1.fastq', fastq_2: '1_2.fastq'],
            [id: 2, fastq_1: '2_1.fastq', fastq_2: '2_2.fastq'],
            [id: 3, fastq_1: '3_1.fastq', fastq_2: null],
        ]

        when:
        new CsvWriter([:]).apply(records, file)
        then:
        file.text == '''\
            1,1_1.fastq,1_2.fastq
            2,2_1.fastq,2_2.fastq
            3,3_1.fastq,
            '''.stripIndent()

        when:
        new CsvWriter([header: true]).apply(records, file)
        then:
        file.text == '''\
            id,fastq_1,fastq_2
            1,1_1.fastq,1_2.fastq
            2,2_1.fastq,2_2.fastq
            3,3_1.fastq,
            '''.stripIndent()
    }

    def 'should quote only values that require it'() {
        given:
        def file = TestHelper.createInMemTempFile()
        and:
        def records = [
            ['plain', 'with,comma', 'with "quote"', 'with\nnewline', 'with\rcarriage', '', null]
        ]

        when:
        new CsvWriter([:]).apply(records, file)
        then:
        file.text == 'plain,"with,comma","with ""quote""","with\nnewline","with\rcarriage",,\n'
    }

    def 'should quote headers and honor a custom separator'() {
        given:
        def file = TestHelper.createInMemTempFile()
        and:
        def records = [
            ['sample': 'S1', 'read|group': 'case|control']
        ]

        when:
        new CsvWriter([header: ['sample', 'read|group'], sep: '|']).apply(records, file)
        then:
        file.text == 'sample|"read|group"\nS1|"case|control"\n'
    }

    def 'should write empty csv file'() {
        given:
        def file = TestHelper.createInMemTempFile()
        and:
        def records = []

        when:
        new CsvWriter([:]).apply(records, file)
        then:
        file.text == ''

        when:
        new CsvWriter([header: true]).apply(records, file)
        then:
        file.text == ''
    }

}
